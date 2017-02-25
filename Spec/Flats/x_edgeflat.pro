;+ 
; NAME:
; x_edgeflat
;     Version 1.1
;
; PURPOSE:
;  To trace the order edges of the individual orders.  And the results are
;  fit (which replaces x_fittflat).  The fit is a low order fit to the slit
;  width, and then a higher order fit to the order centers.
;  1.  Identify the order edges (interactive is recommended).
;  2.  Performs an order by order tracing of the order edges using
;      trace_crude.
;  3.  Perform (iteratively) a PCA analysis on the coefficients of the
;      individual traces.
;  4.  Create and write to disk a structure summarizing the fits.
;  5.  Do a 2-part fit that replaces the fittflat call.
;
;
; CALLING SEQUENCE:
;   
;  x_edgeflat, flat, bin, trc_fil, ordr_fil, /CHK
;
; INPUTS:
;  flat  -- Flat image
;  bin   -- Binning in spatial direction?
;
; RETURNS:
;
; OUTPUTS:
;  trc_fil -- A structure containing the information for the order by order
;  fits.  This structure is then fed to mike_fittflat to create a 2D
;  solution.  Typical name:  'Flats/TStr_B_01.fits'
;  ordr_fil -- Order structure, e.g. 'Flats/OStr_B_01.fits'
;
; OPTIONAL KEYWORDS:
;  /ZERO_OUT  -- 
;  FLATIVAR=  -- Inverse variance array of the flat
;  /CHK -- Check the order edges interactively 
;  /CLOBBER -- Overwrite the output structure
;  NCOEFF  -- Number of Legendre coefficients for tracing [default: 6]
;  SMSHROW=  -- Row at which to 'smash' the spectrum to then look for
;               peaks
;  NOOVERLAP= -- If orders are closer than NOOVERLAP, then they are not
;               included in the PCA fit but are extrapolated from the
;               PCA analysis.
;  KEEP_FRAC -- This value sets the minimum fractional amount that an
;               order can be on the CCD to be included in the PCA
;               [Default: 0.9]
;  MNSEP    -- Minimum separation between orders
;  EXTRAP   -- Number of orders to extrapolate for order definition (e.g. [1,1])
;  SZ_EXTRAP -- Fraction of slit size the extrapolated order must
;               occupy to be included [Default: 0.25]
;  WIDEN    -- Increase the slit width by WIDEN*2. pixels (both edges)
;
; OPTIONAL OUTPUTS:
;  QAFIL -- Filename for QA output
;
; COMMENTS:
;
; EXAMPLES:
;   x_edgeflat, flat, bin, trc_fil, ordr_fil
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   1-Nov-2004  Written by SB, replaces x_trcflat and x_fittflat
;                QA file is still in mike_trcflat
;   1-May-2005  Polished to extrapolate orders, avoid orders with
;                significant overlap
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro x_edgeflat, flat, bin, trc_fil, ordr_fil, CHK=chk, $
                CLOBBER=clobber, DEBUG=debug, TNSIG=tnsig, $
                MNSEP=mnsep, N_X0FIT=n_x0fit, WIDEN=widen, $
                ncoeff=ncoeff, $
                NOSTOP=nostop, SMSHROW=smshrow, ZERO_OUT=zero_out, $
                FLATIVAR=flativar, QAFIL=qafil, SGUESS=sguess, EXTRAP=extrap, $
                POLY_N=poly_n, N_XPFIT=n_xpfit, EORDR=EORDR, $
                NOOVERLAP=nooverlap, LOGFIL=logfil, KEEP_FRAC=keep_frac, $
                _EXTRA=EXTRA

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'x_edgeflat, flat, bin, trc_fil, ordr_fil, ' + $
        'THRESH=, /CLOBBER,'
      print, '           /CHK, /INTER, /DEBUG, NSIG=, MINFLAT=, /ZERO_OUT, NOOVERLAP=,'
      print, '           WIDEN= '
      print, '    /EXTRAP [v1.2]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( POLY_N ) then poly_n=3
  if not keyword_set( WIDEN ) then widen = 0.
  if not keyword_set( OFFSET ) then offset = 0.5
  if not keyword_set( NCOEFF) then ncoeff = 6L
  if not keyword_set( BIN) then bin = 1
  if not keyword_set(NITER) then niter = 2L
  if not keyword_set(NY_DIFF) then ny_diff = 2
  if not keyword_set(NT_DIFF) then nt_diff = 3
  if not keyword_set(KEEP_FRAC) then keep_frac = 0.9
  if not keyword_set(SZ_EXTRAP) then sz_extrap = 0.25
  if not keyword_set(N_X0FIT) then N_X0FIT = 4  ;; Order number for x0 fit

  x_psclose
  !p.multi=[0,1,1]


  if not keyword_set(FLATIVAR) then stop
  flativar = flativar*(smooth(1.0*(flativar EQ 0),3) EQ 0)
  sz = size(flat, /dimensions)

  ntmp = 0L
  
  ;; Sawtooth
  sf =  shift(flat ,1) 
  sb =  shift(flat,-1) 
  
  saw =  -1.*(sf - flat)
  sawivar = flativar/2.0

  ;; Smash
  if not keyword_set( smshrow ) then smshrow = sz[1]/2 else begin
      if abs(smshrow - sz[1]/2) GT 100 THEN begin
          print, 'x_trcflat:  Continue at your own risk!'
          stop
      endif
  endelse
  print, 'x_edgeflat: Using row ',strtrim(smshrow,2), $
    ' to identify orders'
  smsh = djs_median(saw[*,smshrow-15:smshrow+15],2)

  ;; Smooth
  nfilter = long(sz[0]/20.)
  smsh_denom = djs_median(convol(flat[*,smshrow-15:smshrow+15],  $
                                 replicate(1./nfilter,nfilter), $
                                 /edge_truncate) + 1,2) 
  smsh[0] = 0
  smsh[sz[0]-1]= 0
  smsh_search = smooth((smsh/(smsh_denom > 1)),3) 
;  smsh_search = smsh
  
  x_trcflat_edges, smsh_search, qq, bin, ledge, redge, CHK=chk, $
    NSIG=TNSIG, SORDR=sordr, _EXTRA=extra, MNSEP=mnsep, EORDR=eordr

  ledge = ledge[sort(redge)]
  redge = redge[sort(redge)]

  npk = n_elements(redge)
  if redge[0] LT ledge[0] then begin
      width = (redge[1:*] - ledge)
      width_p = ladfit(findgen(npk-1),width)
      
      ledge = [redge[0]-width_p[0] + width_p[1], ledge[0:npk-2]]
  endif

          
  ystart = replicate(smshrow, npk)
  lsaw = saw > 0. 
  lsaw[0,*] = 0
  lsaw[sz[0]-1,*] = 0
  
  rsaw = (-1.*saw) > 0.
  rsaw[0,*] = 0
  rsaw[sz[0]-1,*] = 0
  
  lstart = ledge
  rstart = redge
  for j=0L,9 do begin 
      lstart = trace_fweight(smsh>0, lstart, ystart*0, radius=2.) 
      rstart = trace_fweight((-1.0*smsh)>0,rstart,ystart*0, radius=2.) 
  endfor
  
  lcen = trace_crude(lsaw, yset=ly,XSTART=lstart,$
                     radius=3.5, ystart=smshrow, xerr=lerr , $
                     maxshifte=0.2)
  rcen = trace_crude(rsaw, yset=ry,XSTART=rstart,$
                     radius=3.5, ystart=smshrow, xerr=rerr , $
                     maxshifte=0.2)

  all_lcen = lcen
  all_rcen = rcen

  ;;;;;;;;;;;;;;;;
  ;;  Deal with overlapping orders (do not fit but extrapolate)
  svnpk = npk
  nbad = 0
  if keyword_set(NOOVERLAP) then begin
      sep = djs_median( shift(all_lcen,0,-1) - all_rcen, 1)
      sep[n_elements(sep)-1] = sep[n_elements(sep)-2] 
      nor = n_elements(sep)
      ov_msk = replicate(1B, nor)
      if median(sep[0:nor/2]) GT median(sep[nor/2:*]) then begin
          bad = where( sep[nor/3:*] LT NOOVERLAP, nbad)
          if nbad NE 0 then begin
              ov_msk[nor/3 + bad[0]:*] = 0B
              ;; Actual ignored orders
              bdd = where(ov_msk EQ 0B, nbad,ncomplement=npk)
              all_ordr = [lindgen(npk), npk+lindgen(nbad)]
              ;; flag
              flg_over = 1
              print, 'x_edgeflat:  Order overlap!!  Will extrapolate ', $
                     nbad, ' orders of', nor, ' total'
          endif else all_ordr = lindgen(nor)
      endif else begin
          bad = where( sep[0:2*nor/3] LT NOOVERLAP, nbad) 
          if nbad NE 0 then begin
              ov_msk[0:bad[nbad-1]] = 0B
              bdd = where(ov_msk EQ 0B, nbad,ncomplement=npk)
              all_ordr = [lindgen(nbad)-nbad, lindgen(npk)]
              ;; flag
              flg_over = 2
              print, 'x_edgeflat:  Order overlap!!  Will extrapolate ', $
                     nbad, ' orders of',  nor, ' total'
          endif else all_ordr = lindgen(nor)
      endelse
      lcen = lcen[*,where(ov_msk)]
      rcen = rcen[*,where(ov_msk)]
      lerr = lerr[*,where(ov_msk)]
      rerr = rerr[*,where(ov_msk)]
      ;; Reset npk
      npk = n_elements(where(ov_msk))
  endif else begin
      nbad = 0
      ov_msk = replicate(1B,npk)
      all_ordr = lindgen(npk)
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Only use orders with KEEP_FRAC on the chip
  if KEEP_FRAC GT 0. then begin
      print, 'x_edgeflat: KEEP_FRAC', keep_frac
      cen = (lcen + rcen) / 2.
      szc = size(cen, /dimensions)
      gdo = where(ov_msk,ngdo)
      mskgd = replicate(1B,szc[1])
      for qq=0L,szc[1]-1 do begin
          gd = where(cen[*,qq] GT 5. AND $
                     cen[*,qq] LT sz[0]-5. AND $
                     lerr[*,qq] LT 1. AND $
                     rerr[*,qq] LT 1., kpnum)
          if float(kpnum)/float(sz[1]) LT KEEP_FRAC then begin
              ;; Reject
              print, 'x_edgeflat: Rejecting order ', all_ordr[gdo[qq]], $
                ' because it is partial.'
              print, 'x_edgeflat:     Will extrapolate it..'
              mskgd[qq] = 0B
              ov_msk[gdo[qq]] = 0B
              if size( EDGORD, /type ) EQ 0 then $
                edgord = [all_ordr[gdo[qq]]] else $
                edgord = [edgord,all_ordr[gdo[qq]]]
          endif
      endfor
      ;; Reset lcen, rcen
      lcen = lcen[*,where(mskgd)]
      rcen = rcen[*,where(mskgd)]
      rej = where(ov_msk EQ 0B, ncomplement=npk)
   endif

  FOR iiter=1,niter do begin
      l2cen = x_fweight(lsaw, lcen, ly, invvar=sawivar, $
                        radius=2, xerr=l2err,sig=0.2) 
      
      r2cen = x_fweight(rsaw, rcen, ry, invvar=sawivar, $
                        radius=2, xerr=r2err,sig=0.2) 
      
      
      slit_center = 0.5*[r2cen+l2cen]
      center_ivar = 4.0*(r2err LT 10 AND l2err LT 10)/(l2err^2 + r2err^2)
      
      x_xy2traceset, ly, slit_center, tset_center, invvar=center_ivar, $
        maxdev=0.5, maxrej=1, maxiter=40, /groupbadpix, $
        yfit=center_fit, ncoeff=ncoeff, outmask=outmask2, /silent
      outmask2 = (abs(slit_center - center_fit) LT 0.5) AND (center_ivar GT 0)
      
      nrow = (size(ly))[1] 
      msktrc = where(total(outmask2,1) LT nrow*0.5, complement=goodpks )
      
      
      ;; PCA
      center_fit_rel = x_basis(tset_center, slit_center, outmask2, $
                               ncoeff=ncoeff,eigenvec=eigenvec, poly_n=poly_n, $
                               msktrc=msktrc, outstr=pca_out, /skipx0, n_xpfit=n_xpfit) 
      
      ;; Did PCA work?
      if not keyword_set(PCA_OUT) then begin
          print, 'x_edgeflat:  PCA failed.  Assuming too few orders (<7)'
          print, 'x_edgeflat:  Continue if this is the case and the code'
          print, '               will go ahead without PCA'
          print, 'x_edgeflat:  Note, the code will NOT extrapolate. 
          stop
          flg_pca = 0
          x_xy2traceset, ly, l2cen, tset_center, invvar=center_ivar, $
                         maxdev=0.5, maxrej=1, maxiter=40, /groupbadpix, $
                         yfit=lcen, ncoeff=ncoeff, outmask=outmask2, /silent
          x_xy2traceset, ry, r2cen, tset_center, invvar=center_ivar, $
                         maxdev=0.5, maxrej=1, maxiter=40, /groupbadpix, $
                         yfit=rcen, ncoeff=ncoeff, outmask=outmask2, /silent
          pca_out = 0
          diff_fit = (rcen - lcen) 
          break
      endif else begin
          est_center = pca_out.x0
          x0msk = (est_center NE 0.0)
          flg_pca = 1
      
          ;; Allow for Edges
          xx = all_ordr[where(ov_msk,nxx)]
          if nxx NE npk then stop
          
          for jj=0,npk do begin 
              aset = func_fit(xx, est_center, N_X0FIT, invvar = x0msk, yfit=x0fit, $
                              function_name='poly') 
              qdone = djs_reject(est_center, x0fit, outmask=x0msk, inmask=x0msk,$
                                 maxdev=1.0, maxrej=1)
              if qdone EQ 1 then break 
          endfor
          ;;
          ;;  Make up a center_set and diff_set
          ;; 
          
          avg_fit = x_basis(tset_center, slit_center, outmask2, $
                            ncoeff=ncoeff,eigenvec=eigenvec, poly_n=poly_n, $
                            msktrc=msktrc, outstr=pca_out, x0in=x0fit[*], n_xpfit=n_xpfit) 
          
          
          diff = (r2cen - l2cen) 
          differr = 900.0*(l2err GT 10 OR r2err GT 10)+sqrt(l2err^2 + r2err^2)
          t= 0
          x_trace2d, diff, differr, diff_fit, nycoeff=ny_diff, ntcoeff=nt_diff, $
                     res=diffres, t=t, mask=diffmask, sigrej=10.
          ;; 
          lcen = avg_fit-0.5*diff_fit     
          rcen = avg_fit+0.5*diff_fit     
      endelse
       
    endfor 

  ;; Extrapolate beyond smash row
    if keyword_set( EXTRAP ) and flg_pca then begin
        if EXTRAP[0] GT 0 then begin
            print, 'x_edgeflat:  Extrapolating ', extrap[0], ' orders on low side'
            ov_msk = [replicate(0B,extrap[0]),ov_msk]
            all_ordr = [-1*extrap[0] + lindgen(extrap[0]), all_ordr]
        endif
        if EXTRAP[1] GT 0 then begin
            print, 'x_edgeflat:  Extrapolating ', extrap[1], ' orders on high side'
            ov_msk = [ov_msk,replicate(0B,extrap[1])]
            all_ordr = [all_ordr, lindgen(extrap[1]) + svnpk]
        endif
     endif

    ;;;;;;;;;;;;;;
    ;; Extrapolate as necessary
    nordr = n_elements(all_ordr)
    if nordr GT n_elements(x0fit) AND flg_pca then begin
        i_ex = where(ov_msk EQ 0B, nextrap, complement=gd)
        ;i_ex = 26+lindgen(7)
        ;nextrap=7
        if nextrap EQ 0 then stop

        x0_ex = poly(all_ordr[i_ex], aset)

        ;;  Center of order
        high_matr = fltarr(nextrap,ncoeff-1)
        high_matr[*,0] = poly(x0_ex, pca_out.coeff0)
        high_matr[*,1] = poly(x0_ex, pca_out.coeff1)
        high_matr[*,2:*] = replicate(1.,nextrap) # (pca_out.high_fit[0,2:*])[*]
        extrap_fit = eigenvec # transpose(high_matr) + x0_ex ## replicate(1.,nrow)

        ;; Difference
        nrm0 = npk-1.
        t = (2*all_ordr[i_ex]/nrm0 - 1) ## replicate(1,nrow)
        y = ((2*dindgen(nrow) - nrow) / nrow) # replicate(1,nextrap)
        work2d = dblarr(nrow*nextrap,ny_diff*nt_diff)
        worky = flegendre(y[*], ny_diff)
        workt = flegendre(t[*], nt_diff)
        for i=0,nt_diff-1 do begin
            for j=0,ny_diff-1 do begin
                work2d[*,j*nt_diff+i] = worky[*, j] * workt[*,i]
            endfor
        endfor
        extrap_diff = dblarr(nrow,nextrap)
        extrap_diff[*] = work2d # diffres
        ;; EXTRAP
        extrap_lcen = extrap_fit-0.5*extrap_diff
        extrap_rcen = extrap_fit+0.5*extrap_diff
    endif 
    
    x2cen = [[l2cen],[r2cen]]
    x2fit = [[lcen], [rcen]]
    x2err = [[l2err],[r2err]]
    

    left_chi2 = total((l2cen - lcen)^2/l2err^2 * (l2err LT 10.)) / $
                   total(l2err LT 10)
    right_chi2 = total((r2cen - rcen)^2/r2err^2 * (r2err LT 10.)) / $
                   total(r2err LT 10)
    total_chi2 =     total((x2cen-x2fit)^2/[x2err^2] * (x2err LT 10.)) / $
                   total(x2err LT 10.)

    splog, 'Left edge reduced chi^2 ',  left_chi2, filename=logfil, /appen, $
      /close
    splog, 'Right edge reduced chi^2',  right_chi2, filename=logfil, /appen, $
      /close
    splog, 'Total reduced chi^2     ',  total_chi2, filename=logfil, /appen, $
      /close

    ;; Save for QA
    sv_qa1 = { $
               pca_out: pca_out, $
               x2cen: x2cen, $
               x2fit: x2fit, $
               ycen: ly,  $
               npk: npk,  $
               msk: outmask2 $
             }
    
    
    nneg = npk
    npos = npk
    tmp_xcen = x2cen - 0.5
    tmp_xerr=  [[l2err],[r2err]]
    tmp_xfit = x2fit - 0.5
    ntmp = npos+nneg
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; QA
    if keyword_set( QAFIL ) then begin
        print, 'x_edgeflat: QA file is '+qafil
        x_psopen, qafil, /maxs
        clr = getcolor(/load)
        
        if keyword_set(sv_qa1.pca_out) then begin
            !p.multi=[0,ncoeff/2 + (ncoeff mod 2),2,0,1]
            
            in = lindgen(npk) 
            bad = where(sv_qa1.pca_out.x0msk EQ 0B, nbadd)
            plot, in, sv_qa1.pca_out.x0, psym=1, $
              color=clr.black, background=clr.white, charsize=1.5, /yno
            if nbadd NE 0 then oplot, [in[bad]], [sv_qa1.pca_out.x0[bad]], $
              psym=2, color=clr.red
            oplot, sv_qa1.pca_out.x0fit[*], color=clr.blue
            mn = min(sv_qa1.pca_out.x0fit, max=mx)
            xyouts, 1., mn + (mx-mn)*0.9, 'x0', charsize=1.5
            
            ;; PCA
            for ii=0L,ncoeff-2 do begin
                plot, sv_qa1.pca_out.usetrc, sv_qa1.pca_out.hidden[ii,*], $
                  color=clr.black, background=clr.white, $
                  xrange=[0L,sv_qa1.npk-1], psym=1, charsize=1.5, /yno
                oplot, sv_qa1.pca_out.high_fit[*,ii], color=clr.blue
                ;; Rej points
                case ii of 
                    0L: begin
                        if sv_qa1.pca_out.rejpt.rej0[0] NE -1 then $
                          oplot, [sv_qa1.pca_out.usetrc[sv_qa1.pca_out.rejpt.rej0]] mod npk, $
                          [sv_qa1.pca_out.hidden[ii,sv_qa1.pca_out.rejpt.rej0]], $
                          color=clr.red, psym=2
                    end
                    1L: begin
                        if sv_qa1.pca_out.rejpt.rej1[0] NE -1 then $
                          oplot, [sv_qa1.pca_out.usetrc[sv_qa1.pca_out.rejpt.rej1]] mod npk, $
                          [sv_qa1.pca_out.hidden[ii,sv_qa1.pca_out.rejpt.rej1]], $
                          color=clr.red, psym=2
                    end
                    else: 
                endcase
                ;; Label
                mn = min(sv_qa1.pca_out.hidden[ii,*], max=mx)
                xyouts, 1., mn + (mx-mn)*0.9, 'PCA'+strtrim(ii,2), charsize=1.5
            endfor
        endif


       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Left edge Traces 
        ;;
        !p.multi=[0,1,1]
        tt = tmp_xfit+0.5  ;; Offset
        gd = where(sv_qa1.msk EQ 1B, complement=bad)
        sd = n_elements(sv_qa1.msk)
        
        plot, sv_qa1.ycen[gd], 10*(sv_qa1.x2cen[gd]-tt[gd])+tt[gd], $
              psym=3, color=clr.black, background=clr.white, charsize=1.5, $
              xrange=[0., sz[1]], yrange=[0.,sz[0]], /xs, /ys
        !p.thick=1
        if bad[0] NE -1 then $
          oplot, sv_qa1.ycen[bad], 10*(sv_qa1.x2cen[bad]-tt[bad])+tt[bad], $
                 psym=3, color=clr.red
        for ii=0L,npk-1 do begin
            oplot, findgen(sz[1]), tt[*,ii], color=clr.green
        endfor
        xyouts, 0.5, 0.96, $
                'Left edge Reduced chi^2 = '+string(left_chi2,format='(f8.4)'), $
                color=clr.black, /normal, charsize=2.5, alignment=0.5
        
        ;; Extrapolate 
        if keyword_set( EXTRAP ) AND keyword_set(NEXTRAP) then begin
            for ii=0L,extrap[0]-1 do $
                   oplot, findgen(sz[1]), extrap_lcen[*,extrap[0]-ii-1], $
                          color=clr.blue
            for ii=nextrap-extrap[1],nextrap-1 do $
                   oplot, findgen(sz[1]), extrap_lcen[*,ii], color=clr.blue
        endif
        
        ;; Rej orders
        if size( EDGORD, /type ) NE 0 then begin
            nedg = n_elements(edgord)
            if not keyword_set(EXTRAP) then dex = [0,0L] else dex=extrap
            for ii=0L,nedg-1 do begin
                mtch = where(all_ordr[i_ex] EQ edgord[ii])
                oplot, findgen(sz[1]), extrap_lcen[*,mtch[0]], $
                       color=clr.skyblue
            endfor
        endif
        
        ;; OVERLAP 
        if keyword_set(NOOVERLAP) AND (nbad NE 0) then begin
            if not keyword_set(EXTRAP) then dex = [0,0L] else dex=extrap
            for ii=dex[0],nbad+dex[0]-1 do begin
                ;; Avoid rej
                if size(edgord, /type) NE 0 then begin
                    b = where(ii EQ (dex[0]+edgord), nb)
                    if nb then continue
                endif
                oplot, findgen(sz[1]), extrap_lcen[*,ii], color=clr.orange
            endfor
        endif
        
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Right edge traces 
        ;;
        plot, sv_qa1.ycen[gd], 10*(sv_qa1.x2cen[gd+sd]-tt[gd+sd])+tt[gd+sd], $
              psym=3, color=clr.black, background=clr.white, charsize=1.5, $
              xrange=[0., sz[1]], yrange=[0.,sz[0]], /xs, /ys
        !p.thick=1
        if bad[0] NE -1 then oplot, sv_qa1.ycen[bad], $
          10*(sv_qa1.x2cen[bad+sd]-tt[bad+sd])+tt[bad+sd], $
          psym=3, color=clr.red
        for ii=npk,2*npk-1 do begin
            oplot, findgen(sz[1]), tt[*,ii], color=clr.green
        endfor
        xyouts, 0.5, 0.96, $
                'Right edge reduced chi^2 = '+string(right_chi2,format='(f8.4)'), $
                color=clr.black, /normal, charsize=2.5, alignment=0.5
        
        ;; Extrapolate
        if keyword_set( EXTRAP ) and keyword_set(NEXTRAP) then begin
            for ii=0L,extrap[0]-1 do $
                   oplot, findgen(sz[1]), extrap_rcen[*,extrap[0]-ii-1], $
                          color=clr.blue
            for ii=nextrap-extrap[1],nextrap-1 do $
                   oplot, findgen(sz[1]), extrap_rcen[*,ii], color=clr.blue
        endif
        
        ;; Rej orders
        if size( EDGORD, /type ) NE 0 then begin
            nedg = n_elements(edgord)
            if not keyword_set(EXTRAP) then dex = [0,0L] else dex=extrap
            for ii=0L,nedg-1 do begin
                mtch = where(all_ordr[i_ex] EQ edgord[ii])
                oplot, findgen(sz[1]), extrap_rcen[*,mtch[0]], $
                       color=clr.skyblue
            endfor
        endif
        
        ;; OVERLAP
        if keyword_set(NOOVERLAP) AND (nbad NE 0) then begin
            if not keyword_set(EXTRAP) then dex = [0,0L] else dex=extrap
            for ii=dex[0],nbad+dex[0]-1 do begin
                ;; Avoid rej
                if size(edgord, /type) NE 0 then begin
                    b = where(ii EQ (dex[0]+edgord), nb)
                    if nb then continue
                endif
                oplot, findgen(sz[1]), extrap_rcen[*,ii], color=clr.orange
            endfor
        endif
        
        x_psclose
        !p.multi=[0,1,1]
        replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
        ps_replacetitle, replace_title, qafil

        spawn, 'gzip -f '+qafil
    endif
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ACCOUNTING
    
    if not keyword_set(SGUESS) then sguess = 0L
    
    swidth = median(diff_fit)
    
    print, 'x_edgeflat: Adopting order width of ', swidth, ' binned pixels'
    
    ;; Fill up the output
    trc_str = { $
                xcen: tmp_xcen,$
                xerr: tmp_xerr,$
                xfit: tmp_xfit,$
                flg: intarr(ntmp), $
                nordr: npk, $
                smshrow: smshrow, $
                swidth: swidth, $
                orderstrt: SGUESS $
              }
    
    ;; Zero out bad few columns of red side:: 
    ;; Is this setup dependent?!  (i.e. the dichroic)
    ;; If so, it will need to be a function of order #
    
    if keyword_set(ZERO_OUT) then begin
        bd = where(trc_str.xcen LT 10., nbd)
        if nbd NE 0 then trc_str.xerr[bd] = 99.
    endif
    
    ;; Fill up flg
    trc_str.flg[0:nneg-1] = 1   ; LHS of slit
    trc_str.flg[nneg:ntmp-1] = 2 ; RHS of slit
    
    ;; Write 
    mwrfits, trc_str, trc_fil, /create
    spawn, 'gzip -f '+trc_fil

;
;     Now write out the tflat fit
;
    tmp = { $
            order: 0L, $        ;  Physical order #
            flg_anly: 0, $      ;  1=ok, 0=ng
            xcen: 0., $         ;  Center of order at ycen
            ycen: 0., $         ;  Smash row
            ymin: 0L, $         ;  Minimum used in y direction (-> -1)
            ymax: nrow, $       ;  Maximum used in y direction (-> +1)
            lcen: 0.d, $        ;  Wavelength correspondoning to xcen,ycen
            lhedg: fltarr(nrow), $ ; LH edge
            rhedg: fltarr(nrow),  $ ; RH edge
            arc_m: fltarr(nrow), $
            profile0 : fltarr(251),$ ; Slit profile
            profile1 : fltarr(251),$ ; Slit profile
            lhc : fltarr(ncoeff),$ ; LH coefficients
            rhc : fltarr(ncoeff) $ ; RH coefficients
          }
    
    smrow = trc_str.smshrow
    
    ;; EXTRAPOLATE
    nordr = n_elements(all_ordr)

    ord_str = replicate(tmp, nordr)
    ord_str.ycen = float(smrow)
        
    if flg_pca then begin
        allx = x0fit 
        nfit = n_elements(x0fit)
        if keyword_set(nextrap) then begin
            if nextrap NE 0 then allx = [allx,x0_ex]
        endif
        srt = sort(allx)
    endif else begin
        srt = lindgen(nordr)    ;  JXP -- This might make for a bug
        nfit = nordr
    endelse
        
    msk_ordr = replicate(1B, nordr)
    for ii=0L,nordr-1 do begin
        idx = srt[ii]
        ;; Fill in edges
        if idx LT nfit then begin ;; FIT
            ord_str[ii].lhedg = lcen[*,idx] - offset - WIDEN
            ord_str[ii].rhedg = rcen[*,idx] - offset + WIDEN
            ord_str[ii].xcen = (ord_str[ii].lhedg[smrow]+ord_str[ii].rhedg[smrow])/2.
        endif else begin  ;; EXTRAP
            jj = where(abs(allx[idx]-x0_ex) LT 0.1)
            jj = jj[0]
            ord_str[ii].lhedg = extrap_lcen[*,jj] - offset - WIDEN
            ord_str[ii].rhedg = extrap_rcen[*,jj] - offset + WIDEN 
            ord_str[ii].xcen = (ord_str[ii].lhedg[smrow]+ord_str[ii].rhedg[smrow])/2.
        endelse
        ;; Sanity check
        cen = (ord_str[ii].rhedg+ord_str[ii].lhedg)/2.
        gd = where(cen GT 5. AND cen LT sz[0]-5, ngd)
        if ngd LT sz[1]*SZ_EXTRAP then begin
            msk_ordr[ii] = 0B
            print, 'x_edgeflat:  Dropping order due to short size', ii
        endif
    endfor
    keep = where(msk_ordr)
    ord_str = ord_str[keep]
    
;        ord_str.lhedg = trc_str.xfit[*,0:npk-1]
;        ord_str.rhedg = trc_str.xfit[*,npk:*]
;        ord_str.xcen = (ord_str.lhedg[smrow]+ord_str.rhedg[smrow,*])/2.
        ;;
        ;; just do a fit
    ;; 
    ly = ly[*,0] # replicate(1.,nordr)
    ry = ry[*,0] # replicate(1.,nordr)
    xy2traceset, ly, ord_str.lhedg, lset, ncoeff=ncoeff, /silent
    xy2traceset, ry, ord_str.rhedg, rset, ncoeff=ncoeff, /silent
        
        
    ord_str.lhc = lset.coeff
    ord_str.rhc = rset.coeff
        
    
    ord_str.order = lindgen(n_elements(ord_str)) + 1000
    
    ;; Write to fits
    print, 'x_edgeflat: Writing order structure ', ordr_fil
    mwrfits, ord_str, ordr_fil, /create
    
;; 
    print, 'x_edgeflat: You may now wish to check the results with x_chktrcflat'
    print, 'x_edgeflat: All done!'
    
  return
end
