;+ 
; NAME:
; mike_trcflat
;     Version 1.1
;
; PURPOSE:
;  To trace the order edges of the individual orders.  The results are
;  then fed to mike_fittflat to calculate a 2D solution.  The
;  following steps are performed:
;  1.  Identify the order edges (interactive is recommended).
;  2.  Performs an order by order tracing of the order edges using
;      trace_crude.
;  3.  Perform (iteratively) a PCA analysis on the coefficients of the
;      individual traces.
;  4.  Create and write to disk a structure summarizing the fits.
;
;
; CALLING SEQUENCE:
;   
;  mike_trcflat, mike, setup, [side], /CHK, 
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  A structure containing the information for the order by order
;  fits.  This structure is then fed to mike_fittflat to create a 2D
;  solution.  Typical name:  'Flats/TStr_B_01.fits'
;
; OPTIONAL KEYWORDS:
;  /CHK -- Check the order edges interactively (similar to INTER)
;  /INTER -- Determine the order edges interactively
;  SMSHROW -- Row where order edges are determined (default: 1/2 way
;             up the image)
;  THRESH  -- Threshold for an order edge on the red side
;              (default: 100.)
;  /CLOBBER -- Overwrite the output structure
;  P_NSIG  --  Number of sigma significance that an order edge should
;             have for the red side (default: 50.)
;  NSIG  --  Number of sigma significance that an order edge should
;             have fo the blue side (default: 2.0)
;  MINFLAT -- Mininum counts in flat on blue side to use orders (100)
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_trcflat, mike, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;   mike_trcflat_clean
;   mike_getfil
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;   30-Apr-2003 Updated by JXP
;   04-Aug-2003 Updated by JXP ::  Added robust peak cleaning routine
;   15-Mar-2004 Updated by JXP ::  Added interactive mode
;-
;------------------------------------------------------------------------------

;; Cleans out the bad peaks and replaces orders

pro mike_trcflat_clean, side, smsh, cen, ord_sep, NPK=npk, CHK=chk, CBIN=cbin

  ;; Init
  if not keyword_set( NPK ) then npk = n_elements(cen)
  if not keyword_set( CBIN ) then cbin = 2

  ;; Npix
  npix = n_elements(smsh)

  ;; Bin ratio
  rt_bin = 2./float(cbin)
  ;; Mask
  msk = bytarr(npk) + 1B

  ;; Scaled cen
  scl_cen = cen / rt_bin

  ;; Read in fit
  if side EQ 1 then fil_sep = getenv('MIKE_DIR')+'/pro/Flat/blue_sep2x2.idl' $
  else begin
      if not keyword_set( RED ) then $
        fil_sep = getenv('MIKE_DIR')+'/pro/Flat/red_sep2x2B.idl' $
      else stop
  endelse
  restore, fil_sep

  ;; Loop on centroids
  ordsep10 = alog10(ord_sep > 1)
  all_off = ordsep10 - x_calcfit(scl_cen[1:npk-1],fitstr=fit_sep) $
    + alog10(rt_bin)
  off = median(all_off)

  ;; Find first good pair (correct separation)
  for isv=1L,npk-1 do begin
      diff = ordsep10[isv-1] - x_calcfit(scl_cen[isv], $
                                        fitstr=fit_sep) + alog10(rt_bin) - off
      if abs(diff) LT 0.1 then break
  endfor

  if isv NE 1 then begin 
      print, 'mike_trcflat: WARNING -- Trouble at first peak!  ' + $
        'This is the most difficult to handle'
      msk[0:isv-1] = 0B

      ;; Now go back until hitting edge of CCD
      scen = cen[isv]
      cedg = round(scen - 10^(x_calcfit(scen/rt_bin, $
                                        fitstr=fit_sep)+alog10(rt_bin)))
      while(cedg GT 0.) do begin
          ;; Check for good value already existing
          old_cen = where(abs(cen - cedg) LT 3., nold)
          if nold EQ 1 then begin
              cen = [cen, cen[old_cen]]
              msk = [msk,1B]
              scen = (cen[old_cen])[0]
          endif else begin
              ;; Take max flux near the expected spot
              wind = 5L
              mx = max(smsh[(cedg-wind)>0:cedg+wind],imx)
              ;; Add it in
              cen = [cen, cedg-wind+imx]
              msk = [msk,1B]
              scen = cedg-wind+imx
          endelse
          cedg = round(scen - 10^(x_calcfit(scen/rt_bin, $
                                            fitstr=fit_sep)+alog10(rt_bin)))
      endwhile
  endif
      
  for ii=isv,npk-1 do begin
      ;; Check order sep with fit
      diff = ordsep10[ii-1] - x_calcfit(scl_cen[ii], $
                                        fitstr=fit_sep) + alog10(rt_bin) - off
      if diff LT -0.1 then begin ;; False positive
          if not keyword_set(silent) then $
            print, 'mike_trcflat: Found and removed a false positive!', ii
          mn = min(smsh[round([cen[ii-1],cen[ii]])],imn)
          msk[ii+imn-1] = 0B
      endif
      
      ;; Missed a peak?
      if diff GT 0.1 then begin ;; Missed a peak!
          if not keyword_set(silent) then $
            print, 'mike_trcflat: Missed a peak!  Adding it in'
          ;; Find nearest good peak below the spot
          a = where(msk[lindgen(ii)] EQ 1B,na)
          ;; Expected spot
          exp_cen = round(cen[a[na-1]] + $
                          10^(x_calcfit(scl_cen[a[na-1]], $
                                       fitstr=fit_sep)+alog10(rt_bin)))

          ;; Take max flux near the expected spot
          wind = 5L
          mx = max(smsh[0>(exp_cen-wind)<(npix-1): $
                        (exp_cen+wind)<(npix-1)],imx)
          ;; Add it in
          cen = [cen, exp_cen-wind+imx]
          msk = [msk,1B]
      endif
  endfor
  ;; Reset
  cen = cen[where(msk EQ 1B)]
  cen = cen[sort(cen)]
  npk = n_elements(cen)
  ord_sep = (cen - shift(cen,1))[1:npk-1]

  ;; CHK
  if keyword_set( CHK ) then begin
      x_prspeaks, smsh, cen, outmsk, ytwo=smsh[round(cen)], /block
      cen = cen[where(outmsk EQ 1B)]
      ;; Reset to max
      for nn=0L,n_elements(cen)-1 do begin
          ip = round(cen[nn])
          mx = max(smsh[ip-1:(ip+1)<(npix-1)],imx)
          cen[nn] = ip+imx-1
      endfor
  endif

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_trcflat_edges, smsh, side, cbin, ledge, redge, CHK=chk, THIN=thin

  ;; Slit width
  swdth = 20.5 * 2. / cbin

  ;; Npix
  npix = n_elements(smsh)

  ;; Auto Thin
  if cbin GT 4 then thin = 1

  ;; Bin ratio
  rt_bin = 2./float(cbin)

  ;; Read in info
  if side EQ 1 then fil_sep = getenv('MIKE_DIR')+'/pro/Flat/blue_sep2x2.idl' $
  else begin
      fil_sep = getenv('MIKE_DIR')+'/pro/Flat/red_sep2x2B.idl' 
  endelse
  restore, fil_sep

  ;; Rebin slit_profile
  nslit = n_elements(slit_prof)
  slit_prof = congrid(slit_prof, round(nslit*2./cbin))

  ;; Convolve
  conv_saw = convol(smsh, slit_prof)

  ;; Find peaks
  srt = sort(conv_saw)
  thresh = 0.2 * conv_saw[ srt[npix-5] ]  ;; -5 deals with Cosmic rays
  dumi = lindgen(npix)
  ap = where(conv_saw GT thresh AND $
             dumi GT 3 AND dumi LT (npix-3))

  if not keyword_set( THIN ) then begin
      ppk = where( conv_saw[ap] GT conv_saw[ap+1] AND $
                  conv_saw[ap] GT conv_saw[ap+2] AND $
                  conv_saw[ap] GT conv_saw[ap-1] AND $
                  conv_saw[ap] GT conv_saw[ap-2], npk)
  endif else begin
      ppk = where( conv_saw[ap] GT conv_saw[ap+1] AND $
                  conv_saw[ap] GT conv_saw[ap-1], npk)
  endelse
  peak = ap[ppk]

;  x_splot, conv_saw, xtwo=peak, ytwo=conv_saw[peak], /block, psym2=1
; peak = [peak, 450., 728., 844.]
; peak = peak[sort(peak)]

  npk = n_elements(peak)

  ;; Fit the separation
  stp = peak - shift(peak,1)
  fitstr = x_setfitstrct(MAXREJ=5, NORD=1, /FLGREJ, LSIG=2., HSIG=2., NITER=3)
  fit = x_fitrej(peak[1:*], stp[1:*], FITSTR=fitstr, REJPT=rejpt)

  ;; Start on 5th peak (or 7th if 5th rejected)
  rejpt = rejpt + 1L
  a = where(rejpt EQ 5L, na)
  if na NE 0 then begin
      a = where(rejpt EQ 7L, na)
      if na NE 0 then stop else istrt = 7L
  endif else istrt = 5L
  
  ;; Head left
  flg_f = 1
  fin_pk = [peak[istrt]]
  dx = x_calcfit(fin_pk[0], FITSTR=fitstr)
  expx = fin_pk[0] - dx
  while expx GT 2. do begin
      ;; Set tolerance
      if flg_f EQ 1 then begin
          toler = 5. 
          flg_f = 0
      endif else toler = 3.
      ;; Look for peak
      mtch = where( abs(expx-peak) LT toler, nmtch)
      case nmtch of
          0: begin  ;; Create a peak
              i0 = (round(expx)-1.5) > 0L
              i1 = (round(expx)+1.5) < (npix-1)
              mx = max( conv_saw[i0:i1], imx)
              fin_pk = [float(imx+i0), fin_pk]
          end
          1: fin_pk = [peak[mtch], fin_pk]
          2: begin
              mx = max( conv_saw[round(peak[mtch])], imx)
              fin_pk = [peak[mtch[imx]], fin_pk]
          end
      endcase
      ;; New expected
      dx = x_calcfit(fin_pk[0], FITSTR=fitstr)
      expx = fin_pk[0] - dx
  endwhile

  ;; Head right
  flg_f = 1
  nfin = n_elements(fin_pk)
  dx = x_calcfit(fin_pk[nfin-1], FITSTR=fitstr)
  expx = fin_pk[nfin-1] + dx
  while expx LT (npix-3.) do begin
      ;; Set tolerance
      if flg_f EQ 1 then begin
          toler = 5. 
          flg_f = 0
      endif else toler = 3.
      ;; Look for peak
      mtch = where( abs(expx-peak) LT toler, nmtch)
      case nmtch of
          0: begin  ;; Create a peak
              i0 = (round(expx)-1.5) > 0L
              i1 = (round(expx)+1.5) < (npix-1)
              mx = max( conv_saw[i0:i1], imx)
              fin_pk = [fin_pk,float(imx+i0)]
          end
          1: fin_pk = [fin_pk,peak[mtch]]
          2: begin
              mx = max( conv_saw[round(peak[mtch])], imx)
              fin_pk = [fin_pk,peak[mtch[imx]]]
          end
      endcase
      ;; New expected
      nfin = nfin + 1
      dx = x_calcfit(fin_pk[nfin-1], FITSTR=fitstr)
      expx = fin_pk[nfin-1] + dx
  endwhile

  if keyword_set( CHK ) then begin
      x_splot, conv_saw, xtwo=fin_pk, ytwo=conv_saw[round(fin_pk)], $
        /block, psym2=1
  endif

  ;; Analyse slit_prof
  nprof= n_elements(slit_prof)
  cen_sprof = float(nprof-1)/2.
  xd = findgen(nprof)
  mx = max(slit_prof,imx)
  i0 = (imx-3) > 0
  i1 = (imx+3) < (nprof-1)
  ledg_sprof = x_centspln(xd[i0:i1], slit_prof[i0:i1], /force)
  mn = min(slit_prof,imn)
  i0 = (imn-3) > 0
  i1 = (imn+3) < (nprof-1)
  redg_sprof = x_centspln(xd[i0:i1], -1.*slit_prof[i0:i1], /force)
  dl = cen_sprof - ledg_sprof
  dr = redg_sprof - cen_sprof

  ;; Grab left edges
  tmp_le = fin_pk - dl
  if tmp_le[0] LT 2. then tmp_le = tmp_le[1:*]
  ledge = fltarr(n_elements(tmp_le))
  for jj=0L,n_elements(tmp_le)-1 do begin
      i0 = round(tmp_le[jj]-3) > 0
      i1 = round(tmp_le[jj]+3) < (npix-2)
      ledge[jj] = i0 + x_centspln(findgen(i1-i0+1), smsh[i0:i1]>0., /force)
  endfor

  ;; Grab right edges
  tmp_re = fin_pk + dr
  
  if tmp_re[nfin-1] GT npix-2. then tmp_re = tmp_re[0:nfin-2]
  redge = fltarr(n_elements(tmp_re))
  for jj=0L,n_elements(tmp_re)-1 do begin
      i0 = round(tmp_re[jj]-3) > 0
      i1 = round(tmp_re[jj]+3) < (npix-2)
      redge[jj] = i0 + x_centspln(findgen(i1-i0+1), (-1.*smsh[i0:i1])>0., /force)
  endfor

  ;; Try to pick up LHS if possible
  l0 = redge[0] - swdth
  if l0 GT 0. AND abs(ledge[0]-l0) GT 5. then begin
      i0 = (round(l0)-1.5) > 0L
      i1 = (round(l0)+1.5) < (npix-1)
      mx = max( smsh[i0:i1], imx)
      ledge = [float(imx + i0), ledge]
  endif

  ;; Final check
  if keyword_set( CHK ) then begin
      x_splot, smsh, xtwo=ledge, ytwo=smsh[round(ledge)], $
        /block, psym2=1, xthr=redge, ythr=smsh[round(redge)], psym3=2
  endif

  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_trc2d, xcen, mask 

   if NOT keyword_set(sigrej) then sigrej=4.

   ; xcen and xerr come from trace_crude

   ntrace = (size(xcen))[2]
   npix   = (size(xcen))[1]

   good = where(mask,ngood)

   y = (2.*(lindgen(npix,ntrace) mod npix) - npix) / npix 

   nrm0 = 1.0*(ntrace-1)
   t = 2.*(lindgen(npix,ntrace) / npix)/nrm0 - 1  

   work2d = dblarr(npix*ntrace,13)
   worky = flegendre(y[*], 5)
   workt = flegendre(t[*], 5)

   work2d[*,0:4] = workt[*,0:4]
   work2d[*,5:8] = worky[*,1:4] 
   work2d[*,9:12] = worky[*,1:4] * (workt[*,1] # replicate(1,4))

   work2dgood = work2d[good,*]
   alpha = transpose(work2dgood) # work2dgood
   beta = transpose(work2dgood) # (xcen[good])
   choldc, alpha, p, /double
   res = cholsol(alpha,p,beta)
   surffit = xcen * 0.0
   surffit[*] = work2d # res
 
return, surffit
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro mike_trcflat, mike, setup, side, CHK=chk, SMSHROW=smshrow, $
                  CLOBBER=clobber, P_NSIG=p_nsig, THRESH=thresh, DEBUG=debug, $
                  INTER=inter, NSIG=nsig, MINFLAT=minflat, NOBASIS=nobasis, $
                  NOSTOP=nostop

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_trcflat, mike, setup, [side], P_NSIG=, THRESH=, /CLOBBER,'
      print, '           /CHK, /INTER, /DEBUG, NSIG=, MINFLAT= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  if not keyword_set( P_NSIG) then p_nsig = 50.
  if not keyword_set( NSIG) then nsig = 5.0
  if not keyword_set( THRESH) then thresh = 100.
  if not keyword_set( BGUESS) then BGUESS = 75L
  if not keyword_set( NCOEFF) then ncoeff = 6L
  if not keyword_set( MINFLAT) then minflat = 100.


  x_psclose
  !p.multi=[0,1,1]

; Loop on side

  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      ;; SIDE
      if qq EQ 1 then print, 'mike_trcflat: Tracing BLUE trace flat' $
      else print, 'mike_trcflat: Tracing RED trace flat'

      ;; Check output file
      trc_fil = mike_getfil('tflat_str', setup, SIDE=qq, /name, CHKFIL=chkf)
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'mike_trcflat: Trace exists, moving on..'
          continue
      endif

      ;; Read in flat 
      flat = mike_getfil('tflat_fil', setup, SIDE=qq)
      sz = size(flat, /dimensions)

      ;; Set binning
      cbin = round(2048. / sz[0])
      rbin = round(4096. / sz[1])
      
      ;; Disallow trimming
;;      if 2048L mod sz[0] NE 0 or 4096L mod sz[1] NE 0 then begin
;;          print, 'mike_trcflat:  No trimming allowed!!'
;;          stop
;;      endif

      ;; Create tmp_xcen
      tmp_xcen = fltarr(sz[1],1000L)
      tmp_xerr = fltarr(sz[1],1000L)
      tmp_xfit = fltarr(sz[1],1000L)
      ntmp = 0L
  
      ;; Sawtooth
      sf =  shift(flat ,1) 
      sb =  shift(flat,-1) 

      ;; the denominator was biasing the centroid with this earlier definition:
      ;;    saw =  (sf-sb)/(abs(sf+sb) - 2.0*min(flat) + 10)
      ;; 
      ;; Let's normalize with this function
      ;; sawnorm = djs_median(flat, 1) ## replicate(1,sz[0])
      ;; saw =  (sf - sb)/(sawnorm + (sawnorm EQ 0)) * (sawnorm GT 0)
      ;;
      ;; Let's just subtract the two, otherwise we might bias the centroid

      ;; ipos = 1 are the LHS of each order; ipos=2 are the RHS

      saw =  -1.*(sf - flat)
      sawivar = 1/(flat + sf + 100.0)

      ;; Smash
      if not keyword_set( smshrow ) then smshrow = sz[1]/2 else begin
          if abs(smshrow - sz[1]/2) GT 100 THEN begin
              print, 'mike_trcflat:  Continue at your own risk!'
              stop
          endif
      endelse
      print, 'mike_trcflat: Using row ',strtrim(smshrow,2), $
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

      mike_trcflat_edges, smsh_search, qq, cbin, ledge, redge, CHK=chk

      for ipos=1,2 do begin
          
          case ipos of 
              1: cen = ledge
              2: cen = redge
          endcase

          cen = cen[sort(cen)]
          npk = n_elements(cen)
          
          ;; Check for next neighbors!  (bug:  JXP 3/16/04)
          for i=1L, n_elements(cen)-2 do begin
              if abs(cen[i]-cen[i-1]) LT 5. $
                OR abs(cen[i+1]-cen[i]) LT 5. then begin
                  print, 'mike_trcflat: Problem..'
                  stop
              endif
          endfor
        
          
          ;; Trace from middle
          if ipos EQ 1 then $
            print, 'mike_trcflat: Tracing [negative peaks] from the middle...' $
          else print, 'mike_trcflat: Tracing [positive peaks] from the middle...' 
          if ipos EQ 2 then begin
              ;; Require no RHS exists to the left of lowest LHS !
              if min(cen) LT min(svcen) then begin
                  cen = cen[1:*]
                  print, 'mike_trcflat:  Tossed out RHS of first order'
              endif
          endif
          npk = n_elements(cen)

          xstart = cen
          ystart = replicate(smshrow, n_elements(cen))
          if ipos EQ 1 then tsaw = saw > 0. else tsaw = (-1.*saw) > 0.
          tsaw[0,*] = 0
          tsaw[sz[0]-1,*] = 0
          
          
          for j=0L,9 do $
            xstart = trace_fweight(tsaw, xstart, ystart, radius=2.)
          
          xcen = trace_crude(tsaw, yset=ycen ,XSTART=xstart,$
                             radius=3.5, ystart=smshrow, xerr=xerr , $
                             maxshifte=0.2)
          
          ;; Change to x_fweight in the next release!
;          x2cen = x_fweight(tsaw, xcen, ycen, invvar=sawivar, $
          x2cen = mike_fweight(tsaw, xcen, ycen, invvar=sawivar, $
                               radius=2, xerr=x2err,sig=0.2) 
          
          outmask=0

          ;; JXP -- Kludge to fix the first 2 orders on red side
;;                  (10feb04) [x2cen requirement]
          x2ivar = (x2err LT 1.0 AND x2cen GT 10)/x2err^2
          
          xy2traceset, ycen, x2cen, tset, invvar=x2ivar, $
            yfit=x2fit, ncoeff=ncoeff, maxdev=1.0,outmask=outmask2
          
          ;; Mask out reddest orders on blue side
          if qq EQ 1 then msktrc = [0,1,2] else begin
              if keyword_set( MSKTRC ) then delvarx, msktrc
          endelse

          ;; PCA
          if not keyword_set( NOBASIS ) then $
;          if not keyword_set( NOBASIS ) and qq EQ 1 then $
            ;; Change to x_basis in the next release!
;            x3fit = x_basis(tset, x2cen, outmask2, $
            x3fit = mike_basis(tset, x2cen, outmask2, $
                               ncoeff=ncoeff,eigenvec=eigenvec, $
                               msktrc=msktrc, outstr=pca_out) $
          else begin
              x3fit = x2fit
              pca_out = 0L
          endelse

          ;; Recentroid
          x4cen = x3fit
          for j=1L,4 do $
            x4cen = trace_gweight(tsaw, x4cen, ycen, sigma=1)
          
          ;; Final fit
          x4flux = extract_boxcar(tsaw, x4cen, radius=3)
          x4ivar = (x4flux>0) + 1.
          x4err = 1./sqrt(x4ivar)
          
          xy2traceset, ycen, x4cen, tset4, invvar=x4ivar, $
            yfit=x4fit, ncoeff=ncoeff, maxdev=1.0,outmask=outmask4

          x4err = x4err + ((outmask4 EQ 0) OR (x4err LE 0))*900.

          ;; Save for QA
          sv_qa = { $
                    pca_out: pca_out, $
                    x4cen: x4cen, $
                    xfit: x4fit, $
                    ycen: ycen,  $
                    npk: npk,  $
                    msk: outmask2 $
                  }
          if ipos EQ 1 then sv_qa1 = sv_qa else sv_qa2 = sv_qa

          
          ;; Save into tmp arrays
          if ipos EQ 1 then begin  ; LHS
              ntmp = (size(x4cen))[2]
              nneg = ntmp
              xcen_neg = x4cen - 0.5
              xerr_neg = x4err
              xfit_neg = x4fit - 0.5
              tmp_xcen[*,0:ntmp-1] = xcen_neg
              tmp_xerr[*,0:ntmp-1] = xerr_neg
              tmp_xfit[*,0:ntmp-1] = xfit_neg
              svcen = cen
          endif else begin
              npos = (size(x4cen))[2]
              xcen_pos = x4cen - 0.5
              xerr_pos = x4err
              xfit_pos = x4fit - 0.5
              tmp_xcen[*,ntmp:ntmp+npos-1] = xcen_pos
              tmp_xerr[*,ntmp:ntmp+npos-1] = xerr_pos 
              tmp_xfit[*,ntmp:ntmp+npos-1] = xfit_pos
              ntmp = ntmp + npos
          endelse
      endfor

      ;; Require significant flux in bluest blue orders
      if qq EQ 1 then begin 
          flg_min = 0
          kmin = 5L
          for kk=kmin,nneg-1 do begin
              ;; Grab RHS
              rhs = where(tmp_xcen[sz[1]/2,nneg:ntmp-1] GT tmp_xcen[sz[1]/2,kk],nrhs)
              if nrhs EQ 0 then continue
              ;; Measure flux
              fx_flat = median(flat[tmp_xcen[sz[1]/2,kk]: $
                                    tmp_xcen[sz[1]/2,rhs[0]+nneg],sz[1]/2])
              ;; Measure scattered light in between orders
              flg_sctt = 1
              case kk of
                  kmin: scatt = median(flat[tmp_xcen[sz[1]/2,rhs[0]+nneg]:$
                                         tmp_xcen[sz[1]/2,kk+1],sz[1]/2])
                  (nneg-1): scatt = median(flat[tmp_xcen[sz[1]/2,rhs[0]+nneg-1]:$
                                         tmp_xcen[sz[1]/2,kk],sz[1]/2])
                  else: begin
                      x0 = tmp_xcen[sz[1]/2,rhs[0]+nneg]
                      x1 = tmp_xcen[sz[1]/2,kk+1]
                      if x1 - x0 LT 2 then begin
                          x1 = x0 + 2
                          flg_sctt = 0
                      endif
                      scatt1 = median(flat[x0:x1,sz[1]/2])
                      x0 = tmp_xcen[sz[1]/2,rhs[0]+nneg-1]
                      x1 = tmp_xcen[sz[1]/2,kk]
                      if x1 - x0 LT 2 then begin
                          x1 = x0 + 2
                          flg_sctt = 0
                      endif
                      scatt2 = median(flat[x0:x1,sz[1]/2])
                      scatt = (scatt1+scatt2) / 2.
                  end
              endcase
              if flg_sctt NE 0 then sv_scatt = scatt else begin
                  if not keyword_set(sv_scatt) then begin
                      print, 'mike_trcflat:  Something is probably wrong with'+$
                        ' your order traces.  Not enouch blue counts?'
                      stop
                  endif else scatt = sv_scatt
              endelse
;              print, 'mike_trcflat: Testing counts -- ', kk, fx_flat, scatt
              ;; Test
              if (fx_flat-scatt) LT MINFLAT then begin
                  flg_min = flg_min + 1
                  if flg_min EQ 1 then begin
                      kksv = kk
                      rhsv = rhs[0]
                  endif
              endif
          endfor
          if flg_min GE 2 then begin
              print, 'mike_trcflat:  Counts have dipped below 100.  The '+$
                'code is going to dump those orders!!: '  , kksv, nneg
              print, 'mike_trcflat:  If you wish to recover them, contact the'+$
                ' desginers of the pipeline.'
              print, 'mike_trcflat:  Otherwise, continue on..'
              if not keyword_set( NOSTOP ) then stop
              tmp = nneg
              nneg = kksv
              npos = rhsv
              tmp_xcen[*,nneg:nneg+npos-1] = tmp_xcen[*,tmp:tmp+npos-1]
              tmp_xerr[*,nneg:nneg+npos-1] = tmp_xerr[*,tmp:tmp+npos-1]
              tmp_xfit[*,nneg:nneg+npos-1] = tmp_xfit[*,tmp:tmp+npos-1]
              ntmp = npos+nneg
              ;; QA
              sv_qa1.msk[*,nneg:*] = 0B
              sv_qa2.msk[*,npos:*] = 0B
          endif
      endif
    
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; QA
      qafil = mike_getfil('qa_trcflat', setup, SIDE=qq)
      print, 'mike_trcflat: QA file is '+qafil
      x_psopen, qafil, /maxs
      clr = getcolor(/load)

      if not keyword_set(NOBASIS) AND keyword_set(sv_qa1.pca_out) then begin
          !p.multi=[0,ncoeff/2 + (ncoeff mod 2),2,0,1]
          
          ;; LHS ;;
          bad = where(sv_qa1.pca_out.x0msk EQ 0B, nbad)
          plot, findgen(sv_qa1.npk), sv_qa1.pca_out.x0, psym=1, $
              color=clr.black, background=clr.white, charsize=1.5
          if nbad NE 0 then oplot, [bad], [sv_qa1.pca_out.x0[bad]], $
              psym=2, color=clr.red
             oplot, findgen(sv_qa1.npk), sv_qa1.pca_out.x0fit, color=clr.blue
          mn = min(sv_qa1.pca_out.x0fit, max=mx)
          xyouts, 1., mn + (mx-mn)*0.9, 'x0', charsize=1.5
          
          ;; PCA
          for ii=0L,ncoeff-2 do begin
              if ii GT 1 then begin
                median_hidden = median(sv_qa1.pca_out.hidden[ii,*])
                plot, sv_qa1.pca_out.usetrc, sv_qa1.pca_out.hidden[ii,*], $
                  color=clr.black, background=clr.white, $
                  xrange=[0L,sv_qa1.npk-1], psym=1, charsize=1.5, $
                  yrange=median_hidden+[-1.,1.],/ys
              endif else $
                plot, sv_qa1.pca_out.usetrc, sv_qa1.pca_out.hidden[ii,*], $
                  color=clr.black, background=clr.white, $
                  xrange=[0L,sv_qa1.npk-1], psym=1, charsize=1.5, /ys

              oplot, findgen(sv_qa1.npk), sv_qa1.pca_out.high_fit[*,ii], $
                  color=clr.blue
              ;; Rej points
              case ii of 
                  0L: begin
                   if sv_qa1.pca_out.rejpt.rej0[0] NE -1 then $
                    oplot, [sv_qa1.pca_out.usetrc[sv_qa1.pca_out.rejpt.rej0]], $
                        [sv_qa1.pca_out.hidden[ii,sv_qa1.pca_out.rejpt.rej0]], $
                        color=clr.red, psym=2
                  end
                  1L: begin
                      if sv_qa1.pca_out.rejpt.rej1[0] NE -1 then $
                        oplot, [sv_qa1.pca_out.usetrc[sv_qa1.pca_out.rejpt.rej1]], $
                        [sv_qa1.pca_out.hidden[ii,sv_qa1.pca_out.rejpt.rej1]], $
                        color=clr.red, psym=2
                  end
                  else: 
              endcase
          ;; Label
              mn = min(sv_qa1.pca_out.hidden[ii,*], max=mx)
              xyouts, 1., mn + (mx-mn)*0.9, 'PCA'+strtrim(ii,2), charsize=1.5
          endfor

          ;; RHS ;;
          bad = where(sv_qa2.pca_out.x0msk EQ 0B, nbad)
          plot, findgen(sv_qa2.npk), sv_qa2.pca_out.x0, psym=1, color=clr.black, $
            background=clr.white, charsize=1.5
          if nbad NE 0 then oplot, bad, sv_qa2.pca_out.x0[bad], psym=2, color=clr.red
          oplot, findgen(sv_qa2.npk), sv_qa2.pca_out.x0fit, color=clr.blue
          mn = min(sv_qa2.pca_out.x0fit, max=mx)
          xyouts, 1., mn + (mx-mn)*0.9, 'x0', charsize=1.5
          
          ;; PCA
          for ii=0L,ncoeff-2 do begin
              plot, sv_qa2.pca_out.usetrc, sv_qa2.pca_out.hidden[ii,*], $
                color=clr.black, background=clr.white, $
                xrange=[0L,sv_qa2.npk-1], psym=1, charsize=1.5
              oplot, findgen(sv_qa2.npk), sv_qa2.pca_out.high_fit[*,ii], $
                color=clr.blue
              ;; Rej points
              case ii of 
                  0L: begin
                      if sv_qa2.pca_out.rejpt.rej0[0] NE -1 then $
                        oplot, [sv_qa2.pca_out.usetrc[sv_qa2.pca_out.rejpt.rej0]], $
                        [sv_qa2.pca_out.hidden[ii,sv_qa2.pca_out.rejpt.rej0]], $
                        color=clr.red, psym=2
                  end
                  1L: begin
                      if sv_qa2.pca_out.rejpt.rej1[0] NE -1 then $
                        oplot, [sv_qa2.pca_out.usetrc[sv_qa2.pca_out.rejpt.rej1]], $
                        [sv_qa2.pca_out.hidden[ii,sv_qa2.pca_out.rejpt.rej1]], $
                        color=clr.red, psym=2
                  end
                  else: 
              endcase
              ;; Label
              mn = min(sv_qa2.pca_out.hidden[ii,*], max=mx)
              xyouts, 1., mn + (mx-mn)*0.9, 'PCA'+strtrim(ii,2), charsize=1.5
          endfor
      endif

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Traces 
      ;; LHS
      !p.multi=[0,1,1]
      tt = tmp_xfit[*,0:nneg-1]+0.5  ;; Offset
;      tt = sv_qa1.xfit[*,0:nneg-1]+0.5  ;; Offset
      gd = where(sv_qa1.msk EQ 1B, complement=bad)
      plot, sv_qa1.ycen[gd], 10*(sv_qa1.x4cen[gd]-tt[gd])+tt[gd], $
        psym=3, color=clr.black, $
        background=clr.white, charsize=1.5, xrange=[0., sz[1]], yrange=[0.,sz[0]]
      !p.thick=1
      if bad[0] NE -1 then $
      oplot, sv_qa1.ycen[bad], 10*(sv_qa1.x4cen[bad]-tt[bad])+tt[bad], $
        psym=3, color=clr.red
      for ii=0L,nneg-1 do begin
          oplot, findgen(sz[1]), tt[*,ii], color=clr.green
      endfor
      if not keyword_set(NOBASIS) and keyword_set(sv_qa1.pca_out) then $
        xyouts, 0.5, 0.96, $
        'Reduced chi^2 = '+string(sv_qa1.pca_out.red_chi2,format='(f7.4)'), $
        color=clr.black, /normal, charsize=2.5, alignment=0.5

      ;; RHS
      tt = tmp_xfit[*,nneg:ntmp-1]+0.5  ;; Offset
;      tt = sv_qa2.xfit[*,0:npos-1]+0.5  ;; Offset
      gd = where(sv_qa2.msk EQ 1B, complement=bad)
      plot, sv_qa2.ycen[gd], 10*(sv_qa2.x4cen[gd]-tt[gd])+tt[gd], $
        psym=3, color=clr.black, $
        background=clr.white, charsize=1.5, xrange=[0., sz[1]], yrange=[0.,sz[0]]
      !p.thick=1
      if bad[0] NE -1 then $
      oplot, sv_qa2.ycen[bad], 10*(sv_qa2.x4cen[bad]-tt[bad])+tt[bad], $
        psym=3, color=clr.red
      for ii=0L,npos-1 do begin
          oplot, findgen(sz[1]), tt[*,ii], color=clr.green
      endfor
      if not keyword_set(NOBASIS) and keyword_set(sv_qa1.pca_out) then $
        xyouts, 0.5, 0.96, $
        'Reduced chi^2 = '+string(sv_qa2.pca_out.red_chi2,format='(f7.4)'), $
        color=clr.black, /normal, charsize=2.5, alignment=0.5

      x_psclose
      !p.multi=[0,1,1]
      spawn, 'gzip -f '+qafil

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ACCOUNTING

      ;; Only ready for red side right now
      if qq EQ 2 then begin
          fitordname = getenv('MIKE_DIR')+'/pro/Flat/model_R.fits'
          x_fitstrtofits, fitordstrct, fitordname, /REVERSE
          
          ord_sep = (cen - shift(cen,1))[1:*]
          ord_sep = ord_sep*(cbin)/2.0
          
          ;; 
          orders = x_calcfit(ord_sep,FITSTR=fitordstrct) - 1.
          orders = round(orders)
          
          ;; Find offset
          offset = round(median(orders[0:6] + lindgen(7)))
          orders = offset - lindgen(n_elements(orders))
          
          start_guess_order = orders[0]
;          start_guess_order = round((order_coeff[0] - 1.665)/0.0103 + 76.05)
          print, 'Guess starting order number traced is ', start_guess_order
      endif else start_guess_order = BGUESS
      
      ;; Determine slit length
      npk = nneg
      swid_val = fltarr(npk)
      for q=0L,npk-1 do begin
          diff = xcen_pos[smshrow,*] - xcen_neg[smshrow,q]
          a = where(diff GT 0,na)
          if na NE 0 then mn = diff[a[0]] else mn = 0.
          swid_val[q] = mn
      endfor

      swidth = median(swid_val)
      print, 'mike_trcflat: Adopting order width of ', swidth, ' binned pixels'

      ;; Fill up the output
      trc_str = { $
                  xcen: tmp_xcen[*,0:ntmp-1],$
                  xerr: tmp_xerr[*,0:ntmp-1],$
                  xfit: tmp_xfit[*,0:ntmp-1],$
                  flg: intarr(ntmp), $
                  nordr: npk, $
                  smshrow: smshrow, $
                  swidth: swidth, $
                  orderstrt: start_guess_order $
                }

      ;; Zero out bad few columns of red side:: 
      ;; Is this setup dependent?!  (i.e. the dichroic)
      ;; If so, it will need to be a function of order #
      
      if qq EQ 2 then begin
          bd = where(trc_str.xcen LT 10., nbd)
          if nbd NE 0 then trc_str.xerr[bd] = 99.
      endif

      ;; Zero out orders 75-77 on blue side.  
      ;; They are too tightly packed or it is the dichroic 'echo'
;      if qq EQ 1 then begin
;          stop
;          if trc_str.orderstrt LE 78 then begin
;              nbd = 78-trc_str.orderstrt+1
;              for i=0L,nbd-1 do trc_str.xerr[*,i] = 99.
;              for i=0L,nbd-1 do trc_str.xerr[*,i+nneg] = 99.
;          endif
;      endif
      
      ;; Fill up flg
      trc_str.flg[0:nneg-1] = 1 ; LHS of slit
      trc_str.flg[nneg:ntmp-1] = 2 ; RHS of slit
      
      ;; Write 
      mwrfits, trc_str, trc_fil, /create

  endfor

  ;; 
  print, 'mike_trcflat: You may now wish to check the results with mike_chktrcflat'
  print, 'mike_trcflat: All done!'

  return
end
