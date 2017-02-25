;+ 
; NAME:
; x_fntobj   
;     Version 1.1
;
; PURPOSE:
;    Identify the object within the slit (single object) and trace it
;    along each order.  The code rectifies each order, collapses it
;    and looks for the object.  IF it finds it in >7 orders then it
;    peaks up on the object.  Otherwise, it guesses the object is at
;    the center of the order.  It then uses xy2tracset (trace_crude)
;    to trace the object along each order individually.  Finally, it
;    performs a PCA analysis on the trace_crude coefficients to create
;    a smoothed (pseudo-2D) solution for the trace.  The code outputs
;    an object structure which includes the trace and other info.
;
;
; CALLING SEQUENCE:
; x_fntobj, objstr, ordr_str, img, ivar, qafil
;
; INPUTS:
; objstr -- Object structure
; ordr_str -- Order structure defining the echelle footprint
; img    -- Image array
; ivar   -- Inverse variance array
;
; RETURNS:
;
; OUTPUTS:
;  Object structure including the trace and other inrormation
;  regarding the object.  This structure is then filled up with the 1D
;  spectrum, etc.
;
; OPTIONAL KEYWORDS:
;  /STD     - Find object for standard star
;  /CHK     - Show diagnostics of the code performance
;  /NOCLOB  - Do not clobber existing files (default: clobber)
;  /DEBUG   - Debug
;  OBJAPER  - Set aperture to mask for sky subtraction (default: 20
;             unbinned pixels for obj, 26 for STD)
;  FWIDTH   - Fraction of width used for radius in tracing (default:
;             1/4 slit width)
;  POLY_NCOEFF - Order to fit PCA0 with in x_basis [default: set by number of
;                good traces]
;  MIN_PCA  - Minimum number of orders for PCA trace [default: 10]
;
; OPTIONAL OUTPUTS:
;  QAFIL=  -- Filename of QA
;
; COMMENTS:
;
; EXAMPLES:
;   x_fntobj
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Sep-2003 Written by JXP (combined aspects of fndobj with trcobj)
;   10-Nov-2003 Changed fntobj_off, SMB
;-
;------------------------------------------------------------------------------

function m_fntobj_off, ordr_str, img, qq, DEBUG=debug, xerr=xerr, $
                          CHK=chk

  ;; Rectify the image
  rect_image = x_ordrectify(img, ordr_str[qq].lhedg, ordr_str[qq].rhedg)

  ;; Smooth
  sz = size(img,/dimen)
  vertical_median_smooth = 63*sz[1]/4096L
  half_smooth = (vertical_median_smooth+1)/2
  nf = sz[1] / half_smooth*half_smooth
  rect_rebin = djs_median(reform(rect_image[*,0:nf-1], $
              (size(rect_image))[1], half_smooth , sz[1]/half_smooth), 2)
  rect_smooth = x_medianrow(transpose(rect_image), vertical_median_smooth,$
                   reflect=1)
  ;; Trace down the image
  down1= trace_crude(rect_rebin, xstart=(size(rect_image))[1]/2, $
                maxshift0=5,maxshifte=0.3)
  down = x_fweight(rect_rebin, down1, niter=10)

  smash_image= total(rect_smooth,1) / (size(rect_image))[2] 
  smash_sz = n_elements(smash_image)
                                ;This is only used for error
                                ;estimates...because center is
                                ;redefined on the next step
  center = find_npeaks(smash_image, nfind=1, width=3L, ypeak=peak, $
                 xerr=peak_err)
  center= median(down)
  
  if keyword_set( CHK ) then begin
      plot, smash_image, /yno, psym = 10, ytitle='Smashed Image - Counts'
      wait, 0.2
  endif
  slitfrac = center/(smash_sz-1)

  ;; Fraction
  frac = center - fix(center)
  value = convol(djs_median(smash_image,width=3,boundary='reflect'), $
                 [1.0-frac,frac],/edge_wrap)
  if keyword_set( CHK ) then djs_oplot, value, color = 'red'
  pixl  = fix(center) - 2 > 0
  pixmid = fix(center) + 1 < (smash_sz-1)
  pixr = fix(center) + 4 < (smash_sz-1)
  
  yl    = value[pixl]
  yr    = value[pixr]
  ymid  = value[pixmid]
  
  xerr = 2.0*sqrt(((yl> 0) +  (yr > 0) + 1.))/ ((2.0*ymid - yl -yr) > 0.1) / $
    sqrt((size(img,/dimen))[1]) 
  

  if keyword_set(DEBUG) then print, center, slitfrac, peak, xerr, peak_err
  xerr = xerr * (peak_err GT 0)
  return, slitfrac
  
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro m_fntobj, objstr, ordr_str, img, ivar, qafil, $
              CHK=chk, NOCLOB=noclob, OBJAPER=objaper, FWIDTH=fwidth, $
              DEBUG=debug, STD=std, MINPCA=minpca, FCOEFF=fcoeff, $
              FULLTHR=fullthr, NONLIN=nonlin, POLY_NCOEFF=poly_ncoeff

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'x_fntobj, objstr, ordr_str, img, ivar /CHK, OBJAPER= [v1.1]'
      return
  endif 
  
;  if keyword_set(STD) then minpca=20L

  ;;  Optional Keywords
  if not keyword_set(OBJAPER) then begin
     if not keyword_set(STD) then objaper = [0.75, 0.75] $
     else objaper=[-0.6,0.6]
  endif
;  if not keyword_set(OBJAPER) then objaper = [10.,10.]
  if not keyword_set(NOCEFF) then ncoeff = 7
  if not keyword_set(FCOEFF) then fcoeff = 15
  if not keyword_set(MINPCA) then minpca = 10L
  if not keyword_set(MINERR) then minerr = 0.1
  if not keyword_set(FWIDTH) then fwidth = 0.25
  if not keyword_set(FULLTHR) then fullthr = 1000.
  sz_img = size(img,/dimensions)

  ;; FIND OBJ
  print, 'x_fntobj: Finding the object...'
  nordr = n_elements(ordr_str)
  frac_array = fltarr(nordr)
  xerr_array = fltarr(nordr)
  for qq = 0,nordr-1 do begin  
      frac_array[qq] = m_fntobj_off(ordr_str, img, qq, xerr=xerr_temp, $
                                    debug=debug, CHK=chk) 
      xerr_array[qq] = xerr_temp 
  endfor
  good = where(xerr_array GT 0 AND xerr_array LT 0.4, ngood)
      
  if ngood GT 7 then begin
      lad_coeff = ladfit(good, frac_array[good])
      lad_fit = poly(findgen(nordr),lad_coeff)
      ;; frac_guess = (lad_fit + median(frac_array[good]))/2.
      frac_guess = lad_fit 
  endif else $
    if ngood GE 1 then frac_guess=fltarr(nordr) + median(frac_array[good]) $
  else                 frac_guess=fltarr(nordr) + 0.5
  
  ordr_width = (ordr_str.rhedg-ordr_str.lhedg)
  trace_guess = ordr_width * (frac_guess##replicate(1,sz_img[1])) + $
    ordr_str.lhedg
  
  if keyword_set(CHK) then begin
      plot, frac_array, xtitle='Trace Number', ytitle='slit fraction', $
        xmargin=[13,1], xstyle=1
      clr = getcolor(/load)
      if good[0] NE -1 then oplot, [good], [frac_array[good]],ps=1 
      oplot, frac_guess, color=clr.red
  endif
  frac_good = good              ;; Save for QA
          
  
  ;; TRACE
  print, 'x_fntobj: Tracing...'
  yset = findgen(sz_img[1]) # replicate(1,nordr)
  radius = fix(median(ordr_width)*fwidth) +1

  x2 = x_fweight(img, trace_guess, yset, $; invvar=ivar, $
                 niter=10, xerr=xerr, radius=radius, sig=0.7)

  ;; Min err  -- JXP 10/2/05
;  for qq=0L,nordr-1 do begin
;      gd = where(xerr[*,qq] LT 900 AND xerr[*,qq] GT 0.,ngd)
;      if ngd NE 0 then begin
;          djs_iterstat, xerr[gd,qq], median=mede, sigma=sige
;          xerr[gd,qq] = (mede - sige) > xerr[gd,qq] < (mede+sige)
;      endif
;      xerr[*,qq] = xerr[*,qq] > MINERR  ;; Min error in trace
;  endfor
;  x_splot, xerr[*,1], /blo, psym1=1

  ;; Too many data points with xerr=999 (red side)

  shift = total((x2-trace_guess)/xerr^2,1)/total(1./xerr^2,1)
  xivar = (xerr LT 900)/(xerr^2 + 0.01^2)
  xy2traceset, yset, x2, tset, invvar=xivar, ncoeff=ncoeff, yfit=x2fit, $
    outmask=outmask, upper=10, lower=10

  ;;  Find eigenvectors.....
;  print, 'MIKE_FNTOBJ:  PCA...'
  if ngood GT MINPCA AND not keyword_set(STD) then begin
      x3fit = x_basis(tset, x2, outmask, ncoeff=ncoeff, eigenvec=eigenvec,$
                      msktrc=msktrc, outstr=pca_out, NONLIN=nonlin, $
                      POLY_NCOEFF=poly_ncoeff)
  endif else begin
      print, 'x_fntobj: Not enough Good fits for PCA'
      x3fit = x2fit
  endelse
  
  if keyword_set( DEBUG ) then begin
      tmp = img
      for qq=0L,nordr-1 do begin
          rnd_trc = round(x3fit[*,qq]) 
          plot_this = where(rnd_trc GE 0 AND rnd_trc LT sz_img[0])
          if plot_this[0] NE -1 then $
            tmp[rnd_trc[plot_this], plot_this] = -10000
      endfor
;          print, objstr[mm].xcen, objstr[mm].ycen
      xatv, tmp, /block, /histeq
      stop
  endif
          

  if n_elements(x3fit) EQ 1 then begin
      ;; time for plan B  use order center
      print, 'MIKE_FNTOBJ: Using order centers as guide for trace'
      
      cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2. 
      xy2traceset, yset, cen_ref, ordrcen_set, ncoeff=ncoeff, $
        yfit=x2fit, outmask=outmask
      
      x3fit = x_basis(ordrcen_set, x2, outmask, $
                      ncoeff=ncoeff, eigenvec=eigenvec, $
                      msktrc=msktrc, outstr=pca_out)
      
  endif
;
;     final tweak to traces...
;
  if x3fit[0] NE -1L then begin
      
      ;; KLUDGE!
      x4 = x_fweight(img, x3fit, yset, $; invvar=ivar, $
                     niter=5, xerr=xerr4, radius=radius/1.5, sig=0.7)
      
      x_ivar = 1./(xerr4 + (xerr4 EQ 0))^2 * (xerr4 GT 0 AND xerr4 LT 1.5) * $
        (x4 GT radius) * (x4 LT sz_img[0]-radius)
      ;; rebin into 30 sections x nord
      nbin = long(sz_img[1] / 31)
      lo = long(0.5*nbin)
      x_ivar_check = total(reform(x_ivar[lo:lo+30*nbin-1,*],nbin, 30, nordr),1)
      fullfit_fib = where(total(x_ivar_check GT FULLTHR,1) EQ 30, nfullfit)
      print, 'x_fntobj: There are ', strtrim(nfullfit,2), $
        " fibers which don't require PCA interpolation"
      
      
      if keyword_set( DEBUG ) then begin
          tmp = img
          for qq=0L,nordr-1 do begin
              rnd_trc = round(x4[*,qq]) 
              plot_this = where(rnd_trc GE 0 AND rnd_trc LT sz_img[0])
              if plot_this[0] NE -1 then $
                tmp[rnd_trc[plot_this], plot_this] = -10000
          endfor
;          print, objstr[mm].xcen, objstr[mm].ycen
          xatv, tmp, /block, /histeq
          stop
      endif
          
      xy2traceset, yset, 1.0d*(x4-x3fit), tset4, invvar=(xerr4 LT 1.5), $
        ncoeff=5, yfit=x4fit, outmask=outmask
      xfinal_fit = x3fit 
      addin = where(total(x4fit^2,1) LT 4*sz_img[1])
      xfinal_fit[*,addin] =  x3fit[*,addin] + x4fit[*,addin]
      
      ;; Final full fit on the highest SNR orders
      if nfullfit NE 0 then begin
          xy2traceset, yset, x4, tset4full, $
            ncoeff=FCOEFF, yfit=x4full_fit, outmask=outmask
          xfinal_fit[*,fullfit_fib] = x4full_fit[*,fullfit_fib]
      endif
      
  endif else begin
      xfinal_fit = x2fit
  endelse


  ;; Fill up the Object Structure
  nordr = n_elements(ordr_str)
  objstr[0:nordr-1].trace[0:sz_img[1]-1] = xfinal_fit
  objstr[0:nordr-1].aper = objaper # replicate(1,nordr)
  objstr[0:nordr-1].obj_id = 'a'
  objstr[0:nordr-1].order = ordr_str.order
  medx = sz_img[1]/2.
  objstr[0:nordr-1].xcen = medx
  
  if x3fit[0] EQ -1L then begin
      ;; No basis
      objstr[0:nordr-1].ycen = (x2fit[medx,*] - trace_guess[medx,*])[*]
  endif else begin
      objstr[0:nordr-1].ycen = (x3fit[medx,*] - trace_guess[medx,*])[*]
  endelse

  ;; 
  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      for qq=0L,nordr-1 do begin
          rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
          plot_this = where(rnd_trc GE 0 AND rnd_trc LT sz_img[0])
          if plot_this[0] NE -1 then $
            tmp[rnd_trc[plot_this], plot_this] = -10000
      endfor
;          print, objstr[mm].xcen, objstr[mm].ycen
      xatv, tmp, /block, /histeq, min=-50, max=500.
  endif
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
  ;; QA
  if keyword_set(QAFIL) then begin
      print, 'x_fntobj: QA file in ', qafil+'.gz'
      x_psopen, qafil, /maxs
      clr = getcolor(/load)
      !p.multi=[0,1,1]
      
      ;; Slit fraction
      plot, frac_array, xtitle='Trace Number', ytitle='slit fraction', $
        color=clr.black, background=clr.white, xmargin=[13,1], xstyle=1
      if frac_good[0] NE -1 then oplot, [frac_good], [frac_array[frac_good]],ps=1, color=clr.red
      oplot, frac_guess, linestyle=2
      
      !p.multi=[0,ncoeff/2 + (ncoeff mod 2),2,0,1]
      ;; PCA ;;
      if keyword_set( pca_out ) then begin
          npk = n_elements(pca_out.x0)
          bad = where(pca_out.x0msk EQ 0B, nbad)
          plot, findgen(npk), pca_out.x0, psym=1, color=clr.black, $
            background=clr.white, charsize=1.5
          if nbad NE 0 then oplot, [bad], [pca_out.x0[bad]], psym=2, color=clr.red
          oplot, findgen(npk), pca_out.x0fit, color=clr.blue
          mn = min(pca_out.x0fit, max=mx)
          xyouts, 1., mn + (mx-mn)*0.9, 'x0', charsize=1.5
          
          ;; PCA
          for ii=0L,ncoeff-2 do begin
              plot, pca_out.usetrc, pca_out.hidden[ii,*], $
                color=clr.black, background=clr.white, $
                xrange=[0L,npk-1], psym=1, charsize=1.5, /yno
              oplot, findgen(npk), pca_out.high_fit[*,ii], $
                color=clr.blue
              ;; Rej points
              case ii of 
                  0L: begin
                      if pca_out.rejpt.rej0[0] NE -1 then $
                        oplot, [pca_out.usetrc[pca_out.rejpt.rej0]], $
                        [pca_out.hidden[ii,pca_out.rejpt.rej0]], $
                        color=clr.red, psym=2
                  end
                  1L: begin
                      if pca_out.rejpt.rej1[0] NE -1 then $
                        oplot, [pca_out.usetrc[pca_out.rejpt.rej1]], $
                        [pca_out.hidden[ii,pca_out.rejpt.rej1]], $
                        color=clr.red, psym=2
                  end
                  else: 
              endcase
              ;; Label
              mn = min(pca_out.hidden[ii,*], max=mx)
              xyouts, 1., mn + (mx-mn)*0.9, 'PCA'+strtrim(ii,2), charsize=1.5
          end
      endif

      !p.multi=[0,1,1]
      ;; Points
      if keyword_set(x4) then xcen = x4 else xcen = x2
      ycen = findgen(sz_img[1]) # replicate(1.,nordr)
      gd = where(outmask EQ 1B, complement=bad, ncomplement=nbad)

      plot, ycen[gd], 10*(xcen[gd]-xfinal_fit[gd])+xfinal_fit[gd], $
        psym=3, color=clr.black, thick=5, xstyle=1, ystyle=1, $
        background=clr.white, charsize=1.5, xrange=[0., sz_img[1]], $
        yrange=[0.,sz_img[0]], /nodata
      !p.thick=1
      oplot, ycen[gd], 5*(xcen[gd]-xfinal_fit[gd])+xfinal_fit[gd], $
        color=clr.black, psym=3
      ;; Rejected
      if nbad NE 0 then $
        oplot, ycen[bad], (xcen[bad]-xfinal_fit[bad])+xfinal_fit[bad], $
        psym=3, color=clr.red
      ;; Fit
      for ii=0L,nordr-1 do begin
          oplot, findgen(sz_img[1]), xfinal_fit[*,ii], color=clr.blue
      endfor
      ;; Label
      xyouts, 0.5, 0.96, 'OBJ FIT (5x resid)', color=clr.blue, $
        charsize=1.5, /normal
      
      ;; Close
      x_psclose
      !p.multi=[0,1,1]
      replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
      ps_replacetitle, replace_title, qafil

      spawn, 'gzip -f '+qafil
      
      ;; Clean up
      if keyword_set(x4) then delvarx, x4
      if keyword_set(pca_out) then delvarx, pca_out
      
  endif
  
;  DONE
  print, 'x_fntobj: All done! '
  return
end


