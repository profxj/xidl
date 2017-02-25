;+ 
; NAME:
; x_extobjbox   
;    Version 1.0
;
; PURPOSE:
;    Given the flux and wave images and the starting position of the
;       object to extract (x,y),  do a boxcar extraction
;
; CALLING SEQUENCE:
;   
;  x_extobjbox, fx, wv, xyguess, fin_spec, MSK=, VAR=, WVMNX=
;
; INPUTS:
;  fx  -- Image array
;  wv  -- Wavelength array
;  xyguess -- Guess at the trace
;
; RETURNS:
;
; OUTPUTS:
;   fin_spec -- Final spectrum
;
; OPTIONAL KEYWORDS:
;   WVMNX - Endpoints for extraction (default: [3200., 11000.])
;   COLLMNX - Endpoints for collapsing the spectrum prior to extraction
;              (default: [3400., 8000.])
;   REDBLUE - Spectrum runs from red to blue
;   NEWWV -  Finaly 1D wavelength array desired
;   CRVAL1 -  Starting wavelength of final 1D array
;   CDELT -   Delta lambda (lor or linear)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_extobjbox, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_extobjbox, fx, wv, xyguess, fin_spec, MSK=msk, VAR=var, WVMNX=wvmnx, $
                 NMED=nmed, PIX=pix, FRAC=frac, RADIUS=radius, $
                 APER=aper, TOT_TRC=tot_trc, CRVAL1=crval1, CDELT=cdelt, $
                 NPIX=npix, COLLMNX=collmnx, DEBUG=debug, REDBLUE=redblue, $
                 REJ_CR=rej_cr, MXSHIFT=MXSHIFT, TRADIUS=tradius,$
                 NEWWV=newwv, TRC_ORD=trc_ord, NOCRREJ=nocrrej, REJSIG=rejsig,$
                 REBINC=rebinc, BKAPER=bkaper, SIG_COLL=sig_coll, SKY=sky, $
                 SILENT=silent, CHK=chk


;  Error catching
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_extobjbox, fx, wv, xyguess, fin_spec, MSK=, VAR=, '
    print, '       WVMNX=, NMED=, PIX=, FRAC=j, RADIUS=, APER=, BKAPER='
    print, '       TOT_TRC=, CRVAL1=, CDELT=, /REDBLUE, /CREBIN [v1.0]'
    return
  endif 


;  Optional Keywords

  sz = size(fx, /dimensions)
  if not keyword_set(REJSIG) then rejsig = 7.
  if not keyword_set(WVMNX) then wvmnx = [3000., 12000.]
  if not keyword_set( VAR ) then begin
      var = fltarr(sz[0],sz[1])
      var = (fx>0.) + 25.  ; Assumes readnoise = 25
  endif

  ;; Grab the aperture and profile
  if not keyword_set( APER ) then begin
      ;; Get the aperture
      x_extapprof, fx, wv, xyguess, APSTRCT=apstrct, FLG_APER=flg_aper, MSK=msk,$
        TOT_TRC=tot_trc, VAR=var, COLLMNX=collmnx, NMED=nmed, NAVE=nave, $
        RADIUS=radius, frac=frac, TRADIUS=tradius, SIG_COLL=sig_coll, $
        TRC_ORD=trc_ord, CHK=chk
      if flg_aper EQ -1 then begin
          print, 'x_extobjbox: No object'
          fin_spec = { npix: 0 }
          return
      endif
      if apstrct.flg_smsh NE 1 then begin
          if keyword_set( BKAPER ) then apstrct.aper = bkaper $
          else apstrct.aper = [5.,5.]
          print, 'x_extobjbox: No flux to set aperture for! Using backup!'
      endif
  endif else begin
      print, 'x_extobjbox: Using an aperture ', aper
      apstrct = { $
                aper: aper, $
                gcoeff: fltarr(3), $
                flg_smsh: 0 $
                }
  endelse



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;  EXTRACT ;;;;;;;

  fx_img = fx
  msk = bytarr(sz[0],sz[1])

  ; LOOP
  for i=0L,sz[0]-1 do begin
      ; Check wavelength
      a = where(wv[i,*] GT wvmnx[0] AND wv[i,*] LT wvmnx[1], nwv)
      if nwv EQ 0 then continue
      ; Find aperture edges
      lmin = (sz[1]-1) < (long(tot_trc[i]-apstrct.aper[0]-0.5)+1) > 0L
      lmax = 0L > (long(tot_trc[i]+apstrct.aper[1]+0.5) < (sz[1]-1)) 
      msk[i,lmin:lmax] = 1
      ; Deal with endpoints
      if (lmin-0.5) LT tot_trc[i]-apstrct.aper[0] then $
        fx_img[i,lmin] = fx[i,lmin]*(1.-(tot_trc[i]-apstrct.aper[0]-lmin+0.5))
      if lmax GT tot_trc[i]+apstrct.aper[1] then $
        fx_img[i,lmax] = fx[i,lmax]*(tot_trc[i]+apstrct.aper[1]+0.5-lmax)
      ; Find total 
      tot = total(fx_img[i,lmin:lmax])
      ; Ignore tot LT 0 case
      if tot GT 0. AND not keyword_set( NOCRREJ ) AND $
        apstrct.flg_smsh EQ 1 then begin
          xpix = lmin + findgen(lmax-lmin+1) - tot_trc[i] + apstrct.sv_cen
          ; Get spline values
          spval = spl_interp(apstrct.ax,apstrct.ay,apstrct.splin, $
                             xpix,/double)/apstrct.spln_area
          ; Calc number of sigma [need to do correct variance]
          nsig = abs(fx_img[i,lmin:lmax]/tot - spval) / $
            (sqrt(var[i,lmin:lmax]) / tot)
          ; Reject (require >2 factor difference)
          rej = where(nsig GT rejsig AND $
                      abs(fx_img[i,lmin:lmax]/tot/spval) GT 2, nrej)
          if nrej NE 0 then var[i,lmin+rej] = -1.
      endif
  endfor
  gdpix = where(msk EQ 1) 

;;;;;;;; REBIN ;;;;;;;;

  ; Set wavelength scale (constant velocity (resolution) as default)
  if not keyword_set( NEWWV ) then begin
      if not keyword_set( CRVAL1 ) then $
        crval1 = double(alog10( wv[0,tot_trc[0]] < wv[sz[0]-1,tot_trc[sz[0]-1]] ))
      if not keyword_set( CDELT ) then $
        cdelt = (alog10( wv[0,tot_trc[0]] > wv[sz[0]-1,tot_trc[sz[0]-1]] ) - $
                 crval1)/double(sz[0])
      if not keyword_set( NPIX ) then npix = 2000L
      newwv = 10^( crval1 + dindgen(npix)*cdelt)
  endif else npix = n_elements(newwv)

  ; Check wavelength direction
  if not keyword_set( REDBLUE ) then begin
      if (tot_trc[sz[0]/2] LT 0. OR tot_trc[sz[0]/2] GT (sz[1]-1)) then begin
          mxt = where(tot_trc GT 0. AND tot_trc LT (sz[1]-1))
          mn = min(abs(sz[0]/2 - mxt),rb_idx)
          rb_idx = mxt[rb_idx]
      endif else rb_idx = sz[0]/2
      if wv[rb_idx,tot_trc[rb_idx]] GT wv[rb_idx + 1, tot_trc[rb_idx + 1]] $
        then redblue = 1 else redblue = 0
  endif 
  
  ; Rebin
  if not keyword_set( SILENT ) then $
    print, 'x_extobjbox: Rebinning...'
  x_rebin2dspec, wv, fx, newwv, newfx, GDPIX=gdpix, $
    VAR=var, NWVAR=newvar, SILENT=silent, REDBLUE=redblue, CR=rej_cr, $
    REBINC=rebinc
if keyword_set(SKY) then begin
  x_rebin2dspec, wv, sky, newwv, skyspec, GDPIX=gdpix, $
    VAR=var, NWVAR=newvar, SILENT=silent, REDBLUE=redblue, CR=rej_cr, $
    REBINC=rebinc
  norm = sky
  norm[*] = 1.
  x_rebin2dspec, wv, norm, newwv, nrmspec, GDPIX=gdpix, $
    VAR=var, NWVAR=newvar, SILENT=silent, REDBLUE=redblue, CR=rej_cr, $
    REBINC=rebinc
  ;; Normalize by boxcar aperture -- Gives flux per pixel
  gdn= where(nrmspec > 0.)
  skyspec[gdn] = skyspec[gdn] / nrmspec[gdn]
endif else skyspec = 0*newvar
;;;;;; PASS IT BACK ;;;;;;;
 fin_spec = { $
               npix: npix, $
               wv: newwv, $
               fx: newfx, $
               var: newvar, $
               novar: 0*newvar, $
               sky: skyspec, $
               trc: tot_trc, $
               aper: apstrct.aper, $
               gauss_sig: apstrct.gcoeff[2] $
             }
  
;;   fin_spec = { $
;;                npix: npix, $
;;                wv: newwv, $
;;                fx: newfx, $
;;                var: newvar, $
;;                trc: tot_trc, $
;;                aper: apstrct.aper $
;;              }

  return
end
