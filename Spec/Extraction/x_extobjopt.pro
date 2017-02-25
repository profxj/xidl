;+ 
; NAME:
; x_extobjopt   
;    Version 1.0
;
; PURPOSE:
;    Given the flux and wave images, extract a 1D spectrum using the
;    Horne optimal algorithm.  Assumes a Gaussian or input profile
;    image for now.
;
; CALLING SEQUENCE:
;   
;   x_extobjopt, fx, wv, xyguess, fin_spec
;
; INPUTS:
;    fx -- 2D Image
;    wv -- 2D wavlenegth array
;    xyguess -- One spot on the image that the object is
;
; RETURNS:
;
; OUTPUTS:
;   fin_spec -- Final 1D spectrum
;
; OPTIONAL KEYWORDS:
;   WVMNX - Endpoints for extraction (default: [3200., 11000.])
;   COLLMNX - Endpoints for collapsing the spectrum prior to extraction
;              (default: [3400., 8000.])
;   REDBLUE - Spectrum runs from red to blue
;   NEWWV -  Finaly 1D wavelength array desired
;   CRVAL1 -  Starting wavelength of final 1D array
;   CDELT -   Delta lambda (lor or linear)
;   DEFGAU = Value for Gaussian profile if data is insufficient
;   USEGAU = Use this Gaussian sigma for the profile 
;   NCRITER = Number of iterations of CR rejection (default: 1)
;
; OPTIONAL OUTPUTS:
;  QAPROF -- The object profile for QA purposes
;
; COMMENTS:
;
; EXAMPLES:
;   x_extobjopt, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_extobjopt, fx, wv, var, sky, xyguess, fin_spec, $
                 NOVAR = NOVAR1, $
                 NMED=nmed, PIX=pix, FRAC=frac, RADIUS=radius, $
                 APER=aper, TOT_TRC=tot_trc, CRVAL1=crval1, CDELT=cdelt, $
                 NPIX=npix, COLLMNX=collmnx, DEBUG=debug, REDBLUE=redblue, $
                 MXSHIFT=MXSHIFT, TRADIUS=tradius, READNO=readno, $
                 NEWWV=newwv, TRC_ORD=trc_ord, REJSIG=rejsig,$
                 REBINC=rebinc, SIG_COLL=sig_coll, MAXGAU=maxgau1, $
                 MINGAU = MINGAU1,  $
                 SILENT=silent, DEFGAU=defgau, USEGAU=usegau, CHK=chk,$
                 APSTRCT=apstrct, PROFILE=profile, PROF_LIM=prof_lim, $
                 NCRITER=ncriter, fit_aper = fit_aper, GROW=grow, USE_INPUT_VAR = USE_INPUT_VAR
;                 SCATTIMG=scattimg



;  Error catching
  if  N_params() LT 6  then begin 
    print,'Syntax - ' + $
             'x_extobjopt, fx, wv, var, sky, xyguess, fin_spec, '
    print, '       WVMNX=, NMED=, PIX=, FRAC=j, RADIUS=, APER=, BKAPER='
    print, '       TOT_TRC=, CRVAL1=, CDELT=, /REDBLUE, /CREBIN [v1.0]'
    return
  endif 


;  Optional Keywords

  sz = size(fx, /dimensions)
  if not keyword_set(READNO) then readno = 5.
  if not keyword_set(REJSIG) then rejsig = 7.
;  if not keyword_set(WVMNX) then wvmnx = [3200., 11000.]
  if n_elements(NCRITER) EQ 0 then ncriter = 1L
  if keyword_set(MAXGAU1) THEN maxgau = maxgau1 $
  ELSE maxgau = 9999.0
  if keyword_set(MINGAU1) THEN mingau = mingau1 $
  ELSE mingau = -9999.0
  

;  if not keyword_set(SCATTIMG) then stop

  ;; Ivar
  good = where(var GT 0.,complement=bad,ngood,ncomplement=nbad)
  ivar = fltarr(sz[0],sz[1])
  ivar[good] = 1./var[good]
  
  ;; PROFILE
  ;; Grab the aperture and profile
  if not keyword_set( PROFILE ) then begin
      flg_aper = 1
      if not keyword_set( APER ) AND NOT keyword_set( GCOEFF ) then $
        x_extapprof, fx, wv, xyguess, APSTRCT=apstrct, FLG_APER=flg_aper, MSK=msk,$
        TOT_TRC=tot_trc, VAR=var, COLLMNX=collmnx, NMED=nmed, NAVE=nave, $
        RADIUS=radius, frac=frac, TRADIUS=tradius, SIG_COLL=sig_coll, $
        TRC_ORD = trc_ord, CHK = chk, USEGAU = USEGAU, FIT_APER = FIT_APER $
      else apstrct = { $
                       aper: aper, $
                       gcoeff: gcoeff, $
                       flg_smsh: 1 $
                     }
      
      if flg_aper EQ -1 then begin
          print, 'x_extobjopt: No object'
          fin_spec = { npix: 0 }
          return
      endif
      if keyword_set(USEGAU) then begin
          print, 'x_extobjopt:  Using input gaussian', usegau,$
                 ' instead of fit ', apstrct.gcoeff[2]
          apstrct.gcoeff[2] = usegau 
      endif else begin
          ;; Take calculated one if it exists
          if apstrct.flg_smsh NE 1 then begin
              if keyword_set(DEFGAU) then begin
                  print, 'x_extobjopt:  Profile too weak.  Using input value of', $
                    defgau
                  apstrct.gcoeff[2] = defgau
              endif else begin
                  print, 'x_extobjopt:  Profile too weak.  Returning..'
                  fin_spec = { npix: 0 }
                  return
              endelse
          endif else begin
             IF apstrct.gcoeff[2] GT MAXGAU OR apstrct.gcoeff[2] LT MINGAU $
             then begin
                IF apstrct.gcoeff[2] GT MAXGAU THEN BEGIN
                   print, 'x_extobjopt:  Gaussian fit of sigma = ' $ $
                          + strcompress(string(apstrct.gcoeff[2] $
                                               , FORMAT = '(F6.2)'), /rem) $
                          + '  exceeds maximum allowed value sig_max= ' $
                          + strcompress(string(maxgau $
                                               , FORMAT = '(F6.2)'), /rem) 
                   print, 'x_extobjopt:  Using default gaussian value', $
                          defgau
                ENDIF ELSE IF  apstrct.gcoeff[2] LT MINGAU THEN BEGIN
                   print, 'x_extobjopt:  Gaussian fit of sigma = ' $ $
                          + strcompress(string(apstrct.gcoeff[2] $
                                               , FORMAT = '(F6.2)'), /rem) $
                          + '  is less than minimum allowed value sig_min= ' $
                          + strcompress(string(mingau $
                                               , FORMAT = '(F6.2)'), /rem) 
                   print, 'x_extobjopt:  Using default gaussian value', $
                          defgau
                ENDIF
                apstrct.gcoeff[2] = defgau
                                ;apstrct.gcoeff[2] = maxgau
             ENDIF
             defgau = apstrct.gcoeff[2]
           endelse
      endelse
      
      if not keyword_set(SILENT) then $
        print, 'x_extobjopt:  Using a Gaussian with sigma = ', $
        apstrct.gcoeff[2]       ;,$

      ;; Profile
      xc = replicate(1.,sz[0])#findgen(sz[1])
      ;; Restrict dx according to trace
      dx = (-1.*double(sz[1])) > $
        (xc - tot_trc#replicate(1.d,sz[1])) < double(sz[1])
      P =  exp( -(dx)^2 / (2.d*apstrct.gcoeff[2]^2) )

      nsigma = 2.
      Pgauss = 1.0 / sqrt(2.0 * !pi) / apstrct.gcoeff[2] * $
        exp( (-1.0/2.0) * nsigma^2)
      outcs = where(P LT Pgauss, outcnt)
      if outcnt NE 0 then P[outcs] = 0. 
      ;; Normalize
      Ptot = total(P, 2) > 1e-5
      P = P / (Ptot # (dblarr(sz[1])+1.) )
  endif else begin
      P = profile
      apstrct = { $
                  aper: [0,0.], $
                  gcoeff: fltarr(3), $
                  flg_smsh: 1 $
                }
  endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  msk = bytarr(sz[0],sz[1])
  P2 = P^2

  ;; Optimal with CR rejection
  for qq=0L,NCRITER-1 do begin
      fopt_num = (P * fx * ivar) 
      fopt_den = (P2 * ivar) 
      dd = total(fopt_den,2,/double) 
      gd = where(dd GT 0.)
      fopt_ini = dblarr(sz[0])
      fopt_ini[gd] = (total(fopt_num,2,/double))[gd] / dd[gd]
      fopt_img = fopt_ini # replicate(1., sz[1])

      ;; New Variance
      IF KEYWORD_SET(USE_INPUT_VAR) THEN nwvar = var $
      ELSE nwvar = readno^2 + P * fopt_img + sky
      if nbad NE 0 then nwvar[bad] = 0.  ; Bad pixels
      
      ;; Check for CR (only 1 iteration for now)
      cr = where((fx - P*fopt_img)^2 GT rejsig^2*nwvar, ncr)
  
      ;; Recalculate ivar
      if ncr NE 0 then ivar[cr] = 0.
   endfor

  ;; Recalulate variance
  ;delvarx, var
  fopt_num = (P * fx * ivar) 
  fopt_den = (P2 * ivar) 
  dd = total(fopt_den,2,/double) 
  gd = where(dd GT 0.)
  fopt_ini = dblarr(sz[0])
  fopt_ini[gd] = (total(fopt_num,2,/double))[gd] / dd[gd]
  fopt_img = fopt_ini # replicate(1., sz[1])
  
  ;; New Variance
  IF KEYWORD_SET(USE_INPUT_VAR) THEN nwvar = var $
  ELSE nwvar = readno^2 + P * fopt_img + sky
  good = where(nwvar GT 0. AND ivar GT 0., complement=bad, ncomplement=nbad)
  ivar = fltarr(sz[0],sz[1])
  ivar[good] = 1./nwvar[good]

  ;; Mask
  msk = fltarr(sz[0],sz[1]) + 1.
  if nbad NE 0 then msk[bad] = 0.

  ;; Optimal values now
  fopt_num = float(msk * P * fx * ivar) 
  fopt_den = float(msk * P2 * ivar) 
  fvar_num = float(msk * P ) 
  fvar_den = float(msk * P2 * ivar) 

;; Added by JFH for NOVAR 27-December 2006 
;; This is the noise from sky and read noise alone. Since this variance
;; is not the same as that used for the optimal weights, we don't get the usual
;; cancellation.
  IF KEYWORD_SET(USE_INPUT_VAR) AND KEYWORD_SET(NOVAR1) THEN novar = novar1 $
  ELSE novar = readno^2 + sky
  fnovar_den =  float(msk * P2 ) ; (set ivar = 1 for novar)
  fnovar_num =  float(msk * P2 * novar)
  sopt_num   =  float(msk * P2 *sky)
  
  
;;;;;;;; REBIN ON WAVE SCALE ;;;;;;;;

  ; Set wavelength scale (constant velocity (resolution) as default)
  if not keyword_set( NEWWV ) then begin
      if not keyword_set( CRVAL1 ) then $
        crval1 = double(alog10( wv[0,tot_trc[0]] $
                                < wv[sz[0]-1,tot_trc[sz[0]-1]] ))
      if not keyword_set( CDELT ) then $
        cdelt = (alog10( wv[0,tot_trc[0]] > wv[sz[0]-1,tot_trc[sz[0]-1]] ) - $
                 crval1)/double(sz[0])
      if not keyword_set( NPIX ) then npix = 2000L
      newwv = 10.^( crval1 + dindgen(npix)*cdelt)
  endif else npix = n_elements(newwv)

  ; Check wavelength direction
  if size(REDBLUE,/type) EQ 0 then begin
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
    print, 'x_extobjopt: Rebinning...'
  ;; Flux -- Numerator
  gdpix = where(msk * P GT 0.)
  x_rebin2dspec, wv, fopt_num, newwv, newf_num, GDPIX=gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC
  ;; Flux -- Denom
  x_rebin2dspec, wv, fopt_den, newwv, newf_den, GDPIX=gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC

  ;; Sky -- Numerator
  x_rebin2dspec, wv, sopt_num, newwv, news_num, GDPIX = gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC

  ;; Variance -- Numerator
  x_rebin2dspec, wv, fvar_num, newwv, newv_num, GDPIX=gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC
  ;; Variance -- Denom
  x_rebin2dspec, wv, fvar_den, newwv, newv_den, GDPIX=gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC
  ;; No object variance -- Numerator
  x_rebin2dspec, wv, fnovar_num, newwv, newnov_num, GDPIX = gdpix, $
    SILENT=silent, REDBLUE=redblue, /REBINC
  ;; No object variance -- Denom
  x_rebin2dspec, wv, fnovar_den, newwv, newnov_den, GDPIX = gdpix, $
                 SILENT = silent, REDBLUE = redblue, /REBINC

  ;; Profile
  if keyword_set(PROF_LIM) then begin
     x_rebin2dspec, wv, P, newwv, new_P, GDPIX = gdpix $
                    , SILENT = silent, REDBLUE = redblue, /REBINC
     ;; profile with no masking
     x_rebin2dspec, wv, P, newwv, new_P_nomask $
                    , SILENT = silent, REDBLUE = redblue, /REBINC
     gdp = WHERE(new_P_nomask GT 0.0)
     frac_mask = fltarr(npix)
     frac_mask[gdp] = new_P[gdp]/new_P_nomask[gdp]
     bad_pix = where(frac_mask LT PROF_LIM, nbad)
     if nbad NE 0 then begin
        ;; Grow -- Taking this off as a default!  JXP -- 2012 Oct 15
        if keyword_set(GROW) then $
           bad_pix = 0 > [bad_pix, (bad_pix + 1), (bad_pix-1)] < (npix-1)
        newf_den[bad_pix] = 0.
     endif
  endif

  ;; Final flux array
  gd = where(newf_den GT 0., ngd)
  newfx = fltarr(npix)
  newvar = fltarr(npix)
  newnovar = fltarr(npix)
  newsky   = fltarr(npix)
  if ngd NE 0 then newfx[gd]  = newf_num[gd] / newf_den[gd]
  if ngd NE 0 then newvar[gd]= newv_num[gd] / newv_den[gd]
  if ngd NE 0 then newnovar[gd] = newv_num[gd]*newnov_num[gd]/(newnov_den[gd])^2
  if ngd NE 0 then newsky[gd] = newv_num[gd]*news_num[gd]/(newnov_den[gd])^2
  ;; Extract sky and novar with (ivar=1) so as to not depend on the object. 
  


;;;;;; PASS IT BACK ;;;;;;;
  fin_spec = { $
             npix: npix, $
             wv: newwv, $
             fx: newfx, $
             var: newvar, $
             novar: newnovar, $
             sky: newsky, $
             trc: tot_trc, $
             aper: apstrct.aper, $
             gauss_sig: apstrct.gcoeff[2] $
             }

  return
end
