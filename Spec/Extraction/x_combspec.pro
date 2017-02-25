;+ 
; NAME:
; x_combspec
;    Version 1.2
;
; PURPOSE:
;   Combines multiple exposures of the same slit
;
; CALLING SEQUENCE:
;  x_combspec, flux, var, fflux, fvar, NRMFLUX=, WAVE=, NSIG=,
;               SCALE=, SNR=
;
; INPUTS:
;   flux -- 2D flux array
;   var  -- 2D variance array
;
; RETURNS:
;  fflux -- Combined flux array
;  fvar  -- Combined variance array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NRMFLUX  -  Array of 2 elements (wave min,max) to normalize flux
;  WAVE=    -- Wavelength array for normalizing
;  WVMNX=   -- Wavelength region for normalizing
;  NSIG=    -- Number of sigma to reject on
;  SCALE=   -- Scale the flux arrays by this array
;  SNR=     -- Signal-to-noise ratios to weight by
;  SKY=     -- Array containing the sky level
;
; OPTIONAL OUTPUTS:
;  FSKY=    -- Combined sky spectrum
;  FNOVAR=  -- Variance without object noise
;
; COMMENTS:
;
; EXAMPLES:
;   x_combspec, flux, var, fflux, fvar
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-May-2002 Written by JXP
;   03-Sep-2002 Added SCALE keyword + weighted means for 2 exposures
;   25-Aug-2004 Modified x_combspec_two to scale both spectra and to
;               combine if one spectra has zero flux (Kathy Cooksey)
;   27-Dec-2006 Added novar functionality to combine variances which do 
;               not include object noise (Joe Hennawi). 
;-
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
pro x_combspec_two, flux, var, fflux, fvar, WAVE = wave, WVMNX = wvmnx $
                    , MEDFLUX = medflux, novar = novar, fnovar = fnovar $
                    , sky = sky, FSKY = FSKY

;; SIZE
sz = size(flux, /dimensions)

  ; Normalize the flux
  if not keyword_set(MEDFLUX) then begin
      if keyword_set(WVMNX) then begin
          gdwv = where(wave GT wvmnx[0] AND $
                       wave LT wvmnx[1] AND $
                       var[*,0] GT 0. AND $
                       var[*,1] GT 0., ngd)
          if ngd EQ 0 then begin
              print, 'x_combspec_two: No way to normalize!'
              stop
              return
          endif
          rtio = flux[gdwv,0] / flux[gdwv,1]
          scale = median(rtio)
          medflux = 1.
      endif else begin
          scale = 1.
          medflux = 1.
      endelse
  endif else scale = medflux[1]

  ; Reset the flux of second exposure
  
  fx2 = flux[*,1] * scale
  var2 = var[*,1] * scale^2
  novar2 = novar[*, 1]*scale^2
  sky2 = sky[*, 1]* scale
  ; Reset the flux of first exposure

  fx1 = flux[*,0] * medflux[0]
  var1 = var[*,0] * medflux[0]^2
  novar1 = novar[*, 0] * medflux[0]^2
  sky1 = sky[*, 0]*medflux[0]

;;;;;;;
; Add it up
               
  fflux = fltarr(sz[0])
  fvar = replicate(-1.d, sz[0]) 
  fnovar = replicate(-1.d, sz[0])
  fsky = fltarr(sz[0])

 ; Good pixels in both (weighted mean)
  gdpix = where(var[*,0] GT 0. AND var[*,1] GT 0., ngd)

  IF ngd NE 0 THEN BEGIN
      ;; Weight
      wgt1 = 1./var1[gdpix]
      wgt2 = 1./var2[gdpix]
      
      smm_weight = wgt1 + wgt2
      fflux[gdpix] = (fx1[gdpix]*wgt1 + fx2[gdpix]*wgt2)/smm_weight
      fsky[gdpix]  = (sky1[gdpix]*wgt1 + sky2[gdpix]*wgt2)/smm_weight
      fnovar[gdpix] = (novar1[gdpix]*wgt1^2 + novar2[gdpix]*wgt2^2)/smm_weight^2
      fvar[gdpix] = 1./smm_weight
  ENDIF
  
 ; Good pixel in first
  onepix = where(var[*,0] GT 0. AND var[*,1] LE 0., ngd)
  if ngd GT 0 then begin
      fflux[onepix] = fx1[onepix] 
      fvar[onepix] = var1[onepix]
      fnovar[onepix] = novar1[onepix]
      fsky[onepix] = sky1[onepix]
  endif

 ; Good pixel in the 2nd
  twopix = where(var[*,0] LE 0. AND var[*,1] GT 0., ngd)
  if ngd GT 0 then begin
      fflux[twopix] = fx2[twopix] 
      fvar[twopix] = var2[twopix]
      fnovar[twopix] = novar2[twopix]
      fsky[twopix] = sky2[twopix]
  endif

  return
end

;------------------------------------------------------------------------------
pro x_combspec_all, flux, var, fflux, fvar, WAVE=wave, WVMNX=wvmnx, NSIG=nsig,$
                    NRMFLX=nrmflx, MEDFLUX=medflux, SNR=snr, MASK=mask,$
                    novar = novar, fnovar = fnovar, sky = sky, FSKY = FSKY

  ;; KEYWORD
  if not keyword_set( NSIG ) then nsig = 5.

  ;; SIZE
  sz = size(flux,/dimensions)
  if not keyword_set( SNR ) then begin
      snr = replicate(1., sz[1])
      if keyword_set(WVMNX) then begin
          for q=0L,sz[1]-1 do begin
              gdwv = where(wave GT wvmnx[0] AND $
                           wave LT wvmnx[1] AND $
                           var[*,0] GT 0. AND $
                           var[*,q] GT 0., ngd)
              snr[q] = median(flux[gdwv,q] / sqrt(var[gdwv,q]))
          endfor
      endif
  endif

  ;; Normalize the flux
  nrmflx = fltarr(sz[1])
  if not keyword_set( MEDFLUX ) then begin
      if keyword_set(WVMNX) then begin
          nrmflx[0] = 1.
          ;; Normalize to first spectrum
          for q=1L,sz[1]-1 do begin
              gdwv = where(wave GT wvmnx[0] AND $
                           wave LT wvmnx[1] AND $
                           var[*,0] GT 0. AND $
                           var[*,q] GT 0., ngd)
              if ngd EQ 0 then begin
                  print, 'x_combspec_all: No way to normalize!'
                  stop
                  return
              endif
              rtio = flux[gdwv,0] / flux[gdwv,q]
              nrmflx[q] = median(rtio)
          endfor
      endif else nrmflx[*] = 1.
  endif else nrmflx = medflux

  ;; Final arrays
  fflux = fltarr(sz[0])
  fsky  = fltarr(sz[0])
  fvar = replicate(-1.d,sz[0])
  fnovar = replicate(-1.d, sz[0])
  if not keyword_set( MASK ) then mask = replicate(0,sz[0],sz[1])

  ;; Coadd with rejection
  ;stop
  for i=0L,sz[0]-1 do begin
      gd = where(var[i,*] GT 0., ngd)
      if ngd EQ 0 then begin 
          mask[i,*] = 0
          continue
      endif
      mask[i,gd] = 1

      ;; Median
      mdval = median(flux[i,gd]*nrmflx[gd],/even)

      ;; Reject
      dff = abs(flux[i,gd]*nrmflx[gd] - mdval)/(sqrt(var[i,gd])*nrmflx[gd])
      if nsig GT 0. then begin
          gdgd = where(dff LT nsig, ngdgd,complement=msk,ncomplement=nmsk)
          if nmsk GT 0 then mask[i,gd[msk]] = 0
          if ngdgd EQ 0 then continue
          gd = gd[gdgd]
      endif
      ;; Weighted mean
      fflux[i] = total((snr[gd]^2)*flux[i,gd]*nrmflx[gd])/total(snr[gd]^2)
      fsky[i]  = total((snr[gd]^2)*sky[i, gd]*nrmflx[gd])/total(snr[gd]^2)
      fvar[i] = total((snr[gd]^4)*var[i,gd]*(nrmflx[gd]^2))/$
        (total(snr[gd]^2))^2
      fnovar[i] =  total((snr[gd]^4)*novar[i, gd]*(nrmflx[gd]^2))/$
        (total(snr[gd]^2))^2
   endfor

  return
end

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro x_combspec, flux, var, fflux, fvar, NRMFLUX=nrmflux, WAVE=wave, NSIG=nsig,$
                SCALE = scale, SNR = snr, MASK = mask $
                , NOVAR = NOVAR, FNOVAR = FNOVAR, SKY = SKY, FSKY = FSKY

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
          'x_combspec, flux, var, fflux, fvar, WAVE=, NRMFLUX=, NSIG=, SCALE=, ' + $
          'SKY=, SNR=, FNOVAR= [v1.2]'
    return
  endif 
  IF NOT KEYWORD_SET(NOVAR) THEN NOVAR = var
  IF NOT KEYWORD_SET(SKY) THEN SKY = 0.0*flux

;  Optional Keywords

;  Size
  sz = size(flux, /dimensions)
  if n_elements(sz) NE 2 then begin
      print, 'x_combspec: Only one spectra input! Returning...'
      return
  endif

; 
  case sz[1] of 
      2: x_combspec_two, flux, var, fflux, fvar, WVMNX=nrmflux, WAVE=wave, $
        MEDFLUX=SCALE, NOVAR = NOVAR, FNOVAR = FNOVAR, SKY = SKY, FSKY = FSKY
      else: x_combspec_all, flux, var, fflux, fvar, $
        WVMNX=nrmflux, WAVE=wave, NSIG=nsig, MEDFLUX=SCALE, SNR=snr, MASK=mask $
        , novar = novar, FNOVAR = FNOVAR, SKY = SKY, FSKY = FSKY
  endcase

;  print, 'x_combspec:  All done!'


  return
end
  

