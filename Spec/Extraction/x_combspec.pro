;+ 
; NAME:
; x_combspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same slit
;
; CALLING SEQUENCE:
;   
;   x_combspec, wave, flux, var, fflux, fvar, NRMFLUX=
;
; INPUTS:
;   spec -- Array of spectra
;
; RETURNS:
;
; OUTPUTS:
;   
;
; OPTIONAL KEYWORDS:
;  NRMFLUX  -  Array of 2 elements (wave min,max) to normalize flux
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_combspec, wfccd, mask_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-May-2002 Written by JXP
;   03-Sep-2002 Added SCALE keyword + weighted means for 2 exposures
;-
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
pro x_combspec_two, flux, var, fflux, fvar, WAVE=wave, WVMNX=wvmnx, $
                    MEDFLUX=medflux

  ;; SIZE
  sz = size(flux,/dimensions)

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
          medflux = median(rtio)
      endif else medflux = 1.
  endif else medflux = medflux[1]

  ; Reset the flux of second exposure
  
  fx2 = flux[*,1] * medflux
  var2 = var[*,1] * medflux^2

;;;;;;;
; Add it up
               
  fflux = fltarr(sz[0])
  fvar = replicate(-1.d, sz[0]) 

 ; Good pixels in both (weighted mean)
  gdpix = where(var[*,0] GT 0. AND var[*,1] GT 0., ngd)

  ;; Weight
  wgt1 = 1./var[gdpix,0]
  wgt2 = 1./var2[gdpix]

  smm_weight = wgt1 + wgt2
  fflux[gdpix] = (flux[gdpix,0]*wgt1 + fx2[gdpix]*wgt2)/smm_weight
  fvar[gdpix] = 1./smm_weight

 ; Good pixel in first
  onepix = where(var[*,0] GT 0. AND var[*,1] LE 0., ngd)
  if ngd GT 0 then begin
      fflux[onepix] = flux[onepix,0] 
      fvar[onepix] = var[onepix,0]
  endif

 ; Good pixel in the 2nd
  twopix = where(var[*,0] LE 0. AND var[*,1] GT 0., ngd)
  if ngd GT 0 then begin
      fflux[twopix] = fx2[twopix] 
      fvar[twopix] = var2[twopix] 
  endif

  return
end

;------------------------------------------------------------------------------
pro x_combspec_all, flux, var, fflux, fvar, WAVE=wave, WVMNX=wvmnx, NSIG=nsig,$
                    NRMFLX=nrmflx, MEDFLUX=medflux, SNR=snr

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
  fvar = replicate(-1.d,sz[0])

  ;; Coadd with rejection
  for i=0L,sz[0]-1 do begin
      gd = where(var[i,*] GT 0., ngd)
      if ngd EQ 0 then continue
      ;; Median
      mdval = median(flux[i,gd]*nrmflx[gd],/even)

      ;; Reject
      dff = abs(flux[i,gd]*nrmflx[gd] - mdval)/(sqrt(var[i,gd])*nrmflx[gd])
      gdgd = where(dff LT nsig, ngdgd)

      if ngdgd EQ 0 then continue
      gd = gd[gdgd]

      ;; Weighted mean
      fflux[i] = total((snr[gd]^2)*flux[i,gd]*nrmflx[gd])/total(snr[gd]^2)
      fvar[i] = total((snr[gd]^4)*var[i,gd[gdgd]]*(nrmflx[gd]^2))/$
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
                SCALE=scale, SNR=snr

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'x_combspec, flux, var, fflux, fvar, WAVE=, NRMFLUX= [v1.0]'
    return
  endif 

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
        MEDFLUX=SCALE
      else: x_combspec_all, flux, var, fflux, fvar, $
        WVMNX=nrmflux, WAVE=wave, NSIG=nsig, MEDFLUX=SCALE, SNR=snr
  endcase

;  print, 'x_combspec:  All done!'


  return
end
  

