;+
; NAME:
;   bspline_magfit
;
; PURPOSE:
;   Correct a calibrated quasar with known photometry
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------
FUNCTION bspline_magfit, wave, flux, invvar, flux_std $
                         , wave_min = wave_min, wave_max = wave_max $
                         , _EXTRA = EXTRA, sensfit = sensfit $
                         , outmask = outmask, sensfunc = sensfunc, nointerp=nointerp

nx = n_elements(wave)
pos_error = 1./sqrt((invvar > 0) + (invvar EQ 0))
pos_mask = (flux GT pos_error/10.0) AND (invvar GT 0) AND flux_std GT 0.0
pos = where(pos_mask, npos)

fluxlog = 2.5*alog10(flux > (pos_error/10))
logivar = invvar * flux^2 * pos_mask*1.08574D
magfunc = 2.5*alog10(flux_std > 1.0e-2) - fluxlog
; cap the magfunc so that sensfunc < 1.0e10
magfunc = magfunc <  25.0
sensfunc = 10.0^(0.4*magfunc)*pos_mask

; Interpolate over masked pixels
if (keyword_set(nointerp) eq 0) then begin ; jm11jun10ucsd
   bad_inds = WHERE(logivar EQ 0.0, nbad $
     , COMPLEMENT = good_inds, NCOMPLEMENT = ngood)
   
   IF nbad NE 0 THEN BEGIN
      magfunc[bad_inds] = interpol(magfunc[good_inds], wave[good_inds] $
        , wave[bad_inds])
      logivar[bad_inds] = interpol(logivar[good_inds], wave[good_inds] $
        , wave[bad_inds])>0
   ENDIF
endif

;  First iteration
bset1 = bspline_iterfit(wave, magfunc, _EXTRA = EXTRA, invvar = logivar, $
                        yfit = logfit1)

modelfit1 = 10.0^(0.4*(logfit1))

residual = sensfunc/(modelfit1 + (modelfit1 EQ 0)) - 1.
new_mask = pos_mask AND sensfunc GT 0
residual_ivar = (modelfit1*flux/(sensfunc + (sensfunc EQ 0.0)))^2*invvar
residual_ivar = residual_ivar*new_mask
if (keyword_set(nointerp) eq 0) then begin ; jm11jun10ucsd
   IF nbad NE 0 THEN BEGIN
      residual[bad_inds] = interpol(residual[good_inds], wave[good_inds] $
        , wave[bad_inds])
      residual_ivar[bad_inds] = interpol(residual_ivar[good_inds] $
        , wave[good_inds], wave[bad_inds])>0
   ENDIF
endif

;residual_ivar = invvar * modelfit1^2 * (modelfit1 GT 0) * $
;  (abs(residual) LT 1)*pos_mask

;      Now do one more fit to the ratio of data/model - 1.
bset_residual = bspline_iterfit(wave, residual, _EXTRA = EXTRA  $
                                , invvar = residual_ivar, yfit = residual_fit $
                                , fullbkpt = bset1.fullbkpt, outmask = outmask)


bset_log1 = bset1
bset_log1.coeff = bset_log1.coeff + bset_residual.coeff

sensfit = 10.0^(0.4*(bspline_valu(wave, bset_log1)))

absdev = median(abs(sensfit/modelfit1-1))
splog, 'Difference between fits is ', absdev

; compute maximum and minimum wavelength
;;TESTING
wave_min = min(wave)
wave_max = max(wave)
;wave_min = min(wave[WHERE(outmask)])
;wave_max = max(wave[WHERE(outmask)])

temp = create_struct('WAVE_MIN', wave_min, 'WAVE_MAX', wave_max)

bset_log = struct_addtags(bset_log1, temp)

return, bset_log
end
      
