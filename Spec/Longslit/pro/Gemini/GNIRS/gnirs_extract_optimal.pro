;+
; NAME:
;   gnirs_extract_optimal
;
; PURPOSE:
;   Perform an optical extraction of the spectra using profile fitting
;   techinques.
;
; CALLING SEQUENCE:
;  spec = long_extract_optimal(wave, image, ivar, oprof, mask,
;  skyimage, trace)
;
; INPUTS:
;  wave -- Wavelength image
;  image -- Data image
;  ivar -- Inverse variance
;  oprof -- Object profile
;  mask  -- Image which defines the useful data regions
;
; RETURNS:
;  spec -- A comprehensive structure containing the extracted spectra,
;          trace, wavelength array, etc.
;
; OPTIONAL INPUTS:
;  BOX_RAD=  -- Radius for boxcar extraction
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
; MODELIVAR= -- Model of the inverse variance
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   traceset2xy
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;------------------------------------------------------------------------------

; routine to take an subregion of an image, with inverse variance,
; and spatial profiles of object and sky.  Assumes that sky_profile
; represents spatial sky_illumination at constant flux density
;  

function gnirs_extract_optimal, wave, image, ivar, oprof, mask, skyimage $
                                , trace, BOX_RAD = box_rad 

nc = (size(wave))[1]
nr = (size(image))[2]

nt = n_elements(trace)

struct = create_struct('WAVE_OPT', dblarr(nt) $      ; optimal wavelengths
                       , 'FLUX_OPT', fltarr(nt) $ ; optimal flux
                       , 'IVAR_OPT', fltarr(nt) $ ; optimal inverse var
                       , 'SIG_OPT', fltarr(nt) $  ; optimal sigma
                       , 'SKY_OPT', fltarr(nt) $ ; optimally exttracted sky
                       , 'MASK_OPT', bytarr(nt) $ ; optimal mask
                       , 'FRAC_USE', fltarr(nt) $ ; frac of profile pixels used
                       , 'CHI2', fltarr(nt) $      ; chi2 of optimal model fit
                       , 'WAVE_BOX', dblarr(nt) $  ; boxcar wavelengths
                       , 'FLUX_BOX', fltarr(nt) $  ; boxcar flux
                       , 'IVAR_BOX', fltarr(nt) $  ; boxcar inverse var
                       , 'SIG_BOX', fltarr(nt) $  ; boxcar sigma
                       , 'SKY_BOX', fltarr(nt) $   ; boxcar sky
                       , 'MASK_BOX', bytarr(nt) $    ; optimal mask
                       , 'MINCOL', 0L, 'MAXCOL', 0L) ; minmax of profile fits

mincol = min(where(oprof GT 0) mod nc, max = maxcol)
if (mincol EQ -1 OR maxcol EQ -1) then begin
    splog, 'WARNING: Invalid trace; skipping optimal extraction'
    return, struct
endif
masksub = mask[mincol:maxcol, *]
;; Do not use ivar weights for optimal extraction since this is 
;; likely to be noisy to CRs, defects, and hot pixels. 
mweight = masksub
thisprof = oprof[mincol:maxcol, *]

denom  = total(thisprof*masksub, 1, /dou, /nan)
denom2 = total(thisprof^2 * mweight, 1, /dou, /nan)

varimg = 1. / (ivar + (ivar EQ 0))
var_opt = denom*total(thisprof^2*mweight^2*varimg[mincol:maxcol, *], 1 $
                      , /dou, /nan)/(denom2^2 + (denom2^2 EQ 0))
ivar_opt = 1.0/(var_opt + (var_opt EQ 0))
flux_opt = total(thisprof*mweight*image[mincol:maxcol, *], 1, /dou, /nan) $
  /(denom2 + (denom2 EQ 0))
sky_opt  = total(thisprof*mweight*skyimage[mincol:maxcol, *], 1, /dou, /nan) $
  /(denom2 + (denom2 EQ 0))
mask_opt  = (denom2 GT 0.0 AND denom GT 0)

frac_use = total(thisprof * (mweight GT 0), 1, /dou, /nan)
wave_opt = total(thisprof^2* wave[mincol:maxcol, *] * mweight, 1, /dou, /nan) $
  /(denom2 + (denom2 EQ 0))
;; interpolate wavelengths over masked pixels
badwvs = where(denom2 LE 0 OR finite(wave_opt) NE 1 OR wave_opt LE 0.0, nbad)
;; interpolate wavelengths over masked pixels
IF badwvs[0] NE -1 then begin
    ;; can we use the profile average wavelengths instead?
    oo = WHERE(total(thisprof[*, badwvs], 1, /dou, /nan) GT 0.0, noo)
    IF noo GT 0 THEN $
      wave_opt[badwvs[oo]] = total(thisprof[*, badwvs[oo]]^2* $
                                   wave[mincol:maxcol, badwvs[oo]], 1 $
                                   , /dou, /nan)/ $
      total(thisprof[*, badwvs[oo]]^2, 1, /dou, /nan)
    ;; for pixels with completely bad profile values, interpolate from trace. 
    bb = WHERE(total(thisprof[*, badwvs]^2, 1, /dou, /nan) LE 0.0 OR $
               finite(total(thisprof[*, badwvs], 1, /dou, /nan)) NE 1 OR $
               wave_opt LE 0.0 OR finite(wave_opt) NE 1, nbb)
    IF nbb GT 0 THEN $
      wave_opt[badwvs[bb]] = interpolate(wave, trace[badwvs[bb]], badwvs[bb])
ENDIF




badpix = where(mask_opt EQ 0, nbad)
if badpix[0] NE -1 then begin
    flux_opt[badpix] = 0
    wave_opt = djs_maskinterp(wave_opt, (mask_opt EQ 0))
;    wave_opt[badpix] = total(oprof[*, badpix]* wave[*, badpix], 1) / $
;      total(oprof[*, badpix], 1)
endif

flux_model = flux_opt ## replicate(1., nc) * oprof
chi2 = total((image - flux_model)^2 * mweight, 1, /dou, /nan) / $
  ((total(mweight GT 0, 1, /dou, /nan) - 1) > 1)
;struct.trace = trace
struct.wave_opt = wave_opt
struct.flux_opt = flux_opt
struct.ivar_opt = ivar_opt
struct.sig_opt = (ivar_opt GT 0.0)/sqrt(ivar_opt + (ivar_opt LE 0.0))
struct.mask_opt  = mask_opt
struct.sky_opt = sky_opt
struct.frac_use = frac_use
struct.chi2 = chi2

if NOT keyword_set(BOX_RAD) OR (nt NE nr) then return, struct

;   Do not use any weights or masks in the boxcar extraction,
;   except to estimate the errors.  -DJS
;   We compute "denom" in case the trace goes off the edge of the image.
flux_box = extract_boxcar(image, trace, findgen(nr), radius = box_rad)
;;denom = extract_boxcar(0*image+1.0, trace, findgen(nr), radius = box_rad)
box_denom = extract_boxcar((wave GT 0.0), trace, findgen(nr), radius = box_rad)
wave_box = extract_boxcar(wave, trace, findgen(nr), radius = box_rad) $
  /(box_denom + (box_denom EQ 0))
;wave_box = extract_boxcar(wave, trace, findgen(nr), radius = box_rad) $
;;  /(denom + (denom LE 0))
var_box = extract_boxcar(varimg, trace, findgen(nr), radius = box_rad)
sky_box = extract_boxcar(skyimage, trace, findgen(nr), radius = box_rad)
pixtot = extract_boxcar(float(ivar*0 + 1), trace, findgen(nr) $
                        , radius = box_rad) 
mask_box = extract_boxcar(float((ivar*mask) EQ 0), trace, findgen(nr) $
                          , radius = box_rad) NE pixtot
;   =1 for good, =0 for bad (if every pixel is masked then mask the boxcar)
ivar_box  = 1.0/(var_box  + (var_box EQ 0))
;; interpolate wavelengths over masked pixels
badwvs = where(wave_box LE 0 OR finite(wave_box) EQ 0 OR denom LE 0, nbad)
;; interpolate wavelengths over masked pixels
IF nbad GT 0 THEN wave_box[badwvs] = interpolate(wave, trace[badwvs], badwvs)

struct.wave_box  = wave_box
struct.flux_box  = flux_box
struct.ivar_box  = ivar_box
struct.sig_box = (ivar_box GT 0.0)/sqrt(ivar_box + (ivar_box LE 0.0))
struct.mask_box  = mask_box
struct.sky_box   = sky_box

return, struct
end

