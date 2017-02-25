;+
; NAME:
;   long_extract_optimal
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

function long_extract_optimal, wave, imgminsky, ivar, oprof, mask, skyimage $
                               , trace, BOX_RAD = box_rad $
                               , modelivar = modelivar, varnoobj = varnoobj $
                               , rn_img = rn_img


IF NOT KEYWORD_SET(MODELIVAR) THEN MODELIVAR = ivar
IF NOT KEYWORD_SET(VARNOOBJ) THEN varnoobj = 0*ivar
IF NOT KEYWORD_SET(RN_IMG) THEN RN_IMG = 0*skyimage

nc = (size(wave))[1]
nr = (size(imgminsky))[2]

nt = n_elements(trace)

; SIVAR_OPT = computed from the the sciivar image, i.e.  1/sciimg
; IVAR_OPT  = computed from the modelivar image which has a less noisy
;           model of what the sky and object are at each pixel. 
struct = create_struct('WAVE_OPT', dblarr(nt) $      ; optimal wavelengths
                       , 'FLUX_OPT', fltarr(nt) $ ; optimal flux
                       , 'SIVAR_OPT', fltarr(nt) $ ; optimal inverse var
                       , 'IVAR_OPT', fltarr(nt) $ ; model optimal inverse var
                       , 'SKY_OPT', fltarr(nt) $ ; optimally exttracted sky
                       , 'RN_OPT', fltarr(nt) $ ; sigma from RN in combined 
                       , 'NIVAR_OPT', fltarr(nt) $ ; optimal sky + RN ivar
                       , 'MASK_OPT', bytarr(nt) $ ; optimal mask
                       , 'FRAC_USE', fltarr(nt) $ ; frac of profile pixels used
                       , 'CHI2', fltarr(nt) $      ; chi2 of optimal model fit
                       , 'WAVE_BOX', dblarr(nt) $  ; boxcar wavelengths
                       , 'FLUX_BOX', fltarr(nt) $  ; boxcar flux
                       , 'SIVAR_BOX', fltarr(nt) $  ; boxcar inverse var
                       , 'IVAR_BOX', fltarr(nt) $ ; boxcar model inverse var
                       , 'NIVAR_BOX', fltarr(nt) $ ; boxcar sky + RN ivar
                       , 'SKY_BOX', fltarr(nt) $   ; boxcar sky
                       , 'RN_BOX', fltarr(nt) $    ; sigma from RN in combined
                       , 'MASK_BOX', bytarr(nt) $    ; optimal mask
                       , 'MINCOL', 0L, 'MAXCOL', 0L $ ; minmax of profile fits
                       , 'BOX_RAD', box_rad)          ; boxcar radius


mincol = min(where(oprof GT 0) mod nc, max = maxcol)
nsub = maxcol - mincol + 1L
;stop
if (mincol EQ -1 OR maxcol EQ -1) then begin
    splog, 'WARNING: Invalid trace; skipping optimal extraction'
    ;; JXP -- Set wave to -1 to indicate failures
    struct.wave_opt = -1.
    return, struct
endif
mask_sub  = mask[mincol:maxcol, *]
ivar_sub  = ivar[mincol:maxcol, *]
wave_sub  = wave[mincol:maxcol, *]
img_sub   = imgminsky[mincol:maxcol, *]
sky_sub   = skyimage[mincol:maxcol, *]
rn2_sub   = rn_img[mincol:maxcol, *]^2
;; enforce positivity for all of these quantities, as they are used as weights
mivar_sub = modelivar[mincol:maxcol, *] >  0.0 ; modelivar use for the weights 
vno_sub   = varnoobj[mincol:maxcol, *]  >  0.0
; enforce normalization of profile
norm_oprof = total(oprof[mincol:maxcol, *], 1, /doub, /nan) $ 
  ## replicate(1.0D, nsub)
oprof_sub = oprof[mincol:maxcol, *]/norm_oprof >  0.0

ivar_denom = total(mask_sub*oprof_sub, 1, /doub, /nan)
sivar_num  = total(mask_sub*oprof_sub^2*ivar_sub, 1, /doub, /nan)
sivar_opt  = sivar_num/(ivar_denom + (ivar_denom EQ 0))
mivar_num  = total(mask_sub*oprof_sub^2*mivar_sub, 1, /double, /nan)
mivar_opt  = mivar_num/(ivar_denom + (ivar_denom EQ 0))
flux_opt = total(mask_sub*oprof_sub*mivar_sub*img_sub, 1, /double, /nan) $
  /(mivar_num + (mivar_num EQ 0))
;; Optimally extracted noise variance (sky + read noise) only. Since
;; this variance is not the same as that used for the weights, we
;; don't get the usual cancellation. Additional denom factor is the
;; analog of the numerator in Horne's variance formula. Note that we
;; are only weighting by the profile (mivar_sub=1) because
;; otherwise the result depends on the signal (bad).
nivar_num  = total(mask_sub*oprof_sub^2, 1, /doub, /nan) ;; uses unit weights
nvar_opt = ivar_denom*total(mask_sub*oprof_sub^2*vno_sub, 1, /doub, /nan) $
  /(nivar_num^2 + (nivar_num^2 EQ 0)) ;; unit noise weights
nivar_opt  = 1.0D/(nvar_opt + (nvar_opt EQ 0))
;; this give a nearly identical result although formally they are not
;; the same.  
;nivar_num_test = total(mask_sub*oprof_sub^2/(vno_sub + (vno_sub EQ 0)), 1 $
;                       , /dou, /nan)
;nivar_test = nivar_num_test/(ivar_denom + (ivar_denom EQ 0))

;; Optimally extract sky and (read noise)^2 in a similar way
sky_opt = ivar_denom*total(mask_sub*oprof_sub^2*sky_sub, 1, /dou, /nan) $
  /(nivar_num^2 + (nivar_num^2 EQ 0))
rn2_opt = ivar_denom*total(mask_sub*oprof_sub^2*rn2_sub, 1, /dou, /nan) $
  /(nivar_num^2 + (nivar_num^2 EQ 0))
rn_posind = WHERE(rn2_opt GE 0, npos)
rn_opt = fltarr(nt)
IF npos GT 0 THEN rn_opt[rn_posind] = sqrt(rn2_opt[rn_posind])

tot_weight = total(mask_sub*oprof_sub*mivar_sub, 1, /dou, /nan)
mask_opt  = (tot_weight GT 0.0 AND mivar_num GT 0.0 AND ivar_denom GT 0.0)

frac_use = total(oprof_sub*(mask_sub*mivar_sub GT 0), 1, /dou, /nan)
;; Use the same weights = oprof^2*mivar for the wavelenghts as the flux. 
;; Note that for the flux, one of the oprof factors cancels which does 
;; not for the wavelengths. 
wave_opt = total(mask_sub*oprof_sub^2*mivar_sub*wave_sub, 1, /dou, /nan)/$
  (mivar_num + (mivar_num EQ 0))
;wave_opt = total(oprof_sub^2* wave[mincol:maxcol, *] * mweight, 1) $
;  /(mivar_num + (mivar_num EQ 0))
badwvs = where(mivar_num LE 0 OR finite(wave_opt) NE 1 OR wave_opt LE 0.0, nbad)
;; interpolate wavelengths over masked pixels
IF badwvs[0] NE -1 then begin
    ;; can we use the profile average wavelengths instead?
    oo = WHERE(total(oprof_sub[*, badwvs], 1) GT 0.0, noo)
    IF noo GT 0 THEN $
      wave_opt[badwvs[oo]] = total(oprof_sub[*, badwvs[oo]]^2* $
                                   wave_sub[*, badwvs[oo]], 1, /dou, /nan)/ $
      total(oprof_sub[*, badwvs[oo]]^2, 1, /dou, /nan)
    ;; for pixels with completely bad profile values, interpolate from trace. 
    bb = WHERE(total(oprof_sub[*, badwvs]^2, 1, /dou, /nan) LE 0.0 OR $
               finite(total(oprof_sub[*, badwvs], 1, /dou, /nan)) NE 1 OR $
               wave_opt LE 0.0 OR finite(wave_opt) NE 1, nbb)
    IF nbb GT 0 THEN $
      wave_opt[badwvs[bb]] = interpolate(wave, trace[badwvs[bb]], badwvs[bb])
ENDIF

flux_model = flux_opt ## replicate(1., nsub)*oprof_sub
chi2 = total((img_sub - flux_model)^2*mivar_sub*mask_sub, 1, /dou, /nan)/ $
  ((total(mivar_sub*mask_sub GT 0, 1, /dou, /nan) - 1) > 1)
;struct.trace = trace
struct.wave_opt = wave_opt
struct.flux_opt = flux_opt*mask_opt
struct.sivar_opt = sivar_opt*mask_opt
struct.ivar_opt = mivar_opt*mask_opt
struct.nivar_opt = nivar_opt*mask_opt
struct.mask_opt  = mask_opt
struct.sky_opt = sky_opt
struct.rn_opt  = rn_opt
struct.frac_use = frac_use
struct.chi2 = chi2

IF NOT keyword_set(BOX_RAD) OR (nt NE nr) then return, struct

;; we compute "denom" in case the trace goes off the edge of the image.
flux_box = extract_boxcar(imgminsky, trace, findgen(nr), radius = box_rad)
box_denom = extract_boxcar((wave GT 0.0), trace, findgen(nr), radius = box_rad)
wave_box = extract_boxcar(wave, trace, findgen(nr), radius = box_rad) $
  /(box_denom + (box_denom EQ 0))
varimg = 1.0D/(ivar + (ivar EQ 0))
var_box = extract_boxcar(varimg, trace, findgen(nr), radius = box_rad)
mvarimg = 1.0/(modelivar + (modelivar EQ 0))
mvar_box = extract_boxcar(mvarimg, trace, findgen(nr), radius = box_rad)
nvar_box = extract_boxcar(varnoobj, trace, findgen(nr), radius = box_rad)
sky_box = extract_boxcar(skyimage, trace, findgen(nr), radius = box_rad)
rn2_box  = extract_boxcar(rn_img^2, trace, findgen(nr), radius = box_rad)
rn_posind = WHERE(rn2_box GT 0, npos)
rn_box = fltarr(nt)
IF npos GT 0 THEN  rn_box[rn_posind] = sqrt(rn2_box[rn_posind])

pixtot = extract_boxcar(float(modelivar*0 + 1.0D), trace, findgen(nr) $
                        , radius = box_rad) 
mask_box = extract_boxcar(float((modelivar*mask) EQ 0), trace, findgen(nr) $
                          , radius = box_rad) NE pixtot
;;   =1 for good, =0 for bad (if every pixel is masked then mask the boxcar)
;; interpolate wavelengths over masked pixels
badwvs = where(wave_box LE 0 OR finite(wave_box) EQ 0 OR box_denom EQ 0, nbad)
;; interpolate wavelengths over masked pixels
IF nbad GT 0 THEN wave_box[badwvs] = interpolate(wave, trace[badwvs], badwvs)

;   =1 for good, =0 for bad (if every pixel is masked then mask the boxcar)
sivar_box  = 1.0/(var_box  + (var_box EQ 0))
mivar_box = 1.0/(mvar_box + (mvar_box EQ 0))
nivar_box = 1.0/(nvar_box + (nvar_box EQ 0))

struct.wave_box  = wave_box
struct.flux_box  = flux_box*mask_box
struct.sivar_box  = sivar_box*mask_box
struct.ivar_box = mivar_box*mask_box
struct.nivar_box = nivar_box*mask_box
struct.mask_box  = mask_box
struct.sky_box   = sky_box
struct.rn_box    = rn_box

return, struct
end

