FUNCTION gnirs_gprofile, image, ivar, trace_in, flux, fluxivar, objstruct $
                         , hwidth = hwidth, nccd = nccd   $
                         , fwhmfit = fwhmfit, wave = wave $
                         , thisfwhm = thisfwhm, xnew = xnew, skymask = skymask

if NOT keyword_set(thisfwhm ) then thisfwhm = 4.
if NOT keyword_set(nccd) then nccd = 1L
if NOT keyword_set(hwidth) then hwidth = long(3.*max(thisfwhm)+1)

if NOT keyword_set(objstruct) then title_string = '' else $
  title_string = 'Order# ' + strcompress(string(objstruct.slitid), /rem)  $
              + ' Object#' + string(objstruct.objid, FORMAT = '(i3)') 


xnew = trace_in

ncol = (size(image))[1]
nrow = (size(image))[2]
x = dindgen(nrow)

profile_model = image* 0.
igood = where(ivar GT 0.0, ngd)
IF ngd EQ 0 THEN return, profile_model
top = max(where(total(ivar EQ 0, 2) LT nrow))
bot = min(where(total(ivar EQ 0, 2) LT nrow)) > 0
min_column = long(min(trace_in - hwidth)) >  bot
max_column = long(max(trace_in + hwidth)) <  top

; inconsistency in array sizes sub_obj versus profile_model. FIX!!!!!

; define indices for subimage containing object
isub = (lindgen(ncol, nrow))[min_column:max_column, *]
sub_obj  = image[isub]
sub_ivar = ivar[isub]
sub_trace = trace_in-min_column
n_sub = n_elements(sub_obj[*, 0])
sub_x = dindgen(n_sub)
; compute a bspline fit to the boxcar flux. Iterate twice here to 
; properly mask bad pixels
;f_ivar = 1./(fluxvar + (fluxvar EQ 0)) * (fluxvar GT 0)
b_answer = bspline_iterfit(x, flux, everyn = 1.5, yfit = spline_flux $
                           , invvar = f_ivar, upper = 5 $
                           , lower = 5, /groupbadpix, maxrej = 1 $
                           , outmask = bmask, /silent, /relative)
b_answer = bspline_iterfit(x, flux, everyn = 1.5, yfit = spline_flux  $
                           , invvar = f_ivar*bmask, upper = 5 $
                           , lower = 5, /groupbadpix, maxrej = 1 $
                           , outmask = bmask2, /silent, /relative)
c_answer = bspline_iterfit(x, flux, everyn = 30, yfit = cont_flux  $
                           , invvar = fluxivar*bmask2, upper = 5 $
                           , lower = 5, /groupbadpix, maxrej = 1 $
                           , outmask = cmask, /silent, /relative)

sn2 = ((spline_flux*sqrt(fluxivar)*bmask2) > 0)^2
ind_nonzero = where(sn2 GT 0, nzero)
IF nzero GT 0 THEN djs_iterstat, sn2[ind_nonzero], mean = mean_sn2 $
ELSE mean_sn2 = 0.0
sub_sn2 = djs_median(sn2, width = 10, boundary = 'reflect') $
  ##replicate(1, n_sub)
splog, 'sqrt(med(S/N)^2) is ', sqrt(mean_sn2)
; If SNR < 2 then use smooth model to normalize image
IF mean_sn2 LE 4.0 THEN spline_flux = cont_flux

; Interpolate over points <= zero in  the boxcar flux or masked points
; using a continuum model 
badpix = (spline_flux LE 0.0) OR (bmask2 EQ 0)
indbad1 = WHERE(badpix AND cont_flux GT 0.0, nbad1)
IF nbad1 GT 0 THEN spline_flux[indbad1] = cont_flux[indbad1]
indbad2 = WHERE(badpix AND cont_flux LE 0.0, nbad2)
IF nbad2 GT 0 THEN spline_flux[indbad2] = 10.0D

; create the normalized object image
sub_spline = spline_flux##replicate(1, n_sub)
norm_obj = double(sub_obj/(sub_spline + (sub_spline EQ 0)))
norm_ivar = double(sub_ivar*sub_spline^2) ;*(sub_mask GT 0)) $
; Cap very large inverse variances
ivar_smash = djs_avsigclip(norm_ivar, 1)
ivar_img = ivar_smash ## replicate(1, n_sub)
ivar_mask = (norm_obj GT -0.2 AND norm_obj LT 0.7) $
  AND (norm_ivar LT 7.0*ivar_img)
norm_ivar = norm_ivar*ivar_mask

good = where(norm_ivar GT 0, ngood)

xtemp = total(4. + sqrt(sn2 ## replicate(1, n_sub)), /cumul)
xtemp = xtemp/max(xtemp)

; norm_x is the x position along the image centered on the object trace
norm_x = sub_x#replicate(1.0D, nrow)-sub_trace##replicate(1.0D, n_sub)
x2 = x ## replicate(1.0, n_sub)
sigma = replicate((thisfwhm/2.3548D), nrow)
fwhmfit = sigma*2.3548
trace_corr = replicate(0., nrow) ; no trace corrections at this low SNR

sigma_x = norm_x[*]/(sigma ## replicate(1.0D, n_sub)) - $
  (trace_corr ## replicate(1, n_sub))
profile_model[isub] = exp(-0.5D*sigma_x^2)/sqrt(2.0D*!dpi)

title_string = title_string   $
  + ' FWHM:'  + string(thisfwhm, FORMAT = '(F6.2)') $
  +' S/N:' + string(sqrt(mean_sn2), format = '(f8.3)')

xinf = WHERE(finite(xnew) NE 1, nxinf)
IF nxinf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in trace correction'
    splog, '         Replacing with zeros....'
    xnew[xinf] = 0.0
ENDIF
inf = WHERE(finite(profile_model) NE 1, ninf)
IF ninf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in object profile'
    splog, '         Replacing with zeros....'
    profile_model[inf] = 0.0
 ENDIF
;; Normalize profile
norm = replicate(1.0, ncol) # total(profile_model, 1)
profile_model = profile_model/norm
max_model = max(profile_model) <  2.0D
goodpix = WHERE(profile_model GT 0.01*max_model, ngood)
norm_mask = fltarr(ncol, nrow)
norm_mask[goodpix] = 1
norm_obj_val = replicate(1.0, ncol) # total(norm_mask*norm_obj, 1)

qa_longslit_profile, sigma_x, norm_obj, profile_model[isub] $
                     , title = title_string
wait, 1.0D

RETURN, profile_model
END

