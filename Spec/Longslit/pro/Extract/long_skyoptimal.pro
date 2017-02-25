;+
; NAME:
;   long_skyoptimal
;
; PURPOSE:
; Fit a non-parameteric object profile to the data unless the S/N is
; poor (less than 3) in which case fit a simple Gaussian.
;
; CALLING SEQUENCE:
; struct = long_skyoptimal(wave, data, ivar, oprof, sortpix, sigrej =)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
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
;   11-Mar-2005  Written by JH + SB
;-  
;------------------------------------------------------------------------------
;  routine to take an subregion of an image, with inverse variance,
;  and spatial profiles of object and sky.  Assumes that sky_profile
;  represents spatial sky_illumination at constant flux density

FUNCTION long_skyoptimal, wave, data, ivar, oprof, sortpix, sigrej = sigrej $
                          , obj_bmodel = obj_bmodel, hand_bmodel = hand_bmodel, npoly = npoly $
                          , outmask = outmask, rchi2 = rchi2  $
                          , IGNOREOBJ = IGNOREOBJ $
                          , spatial = spatial, _EXTRA = EXTRA $
                          , pixflat_mode = pixflat_mode, silent = silent

nx = n_elements(data)
if NOT keyword_set(sortpix) then sortpix = lindgen(nx)
if NOT keyword_set(npoly) then npoly = 1L
if NOT keyword_set(sigrej) then sigrej = 3.

nc = (size(oprof))[1]
nobj = n_elements(oprof)/nc

if nc NE nx then begin
   splog, 'Object profile should have 1st dimension of size nx'
   return, 0
endif

if NOT keyword_set(silent) then print, 'Iter  Chi^2  Rejected pts'

xmin = 0.
xmax = 1.

if npoly EQ 1 OR NOT keyword_set(spatial) then $
    profile_basis = [[oprof], [wave[*]*0+1]] $
else begin
   xmin = min(spatial, max=xmax)
   x2 = 2.*(spatial[*]-xmin)/(xmax-xmin) - 1
   poly_basis = flegendre(x2, npoly)
   profile_basis = [[oprof], [poly_basis]]
endelse

if keyword_set(pixflat_mode) then $
   profile_basis = [[profile_basis], [pixflat_mode]]

if nobj EQ 1 then relative_mask = (oprof GT 0) $
else relative_mask = total(oprof,2) GT 0

good = sortpix[where(ivar[sortpix] GT 0, ngood)]
good = good[sort(wave[good])]
relative = where(relative_mask[good])

; RAS MODIFICATION: when inmask is passed to bspline_longslit
; via the _EXTRA structure, we are not selecting out the 
; good pixels.
;if (n_elements(EXTRA.inmask) GT 0) then begin
;   maskgood = EXTRA.inmask[good]
;   EXTRA.inmask[0:n_elements(maskgood)-1] = maskgood
;endif

sset = bspline_longslit(wave[good], data[good], ivar[good] $
                        , profile_basis[good, *] $
                        , yfit = yfit, _EXTRA = EXTRA, upper = sigrej $
                        , lower = sigrej, red_chi = rchi $
                        , relative = relative $
                        , outmask = outmask_good1, /groupbadpix, maxrej = 5)
;stop; -- sset.coeff are all zero (KHRR)
; safe masking of 3-sigma points
;chi = (data[good] - yfit)*sqrt(ivar[good])
;chi_med = median(chi)
;mask1 = abs(chi) LT 1.5*sigrej*(chi_med > 1.)
chi2 = (data[good] - yfit)^2*ivar[good]
chi2_srt = chi2[sort(chi2)]
gauss_prob = 1.0D - 2.0D*gaussint(-double(1.2*sigrej))
sigind = round(gauss_prob*double(ngood)) <  (ngood-1L)
chi2_sigrej = chi2_srt[sigind]
mask1 = chi2 LT chi2_sigrej

if NOT keyword_set(silent) then begin
    print, '2nd round...'
    print, 'Iter  Chi^2  Rejected pts'
endif
sset = bspline_longslit(wave[good], data[good], ivar[good] * mask1 $
                        , profile_basis[good, *], yfit = yfit2, _EXTRA = EXTRA $
                        , upper = sigrej, lower = sigrej $
                        , outmask = outmask_good $
                        , relative = relative $
                        , red_chi = rchi2, /groupbadpix, maxrej = 1 )

if (NOT keyword_set(sset)) then begin
   splog, 'WARNING: B-spline failed!'
   obj_bmodel = 0
   hand_bmodel = 0
   outmask = 0
   rchi2 = 0
   return, 0
endif

  ns = npoly+keyword_set(pixflat_mode)
  ncoeff = ns + nobj
  skyset =  create_bsplineset(sset.fullbkpt, sset.nord, npoly=ns)
  skyset.coeff = sset.coeff[nobj:*, *]  ;; Coefficients for the sky
  skyset.bkmask = sset.bkmask
  skyset.xmin = xmin
  skyset.xmax = xmax

  sky_bmodel = bspline_valu(wave, skyset, x2=spatial) 

  obj_bmodel = sky_bmodel*0.
  hand_bmodel = sky_bmodel*0
  objset =  create_bsplineset(sset.fullbkpt, sset.nord)
  objset.bkmask = sset.bkmask
  for i=0,nobj-1 do begin
      IF KEYWORD_SET(IGNOREOBJ) THEN BEGIN 
         IF IGNOREOBJ[i] EQ 1 THEN BEGIN 
            objset.coeff = sset.coeff[i, *]
            hand_bmodel = hand_bmodel + bspline_valu(wave, objset $
   , action = action, lower = laction, upper = uaction) * profile_basis[*, i]
            CONTINUE
         ENDIF
      ENDIF
      objset.coeff = sset.coeff[i, *]
      obj_bmodel = obj_bmodel + bspline_valu(wave, objset $
   , action = action, lower = laction, upper = uaction) * profile_basis[*, i]
  endfor

  outmask = long(wave*0) 
  IF KEYWORD_SET(OUTMASK_GOOD) THEN outmask[good] = outmask_good

return, sky_bmodel
end

