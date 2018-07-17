;+
; NAME:
;   long_gprofile
;
; PURPOSE:
; Fit a non-parameteric object profile to the data unless the S/N is
; poor (less than 3) in which case fit a simple Gaussian.
;
; CALLING SEQUENCE:
;  profile = long_gprofile(image, ivar, trace_in, flux, fluxivar, objstruct)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  /GAUSS -- Calculate a Gaussian profile
;  PROF_NSIGMA= -- Number of sigma to include in the profile fitting.
;  /NO_DERIV -- Turn off derivative algorithm
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
FUNCTION long_gprofile, image, ivar, waveimg $
                        , trace_in, wave, flux, fluxivar, objstruct $
                        , hwidth = hwidth, nccd = nccd $
                        , fwhmfit = fwhmfit, thisfwhm = thisfwhm $
                        , xnew = xnew, gauss = gauss $
                        , SN_GAUSS = SN_GAUSS, SILENT = SILENT $
                        , MED_SN2 = MED_SN2, PROF_NSIGMA=prof_nsigma $
                        , NO_DERIV=no_deriv, WVMNX=wvmnx

IF NOT KEYWORD_SET(SN_GAUSS) THEN SN_GAUSS = 3.0D
if NOT keyword_set(thisfwhm ) then thisfwhm = 4.
if NOT keyword_set(nccd) then nccd = 1L
if NOT keyword_set(hwidth) then hwidth = long(3.*max(thisfwhm)+1)
IF NOT KEYWORD_SET(MAX_TRACE_CORR) THEN MAX_TRACE_CORR = 2.0D
if not keyword_set(WVMNX) then wvmnx = [2900., 30000]
thisfwhm = thisfwhm > 1.0 ;; require the FWHM to be greater than 1 pixel
;; this prevents large negative FWHM

if NOT keyword_set(objstruct) then title_string = '' else $
  title_string = 'Slit# ' + strcompress(string(objstruct.slitid), /rem)  $
  + ' Object#' + string(objstruct.objid, FORMAT = '(i3)') 

xnew = trace_in

ncol = (size(image))[1] ;; spatial
nrow = (size(image))[2] ;; spectral 

top = max(where(total(ivar EQ 0, 2) LT nrow))
bot = min(where(total(ivar EQ 0, 2) LT nrow)) > 0
min_column = long(min(trace_in - hwidth)) >  bot
max_column = long(max(trace_in + hwidth)) <  top

; define indices for subimage containing object
isub = (lindgen(ncol, nrow))[min_column:max_column, *]
profile_model = image* 0.
sub_obj  = image[isub]
sub_ivar = ivar[isub]
sub_wave = waveimg[isub]
sub_trace = trace_in-min_column
n_sub = n_elements(sub_obj[*, 0])
sub_x = lindgen(n_sub)
;; create some images we will need
dim_sub = size(sub_obj, /dim)
nxsub = dim_sub[0]
nysub = dim_sub[1]
sn2_sub = fltarr(nxsub, nysub)
spline_sub = fltarr(nxsub, nysub)
; compute a bspline fit to the boxcar flux. Iterate twice here to 
; properly mask bad pixels
flux_sm = djs_median(flux, width = 5, boundary = 'reflect')
fluxivar_sm =  djs_median(fluxivar, width = 5, boundary = 'reflect')
fluxivar_sm = fluxivar_sm*(fluxivar GT 0.0)
;indsp = WHERE(wave GT 0.9 AND wave LT 2.5 $  ;; Microns!

indsp = WHERE(wave GT wvmnx[0] AND wave LT wvmnx[1] $  
              AND finite(flux_sm) AND flux_sm LT 5.0d5 $
              AND flux_sm GT -1000.0d $
              AND fluxivar_sm GT 0.0, nsp)
IF nsp GT 10 THEN BEGIN 
   b_answer = bspline_iterfit(wave[indsp], flux_sm[indsp] $
                              , everyn = 1.5, yfit = spline_flux $
                              , invvar = fluxivar_sm[indsp], upper = 5 $
                              , lower = 5, /groupbadpix, maxrej = 1 $
                              , outmask = bmask, /silent, /relative)
   b_answer = bspline_iterfit(wave[indsp], flux_sm[indsp] $
                              , everyn = 1.5, yfit = spline_flux  $
                              , invvar = fluxivar_sm[indsp]*bmask, upper = 5 $
                              , lower = 5, /groupbadpix, maxrej = 1 $
                              , outmask = bmask2, /silent, /relative)
   c_answer = bspline_iterfit(wave[indsp], flux_sm[indsp], everyn = 30 $
                              , yfit = cont_flux  $
                              , invvar = fluxivar_sm[indsp]*bmask2, upper = 5 $
                              , lower = 5, /groupbadpix, maxrej = 1 $
                              , outmask = cmask, /silent, /relative)
   sn2 = ((spline_flux*sqrt(fluxivar_sm[indsp] >  0)*bmask2) > 0)^2
;stop
   ind_nonzero = where(sn2 GT 0, nzero)
   IF nzero GT 0 THEN djs_iterstat, sn2[ind_nonzero], median = med_sn2 $
   ELSE med_sn2 = 0.0

   sn2_med = djs_median(sn2, width = 10, boundary = 'reflect')
   igood = where(sub_ivar GT 0.0, ngd)
   IF ngd GT 0 THEN BEGIN
      isrt = sort(wave[indsp])
      sn2_sub[igood] = interpol(sn2_med[isrt], wave[indsp[isrt]] $
                                , sub_wave[igood])
   ENDIF
   splog, 'sqrt(med(S/N)^2) is ', sqrt(med_sn2)

   min_wave = min(wave[indsp])
   max_wave = max(wave[indsp])
   spline_flux1 = fltarr(nysub)
   cont_flux1 = fltarr(nysub)
   sn2_1 = fltarr(nysub)
   ispline  = where(wave GE min_wave AND wave LE max_wave)
;; don't extrapolate
   spline_flux1[ispline] = bspline_valu(wave[ispline], b_answer)
   cont_flux1[ispline]   = bspline_valu(wave[ispline], c_answer)
   sn2_1[ispline] = interpol(sn2, wave[indsp], wave[ispline])
   bmask = lonarr(nysub) + 1L
   bmask[indsp] = bmask2
   spline_flux1 = djs_maskinterp(spline_flux1, bmask EQ 0)
   cmask2 = lonarr(nysub) + 1L
   cmask2[indsp] = cmask
   cont_flux1 = djs_maskinterp(cont_flux1, cmask2 EQ 0)
;stop
; If SNR^2 < 2.0 then there is not enough object flux so don't normalize
   IF med_sn2 LE 2.0 THEN BEGIN
      ;; not enough object flux to fit a profile, so just use the 
      ;; 1-sigma fluctuation of the sky to normalize things
      djs_iterstat, flux[indsp], sigma = sigma1
      spline_sub[igood] = sigma1 > 1.0d
   ENDIF ELSE BEGIN 
      IF med_sn2 LE 5.0 AND med_sn2 GT 2.0 THEN spline_flux1 = cont_flux1   
      ;; Interp over points <= 0 in boxcar flux or masked points using cont model
      badpix = (spline_flux1 LE 0.5) OR (bmask EQ 0)
      indgd = WHERE(badpix EQ 0, ngood0)
      indbad1 = WHERE(badpix AND cont_flux1 GT 0.0 AND cont_flux1 LT 5.0d5, nbad1)
      IF nbad1 GT 0 THEN spline_flux1[indbad1] = cont_flux1[indbad1]
      indbad2 = WHERE(badpix AND cont_flux1 LE 0.0 OR cont_flux1 GT 5.0d5, nbad2)
      IF nbad2 GT 0 AND ngood0 GT 0 THEN $
         spline_flux1[indbad2] = djs_median(spline_flux1[indgd])
;; take a 5-pixel median to filter out some hot pixels
      spline_flux1 = djs_median(spline_flux1, width = 5, boundary = 'reflect')
; create the normalized object image
      IF ngd GT 0 THEN BEGIN 
         isrt = sort(wave)
         spline_sub[igood] = interpol(spline_flux1[isrt], wave[isrt] $
                                      , sub_wave[igood])
      ENDIF
   ENDELSE


   norm_obj = (spline_sub NE 0.0)*float(sub_obj/(spline_sub + (spline_sub EQ 0.0)))
   norm_ivar = float(sub_ivar*spline_sub^2)
;; Cap very large inverse variances
   ivar_mask = (norm_obj GT -0.2 AND norm_obj LT 0.7) $
               AND (sub_ivar GT 0.0) $
               AND finite(norm_obj) EQ 1 $
               AND finite(norm_ivar) EQ 1
   norm_ivar = norm_ivar*ivar_mask
   good = where(norm_ivar GT 0, ngood)
   
   xtemp = total(4.d + sqrt((sn2_1 > 0.0) ## replicate(1.0d, n_sub)), /cumul)
   xtemp = xtemp/max(xtemp)

;x = findgen(nrow)
ENDIF ELSE BEGIN
   ;; What to do if all pixels are bad, for example near negative
   ;; objects in the IR. 
   ngood=0
   med_sn2=0
   djs_iterstat, flux, sigma = sigma1
   spline_sub[*] = sigma1 > 1.0d
   norm_obj = (spline_sub NE 0.0)*float(sub_obj/(spline_sub + (spline_sub EQ 0.0)))
ENDELSE
; norm_x is the x position along the image centered on the object  trace
norm_x = sub_x#replicate(1.0D, nrow)-sub_trace##replicate(1.0D, n_sub)


sigma = replicate((thisfwhm/2.3548), nrow)
fwhmfit = sigma*2.3548
trace_corr = replicate(0., nrow)
;stop
if ngood LT 10  OR med_sn2 LT SN_GAUSS^2 OR KEYWORD_SET(GAUSS) then begin
    splog, 'Too few good pixels or S/N <' + string(sn_gauss) +  $
           '  or GAUSS flag set'
    splog, 'Returning gaussian profile'
    sigma_x = norm_x[*]/(sigma ## replicate(1, n_sub)) - $
              (trace_corr ## replicate(1, n_sub))
    profile_model[isub] = exp(-0.5*sigma_x^2)/sqrt(2.0d*!dpi)*(sigma_x^2 LT 25.)
    title_string = title_string   $
      + ' FWHM:'  + string(thisfwhm, FORMAT = '(F6.2)') $
      +' S/N:' + string(sqrt(med_sn2), format = '(f8.3)')
    print, title_string
    xinf = WHERE(finite(xnew) NE 1, nxinf)
    IF nxinf NE 0 THEN BEGIN
        splog, 'WARNING: Nan pixel values in trace correction'
        splog, '         Returning original trace....'
        xnew = trace_in
    ENDIF
    inf = WHERE(finite(profile_model) NE 1, ninf)
    IF ninf NE 0 THEN BEGIN
        splog, 'WARNING: Nan pixel values in object profile'
        splog, '         Returning original trace....'
        profile_model[inf] = 0.0
    ENDIF
    ;; Normalize profile
    norm=replicate(1.0,ncol) # total(profile_model,1,/double)
    IF total(norm) GT 0.0 THEN profile_model = profile_model/norm
    IF ngood GT 0 THEN $
       qa_longslit_profile, sigma_x, norm_obj, profile_model[isub] $
                            , title = title_string, ind = good $
                            , XTRUNC = 7.0 $
    ELSE qa_longslit_profile, sigma_x, norm_obj, profile_model[isub] $
                              , title = title_string + $
                              '  no good pix showing all' $
                              , XTRUNC = 7.0
    return, profile_model
endif

sigma_iter = 3L
splog, 'Gaussian vs b-spline of width ', thisfwhm, ' pixels'

area = 1.0
sigma_x = norm_x[*]/(sigma ## replicate(1, n_sub)) - $
  (trace_corr ## replicate(1, n_sub))

;mask    = lonarr(n_sub, nrow) + 1L
;; JXP -- March 2009 fix?
mask    = lonarr(n_sub, nrow) 
skymask = lonarr(n_sub, nrow) + 1L

;;
;;   bkpt choice
;;

;; The following lines set the limits for the b-spline fit
limit     = dierfc(0.1/sqrt(med_sn2)) * sqrt(2.0)
if not keyword_set(PROF_NSIGMA) then begin
    sinh_space = 0.25 * alog10((1000./sqrt(med_sn2)) > 10.)
    abs_sigma = max(abs(sigma_x[good])) < (2. * limit)
    min_sigma = min(sigma_x[good]) > (-abs_sigma)
    max_sigma = max(sigma_x[good]) < (abs_sigma)
    nb = long(asinh(abs_sigma)/sinh_space) + 1  

    rb = sinh((findgen(nb)+0.5)*sinh_space)
    bkpt = [reverse(-rb), rb]
    keep = where(bkpt GE min_sigma AND bkpt LE max_sigma)
    bkpt = bkpt[keep]
endif else begin  ;; This is for extended or very bright objects
    splog, 'Using PROF_NSIGMA=', prof_nsigma, ' for extended/bright objects'
    nb = round(PROF_NSIGMA > 10)
    max_sigma = PROF_NSIGMA
    min_sigma = -1*PROF_NSIGMA
    sinh_space = asinh(PROF_NSIGMA)/nb
    rb = sinh((findgen(nb)+0.5)*sinh_space)
    bkpt = [reverse(-rb), rb]
    keep = where(bkpt GE min_sigma AND bkpt LE max_sigma)
    bkpt = bkpt[keep]
 endelse


;;  attempt b-spline fit first...
GOOD_PIX = sn2_sub GT SN_GAUSS^2 AND norm_ivar GT 0
IN_PIX   = sigma_x GE min_sigma AND sigma_x LE max_sigma AND norm_ivar GT 0
goodp1   = where(GOOD_PIX, ngoodpix)
inp1     = where(IN_PIX, ninpix)
IF ngoodpix GE 0.2*ninpix THEN inside = WHERE(GOOD_PIX AND IN_PIX) $
ELSE inside = WHERE(IN_PIX)
si = inside[sort(sigma_x[inside])]
sr = reverse(si)
bset = bspline_iterfit(sigma_x[si], norm_obj[si], invvar = norm_ivar[si] $
                       , nord = 4, bkpt = bkpt, maxiter = 15, yfit = mode_fit $
                       , /silent, upper = 1, lower = 1)

median_fit = median(norm_obj[where(norm_ivar GT 0)])
if abs(median_fit) GT 0.01 then $
  splog, 'Median flux level in profile is not zero', median_fit $
else median_fit = 0.

;; find peak and FWHM


long_findfwhm, mode_fit - median_fit, sigma_x[si], peak, peak_x, lwhm, rwhm
trace_corr = replicate(peak_x, nrow)
min_level = peak*exp(-0.5*limit^2)

bspline_fwhm = (rwhm - lwhm) * thisfwhm/2.3548
splog, 'Bspline FWHM: ', bspline_fwhm, $
       ' compared to initial FWHM: ', thisfwhm
sigma = sigma * (rwhm-lwhm)/2.3548

rev_fit = reverse(mode_fit)
limit = limit * (rwhm-lwhm)/2.3548
lp = min(where((rev_fit LT (min_level+median_fit) AND $
                sigma_x[sr] LT peak_x) OR sigma_x[sr] LT peak_x-limit)) 
if lp[0] NE -1 then l_limit = sigma_x[sr[lp]] $
else l_limit = min_sigma

rp = min(where((mode_fit LT (min_level+median_fit) AND $
                sigma_x[si] GT peak_x) OR sigma_x[si] GT peak_x+limit)) 
if rp[0] NE -1 then r_limit = sigma_x[si[rp]] $
else r_limit = max_sigma
;<<<<<<< long_gprofile.pro
;splog, limit, min_level, l_limit, r_limit
;stop
;=======
splog, "Limits:",limit, min_level, l_limit, r_limit

;>>>>>>> 1.21
; just grab data points inside limits
mask[si] = (norm_ivar[si] GT 0 AND abs(norm_obj[si] - mode_fit) LT 0.1)
inside = where(sigma_x[si] GT l_limit AND sigma_x[si] LT r_limit  $
               AND mask[si], ninside)
if ninside LT 10 then begin
    splog, 'Too few pixels inside l_limit and r_limit'
    profile_model[isub] = exp(-0.5*sigma_x^2)/sqrt(2*!Pi) * (sigma_x^2 LT 25.)
    title_string = title_string  + ' FWHM:'  + string(bspline_fwhm, $
                                                      FORMAT = '(F6.2)') $
      +' S/N:' + string(sqrt(med_sn2), format = '(f8.3)')
    xinf = WHERE(finite(xnew) NE 1, nxinf)
    IF nxinf NE 0 THEN BEGIN
        splog, 'WARNING: Nan pixel values in trace correction'
        splog, '         Replacing with original trace....'
        xnew = trace_in
    ENDIF 
    inf = WHERE(finite(profile_model) NE 1, ninf)
    IF ninf NE 0 THEN BEGIN
        splog, 'WARNING: Nan pixel values in object profile'
        splog, '         Replacing with zeros....'
        profile_model[inf] = 0.0
     ENDIF
    norm=replicate(1.0,ncol) # total(profile_model,1,/double)
    IF total(norm) GT 0.0 THEN profile_model = profile_model/norm
;;    norm = total(profile_model)/float(nrow)
;    IF norm GT 0.0 THEN profile_model = profile_model/norm 
;    BPH this makes no sense with the previous line commented out
    qa_longslit_profile, sigma_x, norm_obj, profile_model[isub] $
                         , l_limit, r_limit, ind = good $
                         , title = title_string, xrange = 7.0
    return, profile_model
endif

sigma_iter = 3L
splog, 'Gaussian vs b-spline of width ', thisfwhm, ' pixels'

inside = si[inside[sort(xtemp[si[inside]])]]
pb = inside*0. + 1.

for iiter = 1, sigma_iter do begin
    
    mode_zero = bspline_valu(sigma_x[inside], bset) * pb
    mode_shift = (bspline_valu(sigma_x[inside]-0.5, bset) - $
                  bspline_valu(sigma_x[inside]+0.5, bset))* pb  * $
      (sigma_x[inside] GT (l_limit + 0.5) AND $
       sigma_x[inside] LT (r_limit - 0.5))
    mode_stretch = bspline_valu(sigma_x[inside]/1.3, bset) * pb /1.3 - $
      mode_zero
    
    if nccd EQ 1 then nbkpt = long(alog10((med_sn2[0]) > 11)) $
    else nbkpt = 1L
    fullbkpt0 = bspline_bkpts(xtemp[inside], nord = 4, nbkpt = nbkpt)
    
    if nccd GT 1 then begin
        fullbkpt = fullbkpt0
        nx = n_elements(xtemp)
        for i = 1L, nccd-1 do $
          fullbkpt = [fullbkpt, replicate(xtemp[nx*i/nccd], 4)]
        fullbkpt = fullbkpt[sort(fullbkpt)]
    endif else if nccd EQ 1 THEN fullbkpt = fullbkpt0
    
    xx = total(xtemp, 1)/n_sub
    profile_basis = [[mode_zero], [mode_shift]]

    mode_shift_set = bspline_longslit(xtemp[inside], norm_obj[inside] $
                                      , norm_ivar[inside], profile_basis $
                                      , fullbkpt = fullbkpt, maxiter = 1 $
                                      , yfit = mode_shift_fit, /silent)
    
    
    temp_set = create_bsplineset(mode_shift_set.fullbkpt, mode_shift_set.nord)
    temp_set.coeff = mode_shift_set.coeff[0, *]
    h0 = bspline_valu(xx, temp_set) 
    temp_set.coeff = mode_shift_set.coeff[1, *]
    h1 = bspline_valu(xx, temp_set)
    ratio_10 = (h1/(h0 + (h0 EQ 0.0)))
    trace_corr = trace_corr + ratio_10/(1.0 + abs(ratio_10)/0.1)
    
    profile_basis = [[mode_zero], [mode_stretch]]
    mode_stretch_set = bspline_longslit(xtemp[inside], norm_obj[inside] $
                                        , norm_ivar[inside], profile_basis $
                                        , fullbkpt = fullbkpt0 $
                                        , maxiter = 1, yfit = mode_stretch_fit $
                                        , /silent)
    temp_set = create_bsplineset(mode_stretch_set.fullbkpt $
                                 , mode_stretch_set.nord)
    temp_set.coeff = mode_stretch_set.coeff[0, *]
    h0 = bspline_valu(xx, temp_set) 
    temp_set.coeff = mode_stretch_set.coeff[1, *]
    h2 = bspline_valu(xx, temp_set)
    h0 = (h0 + h2 * total(mode_stretch)/total(mode_zero)) > 0.1
    ratio_20 = (h2/(h0 + (h0 EQ 0.0)))
    sigma_factor = 0.3 *ratio_20/(1.0 + abs(ratio_20))
    IF NOT KEYWORD_SET(SILENT) THEN $
      splog, '#', iiter, ': Median trace correction ', $
           median(abs(h1/(h0 + (h0 EQ 0.0))))
    IF NOT KEYWORd_SET(SILENT) THEN $
      splog, '#', iiter, ': Median width correction ', $
           median(abs(sigma_factor))
    
    sigma = sigma * (1 + sigma_factor)
    area = area * h0 / (1 + sigma_factor)
    
    sigma_x = norm_x[*]/(sigma ## replicate(1, n_sub)) - $
      (trace_corr ## replicate(1, n_sub))
    
    if iiter LT sigma_iter-1 then begin
        
        ss = sort(sigma_x[inside])
        pb = (area ## replicate(1, n_sub))[inside]
        keep = where(bkpt GE min(sigma_x[inside]) $
                     AND bkpt LE max(sigma_x[inside]), nkeep)
        IF nkeep EQ 0 THEN keep = lindgen(n_elements(bkpt))
        bset = bspline_longslit(sigma_x[inside[ss]], norm_obj[inside[ss]] $
                                , norm_ivar[inside[ss]], pb[ss] $
                                , nord = 4, bkpt = bkpt[keep] $
                                , maxiter = 2, yfit = mode_fit, /silent)
    endif
endfor 
;stop
; Apply trace corrections only if they are small added by JFH
IF abs(djs_median(trace_corr*sigma)) LT MAX_TRACE_CORR THEN $
  xnew = trace_corr * sigma + trace_in $
ELSE xnew = trace_in

fwhmfit = sigma*2.3548
ss = sort(sigma_x)
inside = where(sigma_x[ss] GE min_sigma $
               AND sigma_x[ss] LT max_sigma $
               AND mask[ss] $
               AND finite(norm_obj[ss]) EQ 1 $
               AND finite(norm_ivar[ss]) EQ 1)
pb = area ## replicate(1, n_sub)
bset = bspline_longslit(sigma_x[ss[inside]],  norm_obj[ss[inside]], $
                        norm_ivar[ss[inside]], pb[ss[inside]], $
                        nord = 4, bkpt = bkpt, outmask = outmask, $
                        upper = 10, lower = 10, yfit = profile_fit, /silent)

skymask[*] = 1 - (sigma_x GT min_sigma AND sigma_x LT max_sigma)
;   testing
full_bsp = 0.0*skymask
bsp_pix = WHERE(skymask EQ 0, nbsp)
;    full_bsp[bsp_pix] = bspline_valu(sigma_x[bsp_pix], bset) 
full_bsp = bspline_valu(sigma_x[*], bset) * (skymask EQ 0)
;<<<<<<< long_gprofile.pro
;;THIS IS IT!!!!! -- KHRR
;=======

;>>>>>>> 1.21
long_findfwhm, full_bsp[ss] - median_fit, sigma_x[ss], peak, peak_x, lwhm, rwhm

lp = min(where(reverse((full_bsp[ss] LT (min_level+median_fit) AND $
                        sigma_x[ss] LT peak_x) $
                       OR sigma_x[ss] LT peak_x-limit))) > 0
rp = min(where((full_bsp[ss] LT (min_level+median_fit) AND $
                sigma_x[ss] GT peak_x) OR sigma_x[ss] GT peak_x+limit)) > 0
l_limit = (reverse(sigma_x[ss]))[lp] - 0.1
r_limit = sigma_x[ss[rp]] + 0.1 

repeat begin
    l_limit += 0.1
    l_fit = bspline_valu(l_limit, bset)
    l2    = bspline_valu(l_limit*0.9, bset)
    l_deriv = (alog(l2) - alog(l_fit))/(0.1*l_limit)
endrep until (l_deriv LT -1.0 OR l_limit GE -1.0)

repeat begin
    r_limit -= 0.1
    r_fit = bspline_valu(r_limit, bset)
    r2    = bspline_valu(r_limit*0.9, bset)
    r_deriv = (alog(r2) - alog(r_fit))/(0.1*r_limit)
endrep until (r_deriv GT 1.0 OR r_limit LE 1.0)

;; JXP kludge
if keyword_set(PROF_NSIGMA) then begin
   ;; By setting them to zero we ensure QA won't plot
   ;; them in the profile QA. 
   l_limit = 0.0
   r_limit = 0.0
   NO_DERIV = 1
endif

;   Hack to fix degenerate profiles which have a positive derivative
IF l_deriv LT 0 and r_deriv GT 0 and not keyword_set(NO_DERIV) THEN BEGIN
    left = where(sigma_x LT l_limit)
    if left[0] NE -1 then full_bsp[left] = $
      exp(-(sigma_x[left]-l_limit)*l_deriv) * l_fit
    right = where(sigma_x GT r_limit)
    if right[0] NE -1 then full_bsp[right] = $
      exp(-(sigma_x[right]-r_limit)*r_deriv) * r_fit
    internal = where(sigma_x GE l_limit AND sigma_x LE r_limit, nint)
    IF nint GT 0 THEN skymask[internal] = 1
ENDIF
profile_model[isub] = full_bsp * pb

res_mode = (norm_obj[ss[inside]] - profile_model[isub[ss[inside]]])* $
  sqrt(norm_ivar[ss[inside]])
good = where(outmask EQ 1 AND norm_ivar[ss[inside]] GT 0, ngood)
chi_med = median(res_mode[good]^2)
chi_zero = median(norm_obj[ss[inside]]^2 * norm_ivar[ss[inside]])

splog, '1-d bspline is fine', $
       min(fwhmfit), max(fwhmfit), chi_med, n_elements(bkpt), $
       format = '(a,f6.2, f6.2, f7.3, i6)'

;;in = where(abs(sigma_x) LT 11 AND mask)

xinf = WHERE(finite(xnew) NE 1, nxinf)
IF nxinf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in trace correction'
    splog, '         Replacing with zeros....'
    xnew = trace_in
ENDIF
inf = WHERE(finite(profile_model) NE 1, ninf)
IF ninf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in object profile'
    splog, '         Replacing with zeros....'
    profile_model[inf] = 0.0
 ENDIF
norm=replicate(1.0,ncol) # total(profile_model,1,/double)
;;norm = total(profile_model)/float(nrow)
IF total(norm) GT 0.0 THEN profile_model = profile_model/norm

title_string = title_string  $
  + ' FWHM range:'  + string(min(fwhmfit), max(fwhmfit), $
                             FORMAT = '(F6.2, F6.2)') $
  + ' S/N:'+ string(sqrt(med_sn2), FORMAT = '(F8.3)') $
  + ' Chi^2'+ string(chi_med, chi_zero, FORMAT = '(F8.3, F8.3)') 
;x1 = sigma_x[in]
;y1 = norm_obj[in]/(pb[in] + (pb[in] EQ 0))
;m1 = full_bsp[in]
;;gd1 = WHERE(norm_ivar[in] GT 0.0, ngd1)
qa_longslit_profile, sigma_x, norm_obj/(pb + (pb EQ 0)), full_bsp $
                     , l_limit, r_limit, ind = ss[inside] $
                     , title = title_string, xrange = PROF_NSIGMA
;;qa_longslit_profile, x1[gd1], y1[gd1], m1[gd1] $
;;                     , l_limit, r_limit, title = title_string
;qa_longslit_profile, sigma_x[in], norm_obj[in]/(pb[in] + (pb[in] EQ 0)), $
;                     full_bsp[in], l_limit, r_limit, title = title_string
;stop
return, profile_model

splog, 'Trying 2-d bspline, 1-d gave: ', chi_med

sf = si[where(abs(res_mode) LT 5.*chi_med AND $
              abs(norm_obj[si] - profile_model[isub[si]]) LT 0.4 AND $
              sigma_x[si] GT l_limit AND sigma_x[si] LT r_limit, nsf)]

profile_2d = bspline_iterfit(sigma_x[sf], norm_obj[sf] $
                             , invvar = norm_ivar[sf], bkspace = 0.4 $
                             , outmask = outmask, yfit = profile_fit2d $
                             , npoly = 5, x2 = xtemp[sf] $
                             , /groupbadpix, maxrej = 5 $
                             , upper = 10, lower = 10, /silent)

;    lh = max(where(profile_fit2d LT 0 AND sigma_x[si] LT 0.)) > 0
;    rh = max(where(reverse(profile_fit2d) LT 0 AND sigma_x[sr] GT 0.)) > 0
;    profile_fit2d = profile_fit2d * (sigma_x[si] GE sigma_x[si[lh]] AND $
;                                 sigma_x[si] LE sigma_x[sr[rh]])


profile_model[*] = 0.
profile_model[isub[sf]] = profile_fit2d
res_2d = (norm_obj - profile_model[isub])*sqrt(norm_ivar)
mask = long(isub*0 + 1)
mask[sf] = mask[sf]*outmask
good2d = where((mask[si] EQ 1) AND (norm_ivar[si] GT 0), ngood2d)
chi_med2d = median(res_2d[si[good2d]]^2)

                                
;;  Need a check to make sure 2-d profile is OK
;;
splog, 'Made it to 2d bspline', chi_med2d
if chi_med2d GT 0.8*chi_med then begin
    splog, '1-d bspline was better?'
    profile_model[*] = 0
    profile_model[isub[si]] = profile_fit
endif

xinf = WHERE(finite(xnew) NE 1, nxinf)
IF nxinf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in trace correction'
    splog, '         Replacing with zeros....'
    xnew = trace_in
ENDIF
inf = WHERE(finite(profile_model) NE 1, ninf)
IF ninf NE 0 THEN BEGIN
    splog, 'WARNING: Nan pixel values in object profile'
    splog, '         Replacing with zeros....'
    profile_model[inf] = 0.0
ENDIF
norm=replicate(1.0,ncol) # total(profile_model,1,/double)
;;norm = total(profile_model)/float(nrow)
IF total(norm) GT 0.0 THEN profile_model = profile_model/norm
return, profile_model 

END
