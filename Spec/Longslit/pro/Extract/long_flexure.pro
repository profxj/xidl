;+
; NAME:
;   long_flexure
;
; PURPOSE:
;   Calculate the flexure shift from a sky-line spectrum and then
;   shift the spectrum onto the sky spectrum.
;
; CALLING SEQUENCE:
; long_flexure, struct, skyfile, QAFILE = QAFILE
;
; INPUTS:
;  struct -- Structure containing the extracted spectrum
;  skyfile -- File containing the sky spectrum
;  flg_skyfile --  Flag indicating whether multiple skyfiles exist
;                 (0=No, 1=Yes)
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
; QAFILE= -- Filename for QA output
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
;   11-Mar-2005  Written by Joe Hennawi
;   10-Jun-2011  Modified by Gabor Worseck, rejection of pixels close
;                to masked pixels to correct for partially masked sky lines
;  ----------------------------------------------------------------------------
pro long_flexure1, struct, skyfile, flg_skyfile, QAFILE = QAFILE $
                   , MAXFLEX = MAXFLEX $
                   , FITWIDTH = FITWIDTH, FLG_FLEX=flg_flex $
                   , AVG_SHIFT = AVG_SHIFT, SHIFT_TOL = SHIFT_TOL

; restore archived sky
IF KEYWORD_SET(QAFILE) THEN BEGIN
    dfpsplot, qafile, /color
    sci_temp = strcompress(strsplit(qafile, '/', /extract), /rem)
    sci_name = repstr(sci_temp[n_elements(sci_temp)-1L], '-flex.ps', '')
ENDIF ELSE sci_name = 'sci-QA'
IF NOT KEYWORD_SET(MAXFLEX) THEN MAXFLEX = 10L
IF NOT KEYWORD_SET(FITWIDTH) THEN FITWIDTH = 3L


flg_flex = 1
a = findfile(skyfile+'*', count=nfil)
if nfil EQ 0 then begin
    flg_flex = 0
    return
endif

if flg_skyfile EQ 0 then begin
   restore, skyfile
   wave_ref1 = wave_calib
   sky_ref1  = sky_calib
endif

nobj = n_elements(struct)
ny = n_elements(struct[0].wave_box)

FOR iobj = 0L, nobj-1L DO BEGIN
    wave_obj1 = struct[iobj].wave_box
    if max(wave_obj1) LE 1e-5 then continue ;; Bad extraction

    ;; Multiple sky files?
    if flg_skyfile EQ 1 then begin
       gdwv = where(wave_obj1 GT 0.)
       specific_skyfile = long_grab_skyfile(skyfile, [min(wave_obj1[gdwv],max=mx),mx])
       restore, specific_skyfile
       wave_ref1 = wave_calib
       sky_ref1  = sky_calib
    endif

    ;; mask sky so that bad pixels dont break cross correlation
    ;; GW: Mask pixels around bad pixels (e.g., CCD gaps) to correct
    ;;     for partial coverage of sky lines
    maxfwhm = max(struct[iobj].FWHMFIT)
    growmask = smooth(struct[iobj].mask_box,5*maxfwhm)
    ;;MF 2015 - growmask may clip lines... be careful with this!
    ;;Problems found with LRISBLUE 400 - if so, use line below
    ;;sky_obj1  = struct[iobj].sky_box
    sky_obj1  = struct[iobj].sky_box*growmask
    sky_obj_ivar1 = (sky_obj1 GT 0.0)/(sky_obj1 + (sky_obj1 LE 0.0))

;   define wavelength limits
    min_wave = min(wave_ref1) > min(wave_obj1)
    max_wave = max(wave_ref1) < max(wave_obj1)
    keep_obj = WHERE(wave_obj1 GE min_wave AND wave_obj1 LE max_wave, nkeep)
    IF nkeep EQ 0 THEN BEGIN
       struct[iobj].FLX_SHFT_WAV = 9999.0 
       splog, 'No overlap with sky spectrum. Skipping this slit'
       CONTINUE
    ENDIF
    wave = wave_obj1[keep_obj]
    sky_obj = sky_obj1[keep_obj]
    sky_obj_ivar = sky_obj_ivar1[keep_obj]
;   interpolate the reference sky spectrum onto object wavelengths 
    sky_ref = interpol(sky_ref1, wave_ref1, wave, /quad)
;   Normalize spectra to unit average sky count
    norm = (total(sky_obj)/double(nkeep))
    sky_obj = sky_obj/norm
    sky_obj_ivar = sky_obj_ivar*norm^2
    sky_ref = sky_ref/(total(sky_ref)/double(nkeep))

    diff_wave = max(wave)-min(wave)
    xvector = (wave-min(wave))/diff_wave
    
;; I think the scaling and may be overkill since these are spectra from the 
;; same instrument/setup. 

;   Scale the object spectrum to be the same as the reference spectrum
;    solve_poly_ratio, xvector, sky_obj, sky_ref, replicate(1.0, nkeep) $
;                      , npoly = 3, nback = 2 $
;                      , yfit = sky_obj_scale, ymult = ymult, yadd = yadd
    
    ;; took out solve_poly_ratio scaling
    sky_obj_scale = sky_obj
;   Subtract off bspline spectrum from each 
    obj_set = bspline_iterfit(wave, djs_median(sky_obj_scale, width = 5.0 $
                                               , boundary = 'reflect') $
                              , upper = 3.0, lower = 3.0 $
                              , nbkpts = 20, yfit = sky_obj_cont $
                              , /silent)
    ref_set = bspline_iterfit(wave, djs_median(sky_ref, width = 5.0D $
                                               , boundary = 'reflect') $
                              , upper = 3.0, lower = 3.0 $
                              , nbkpts = 20, yfit = sky_ref_cont, /silent)
    sky_ref_corr = (sky_ref - sky_ref_cont) ;; > (-100.0D)
    sky_obj_corr  = (sky_obj_scale  - sky_obj_cont) ;; > (-100.0D)
    ;;outmask = 0L
;   mask severe outlier pixels since these will break the cross-correlation
    ;qdone = djs_reject(sky_ref_corr, sky_obj_corr, outmask = outmask $
    ;                   , sigma = sky_ref, upper = 20.0, lower = 20.0)
    ;IF total(outmask) NE n_elements(outmask) THEN $
    ;  splog, long(n_elements(outmask) - total(outmask)) $
    ;  , FORMAT = '(%"WARNING: Masked %d bad sky pixels")' 
    
    ;;MF 2015 - to see what's going on with flex use code below
    ;plot, wave, sky_ref_corr, title=string(iobj+1), xrange=[5000,6000]
    ;oplot, wave, sky_obj_corr, line=2, color=fsc_color('red')
    ;oplot, wave, sky_obj, line=2, color=fsc_color('blue')
    ;plot, struct[iobj].wave_box,  struct[iobj].sky_box
    ;oplot, struct[iobj].wave_box,  struct[iobj].sky_box*growmask, color=255
    ;stop

;   cross-correlate spectra
    nsamp = 50
    step = lindgen(2*MAXFLEX*nsamp) - MAXFLEX*nsamp 
    ;; sharpness filter the spectrum to mask any large pixel values
    ;; that will break the cross-correlation. 
    sharpchi =  (sky_obj_scale-djs_median(sky_obj_scale, width = 3 $
                                          , bound = 'reflect'))  $
      *sqrt(sky_obj_ivar >  0.0)
    mask_obj = sky_obj_scale LE 1.0d6 AND sky_obj_scale GT -100.0 AND $
      finite(sky_obj_scale)
    spix1 = WHERE(mask_obj, ns1)
    IF ns1 GT 0 THEN djs_iterstat, sharpchi[spix1], mean = mean1 $
      , median = median1, sigma = sigma1 $
    ELSE sigma1 = 0.0
    spix2 = WHERE(mask_obj AND abs(sharpchi) GT ((0.2*sigma1) > 0.1D), ns2)
    IF ns2 GT 0 THEN djs_iterstat, sharpchi[spix2] $
      , mean = mean2, median = median2, sigma = sigma2
    IF ns2 EQ 0 OR ns1 EQ 0 THEN BEGIN
        struct[iobj].FLX_SHFT_WAV = 9999.0 
        splog, 'Problem with sharpness filtering or sky masking. Skipping this object.'
    ENDIF ELSE BEGIN
        ;Special caveat for the 5577 sky line; gets masked out by sharpchi
        wave5577 = wave GT 5560 AND wave LT 5590
        sharpmask = sharpchi LT (median2 + 20.0*sigma2) OR wave5577
        sharpmask = long_grow_mask(sharpmask, 3)
        objmask = mask_obj AND sharpmask
        refmask = sky_ref LT 1.0d6 AND sky_ref GT -100.0 AND finite(sky_ref)
        bothmask = objmask AND refmask
        goodpix = WHERE(bothmask, ngood)
        ;; smooth both with the same masking of pixels 
        if ngood gt 0 then begin ;;;;KHRR ADDED
           incorr1 = ivarsmooth(sky_obj_corr, double(bothmask), 3)
           incorr2 = ivarsmooth(sky_ref_corr, double(bothmask), 3)
           corr = c2_correlate(rebin(incorr1[goodpix], ngood*nsamp) $
                               , rebin(incorr2[goodpix], ngood*nsamp) $
                               , step, /double)
           max_corr = max(corr, jmax)
           xpeak = step[jmax]/double(nsamp)
           struct[iobj].FLX_SHFT_WAV = xpeak
        ;; On the second pass once the average shift is computed, apply 
        ;; shifts and create QA file
        endif else begin  ;;KHRR ADDED THIS WHOLE PART
           splog, 'Problem with sharpness filtering or sky masking. KHRR!'
           struct[iobj].FLX_SHFT_WAV = 9999.0 
           xpeak = 9999.0
        endelse 
    ENDELSE
    IF N_ELEMENTS(AVG_SHIFT) GT 0 THEN BEGIN
        IF n_elements(struct) GT 1 THEN BEGIN
            fsuffix = '-' + strtrim(struct[iobj].slitid, 2) $
              + '-' + strtrim(struct[iobj].objid, 2)
            title_string = sci_name + fsuffix
        ENDIF ELSE title_string = sci_name
        IF (abs(xpeak-avg_shift) GT SHIFT_TOL) THEN BEGIN
            xshift = avg_shift
            errcode = 1
            struct[iobj].FLX_SHFT_WAV = xshift
        ENDIF ELSE BEGIN
            xshift = xpeak
            errcode = 0
        ENDELSE
        print, title_string, xshift $
               , FORMAT = '(%"%s flexure shift: %6.2f pixels")' 
        wave_box_old = struct[iobj].WAVE_BOX
        wave_opt_old = struct[iobj].WAVE_OPT
        wave_box_new = interpol(wave_box_old, dindgen(ny) $
                                , dindgen(ny) + xshift)
        wave_opt_new = interpol(wave_opt_old, dindgen(ny) $
                                , dindgen(ny) + xshift)
        struct[iobj].WAVE_BOX   = wave_box_new
        struct[iobj].WAVE_OPT   = wave_opt_new
        ;stop
        IF KEYWORD_SET(QAFILE) THEN $
           long_flex_qa, wave_box_new[keep_obj], wave, sky_obj $
                         , sky_ref, step/double(nsamp), corr, xshift, xpeak, max_corr $
                         , title_string $
                         , maxflex = maxflex, errcode = errcode
     ENDIF
ENDFOR

IF KEYWORD_SET(QAFILE) THEN dfpsclose

return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO long_flexure, struct, skyfile, flg_skyfile, QAFILE = QAFILE, MAXFLEX = MAXFLEX $
                  , FITWIDTH = FITWIDTH, MAXGOOD=MAXGOOD

if  N_params() LT 3  then begin 
   print,'Syntax - ' + $
         'long_flexure, struct, skyfile, flg_skyfile [v1.1]'
   return
endif 

IF NOT KEYWORD_SET(MAXFLEX) THEN MAXFLEX = 10L
if not keyword_set(MAXGOOD) then MAXGOOD = (12. > MAXFLEX) ;; JXP:  2010 Aug 4
;; First pass used to compute average shift. Do not apply shifts
long_flexure1, struct, skyfile, flg_skyfile, maxflex = maxflex, fitwidth = fitwidth, $
               FLG_FLEX = flg_flex 
if flg_flex EQ 0 then return
xpeak = struct.flx_shft_wav 
igood = WHERE(abs(xpeak) LE MAXGOOD, npeak)

IF (npeak LE 2) THEN sigrej = 1.0 $ 
;; Irrelevant for only 1 or 2 files
else if (npeak EQ 3) then sigrej = 1.1 $
else if (npeak EQ 4) then sigrej = 1.3 $
else if (npeak EQ 5) then sigrej = 1.6 $
else if (npeak EQ 6) then sigrej = 1.9 $
else sigrej = 2.0

IF npeak EQ 0 THEN BEGIN
    avg_shift = 0.0D
    shift_tol = 2.0D
ENDIF ELSE IF npeak LE 2 THEN BEGIN
    avg_shift = total(xpeak[igood])/double(npeak)
    shift_tol = 2.0D 
ENDIF ELSE BEGIN
    djs_iterstat, xpeak[igood], mean = mean, median = median, sigma = sigma $
                  , sigrej = sigrej
    avg_shift = median
    shift_tol = 10.0D*sigma >  2.0D
ENDELSE

;; Second pass, apply shifts. If they are inconsisent with the 
;; average, then apply the average shift instead. 
long_flexure1, struct, skyfile, flg_skyfile, maxflex = maxflex, fitwidth = fitwidth $
               , avg_shift = avg_shift, shift_tol = shift_tol, qafile = qafile
        
RETURN
END


