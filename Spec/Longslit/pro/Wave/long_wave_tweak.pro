;+
; NAME:
;  long_reidentify
;  Version 1.1
;
; PURPOSE:
;  Given an archived set of wavelength solution cross correlate to find the 
;  the closest matching arc in the archive and use this as the reference 
;  spectrum, and refernce shift in reidentification. 
;
; CALLING SEQUENCE:
;  .run LONG_WAVE_TWEAK
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
; PROCEDURES CALLED:
;
;
; REVISION HISTORY:
;   25-May-2006  Written by J. Hennawi 
;-
;------------------------------------------------------------------------------

;; Restore the arc from the dump
restore, '/Volumes/scr0/PAIRS/REDUX/lris_redux/Jan_2008/long_1.0/Red600/wave-lred0021.sav'
badslit = 1 ;; Slit number that failed
arc_obj = arc1d[*, badslit-1L]
;; Restore the archive file
restore,  getenv('LONGSLIT_DIR') + $
          '/calib/linelists/lris_red_600-10000_d680_mswave_8051.sav'

ny = n_elements(arc_obj)
nyby2 = ny/2L

bin_ratio = 1
ncalib = n_elements(calib)
;; If binning is different than archive, rebin the archived arc
IF bin_ratio NE 1 THEN archive_arc = rebin(archive_arc, ny, ncalib)
shift_vec = findgen(ncalib)
corr_vec = fltarr(ncalib)
MXSHIFT = 300
step = lindgen(2*MXSHIFT) - MXSHIFT

;; gaussian kernel for smoothing residual arc
sig_res = 4.0D
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])

;; Loop over each member of the calibration vector to find the arc with the
;; largest cross-correlation coefficient. 
FOR j = 0L, ncalib-1L DO BEGIN
    arc_ref = archive_arc[*, j]
    arc_obj_sm = djs_median(arc_obj, width = 30, boundary = 'reflect')
    arc_ref_sm = djs_median(archive_arc[*, j], width = 30, boundary = 'reflect')
    arc_obj_sm = convol(arc_obj_sm, kernel, /edge_truncate)
    arc_ref_sm = convol(arc_ref_sm, kernel, /edge_truncate)
    arc_obj_corr = (arc_obj - arc_obj_sm) <  1.0d6
    arc_ref_corr = (arc_ref - arc_ref_sm) <  1.0d6
    corr = c2_correlate(arc_obj_corr, arc_ref_corr, step, /doub)
    xpeak = long_find_nminima(-corr, step, nfind = 1, minsep = 1 $
                              , ypeak = ypeak $
                              , npeak = npeak, errcode = errcode $
                              , width = 10)
    IF xpeak GT max(step) OR xpeak LT min(step) OR YPEAK GT 1.0 THEN BEGIN
        corr_vec[j] = 0.0
        shift_vec[j] = 0.0
    ENDIF ELSE BEGIN
        corr_vec[j] = -ypeak
        shift_vec[j]  =  -xpeak
    ENDELSE
ENDFOR
best_ind = WHERE(corr_vec GT 0.8 AND abs(shift_vec) LE 30.0, nbest)
IF nbest GT 0 THEN BEGIN
    max_corr = max(corr_vec[best_ind], jmax)
    shift = shift_vec[best_ind[jmax]]
    fit_ref = calib[best_ind[jmax]]
    arc_ref = archive_arc[*, best_ind[jmax]]
ENDIF ELSE BEGIN
    splog, 'No good correlated ref arc within 30 pixels'
    splog, 'Using best correlation'
    max_corr = max(corr_vec, mind)
    shift = shift_vec[mind]
    fit_ref      = calib[mind]
    arc_ref      = archive_arc[*, mind]
ENDELSE
FORMAT = '(%"BEST XCORR: %7.2f  SHIFT=%7.2f pix")'
splog, FORMAT = FORMAT, max_corr, shift

IF bin_ratio NE 1 THEN BEGIN
    ffit = *(fit_ref.ffit)
    fit_ref.nrm[0] = fit_ref.nrm[0]/bin_ratio
    fit_ref.nrm[1] = fit_ref.nrm[1]/bin_ratio
ENDIF

x_identify, arc_obj, calib_temp $
            , mfitstr = fit_ref, mspec = arc_ref, mshift = shift $
            , Xsize = 1200, ysize = 600 $
            , linelist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/lris_red_600.lst'

wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
wave_cen = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]

IF TAG_EXIST(calib_temp, 'WAVE_CEN') EQ 0 THEN $
  calib_temp = struct_addtags(calib_temp, create_struct('WAVE_CEN', wave_cen)) $
ELSE calib_temp.WAVE_CEN = wave_cen 
IF TAG_EXIST(calib_temp, 'DISP_CEN') EQ 0 THEN $
  calib_temp = struct_addtags(calib_temp, create_struct('DISP_CEN', disp_cen)) $
ELSE calib_temp.DISP_CEN = disp_cen

stop

archive_arc = transpose([transpose(archive_arc), transpose(arc_obj)])
calib = [calib, calib_temp]


spawn, 'mv -f ' + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600-10000_d680_mswave_8051.sav   ' $
       +  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600-10000_d680_mswave_8051_old.sav' 

save, archive_arc, calib $
      , file = getenv('LONGSLIT_DIR') + $
      '/calib/linelists/lris_red_600-10000_d680_mswave_8051.sav'
;; 


END


