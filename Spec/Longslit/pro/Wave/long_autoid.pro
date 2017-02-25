;+
; NAME:
;   long_autoid
;
; PURPOSE:
;   Use z_arcpairs to auto-id wavelength solution
;
; CALLING SEQUENCE:
;
; INPUTS:
;    arc_obj  -- 1d arc
;    lines    -- arc line structure
;    wstruct  -- wavelength solution parameter structure
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;  Returns a the fit structure for the wavelength solution. 
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
;   17-July-2007  Written by Joe Hennawi
;-
FUNCTION LONG_AUTOID, arc_obj, lines, wstruct, gdfit = gdfit, rejpt = rejpt $
                      , fit_flag = fit_flag, arc_inter=arc_inter


IF NOT KEYWORD_SET(FORDR) THEN FORDR = 9L
ny = n_elements(arc_obj)
nyby2 = ny/2L

;;---------
;; set parameters from wavestruct
nfind = wstruct.NFIND
TOLER = wstruct.toler
pkwdth = wstruct.pkwdth
psig_crude = wstruct.PSIG[0]
disp_guess = wstruct.DISP_GUESS
dr_wave = wstruct.dr_wave
nord_crude = wstruct.NORD[0]
sigrej_crude = wstruct.SIGREJ[0]

if keyword_set(AUTOFIT) then delvarx, autofit  ;; JXP 9/16/2011

x_fndpeaks, arc_obj, line_x1, NSIG = psig_crude, /silent $
            , PKWDTH = pkwdth $
            , MSK = msk, NORDB = fordr, TOLER = TOLER $
            , autofit = autofit, rms = rms, /THIN, /FWEIGHT
sig_lines   = interpolate((arc_obj-autofit)/rms, line_x1)
sig_sort    = reverse(sort(sig_lines))
ap_lines    = sig_sort[0:(nfind-1L) <  (n_elements(line_x1)-1L)]
line_x = line_x1[ap_lines]
line_x = line_x[sort(line_x)]
iqual = WHERE(lines.FLG_QUAL EQ 1)
line_pix = z_arcpairs(line_x, lines[iqual].WAVE, disp_guess $
                      , dr1 = dr_wave)
igood =  where(line_pix GE 0, ngood)
xpix = line_x[igood]
wave = lines[iqual[line_pix[igood]]].WAVE
line_x = 0
autofit = 0
;; Zeroth order fit based on arc_pairs identifications
lines.flg_plt = 0
fit0 = {fitstrct}
fit0.flg_rej = 1 
fit0.niter = 5
fit0.maxrej = ngood/2
fit0.minpt = ngood/2
fit0.hsig = sigrej_crude
fit0.lsig = sigrej_crude
fit0.nord = nord_crude
fit0.FUNC = wstruct.FUNC 
wave_fit0 = x_fitrej(xpix, wave, FITSTR = fit0, REJPT = rejpt0)
wave_vec0 = x_calcfit(dindgen(ny), fitstr = fit0)
wave_cen = wave_vec0[nyby2]
disp_cen = wave_vec0[nyby2]-wave_vec0[nyby2-1L]
fit0 = struct_addtags(fit0, create_struct('WAVE_CEN', wave_cen $
                                          , 'DISP_CEN', disp_cen))

fin_fit = long_witerfit(arc_obj, lines, fit0, wstruct $
                        , gdfit = gdfit, rejpt = rejpt, fit_flag = fit_flag)

RETURN, fin_fit
END
