;+
; NAME:
;   long_waveqa
;
; PURPOSE:
;   Generate QA for the wavelength solution.
;
; CALLING SEQUENCE:
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
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and 
;                Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
PRO LONG_WAVEQA, wlines, fit, arc1d, rejpt, islit $
                 , fwhm, fwhmset, fwhmmask, wstruct, qafile, panic = panic $
                 , BADFIT = BADFIT


IF n_elements(wlines.WAVE) GT 1 THEN BEGIN 
   IF median(wlines.WAVE) GE 1d4 OR median(wlines.WAVE) LE 2.6 THEN IR = 1
ENDIF ELSE IR = 1 ;; this is a hack for now to fix bad IR fits
;;IF strmatch(wstruct.instrument, '*SOFI*') THEN IR = 0
micron = '\mu' + 'm'
;;micron = cgGreek('mu') + 'm'
;angstrom = STRING(197B)
;angstrom = cgSymbol("angstrom", /ps)
;;angstrom = STRING(197B)
angstrom = '\AA'

IF KEYWORD_SET(IR) THEN units = micron ELSE units = angstrom
IF NOT KEYWORD_SET(IR) OR (strmatch(wstruct.instrument, '*SOFI*') AND wstruct.BAND EQ 'H+K') $
THEN ylog = 1 $
ELSE ylog = 0


charsize = 1.1D
!p.multi = [0, 1, 2]
title_string = 'Slit#' + strcompress(string(islit+1), /rem)
IF KEYWORD_SET(PANIC) THEN title_string = title_string + '  PANIC'
IF KEYWORD_SET(BADFIT) THEN title_string = title_string + '  BADFIT'

ny = n_elements(arc1d)
nlines = n_elements(wlines)

wvfit  = x_calcfit(dindgen(ny), fitstr = fit)
linfit = x_calcfit(wlines.PIX, fitstr = fit)

ffit = 1.0*(*fit.ffit)
nrm = fit.nrm
func = fit.FUNC
nord = fit.NORD
; zero non-linear terms
IF func EQ 'POLY' AND nord GT 1 THEN ffit[0, 2:*] = 0.0 $
ELSE IF func EQ 'CHEBY' AND nord GT 2 THEN ffit[2:*] = 0.0

wvfit_lin = x_calcfit(dindgen(ny), func, ffit $
                      , nrm = nrm, nord = nord)
linfit_lin    = x_calcfit(wlines.PIX, func, ffit $
                          , nrm = nrm, nord = nord)
dwv = wvfit - shift(wvfit, 1)
dwv[0] = dwv[1]
dwv_med = djs_median(dwv)
rms = fit.rms/abs(dwv_med)

fwhmcolor = replicate('blue', nlines)
IF (min(fwhmmask) EQ 0) THEN fwhmcolor[where(fwhmmask EQ 0)] = 'red'
traceset2xy, fwhmset, findgen(ny), fwhmfit


nl_wdev = wlines.WAVE-linfit_lin
all_nl_wdev = wvfit - wvfit_lin
IF KEYWORD_SET(YLOG) THEN yrange = [(0.05*djs_median(arc1d)) > 1.0 $
               , 10.0*max(djs_median(arc1d, width = 5))] $
ELSE yrange = [(0.1*djs_median(arc1d)) > 1.0, 1.3*max(djs_median(arc1d, width = 5))] 
wrange = [min(wvfit), max(wvfit)]

stop
djs_plot, wvfit, abs(arc1d), xrange = wrange $
          , yrange = yrange $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + units + ')' $
          , ytitle = 'Arc Spectrum' $
          , title = title_string $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize, ylog = ylog ;;, ylog = (KEYWORD_SET(IR) EQ 0)
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'blue' $
    ELSE this_color = 'red'
    IF KEYWORD_SET(YLOG) THEN BEGIN 
       xy = 0.5d
       xy2 = 0.1
       xy3 = 0.03
    ENDIF ELSE BEGIN
       xy = 0.95
       xy2 = 0.75
       xy3 = 0.7
    ENDELSE
    IF KEYWORD_SET(IR) THEN w_string = strcompress(string(wlines[k].WAVE, FORMAT = '(F7.5)') $
                                                   , /rem)  $
    ELSE  w_string = strcompress(string(wlines[k].WAVE, FORMAT = '(F8.2)'), /rem)
    djs_xyouts, wlines[k].WAVE, xy*yrange[1] $
                , w_string, color = this_color $
                , charthick = 3.0, charsize = 0.7 $
                , orientation = 270.0
    djs_arrow, wlines[k].WAVE, xy2*yrange[1], wlines[k].WAVE $
               , xy3*yrange[1], /DATA, color = this_color, thick = 3.0 $
               , hsize = 0.0, /solid
ENDFOR

IF KEYWORD_SET(IR) THEN dwv_str = $
   '\Delta\lambda = ' + strcompress(string(dwv_med*1d4, format = '(F7.4)') $
                                    , /rem) +  ' '  + angstrom $
ELSE dwv_str = $
   '\Delta\lambda = ' + strcompress(string(dwv_med, format = '(F7.4)') $
                                    , /rem) +  ' '  + units
rms_str = 'RMS = ' +  $
  strcompress(string(rms, format = '(F7.4)'), /rem) + ' (pix)'

IF KEYWORD_SET(IR) THEN $
   cen_str = '\lambda_{cen} = ' + string(wvfit[ny/2], FORMAT = '(F7.5)') $
             +  ' '  + units $
ELSE cen_str = '\lambda_{cen} = ' + string(wvfit[ny/2], FORMAT = '(F7.1)') $
  +  ' '  + units
  
IF dwv_med LT 0 THEN xrange = [ny, 0] $
ELSE xrange = [0, ny]

djs_plot, findgen(ny), wvfit, xrange = xrange $
          , yrange = [0.98*min(wvfit), 1.02*max(wvfit)] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'pixel' $
          , ytitle =  '\lambda ' + '(' + units + ')' $
          , title = 'Fit' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize

IF dwv_med LT 0 THEN xy_pos = 0.9 $
ELSE xy_pos = 0.1
djs_xyouts, xy_pos*ny, wrange[0] + (wrange[1]-wrange[0])*0.85 $
            ,  dwv_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2
djs_xyouts, xy_pos*ny,  wrange[0] + (wrange[1]-wrange[0])*0.75 $
            , rms_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2
djs_xyouts, xy_pos*ny,  wrange[0] + (wrange[1]-wrange[0])*0.65 $
            , cen_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2

FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'blue' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].PIX], [wlines[k].WAVE] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

djs_plot, wvfit, 0.0*findgen(ny), xrange = wrange $
          , yrange = [-6.0*rms, 6.0*rms] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + units + ')' $
          , ytitle =  'Deviation [pix]' $ 
          , title = 'Residuals' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize
dwv_line = interpol(dwv, wvfit, wlines.WAVE)
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'blue' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].WAVE], [(linfit[k]-wlines[k].WAVE)/dwv_line[k]] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

;     Make plot of nonlinear part of wavelength solution
djs_plot, wvfit, all_nl_wdev/dwv $
          , xrange = wrange $
          , yrange = [min(all_nl_wdev/dwv)-2, max(all_nl_wdev/dwv) + 2] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + units + ')' $
          , ytitle = 'x - x_{lin}' $
          , color = 'default', title = title_string $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.1
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'blue' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].WAVE], [nl_wdev[k]/dwv_line[k]] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

;     Plot fits for FWHM of arc lines
djs_plot, [wlines.PIX], [fwhm] $
          , xrange = xrange, yrange = [-0.1, 12.0] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'Pixel', ytitle = 'FWHM' $
          , title = title_string, psym = 6 $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.0, color = fwhmcolor
djs_oplot, findgen(ny), fwhmfit, thick = 3.0, color = 'black'

!p.multi = 0
RETURN
END
