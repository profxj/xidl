PRO NIRI_WAVEQA, wlines, fit, wvfit, arc1d, rejpt, rms, qafile, title_string $
                 , LOG = LOG

IF KEYWORD_SET(LOG) THEN XYFACT = 0.8 $
ELSE XYFACT = 0.95
!p.multi = [0, 1, 3]
dfpsplot, qafile, /color
yrange = [0.8*min(arc1d) >  1.0, 1.2*max(arc1d)]
wrange = [min(wvfit/1.0d4), max(wvfit/1.0d4)]
djs_plot, wvfit/1.0d4, arc1d, xrange = wrange $
          , yrange = yrange $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(\mu' + 'm)' $
          , ytitle = 'Arc spectrum' $
          , title = title_string $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.5, xticks = 10, ylog = log

nlines = n_elements(wlines)
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    w_string = strcompress(string(wlines[k].WAVE/1.0d4, FORMAT = '(F7.5)') $
                           , /rem)
    djs_xyouts, wlines[k].WAVE/1.0d4, XYFACT*yrange[1] $
                , w_string, color = this_color $
                , charthick = 3.0, charsize = 0.7 $
                , orientation = 270.0
    djs_arrow, wlines[k].WAVE/1.0d4, 0.75*yrange[1], wlines[k].WAVE/1.0d4 $
               , 0.7*yrange[1], /DATA, color = this_color, thick = 3.0 $
               , hsize = 0.0, /solid
ENDFOR
angstrom = STRING(197B)
dwv = abs(wvfit - shift(wvfit, 1))
dwv[0] = dwv[1]
dwv_med = djs_median(dwv)
dwv_str = '\Delta\lambda=' + strcompress(string(dwv_med, format = '(F7.4)') $
                                         , /rem) +  ' '  + angstrom
rms_str = 'RMS=' +  strcompress(string(rms, format = '(F7.4)'), /rem) + ' (pix)'

ny = n_elements(wvfit)
djs_plot, findgen(ny), wvfit/1.0d4, xrange = [0, ny], yrange = wrange $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'pixel' $
          , ytitle =  '\lambda ' + '(\mu' + 'm)' $
          , title = 'Fit' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.5
djs_xyouts, 0.1*ny, wrange[0] + (wrange[1]-wrange[0])*0.85 $
            ,  dwv_str, color = 'default' $
            , charthick = 3.0, charsize = 1.5
djs_xyouts, 0.1*ny,  wrange[0] + (wrange[1]-wrange[0])*0.75 $
            , rms_str, color = 'default' $
            , charthick = 3.0, charsize = 1.5

FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].PIX], [wlines[k].WAVE/1.0d4] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = 1.5, color = this_color, psym = 6
ENDFOR

djs_plot, findgen(ny), 0.0*findgen(ny), xrange = wrange $
          , yrange = [-6.0*rms, 6.0*rms] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(\mu' + 'm)' $
          , ytitle =  'Deviation [pix]' $ 
          , title = 'Residuals' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.5

dwv_line = interpol(dwv, wvfit, wlines.WAVE)
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].WAVE]/1.0d4, [(fit[k]-wlines[k].WAVE)/dwv_line[k]] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = 1.5, color = this_color, psym = 6
 ENDFOR

dfpsclose

RETURN
END
