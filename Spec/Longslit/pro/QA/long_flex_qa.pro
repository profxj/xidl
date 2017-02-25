;+
; NAME:
;   long_flex_qa
;
; PURPOSE:
;   Generate QA for the flexure correction
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
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
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------
PRO LONG_FLEX_QA, wave_obj, wave_ref, sky_obj, sky_ref, step, corr $
                  , xshift, xpeak, ypeak, title_string $
                  , MAXFLEX = MAXFLEX, ERRCODE = ERRCODE

cleanplot,/silent
angstrom = STRING(197B)
charsize = 2D
!p.multi = [0, 1, 2]

IF NOT KEYWORD_SET(MAXFLEX) THEN MAXFLEX = 10
IF n_elements(errcode) EQ 0 THEN errcode = 0

green = 'green'
red = 'red'
black = 'default'
if (errcode EQ 0) then color = green else color = red


; sky lines of interest
wave_line  = [3370.0D, 3914.0D, 4046.56, 4358.34, 5577.338, 6300.304 $
              , 7340.885, 7993.332, 8430.174, 8919.610, 9439.660 $
              , 10013.99, 10372.88]
dw         = 20.0D
nline = n_elements(wave_line)
lineid = lonarr(nline)

wave_min = min(wave_obj)
wave_max = max(wave_obj)
FOR j = 0L, nline-1L DO IF $
  (wave_line[j] GT wave_min AND wave_line[j] LT wave_max) THEN lineid[j] = 1L

gdline = WHERE(lineid, ngood)
IF ngood EQ 0 THEN BEGIN
    long_multipanel, ROWS = 2, COLS = 1, /landscape $
                     , position = p $
                     , MARGIN = [0.02, 0.03, 0.02, 0.03] $
                     , OMARGIN = [0.1, 0.06, 0.1, 0.02] 
    djs_plot, [0.0, 1.0], [0.0, 1.0], /nodata $
              , xstyle = 1, ystyle = 1 $
              , xtitle = '\lambda ' + '(' + angstrom + ')' $
              , title = 'NO LINES', color = black $
              , thick = 5.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
              , charsize = 1.5, /noerase, position = p
ENDIF ELSE BEGIN
    FOR j = 0L, ngood -1L DO BEGIN
        iline = gdline[j]
        IF j EQ 0 THEN BEGIN
            IF ngood GT 4 THEN $
              long_multipanel, ROWS = 4, COLS = ceil(ngood/2.0), /landscape $
              , position = p $
              , MARGIN = [0.02, 0.03, 0.02, 0.03] $
              , OMARGIN = [0.04, 0.04, 0.03, 0.02] $
; [left, bottom, right, top]
            ELSE  long_multipanel, ROWS = 2, COLS = ngood, /landscape $
              , position = p $
              , MARGIN = [0.02, 0.03, 0.02, 0.03] $
              , OMARGIN = [0.1, 0.06, 0.1, 0.02] 
        ENDIF ELSE long_multipanel, /advance, position = p $
          , MARGIN = [0.04, 0.03, 0.02, 0.03] $
          , OMARGIN = [0.05, 0.06, 0.05, 0.02] 
        wrange = wave_line[iline] + [-dw, dw]
        pix = WHERE(wave_obj GT wrange[0] AND wave_obj LT wrange[1], npix)
        IF npix GT 0 THEN  $
          yrange = [0.0, 1.2*max([sky_ref[pix], sky_obj[pix]])] $
        ELSE yrange = [0.0, max([sky_ref, sky_obj])]
;   Plot sky spectra 
        djs_plot, wrange, yrange, /nodata $
                  , xstyle = 1, ystyle = 1 $
                  , xtitle = '\lambda ' + '(' + angstrom + ')' $
                  , title = '\lambda = ' +  $
                  string(wave_line[iline], format = '(F7.1)') $
                  , color = black $
                  , thick = 5.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
                  , charsize = 1.5, /noerase, position = p
        djs_oplot, wave_obj, sky_obj, color = color, thick = 3.0 $
                   , psym = 10
        djs_oplot, wave_ref, sky_ref, color = black, thick = 2.0 $
                   , psym = 10
    ENDFOR
ENDELSE


cleanplot, /silent
pos = [0.11, 0.07, 0.93000, 0.4500000]
xrange = [xpeak-3.0D, xpeak+3.0D]
yrange = [0.2D, 1.0D]
    
; Plot correlation and fit 
djs_plot, step, corr, xrange = xrange $
          , yrange = yrange $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'Lag' $
          , ytitle = 'Correlation' $
          , title = title_string $
          , color = black, psym = 10 $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 0.5*charsize, position = pos, /noerase

djs_oplot, [xpeak], [ypeak], psym = 1, symsize = 2.2, color = green $
           , thick = 5.0
;djs_oplot, [xfit], [yfit], thick = 3.0, color = color

corr_str = '\Delta x = ' + strcompress(string(xpeak, format = '(F7.3)') $
                                       , /rem) +  '  pixels' 
djs_xyouts, 0.15, 0.4, corr_str, color = 'default' $
            , charthick = 3.0, charsize = 0.8*charsize, /norm
IF xshift NE xpeak THEN BEGIN
    shift_str = '\Delta x = ' + strcompress(string(xshift, format = '(F7.3)') $
                                            , /rem) +  '  pixels' 
    djs_xyouts, 0.15, 0.3, shift_str, color = 'red' $
                , charthick = 3.0, charsize = 0.8*charsize, /norm
ENDIF

cleanplot,/silent

RETURN
END
