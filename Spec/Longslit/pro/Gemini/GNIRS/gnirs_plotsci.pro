PRO GNIRS_PLOTSCI, scifile, NSMOOTH = NSMOOTH, HARD_PS = HARD_PS, BOX = BOX

!p.multi = 0
THICK = 2.0
IF N_ELEMENTS(FONT) EQ 0 THEN FONT = -1
cancelled = 0
IF KEYWORD_SET(HARD_PS) THEN BEGIN
    thick = 4.0
    FONT = 0
    !P.FONT = 0
    keywords = create_struct( 'BITS_PER_PIXEL', 8, $
                              'COLOR', 1, $
                              'ENCAPSULATED', 1, $
                              'FILENAME', HARD_PS, $
                              'FONT_SIZE', 12, $
                              'INCHES', 1, $
                              'ISOLATIN1', 1, $
                              'PREVIEW', 0, $
                              'TT_FONT', 0, $
                              'XOFFSET', 1.0, $
                              'XSIZE', 13.00, $
                              'YOFFSET', 1.0, $
                              'YSIZE', 9.00, $
                              'PORTRAIT', 1, $
                              'LANDSCAPE', 0, $
                              'HELVETICA', 1, $
                              'BOLD', 0, $
                              'BOOK', 0, $
                              'DEMI', 0, $
                              'ITALIC', 0, $
                              'LIGHT', 0, $
                              'MEDIUM', 0, $
                              'NARROW', 0, $
                              'OBLIQUE', 0 )
;    keywords = PSConfig(Cancel = cancelled)
;    IF NOT cancelled THEN BEGIN
    thisDevice = !D.Name
;    Set_Plot, "PS"
;    Device, _Extra = keywords
    dfpsplot, hard_ps, /color, /landscape ;_EXTRA = keywords
ENDIF

atm_file =  getenv('LONGSLIT_DIR') + '/calib/extinction/atm_trans_am1.0.dat'
rdfloat, atm_file, lam_atm, trans, skip = 2

xlab = Textoidl('\lambda')
ylab = 'Counts  (e-)'

final_struct = mrdfits(scifile, 4)
nsci = n_elements(final_struct)

IF NOT KEYWORD_SET(NSMOOTH) THEN NSMOOTH = 3L

sci_temp = strcompress(strsplit(scifile, '/', /extract), /rem)
sci_name = repstr(sci_temp[n_elements(sci_temp)-1L], '.fits', '')

order_vec = [3, 4, 5, 6, 7, 8]
norders = 6

nobj = 2L
suffix = ['pos', 'neg']

FOR iobj = 0L, nobj-1L DO BEGIN
    IF iobj NE 0 THEN erase
    long_multipanel, nplots = norders, /portrait 
    FOR iorder = 0L, norders-1L DO BEGIN
        IF KEYWORD_SET(BOX) THEN BEGIN
            wave = final_struct[2*iorder + iobj].WAVE_BOX/1.d4
            flux = final_struct[2*iorder + iobj].FLUX_BOX
            ivar = final_struct[2*iorder + iobj].IVAR_BOX
        ENDIF ELSE BEGIN
            wave = final_struct[2*iorder + iobj].WAVE_OPT/1.d4
            flux = final_struct[2*iorder + iobj].FLUX_OPT
            ivar = final_struct[2*iorder + iobj].IVAR_OPT
        ENDELSE
        sig = (ivar NE 0.0)/sqrt(ivar + (ivar LE 0.0))
        djs_iterstat, flux, sigrej = 3.0, sigma = sigma, median = median
        xrange = [min(wave[WHERE(wave GT 0.0)]), max(wave)]
        yrange = [-1.0*sigma, median + 7.0*sigma]
        IF iorder NE 0 THEN  long_multipanel, /advance, position = p
        fsuffix = '-' + strcompress(string(order_vec[iorder]), /rem) + $
          '-' + suffix[iobj]
        title = sci_name + fsuffix
        plot, XRANGE, YRANGE, XSTYLE = 1, YSTYLE = 1 $
              , /NODATA, xtitle = xlab $
              , ytitle = ylab, title = title $
              , xthick = 3.0, ythick = 3.0, charthick = 3.0, yminor = 2 $
              , charsize = 1.5, FONT = FONT, position = position $
              , xtickformat = '(F6.4)'
        IF KEYWORD_SET(NSMOOTH) THEN flux = smooth(flux, NSMOOTH)
        oplot, wave, flux, psym = 10, THICK = THICK $
               , color = fsc_color('black', 255)
        oplot, wave, sig, psym = 10, col = fsc_color('red', 254) $
               , thick = thick
        scale = 0.40
        oplot, lam_atm, (1.0-scale)*yrange[1] + $
               0.90*scale*yrange[1]*trans, thick = thick $
               , color = fsc_color('blue', 252), psym = 10
    ENDFOR
ENDFOR

IF KEYWORD_SET(HARD_PS) THEN dfpsclose

RETURN
END




