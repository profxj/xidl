PRO NIRI_PLOTSCI, scifile, NSMOOTH = NSMOOTH, HARD_PS = HARD_PS, BOX = BOX $
                  , LUCI = LUCI, SINFONI = SINFONI, SOFI = SOFI

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
    Set_Plot, "PS"
    Device, _Extra = keywords
ENDIF

atm_file =  getenv('LONGSLIT_DIR') + '/calib/extinction/atm_trans_am1.0.dat'
rdfloat, atm_file, lam_atm, trans, skip = 2

micron =  cgGreek('mu') + 'm'
FONT = 0
!P.FONT = 0
xlab = Textoidl('\lambda', FONT = FONT) + ' (' + micron + ')'
ylab = 'Counts  (e-)'

IF KEYWORD_SET(LUCI) OR KEYWORD_SET(SOFI) THEN INDX = 3 $
ELSE INDX = 4 
final_struct = mrdfits(scifile, INDX)
nsci = n_elements(final_struct)

IF NOT KEYWORD_SET(NSMOOTH) THEN NSMOOTH = 3L

IF KEYWORD_SET(BOX) THEN BEGIN
    wave_opt = final_struct.WAVE_OPT
    flux_opt = final_struct.FLUX_OPT
    ivar_opt = final_struct.IVAR_OPT
ENDIF ELSE BEGIN
    wave_opt = final_struct.WAVE_BOX
    flux_opt = final_struct.FLUX_BOX
    ivar_opt = final_struct.IVAR_BOX
ENDELSE

sci_temp = strcompress(strsplit(scifile, '/', /extract), /rem)
sci_name = repstr(sci_temp[n_elements(sci_temp)-1L], '.fits', '')

IF median(wave_opt[*, 0]) GT 1d4 THEN wave_scale = 1d4 ELSE wave_scale = 1.0d

FOR ifile = 0L, nsci-1 DO BEGIN
    wave = wave_opt[*, ifile]/wave_scale
    flux = flux_opt[*, ifile]
    ivar = ivar_opt[*, ifile]
    sig = (ivar NE 0.0)/sqrt(ivar + (ivar LE 0.0))
    djs_iterstat, flux, sigrej = 3.0, sigma = sigma, median = median
    igood = where(wave GT 0.0, ngd)
    IF ngd GT 0 THEN xrange = [min(wave[WHERE(wave GT 0.0)]), max(wave)] $
    ELSE xrange = [0.0, 1.0]
    yrange = [-1.0*sigma, median + 5.0*sigma]
    IF ifile EQ 0 THEN  long_multipanel, nplots = nsci, /portrait $
    ELSE long_multipanel, /advance, position = p
    fsuffix = '-' + strtrim(final_struct[ifile].objid, 2)
    title = sci_name + fsuffix
    plot, XRANGE, YRANGE, XSTYLE = 1, YSTYLE = 1 $
          , /NODATA, xtitle = xlab $
          , ytitle = ylab, title = title $
          , xthick = 3.0, ythick = 3.0, charthick = 3.0, yminor = 2 $
          , charsize = 1.5, FONT = FONT, position = position $
          , xtickformat = '(F6.4)'
    IF KEYWORD_SET(NSMOOTH) THEN flux = smooth(flux, NSMOOTH)
    IF ngd EQ 0 THEN CONTINUE
    oplot, wave, flux, psym = 10, THICK = THICK $
           , color = djs_icolor('black')
    oplot, wave, sig, psym = 10, col = djs_icolor('red') $
           , thick = thick
    scale = 0.40
    oplot, lam_atm, (1.0-scale)*yrange[1] + $
           0.90*scale*yrange[1]*trans, thick = thick $
           , color = djs_icolor('blue'), psym = 10
ENDFOR

IF KEYWORD_SET(HARD_PS) AND NOT CANCELLED THEN BEGIN
    Device, /Close_File
    Set_Plot, thisDevice
ENDIF

RETURN
END




