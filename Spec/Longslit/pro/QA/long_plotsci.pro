;+
; NAME:
;   long_plotsci
;
; PURPOSE:
;   Generate QA for the science extractions
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
PRO LONG_PLOTSCI, scifile, index = index, NSMOOTH = NSMOOTH, HARD_PS = HARD_PS $
                  , wave_opt = wave_opt2, flux_opt = flux_opt2 $
                  , ivar_opt = ivar_opt2, flux_box = flux_box2 $
                  , ivar_box = ivar_box2, xpos = xpos2 $
                  , fwhmfit = fwhmfit2, arc_fwhm = arc_fwhm2 $
                  , NOPLOT = NOPLOT $
                  , final_struct = final_struct

cleanplot,/silent
!p.multi = 0
THICK = 2.0
IF N_ELEMENTS(FONT) EQ 0 THEN FONT = -1
cancelled = 0
IF KEYWORD_SET(HARD_PS) THEN BEGIN
    thick = 2.5
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


xlab = Textoidl('\lambda')
ylab = 'S/N Ratio'

final_struct = xmrdfits(scifile, 5)
nsci = n_elements(final_struct)

wave_opt = final_struct.WAVE_OPT
flux_opt = final_struct.FLUX_OPT
ivar_opt = final_struct.IVAR_OPT
flux_box = final_struct.FLUX_BOX
ivar_box = final_struct.IVAR_BOX
xpos     = final_struct.XPOS
fwhmfit  = final_struct.FWHMFIT
arc_fwhm_fit = final_struct.ARC_FWHM_FIT

sci_temp = strcompress(strsplit(scifile, '/', /extract), /rem)
sci_name = repstr(sci_temp[n_elements(sci_temp)-1L], '.fits', '')

FOR ifile = 0L, nsci-1 DO BEGIN
    wave_opt1 = wave_opt[*, ifile]
    flux_opt1 = flux_opt[*, ifile]
    ivar_opt1 = ivar_opt[*, ifile]
    flux_box1 = flux_box[*,ifile]
    ivar_box1 = ivar_box[*,ifile]
    
    ind = WHERE(wave_opt1 GT 0.0 AND finite(wave_opt1) EQ 1 and $ ; jm11jun08ucsd
      finite(flux_opt1) eq 1 and finite(ivar_opt1) eq 1, ngood)
    IF ngood EQ 0 THEN CONTINUE  
    xrange = [min(wave_opt1[ind]), max(wave_opt1[ind])]
    snr1 = flux_opt1*sqrt(ivar_opt1)
    djs_iterstat, snr1[ind], sigrej=8.0, mask=msk
    gd = where(msk,ngd)
    if (ngd eq 0L) then message, 'Problem making the QAplot!'
    yrange = [0.0, 1.15*max(snr1[gd])]
;   yrange = [0.0, min([1.15*max(flux_opt1*sqrt(ivar_opt1)), 300.0])] ; jm11jun08ucsd
    IF ifile EQ 0 THEN  long_multipanel, nplots = nsci, /portrait $
    ELSE long_multipanel, /advance, position = p
    fsuffix = '-' + strtrim(final_struct[ifile].slitid, 2) $
      + '-' + strtrim(final_struct[ifile].objid, 2)
    title = sci_name + fsuffix
    plot, XRANGE, YRANGE, XSTYLE = 1, YSTYLE = 1, $
          /NODATA, xtitle = xlab $
          , ytitle = ylab, title = title $
          , xthick = 3.0, ythick = 3.0, charthick = 3.0, yminor = 2 $
          , charsize = 1.5, FONT = FONT, position = position $
          , xtickformat = '(I5)'
    IF KEYWORD_SET(NSMOOTH) THEN flux_opt1 = smooth(flux_opt1, NSMOOTH)
    oplot, wave_opt1, flux_opt1*sqrt(ivar_opt1), psym = 10, $
           THICK = THICK
;   oplot, wave_opt1, flux_box1*sqrt(ivar_box1), psym=10, color=djs_icolor('blue'), thick=thick
ENDFOR

IF KEYWORD_SET(INDEX) THEN BEGIN
    wave_opt2 = wave_opt[*, index] 
    flux_opt2 = flux_opt[*, index]
    ivar_opt2 = ivar_opt[*, index]
    flux_box2 = flux_box[*, index]
    ivar_box2 = ivar_box[*, index]
    xpos2     = xpos[*, index]
    fwhmfit2  = fwhmfit[*, index]
    arc_fwhm2 = arc_fwhm_fit[*, index]
ENDIF

IF KEYWORD_SET(HARD_PS) AND NOT CANCELLED THEN BEGIN
    Device, /Close_File
    Set_Plot, thisDevice
 ENDIF

;long_multipanel, /OFF
cleanplot,/silent

RETURN
END




