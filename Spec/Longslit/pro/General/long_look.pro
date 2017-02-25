;+
; NAME:
;   long_look
;
; PURPOSE:
;  Script for viewing data
;
; CALLING SEQUENCE:
;   IDL> cd, 'Science' 
;   IDL> .run long_look
;   After spectrum is plotted, you can plot other objid's with
;   sig =(objstruct[objid-1].IVAR_OPT GT 0.0)/(sqrt(objstruct[objid-1].IVAR_OPT)+ (objstruct[objid-1].IVAR_OPT LE 0.0))
;   wave = objstruct[objid-1].WAVE_OPT
;   flux = objstruct[objid-1].FLUX_OPT
;   x_specplot, flux, sig, wave = wave, inflg = 4, xsize = 1000 $
;            , ysize = 700 
; INPUTS:
;
; OPTIONAL INPUTS:
;    TWO     -- set TWO=1 to view two spectra. 
;
;    FIND_FEATURE   -- set FIND_FEATURE=1 to overplot the positions of a desired at specified redshifts.
;                      Default is to look for [OII] 3727.
;    REDSHIFT   -- set to an array of redshifts.
;                  Intended to be used with FIND_FEATURE flag.
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
;   11-Mar-2005  Written by JH
;   04-Oct-2007  FIND_FEATURE ability added by LKP
;-  
;------------------------------------------------------------------------------

;RESID=1
;FIND_FEATURE=1
FEATURE_WAVE=[3727.5, 3933.68, 3968.49, 4101.74, 4340.47, 4861.34, 4958.91, 5006.84, 6562.82, 6583.0]
FEATURE_NAME=['~[OII]', 'CaII K', 'CaII H', 'Hd', 'Hg', 'Hbeta', '[OIII]', '[OIII]', 'Halpha', '[NII]']
redshift=[0.47285]

fils = findfile('sc*.fits*', count = nfil)
stdfil = findfile('std*.fits*', count = nstd)
IF nstd GT 0 THEN BEGIN 
    fils = [fils, stdfil]
    nfil = nfil + nstd
ENDIF
IF nfil GT 0 THEN path = '../'
IF nfil EQ 0 THEN BEGIN
    fils = findfile('Science/s*.fits*', count = nfil) 
    path = './'
ENDIF
IF nfil EQ 0 THEN message, 'Cannot find files'

file_fits = x_guilist(fils)
psfile = strsplit(file_fits, 'f', /extract) 
spawn, 'gv  ' + psfile[0] + 'ps &'
print, 'Doing ', file_fits

sciimg = xmrdfits(file_fits, 0, scihdr,/silen)
sciivar= xmrdfits(file_fits,1,/silen)
sky_model = xmrdfits(file_fits, 2,/silen)
obj       = xmrdfits(file_fits, 3,/silen)
outmask   = xmrdfits(file_fits, 4,/silen)
objstruct = xmrdfits(file_fits, 5,/silen)
dims = size(sciimg, /dim)
nx = dims[0]
ny = dims[1]

;; RA/DEC
print, 'RA = ', sxpar(scihdr, 'RA')
print, 'DEC = ', sxpar(scihdr, 'DEC')

If TAG_EXIST(objstruct, 'FLX_SHFT_SPA') THEN $
  xshift = objstruct[0].FLX_SHFT_SPA $
ELSE xshift = 0.0

If xshift NE 0.0 AND NOT KEYWORD_SET(NOSHIFT) THEN BEGIN
    slitfile =    findfile(path + 'slits*.fits*', count = nfil)
    IF nfil GT 1 THEN BEGIN
        id = x_guilist(slitfile, indx = indslit)
        slitfile = slitfile[indslit]
    ENDIF
    tset_slits = xmrdfits(slitfile[0], 1)
    tset_slits = long_shiftslits(tset_slits, xshift)
    slitmask = long_slits2mask(tset_slits)
    wavefile = findfile(path + 'wave*.fits*', count = nfil)
    IF nfil GT 1 THEN BEGIN
        id = x_guilist(wavefile, indx = indwv)
        wavefile = wavefile[indwv]
    ENDIF
    pixset  = xmrdfits(wavefile[0], 1)
    wavesvfile = repstr(wavefile[0], '.gz', '')
    wavesvfile =  repstr(wavesvfile, '.fits', '.sav')
    restore, wavesvfile
    piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                             , waveimg = waveimg)
ENDIF ELSE BEGIN
    slitfile =    findfile(path + 'slits*.fits*', count = nfil)
    IF nfil GT 1 THEN BEGIN
        id = x_guilist(slitfile, indx = indslit)
        slitfile = slitfile[indslit]
    ENDIF
    IF nfil GT 0 THEN BEGIN
        slitmask = xmrdfits(slitfile[0], 0) 
        tset_slits = xmrdfits(slitfile[0], 1)
    ENDIF ELSE slitmask = fltarr(nx, ny) + 1.0D
    wavefile = findfile(path + 'wave*.fits*', count = nfil)
    IF nfil GT 1 THEN BEGIN
        id = x_guilist(wavefile, indx = indwv)
        wavefile = wavefile[indwv]
    ENDIF
    IF nfil GT 0 THEN waveimg = xmrdfits(wavefile[0], 0) $
    ELSE waveimg = fltarr(nx, ny)
ENDELSE

IF KEYWORD_SET(RESID) THEN $
  xatv, sqrt(sciivar)*(sciimg-sky_model - obj)*(slitmask GT 0.0) $
  , wvimg = waveimg, min = -6.0, max = 6.0, sig = slitmask $                  
ELSE BEGIN
   inds = WHERE(slitmask GT 0.0 AND outmask GT 0.0, ngood)
   sky, (sciimg-sky_model)[inds], skymode, skysig $
        , /silent
   ;stop
   xatv, (sciimg-sky_model)*(slitmask GT 0.0), wvimg = waveimg $
         , min = -2.5*skysig, max = 25.0*skysig, sig = slitmask
;      , min = -50.0, max = 300.0
ENDELSE

IF KEYWORD_SET(PLOTMASK) THEN BEGIN
    ximg = findgen(nx) # replicate(1.0, ny)
    yimg = replicate(1.0, nx) # findgen(ny)
    maskpix = WHERE(slitmask GT 0 AND outmask EQ 0, nmask)
    PLOTSYM, 0, 1.0, THICK = 1.0
    IF nmask GT 0 THEN xatvplot, ximg[maskpix], yimg[maskpix], psym = 8 $
      , color = 4
ENDIF


IF KEYWORD_SET(tset_slits) THEN BEGIN
; ------
;; Expand slit set to get left and right edge
    traceset2xy, tset_slits[0], rows, left_edge
    traceset2xy, tset_slits[1], rows, right_edge
    xatvplot, left_edge, rows, psym = 3, color = 2
    xatvplot, right_edge, rows, psym = 3, color = 2
ENDIF
xatvplot, objstruct.xpos, objstruct.ypos, psym = 3


objid = strtrim(objstruct.slitid, 2) + '-' + strtrim(objstruct.objid, 2)
nobj = n_elements(objstruct)
FOR j = 0L, nobj-1L DO BEGIN
    xpix = djs_median(objstruct[j].XPOS) - 25.0
    ypix = 1400
    xatvxyouts, xpix, ypix, strcompress(objid[j], /rem) $
               , color = 'red', charsize = 3
ENDFOR

IF NOT KEYWORD_SET(TWO) THEN BEGIN
   splog,'Select object'
   id = x_guilist(objid, indx = ind)
   sig = (objstruct[ind].IVAR_OPT GT 0.0)/(sqrt(objstruct[ind].IVAR_OPT) $
                                           + (objstruct[ind].IVAR_OPT LE 0.0))
   wave = objstruct[ind].WAVE_OPT
   flux = objstruct[ind].FLUX_OPT
   x_specplot, flux, sig, wave = wave, inflg = 4, xsize = 1000 $
            , ysize = 700 
ENDIF ELSE BEGIN 
   splog,'Select first object'
   id = x_guilist(objid, indx = ind)
   sig = (objstruct[ind].IVAR_OPT GT 0.0)/(sqrt(objstruct[ind].IVAR_OPT) $
                                           + (objstruct[ind].IVAR_OPT LE 0.0))
   wave = objstruct[ind].WAVE_OPT
   flux = objstruct[ind].FLUX_OPT
   splog,'Select second object'
   id2 = x_guilist(objid, indx = ind2)
   sig2 = (objstruct[ind2].IVAR_OPT GT 0.0)/(sqrt(objstruct[ind2].IVAR_OPT) $
                                           + (objstruct[ind2].IVAR_OPT LE 0.0))
   wave2 = objstruct[ind2].WAVE_OPT
   flux2 = objstruct[ind2].FLUX_OPT
   x_specplot, flux, sig, wave = wave, inflg = 4, xsize = 1000 $
               , ysize = 700,ytwo=flux2,two_wave=wave2
ENDELSE

IF KEYWORD_SET(FIND_FEATURE) THEN BEGIN
;default is to look for OII
   if not(keyword_set(FEATURE_WAVE)) then feature_wave=[3727.]
   if not(keyword_set(FEATURE_NAME)) then feature_name=['[OII]']
   if not(keyword_set(REDSHIFT)) then redshift=[0.]
   ;loop over each redshift
   for z=0, n_elements(redshift)-1 do begin
      ;loop over each feature
      for wave_i=0, n_elements(feature_wave)-1 do begin
         feature_wavei=feature_wave[wave_i]
         wave_feat = where(waveimg  gt $
                           ((1.+redshift[z])*feature_wavei - 5.0) and $
                           waveimg lt ((1.+redshift[z])*feature_wavei + 5.0))
         if wave_feat[0] ne -1 then begin
            xsizer = (size(waveimg))[1]
            wave_feat_x = wave_feat mod xsizer
            wave_feat_y = wave_feat / xsizer
            xatvplot, wave_feat_x, wave_feat_y, psym=3
            redshift_print_pos=strpos(strcompress(string(redshift[z]), /rem), '.')
            redshift_print=strmid(strcompress(string(redshift[z]), /rem), 0, redshift_print_pos+4)
            xatvxyouts, floor(min(left_edge[*,0]))-120, floor(max(wave_feat_y))+2, charsize=1.8, feature_name[wave_i]+', '+redshift_print
         endif
      endfor
   endfor
ENDIF


END
