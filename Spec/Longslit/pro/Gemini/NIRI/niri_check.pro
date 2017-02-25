;PRO NIRI_CHECK

;z = 3.251
;FEATURE_WAVE = [3727.26, 4341.69, 4862.70, 4960.29, 5008.24]
;FEATURE_NAME = ['~[OII]3727', 'Hgamma4341', 'Hbeta4862', '[OIII]4960' $
;                , '[OIII]5008']

;;path = './'
path = 'Science/*/'
;;path = '//Users/joe/NTT_data/redux/'
;;path = './'
fils = findfile(path + 'sci*.fits*', count = nfil)
tels = findfile(path + 'tel*.fits*', count = ntel)
IF nfil GT 0 AND ntel GT 0 THEN  BEGIN
   fils = [fils, tels]
   ind_str = [lindgen(nfil), lindgen(ntel)] 
ENDIF ELSE IF nfil GT 0 AND ntel EQ 0 THEN BEGIN
   fils = fils
   ind_str = lindgen(nfil)
ENDIF ELSE IF nfil EQ 0 And ntel GT 0 THEN BEGIN
   fils = tels
   ind_str = lindgen(ntel)
ENDIF
if n_elements(fils) EQ 0 then message, 'Cannot find files'
gui_in = string(ind_str, format = '(I2.2)')
gui_out = x_guilist(gui_in + '-' + repstr(fils, 'Science/', ''), indx = indx)
file_fits = fils[indx]
psfile = repstr(file_fits, 'fits','ps')

;spawn, 'gv  ' + psfile
hdr = xheadfits(file_fits)
instrument = strtrim(sxpar(hdr, 'INSTRUME'))
IF strmatch(instrument, 'SINFONI') THEN BEGIN
   diff        = mrdfits(file_fits, 0, scihdr)
   ivar        = mrdfits(file_fits, 1)
   waveimg     = mrdfits(file_fits, 2)
   img_profile = mrdfits(file_fits, 3)
   objstruct   = mrdfits(file_fits, 4)
ENDIF ELSE BEGIN
   sciimg = mrdfits(file_fits, 0, scihdr)
   sky_model = mrdfits(file_fits, 1)
   ivar      = mrdfits(file_fits, 2)
   waveimg   = mrdfits(file_fits, 3)
   objstruct = mrdfits(file_fits, 4)
   diff = sciimg - sky_model
ENDELSE
xatv, diff*sqrt(ivar)*(waveimg GT 0.0) $
      , wvimg = waveimg, min = -7.0, max = 7.0
xatvplot, objstruct.xpos, objstruct.ypos, psym = 3

nobj = n_elements(objstruct)
FOR j = 0L, nobj-1L DO BEGIN
    xpix = djs_median(objstruct[j].XPOS) - 25.0
    ypix = 512
    xatvxyouts, xpix, ypix, strcompress(string(j+1), /rem) $
               , color = 'red', charsize = 3
 ENDFOR

objid = strtrim(objstruct.objid, 2)
splog, 'Select object'
id = x_guilist(objid, indx = ind)
ivar  = objstruct[ind].IVAR_OPT
sig = (ivar GT 0.0)/sqrt(ivar + (ivar LE 0.0))
wave = objstruct[ind].WAVE_OPT
flux = objstruct[ind].FLUX_OPT

instrument = strcompress(sxpar(scihdr, 'INSTRUME'), /rem)
IF instrument EQ 'ISAAC' THEN wave = wave*1d4
x_specplot, flux, sig, wave = wave, inflg = 4, xsize = 1000, ysize = 700 

IF KEYWORD_SET(Z) THEN BEGIN
      ;loop over each feature
   FOR wave_i = 0, n_elements(feature_wave)-1 do begin
      feature_wavei = feature_wave[wave_i]
      wave_feat = where(waveimg  gt ((1.0D + z)*feature_wavei - 3.0) and $
                        waveimg  lt ((1.0D + z)*feature_wavei + 3.0))
      if wave_feat[0] ne -1 then begin
         xsizer = (size(waveimg))[1]
         wave_feat_x = wave_feat mod xsizer
         wave_feat_y = wave_feat / xsizer
         xatvplot, wave_feat_x, wave_feat_y, psym = 3
         redshift_print_pos = strpos(strcompress(string(z), /rem), '.')
         redshift_print = strmid(strcompress(string(z), /rem) $
                                 , 0, redshift_print_pos+4)
         xatvxyouts, 20, floor(max(wave_feat_y))+2, charsize = 1.8 $
                     , feature_name[wave_i]+', '+redshift_print
      endif
   ENDFOR
ENDIF



END
