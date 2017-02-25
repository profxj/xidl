
;; Read in the reference solution
;archive_file = getenv('LONGSLIT_DIR') + $
;               '/calib/linelists/lris_red_400_8500_d560.sav'
;iref = 1
;restore, archive_file
;calib_ref = calib[2]
;arc_ref = archive_arc[*, iref]
calib_temp = 0
;; read in the current arc array
restore, '/Users/joe/lris_mar_2014/redux/cpilot09/Red600/wave-r140326_0014.sav'
pkwdth = wstruct.pkwdth
toler = wstruct.toler
thin = wstruct.thin
psig = wstruct.psig[5]
mxoff = wstruct.mxoff[5]
fweight = wstruct.fweight

dims = size(arc1d, /dim)
ny = dims[0]
nyby2 = ny/2L
nslits = dims[1]
; read in the solution for the 600 line grism with the longslit
; this has the appropriate fit paramters, but not the right fit. 
archive = 1
;;restore, 'calib4.sav'
IF NOT KEYWORD_SET(ARCHIVE) THEN BEGIN
;restore, getenv('LONGSLIT_DIR') + '/calib/linelists/lris_red_600_7500.sav'
;; Add central wavelenght and central dispersion tags, which will 
;; be needed later
;; now run x_identify for each slit. It is redundant to calibrate
;; multiple slits at the same mask x-position, so to save time you might
;; want to choose them at different x values. 
   slit = 9
    x_identify, arc1d[*, slit-1L], calib_temp $
;;                , incalib = calib $
                , xsize = 1200, ysize = 600 $
                , linelist = getenv('XIDL_DIR')+ $
                '/Spec/Arcs/Lists/lris_red_600.lst' $
                , pkwdth = pkwdth $
                , TOLER = TOLER, pksig = psig, THIN = THIN, FWEIGHT = FWEIGHT $
                , mxoff = mxoff, maxqual = maxqual 
;                , mfitstr = calib_ref $
;                , mspec = arc_ref $
    IF TAG_EXIST(calib_temp, 'WAVE_CEN') EQ 0 THEN $
       calib_temp = struct_addtags(calib_temp $
                                   , create_struct('WAVE_CEN', 0.0D $
                                                   , 'DISP_CEN', 0.0D))
    wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
    wave_cen = wave_vec[nyby2]
    disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
    calib_temp.WAVE_CEN = wave_cen
    calib_temp.DISP_CEN = DISP_CEN
;; Add central wavelength and central dispersion to
    archive_arc = arc1d[*, slit-1L]
    calib = calib_temp
    save, archive_arc, calib $
          , file = 'calib' + strcompress(string(slit), /rem) + '.sav'
 ENDIF
;; After you have fit each arc, save them in a file named 
; getenv('LONGSLIT_DIR') + '/calib/linelists/lris_blue_600_5000.sav'
slits = [1, 4, 5, 6, 9]
nslits = n_elements(slits)
FOR j = 0L, nslits-1L DO BEGIN
;; where calib is an array of 25 structures
   restore, file = 'calib' + strcompress(string(slits[j]), /rem) + '.sav'
   IF j EQ 0 THEN BEGIN
      calib_out = replicate(calib, nslits)
      archive_arc_out = fltarr(n_elements(archive_arc), nslits) 
   ENDIF
   calib_out[j] = calib
   archive_arc_out[*, j] = archive_arc
ENDFOR
;archive_arc = arc1d[*, slits-1L]
;stop
savefile = getenv('LONGSLIT_DIR') + '/calib/linelists/lris_red_600_5000.sav'
archive_arc = archive_arc_out
calib = calib_out
save, archive_arc, calib, file = savefile
;restore, archive_file
;archive_arc = [[archive_arc], [arc1d[*, slit-1L]]]
;calib = [calib, calib_temp]
;spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560.sav  ' $
;       +  getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560_old.sav'  
;save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
;      '/calib/linelists/lris_red_400_8500_d560.sav'


restore, '/Users/joe/lris_mar_2014/redux/cpilot09/Red600/wave-r140326_0014.sav'
archive_arc = arc1d
calib = xfit
spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600_5000.sav   ' $
       +  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600_5000_old.sav'
save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600_5000.sav'



END
