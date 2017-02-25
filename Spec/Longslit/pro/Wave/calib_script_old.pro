
;; Read in the reference solution
archive_file = getenv('LONGSLIT_DIR') + $
               '/calib/linelists/lris_red_400_8500_d560.sav'
iref = 1
restore, archive_file
calib_ref = calib[2]
arc_ref = archive_arc[*, iref]
calib_temp = 0
;; read in the current arc array
restore, '/Users/joe/lris_data_2013/redux/UDS_0217/Red400/wave-r131230_0039.sav'
pkwdth = wstruct.pkwdth
toler = wstruct.toler
thin = wstruct.thin
fweight = wstruct.fweight
psig = wstruct.psig[5]

dims = size(arc1d, /dim)
ny = dims[0]
nyby2 = ny/2L
nslits = dims[1]
; read in the solution for the 600 line grism with the longslit
; this has the appropriate fit paramters, but not the right fit. 
archive = 0
IF NOT KEYWORD_SET(ARCHIVE) THEN BEGIN
;restore, getenv('LONGSLIT_DIR') + '/calib/linelists/lris_red_600_7500.sav'
;; Add central wavelenght and central dispersion tags, which will 
;; be needed later
;; now run x_identify for each slit. It is redundant to calibrate
;; multiple slits at the same mask x-position, so to save time you might
;; want to choose them at different x values. 
   slit = 20
    x_identify, arc1d[*, slit-1L], calib_temp $
                , mfitstr = calib_ref $
                , mspec = arc_ref $
                , xsize = 1200, ysize = 600 $
                , linelist = getenv('XIDL_DIR')+ $
                '/Spec/Arcs/Lists/lris_red_300.lst' $
                , pkwdth = pkwdth $
                , TOLER = TOLER, pksig = nsig, THIN = THIN, FWEIGHT = FWEIGHT
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
    save, calib_temp, file = 'calib' + strcompress(string(slit), /rem) + '.sav'
 ENDIF
;; After you have fit each arc, save them in a file named 
;; getenv('LONGSLIT_DIR') + '/calib/linelists/lris_blue_600.sav'
;nslits = 4
;slits = [1, 14, 15, 24]
;FOR j = 0L, nslits-1L DO BEGIN
;;; where calib is an array of 25 structures
;   restore, file = 'calib' + strcompress(string(slits[j]), /rem) + '.sav'
;   IF j EQ 0 THEN calib = replicate(calib_temp, nslits)
;   calib[j] = calib_temp
;ENDFOR
;archive_arc = arc1d[*, slits-1L]
;stop
stop
;restore, archive_file
;archive_arc = [[archive_arc], [arc1d[*, slit-1L]]]
;calib = [calib, calib_temp]
;spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560.sav  ' $
;       +  getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560_old.sav'  
;save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
;      '/calib/linelists/lris_red_400_8500_d560.sav'


restore, '/Users/joe/lris_data_2013/redux/UDS_0217/Red400/wave-r131230_0039.sav'
archive_arc = arc1d
calib = xfit
spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_400_8500_d560.sav  ' $
       +  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_400_8500_d560_old.sav'  
save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
      '/calib/linelists/lris_red_400_8500_d560.sav'
      


END
