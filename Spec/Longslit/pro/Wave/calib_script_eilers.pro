
calib_temp = 0
;; read in the output of the long_reduce with wave_calib=1
restore, '//Users/joe/2016sep26/redux/wave-LR.20160926.11323.sav'
;; This is setting the line identification params to be the same as
;; the code defaults.  
pkwdth = wstruct.pkwdth
toler = wstruct.toler
thin = wstruct.thin
fweight = wstruct.fweight
psig = wstruct.psig[5]

dims = size(arc1d, /dim)
ny = dims[0]
nyby2 = ny/2L
nslits = dims[1]
; read in the solution for the 1200 line grism with the longslit
; this has the appropriate fit paramters, but not the right fit. 
slit = 2 ;; This is the slit that is busted
x_identify, arc1d[*, slit-1L], calib_temp $
            , xsize = 1200, ysize = 600 $
            , linelist = getenv('XIDL_DIR') + $
            '/Spec/Arcs/Lists/lris_red_600.lst' $
            , pkwdth = pkwdth $
            , TOLER = TOLER, pksig = nsig, THIN = THIN, FWEIGHT = FWEIGHT
            ;;, mfitstr = calib_ref $
            ;;, mspec = arc_ref $
;; Add central wavelength and central dispersion to
IF TAG_EXIST(calib_temp, 'WAVE_CEN') EQ 0 THEN $
   calib_temp = struct_addtags(calib_temp $
                               , create_struct('WAVE_CEN', 0.0D $
                                               , 'DISP_CEN', 0.0D))
wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
wave_cen = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
calib_temp.WAVE_CEN = wave_cen
calib_temp.DISP_CEN = DISP_CEN

stop ;; Hit go if everything ran okay, because at this step you are overwriting
archive_file =  getenv('LONGSLIT_DIR') +  '/calib/linelists/lris_red_600_7500.sav'
restore, archive_file 
archive_arc = [[archive_arc], [arc1d[*, slit-1L]]]
calib = [calib, calib_temp]
spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600_7500.sav ' $
       +  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_600_7500_old.sav'
save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
      '/calib/linelists/lris_red_600_7500.sav'


;restore, '/Users/joe/lris_data_2013/redux/UDS_0217/Red400/wave-r131230_0039.sav'
;archive_arc = arc1d
;calib = xfit
;spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560.sav  ' $
;       +  getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/lris_red_400_8500_d560_old.sav'  
;save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
;      '/calib/linelists/lris_red_400_8500_d560.sav'
      


END
