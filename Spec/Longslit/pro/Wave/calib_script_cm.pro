
;; Read in the reference solution that didn't run through
calib_file = getenv('LONGSLIT_DIR') + '/calib/linelists/lris_red_1200_clm.sav'
iref = 0
restore, calib_file
calib_ref = calib[iref-1L]
arc_ref = archive_arc[*, iref-1L]
;; read in the output of the long_reduce with wave_calib=1
;restore, '//Users/cleibler/Documents/DataReduction/lris/2015apr/mask_J1211_c4/Blue/wave-b150420_4026.sav'
;; This is setting the line identification params to be the same as
;; the code defaults.
wave_savefile='/Users/joe/cmartin_LRIS/red/E570/wave-r170816_0009.sav'
slit = 16
restore,wave_savefile
calib_temp = 0

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
;;slit = 12 ;; This is the slit that is busted
x_identify, arc1d[*, slit-1L], calib_temp $
            , mfitstr = calib_ref $
            , mspec = arc_ref $
            , xsize = 1200, ysize = 600 $
            , linelist = getenv('XIDL_DIR')+ $
            '/Spec/Arcs/Lists/lris_red_600.lst' $
            , pkwdth = pkwdth $
            , TOLER = TOLER, pksig = nsig, THIN = THIN, FWEIGHT = FWEIGHT
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
;restore, archive_file
archive_arc = [[archive_arc], [arc1d[*, slit-1L]]]
calib = [calib, calib_temp]
;archive_arc = arc1d[*,slit-1L]
;calib = calib_temp
spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_1200_clm.sav  ' $
       +  getenv('LONGSLIT_DIR') + $
       '/calib/linelists/lris_red_1200_clm_old.sav'  
save, archive_arc, calib, file =  getenv('LONGSLIT_DIR') + $
      '/calib/linelists/lris_red_1200_clm.sav'


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
