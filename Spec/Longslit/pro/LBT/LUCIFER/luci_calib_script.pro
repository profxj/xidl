
;; Read in saved arc
restore, '/Users/joe/lucifer_data_2011/redux/luci_200HK_arc.sav'
dims=size(arc1d,/dim)
ny = dims[0]
nyby2 = ny/2L
linelist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/lowd_ir_ohlines.lst' 

;linelist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 

psig = 4.0d
slit_pix = 8.0d ;; 2.0" slit is 8 pixels
pkwdth = slit_pix
toler = slit_pix/2.0d

x_identify, arc1d, calib_temp $
            , xsize = 1200, ysize = 600 $
            , linelist = linelist $
            , PKWDTH = PKWDTH, TOLER = TOLER,  PKSIG = PSIG $
            , /THIN, /FWEIGHT
calib_temp = struct_addtags(calib_temp $
                            , create_struct('WAVE_CEN', 0.0D $
                                            , 'DISP_CEN', 0.0D))
wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
wave_cen = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
calib_temp.WAVE_CEN = wave_cen
calib_temp.DISP_CEN = DISP_CEN
;; Add central wavelength and central dispersion to
save, calib_temp, file = 'calib_200H+K.sav'

END
