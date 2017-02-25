plate_scale = 0.273d ;; Is this correct?????
slit_width = 0.60d
slit = slit_width/plate_scale ;; slit width in pixels
pkwdth = slit
TOLER = 2.0 ;;slit/2.0D ;; tolerance for the line centers, use this value???
FWHM = slit 
psig = 4.0d ;; Sigma threshold to call a peak a line

FORDR = 9L ;; Order of the fit to the arc continuum which is removed
THIN = 1   ;; Use thin lines
fweight = 1 ;; Use flux weighted centroiding 
restore, '/Users/joe/DATA/SOFI_DATA/2011-09-22/sofi_arcH.sav'

dims=size(arc1d,/dim)
ny = dims[0]
nyby2 = ny/2L

linelist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 
x_identify, arc1d, calib_temp, xsize = 1500, ysize = 900 $
            , linelist = linelist, pkwdth = pkwdth $
            , TOLER = TOLER, pksig = psig, /THIN, /FWEIGHT
;;, FORDR = FORDR


calib_temp = struct_addtags(calib_temp $
                            , create_struct('WAVE_CEN', 0.0D $
                                            , 'DISP_CEN', 0.0D))
wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
wave_cen = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
calib_temp.WAVE_CEN = wave_cen
calib_temp.DISP_CEN = DISP_CEN

END
