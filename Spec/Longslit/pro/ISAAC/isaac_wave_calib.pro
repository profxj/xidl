
  
;; wcen = 1.62217, was this a telluric??
;savefile = '/Users/joe/ISAAC_redux/nite1/1007+0337/Science/Hip041640_0/tel-ISAAC_SWS_STD_0001-ISAAC_SWS_STD_0002.sav'

;; wcen = 1.59674
;savefile = '/Users/joe/ISAAC_redux/nite1/1235-0108/Science/Hip091038_1/wave-ISAAC_SWS_STD_0015-ISAAC_SWS_STD_0016.sav'

;; 1.65177
;savefile = '/Users/joe/ISAAC_redux/nite1/1442-0242/Science/Hip105164_4/wave-ISAAC_SWS_STD_0023-ISAAC_SWS_STD_0024.sav'

;;1.69492
;savefile = '/Users/joe/ISAAC_redux/nite1/1610+0442/Science/Hip082658_5/wave-ISAAC_SWS_STD_0028-ISAAC_SWS_STD_0027.sav'

;; 2.05581
;savefile = '/Users/joe/ISAAC_redux/nite1/PKS1935-962/Science/Hip089022_7/wave-ISAAC_SWS_STD_0031-ISAAC_SWS_STD_0032.sav'

;; 1.53350
;savefile = '/Users/joe/ISAAC_redux/nite2/1348+0531/Science/Hip078530_0/wave-ISAAC_SWS_STD_0014-ISAAC_SWS_STD_0015.sav'

;;1.56741
;savefile = '/Users/joe/ISAAC_redux/nite2/2300+0031/Science/Hip080804_1/wave-ISA;AC_SWS_STD_0032-ISAAC_SWS_STD_0033.sav'

;1.72847
;savefile = '/Users/joe/ISAAC_redux/nite3/1357+0117/Science/Hip076069_3/wave-ISAAC_SWS_STD_0009-ISAAC_SWS_STD_0010.sav'

; wcen = 2.370
savefile = '/Users/joe/isaac_redux/1630+0435/Science/Hip090271_3/wave-ISAAC.2013-04-12T08_53_18.660-ISAAC.2013-04-12T08_54_35.245.sav'





restore, savefile 

band = wstruct.band
mode = wstruct.mode
plate_scale = 0.147d
pkwdth = wstruct.pkwdth
toler = wstruct.toler
thin = wstruct.thin
fweight = wstruct.fweight
psig = wstruct.psig
wave_cen = wstruct.wave_cen
dims = size(arc1d, /dim)
ny = dims[0]
nyby2 = ny/2L

linelist = wstruct.linelist
x_arclist, linelist, lines
modelsky = repstr(linelist, 'OH_linelist_', '')
modelsky = repstr(modelsky, '.lst', '.fits')

;; Read in the sky model as help
CASE band OF 
   'H': delta_lam = 0.079d
   'K': delta_lam = 0.122d
   else: message, 'Unrecognized band'
ENDCASE
dlam = delta_lam/1024.0d
wave_min = wave_cen - delta_lam/2.0d
wave_max = wave_cen + delta_lam/2.0d
model = x_readspec(modelsky, wav = wave_model, inflg = 0)
iwave = where(wave_model GE wave_min AND wave_model LE wave_max)
yrange = [0.9*min(model[iwave]), 1.1*max(model[iwave])]
x_specplot, model, wav = wave_model $
            , xr = [wave_min, wave_max], yrang = yrange $
            , xsize = 600, ysize = 300 $
            , ytwo = yrange[1]*0.7 + 0*lines.wave $
            , two_wave = lines.wave, psym2 = 2

x_identify, arc1d[*, 0], calib_temp, xsize = 1200, ysize = 600 $
            , linelist = linelist, pkwdth = pkwdth $
            , TOLER = TOLER, pksig = nsig, THIN = THIN, FWEIGHT = FWEIGHT
calib_temp = struct_addtags(calib_temp $
                            , create_struct('WAVE_CEN', 0.0D $
                                            , 'DISP_CEN', 0.0D))
wave_vec = x_calcfit(dindgen(ny), fitstr = calib_temp)
wave_cen_fit = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
calib_temp.WAVE_CEN = wave_cen_fit
calib_temp.DISP_CEN = DISP_CEN

stop
restore, getenv('LONGSLIT_DIR') + $
         '/calib/linelists/ISAAC_' + mode + '_' + band + '.sav' 
archive_arc = [[archive_arc], [arc1d[*, 0]]]
calib = [calib, calib_temp]
spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
         '/calib/linelists/ISAAC_' + mode + '_' + band + '.sav   ' $
       + getenv('LONGSLIT_DIR') + $
       '/calib/linelists/ISAAC_' + mode + '_' + band + '_old.sav' 
save, archive_arc, calib $
      , file = getenv('LONGSLIT_DIR') + $
      '/calib/linelists/ISAAC_' + mode + '_' + band + '.sav' 
;; 



;archive_arc = arc1d[*, 0]
;calib = calib_temp
;calibpath = '/Users/joe/IDL/xidl/Spec/Longslit/calib/linelists/'
;save, archive_arc, calib, file = calib_path + 'ISAAC_calib_K_2.060.sav'

END
