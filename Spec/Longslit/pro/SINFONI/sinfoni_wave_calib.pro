

savefile = '/Users/joe/SINFONI_082015/redux/wave-SINFONI_IFS_OBS215_0006.sav'
restore, savefile 

pkwdth = wstruct.pkwdth
toler = wstruct.toler
thin = wstruct.thin
fweight = wstruct.fweight
psig = wstruct.psig
;wave_cen = wstruct.wave_cen
dims = size(arc1d, /dim)
ny = dims[0]
nyby2 = ny/2L

linelist = wstruct.linelist
x_arclist, linelist, lines
modelsky = repstr(linelist, '.lst', '.fits')
wave_min = 1.95
wave_max = 2.45
delta_lam = wave_max - wave_min
dlam = delta_lam/double(ny)


model = x_readspec(modelsky, wav = wave_model, inflg = 0)
iwave = where(wave_model GE wave_min AND wave_model LE wave_max)
yrange = [0.9*min(model[iwave]), 1.1*max(model[iwave])]
x_specplot, model, wav = wave_model $
            , xr = [wave_min, wave_max], yrang = yrange $
            , xsize = 900, ysize = 450 $
            , ytwo = yrange[1]*0.05 + 0*lines.wave $
            , two_wave = lines.wave, psym2 = 2
stop

x_identify, arc1d[*, 5], calib_temp, xsize = 1200, ysize = 600 $
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
;; This is for multislit
;restore, getenv('LONGSLIT_DIR') + $
;         '/calib/linelists/LUCI_200H+K.sav'
;archive_arc = [[archive_arc], [arc1d[*, 0]]]
;calib = [calib, calib_temp]
;spawn, 'mv ' + getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/LUCI_200H+K.sav  ' $
;       + getenv('LONGSLIT_DIR') + $
;       '/calib/linelists/LUCI_200H+K_old.sav' 
;save, archive_arc, calib $
;      , file = getenv('LONGSLIT_DIR') + $
;      '/calib/linelists/LUCI_200H+K.sav'
;; 


archive_arc = arc1d[*, 5]
calib = calib_temp
calibpath = '/Users/joe/IDL/xidl/Spec/Longslit/calib/linelists/'
save, archive_arc, calib, file = calibpath + 'SINFONI_0.25_K.sav'

END
