;; 
pro esi_echtstlai, OUTDIR, SKIP_SNGL=skip_sngl, SKIP_MULTI=skip_multi, SKIP_ARC=skip_arc, $
             SKIP_FLAT=skip_flat, SKIP_STD=skip_std, SKIP_NOFLUX=skip_noflux

 if  N_params() LT 1  then return
 ;; Create the lists
 obj_fil = ['Raw/esi0027.fits', 'Raw/esi0028.fits', 'Raw/esi0029.fits'] 
 arc_fil = ['Raw/esi0033.fits', 'Raw/esi0048.fits']
 flt_fil = ['Raw/esi0038.fits', 'Raw/esi0039.fits','Raw/esi0040.fits']
 std_fil = ['Raw/esi0030.fits']


 rmstuff = '\rm -fr Arcs Bias Extract FSpec Final Flats Logs Maps OV Sky s2n*'
 spawn, rmstuff

 ;; Test 1 -- Single object and no Calibs
 if not keyword_set(SKIP_SNGL) then begin
     file_mkdir, OUTDIR+'/SNGL_OBJ'
     esi_echpipe, obj_fil[0]
     
     spawn, 'cp FSpec/* '+OUTDIR+'/SNGL_OBJ/'
     spawn, 'cp s2n_* '+OUTDIR+'/SNGL_OBJ/'
     
     spawn, rmstuff
 endif

 ;; Test 2 -- Several objects and no Calibs
 if not keyword_set(SKIP_MULTI) then begin
     file_mkdir, OUTDIR+'/MULTI_OBJ'
     esi_echpipe, obj_fil
     spawn, 'cp FSpec/* '+OUTDIR+'/MULTI_OBJ/'
     spawn, 'cp s2n_* '+OUTDIR+'/MULTI_OBJ/'
     
     spawn, rmstuff
 endif

 ;; Test 3 -- Single object with Arcs
 if not keyword_set(SKIP_ARC) then begin
     file_mkdir, OUTDIR+'/ARC_ONLY'
     esi_echpipe, obj_fil[0], arc_fil
     spawn, 'cp FSpec/* '+OUTDIR+'/ARC_ONLY/'
     spawn, 'cp s2n_* '+OUTDIR+'/ARC_ONLY/'
     
     spawn, rmstuff
 endif

 ;; Test 4 -- Single object with Flats
 if not keyword_set(SKIP_FLAT) then begin
     file_mkdir, OUTDIR+'/FLAT_ONLY'
     esi_echpipe, obj_fil[0], 0, flt_fil
     spawn, 'cp FSpec/* '+OUTDIR+'/FLAT_ONLY/'
     spawn, 'cp s2n_* '+OUTDIR+'/FLAT_ONLY/'
     
     spawn, rmstuff
 endif

 ;; Test 5 -- Single object with Standard
 if not keyword_set(SKIP_STD) then begin
     file_mkdir, OUTDIR+'/STD'
     esi_echpipe, obj_fil, STANDARD=std_fil
     spawn, 'cp FSpec/* '+OUTDIR+'/STD/'
     spawn, 'cp s2n_* '+OUTDIR+'/STD/'
     
     spawn, rmstuff
 endif

 ;; Test 6 -- Full deal minus fluxing standard
 if not keyword_set(SKIP_NOFLUX) then begin
     file_mkdir, OUTDIR+'/NOFLUX'
     esi_echpipe, obj_fil, arc_fil, flt_fil, STANDARD=std_fil
     spawn, 'cp FSpec/* '+OUTDIR+'/NOFLUX/'
     spawn, 'cp s2n_* '+OUTDIR+'/NOFLUX/'
     
     spawn, rmstuff
 endif

 return
end


 
