pro extrct_obj, wfccd, maskid, exp_id, obj_nm, MKOBJS=mkobjs, SKYSUB=skysub, $
                   EXTRCT=extrct, FLUX=flux, COMBINE=combine, MKALL=mkall, $
                    ZFIND=zfind, IMG=img, CHK=chk, NOSKY=nosky, PLOT=plot

; extrct_obj, wfccd, 0L, [0L,1L], '2154a', /MKALL

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'extrct_obj, wfccd, mask_id, exp_id, obj_nm, /MKOBJS, /SKYSUB, /EXTRCT'
    print, '    /FLUX, /COMBINE, /MKALL, /NOSKY'
    return
  endif 

  if keyword_set(MKALL) then begin
      skysub = 1
      extrct = 1
      flux = 1
      combine = 1
      zfind = 1
      plot = 1
  endif

  ; Set exp
  exp = [exp_id]
  nexp = n_elements(exp)

  len = strlen(obj_nm)
  slit = long(strmid(obj_nm,0,len-1))

  ;;;;;;;;;;;;;;
  ; Sky Subtract
  if keyword_set( SKYSUB ) then $
    for q=0L,nexp-1 do wfccd_skysub, wfccd, maskid, exp[q], SLIT_ID=slit

  ;;;;;;;;;;;;;;
  ; Extract
  if keyword_set( EXTRCT ) then begin
    for q=0L,nexp-1 do begin
        ; Extract
        wfccd_extobjbox, wfccd, maskid, exp[q], SLIT_ID=slit
        wfccd_objsumm, wfccd, maskid, exp[q]
        ; Check the spectra
        if keyword_set( CHK) then wfccd_chkspec, wfccd, maskid, exp[q]
    endfor
  endif

  ;;;;;;;;;;;;;;
  ; Flux
  if keyword_set( FLUX ) then begin  
      for q=0L,nexp-1 do begin
          wfccd_fluxspec, wfccd, maskid, exp[q], /FORCE
          wfccd_objsumm, wfccd, maskid, exp[q]
      endfor
  endif

  ;;;;;;;;;;;;;;
  ; COMBINE
  if keyword_set( COMBINE ) then $
    wfccd_combspec, wfccd, maskid, OBJ_NM=obj_nm, /PATCH  

  ;;;;;;;;;;;;;;
  ; ZFIND
  if keyword_set( ZFIND ) then begin
      if maskid LT 10 then $
        fspec_fil = 'Extract/Fspec_0'+strtrim(maskid,2)+'.fits' $
      else $
        fspec_fil = 'Extract/Fspec_'+strtrim(maskid,2)+'.fits'
      wfccd_zfind, fspec_fil, obj_nm, PLOT=plot
      if keyword_set(CHK) then wfccd_chkfspec, fspec_fil, /zfind
      wfccd_fspecsumm, fspec_fil
  endif

  return
end
