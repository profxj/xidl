pro extrct_he1029, wfccd, maskid, exp_id, MKOBJS=mkobjs, SKYSUB=skysub, $
                   EXTRCT=extrct, FLUX=flux, COMBINE=combine, MKALL=mkall

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'extrct_he1029, wfccd, mask_id, [exp_id], /MKOBJS, /SKYSUB, /EXTRCT'
    print, '    /FLUX, /COMBINE, /MKALL'
    return
  endif 

  if not keyword_set( exp_id ) then exspr = 0L

  if keyword_set( MKALL ) then begin
      mkobjs = 1
      skysub = 1
      extrct = 1
      flux = 1
      combine = 1
  endif

  ;;;;;;;;;;;;;;
  ; Create Obj Structure
  if keyword_set( MKOBJS ) then begin
      wfccd_mkobjstr, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_17apr01.fits'
      wfccd_objsumm, wfccd, maskid, exp_id
  endif

  ;;;;;;;;;;;;;;
  ; Sky Subtract
  if keyword_set( SKYSUB ) then begin
      wfccd_skysub, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_17apr01.fits'
  endif

  ;;;;;;;;;;;;;;
  ; Extract
  if keyword_set( EXTRCT ) then begin
      wfccd_extobjbox, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_17apr01.fits'
      wfccd_objsumm, wfccd, maskid, exp_id
      ; Check the spectra
      wfccd_chkspec, wfccd, maskid, exp_id
  endif

  ;;;;;;;;;;;;;;
  ; Flux
  if keyword_set( FLUX ) then begin
      wfccd_fluxspec, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_17apr01.fits'
      wfccd_objsumm, wfccd, maskid, exp_id
  endif

  ;;;;;;;;;;;;;;
  ; COMBINE
  if keyword_set( COMBINE ) then begin
      wfccd_combspec, wfccd, maskid, /PATCH  ; Combine all of the frames
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_17apr01.fits'
  endif

  return
end
