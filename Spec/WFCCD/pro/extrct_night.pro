pro extrct_night, wfccd, maskid, exp_id, MKOBJS=mkobjs, SKYSUB=skysub, $
                   EXTRCT=extrct, FLUX=flux, COMBINE=combine, MKALL=mkall, $
                    ZFIND=zfind, IMG=img, NOCHK=nochk, NOSKY=nosky

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'extrct_night, wfccd, mask_id, [exp_id], /MKOBJS, /SKYSUB, /EXTRCT'
    print, '    /FLUX, /COMBINE, /MKALL, /NOSKY'
    return
  endif 

  ; Set fspec_fil
  if maskid LT 10 then $
    fspec_fil = 'Extract/Fspec_0'+strtrim(maskid,2)+'.fits' $
  else $
    fspec_fil = 'Extract/Fspec_'+strtrim(maskid,2)+'.fits'

  if keyword_set( MKALL ) then begin
      mkobjs = 1
      skysub = 1
      extrct = 1
      flux = 1
      clobber=1
  endif


  ;;;;;;;;;;;;;;
  ; Create Obj Structure
  if keyword_set( MKOBJS ) then begin
      wfccd_mkobjstr, wfccd, maskid, exp_id, NOSKY=nosky
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
      wfccd_objsumm, wfccd, maskid, exp_id
  endif

  ;;;;;;;;;;;;;;
  ; Sky Subtract
  if keyword_set( SKYSUB ) then begin
      wfccd_skysub, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;;;
  ; Extract
  if keyword_set( EXTRCT ) then begin
      wfccd_extobjbox, wfccd, maskid, exp_id
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
      wfccd_objsumm, wfccd, maskid, exp_id
      ; Check the spectra
      if not keyword_set( NOCHK) then wfccd_chkspec, wfccd, maskid, exp_id
  endif

  ;;;;;
  ;  EDIT
  if keyword_set( EDIT ) then begin
      wfccd_editspec, wfccd, maskid, exp_id
      wfccd_objsumm, wfccd, maskid, exp_id
  endif

  ;;;;;;;;;;;;;;
  ; Flux
  if keyword_set( FLUX ) then begin
      ; Set Exposures
      if not keyword_set( EXP_ID ) then begin
          expsr = where(wfccd.mask_id EQ maskid AND $
                        wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0, nexp) 
          exp = lindgen(nexp)
      endif else begin
          exp = [exp_id]
          nexp = n_elements(exp)
      endelse
      ; Loop
      for q=0L,nexp-1 do begin
          print, 'extrct_night: Fluxing exposure ', exp[q]
          wfccd_fluxspec, wfccd, maskid, exp[q], /FORCE
          wfccd_objsumm, wfccd, maskid, exp[q]
      endfor
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;;;
  ; COMBINE
  if keyword_set( COMBINE ) then begin
      wfccd_combspec, wfccd, maskid, /PATCH  ; Combine all of the frames
      wfccd_fspecsumm, fspec_fil
      ; Update
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;;;
  ; ZFIND
  if keyword_set( ZFIND ) then begin
      wfccd_zfind, fspec_fil
      if not keyword_set(NOCHK) then wfccd_chkfspec, fspec_fil, /zfind
      wfccd_fspecsumm, fspec_fil
  endif

  ;;;;;;;;;;;;;;
  ; IMG STUFF
  if keyword_set( IMG ) then begin
      case maskid of
          0L: begin
              wfccd_fspecxypix, fspec_fil, $
                '/home/xavier/LCO/data/OVI/FUSE/HE1029-140/HE1029-140_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/HE1029-140/HE1029-140XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end 
          1L : begin
              wfccd_fspecxypix, fspec_fil, $
                '/home/xavier/LCO/data/OVI/FUSE/PG1116+215/PG1116+215_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/PG1116+215/PG1116+215XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          2L : begin
              wfccd_fspecxypix, 'Extract/Fspec_02.fits', $
                '/home/xavier/LCO/data/OVI/FUSE/PG1116+215/PG1116+215_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/PG1116+215/PG1116+215XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          3L : begin
              wfccd_fspecxypix, fspec_fil, $
                '/home/xavier/LCO/data/OVI/FUSE/3C273/3C273_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/3C273/3C273XR/fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          4L : begin
              wfccd_fspecxypix, fspec_fil, $
            '/home/xavier/LCO/data/OVI/FUSE/PKS1302-102/PKS1302-102_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/PKS1302-102/PKS1302-102XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          5L : begin
              wfccd_fspecxypix, fspec_fil, $
            '/home/xavier/LCO/data/OVI/FUSE/MRK1383/MRK1383_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/MRK1383/MRK1383XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          6L : begin
              wfccd_fspecxypix, fspec_fil, $
            '/home/xavier/LCO/data/OVI/FUSE/Q1553+113/Q1553+113_BRphot.dat'
              wfccd_wrfspec, tmpfspec, fspec_fil, /read
              tmpfspec.img_fil = $
                '/home/xavier/LCO/data/OVI/FUSE/Q1553+113/Q1553+113XR.fits'
              wfccd_wrfspec, tmpfspec, fspec_fil
          end
          else :
      endcase
      ; Update Summary
      wfccd_fspecsumm, fspec_fil
  endif

  return
end
