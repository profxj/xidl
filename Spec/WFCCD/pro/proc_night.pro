pro proc_night, wfccd, maskid, expsr, MKFLAT=mkflat, MKBIAS=mkbias, $
                 MKMAP=mkmap, MKSTRCT=mkstrct, MKALL=mkall, $
                 MKSLITS=mkslits, MKAIMG=mkaimg, MKNFLAT=mknflat, $
                 MKOBJ=mkobj, CLOBBER=clobber,CLOBAIMG=clobaimg, FLGCR=flgcr, $
                  NOCHK=nochk
  ; MKALL
  if keyword_set( MKALL ) then begin
      MKFLAT = 1
      MKMAP = 1
      MKSLITS = 1
      MKAIMG = 1
      MKNFLAT = 1
      MKOBJ = 1
      FLGCR = 1
      CLOBBER = 1
      clobaimg = 1
  endif

  ;;;;;;;;;;;;;;
  ; Make the structure
  if keyword_set( MKSTRCT ) then begin
      wfccd_strct, wfccd, /mkdir
      ; Adjust the structure and write fits file
      setcrds_night, wfccd, FITS='wfccd_night.fits'
  endif
                 
  ;;;;;;;;;;;;;;
  ; Make the Bias
  if keyword_set( MKBIAS ) then wfccd_bias, wfccd

;;;;;;; MASK DEPENDENT ;;;;;;;;
  ;;;;;;;;;;;;;;
  ; Make the Flat
  if keyword_set( MKFLAT ) then begin
      wfccd_mkflat, wfccd, maskid
      write_wfccdstr, wfccd, FITS='wfccd_night.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Make the Map
  if keyword_set( MKMAP ) then begin
      wfccd_mkmap, wfccd, maskid
      write_wfccdstr, wfccd, FITS='wfccd_night.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;  
  ; Make the Slit Structure
  if keyword_set( MKSLITS ) then begin
      wfccd_setslits, wfccd, maskid
      write_wfccdstr, wfccd, FITS='wfccd_night.fits', /ANONLY
  endif

  ;;;;;;;;;;;
  ; Normalize the Flat
  if keyword_set( MKNFLAT ) then begin
      wfccd_normflat, wfccd, maskid
      write_wfccdstr, wfccd, FITS='wfccd_night.fits', /ANONLY
  endif

;;;;;   MASKID + EXPOSURE ;;;;;;;
  ;;;;;;;;;;;;;;
  ; proc_night, wfccd, 0L, [0L,1L], /MKAIMG, /CLOBBER
  ; Make the Arc Image(s)
  if keyword_set( MKAIMG ) then begin
      ; Process the Arcs
      wfccd_procarc, wfccd, maskid, CLOBBER=clobber
      ; Arc solution
      wfccd_allarcsol, wfccd, maskid, CLOBBER=clobber
      ; Create Arc Image
      if not keyword_set( expsr ) then begin
          expsr = where(wfccd.mask_id EQ maskid AND $
                        wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0, nexp) 
          for q=0L,nexp-1 do begin
              wfccd_arcimg, wfccd, maskid, q, CLOBBER=clobaimg
              if not keyword_set(NOCHK) then wfccd_chkarc, wfccd, maskid, q
              wfccd_cleanarc, wfccd, maskid, q
          endfor
      endif else begin
          nexp = n_elements(expsr)
          for q=0L,nexp-1 do begin
              wfccd_arcimg, wfccd, maskid, expsr[q], CLOBBER=clobaimg
              if not keyword_set(NOCHK) then wfccd_chkarc, wfccd, maskid, $
                expsr[q] 
              wfccd_cleanarc, wfccd, maskid, expsr[q] 
          endfor
      endelse
      ; All done
      print, 'proc_night: All done with MKAIMG!'
      ; Write
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;;;
  ; Make the Object Images
  if keyword_set( MKOBJ ) then begin
      ; Process the Arcs
      wfccd_procobj, wfccd, maskid, CLOBBER=clobber
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;
  ;  CRUDE CR
  if keyword_set( FLGCR ) then begin
      ; Identify CR in each image and set VAR=-1
      wfccd_objcr, wfccd, maskid
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  ;;;;;;;;;;;;;;
  ; Cleanup
  if keyword_set( CLEANUP ) then begin
      ; Delete unnecessary files
      wfccd_cleanup, wfccd, maskid, 3
      write_wfccdstr, wfccd, FITS='wfccd_night.fits'
  endif

  return
end
