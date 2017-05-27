;+
; NAME:
; get_slitpos
;
; PURPOSE:
; finds position of slit on a mask
;
; CALLING SEQUENCE:
; get_slitpos, datadir, slitn, slitx, slity
; 
; INPUTS:
; 
; datadir - the data directory containing the spectral data for the
;           mask. Empty string (or './') for the current directory.
; slitn   - slit number whose position we want to find.
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; slitx - x coordinate of slit.
; slity - what you think.
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;
; Based on code by MCC, written Sep02.
; hacked by BFG Sep02.
;
;----------------------------------------------------------------------


pro get_slitpos,  datadir, slitn, slitx,  slity

;;; datadir: a string containing the data directory of interest
;;; (including the terminal slash).

;;; FIND THE FILE xxxx.bintabs.fits. EXTRACT THE OBJECT CATALOG
;;; EXTENSION. 
  binfile = FINDFILE(datadir+'*bintab*', COUNT=nfile)
  IF nfile EQ 0 THEN $
    MESSAGE, 'Bin Table file not found in current directory!' $
  ELSE BEGIN 
      binfile = binfile[0]
      len = STRLEN(binfile)
      IF STRMID(binfile, len-3, 3) EQ '.gz' THEN BEGIN
          SPAWN, 'gunzip -v ' + binfile
          binfile = STRMID(binfile, 0, len-3)
      ENDIF 
  ENDELSE

;;; OPEN THEN xxxx.bintabs.fits FILE.
  FITS_OPEN, binfile, fcb
;;; DETERMINE WHICH EXTENSION IS THE Design Slits TABLE.
  extnum = WHERE(fcb.extname EQ 'DesiSlits', cnt)
  IF cnt GT 0 THEN design = MRDFITS(binfile, extnum[0],/SILENT) $
  ELSE BEGIN
;;; IF NO Design Slits TABLE IS FOUND THEN SET eflag AND FOLLOW THAT
;;; PORTION OF THE IF/THEN/ESLE STATEMENT.
      MESSAGE, 'ERROR: No DesiSlits table found in '+ $
        'xxxx.bintabs.fits file!'
  ENDELSE

 ;;; DETERMINE WHICH EXTENSION IS THE BluSlits TABLE.
  extnum = WHERE(fcb.extname EQ 'BluSlits', cnt)
  IF cnt GT 0 THEN bluslits = MRDFITS(binfile, extnum[0],/SILENT) $
  ELSE BEGIN
;;; IF NO BluSlits TABLE IS FOUND THEN SET eflag AND FOLLOW THAT
;;; PORTION OF THE IF/THEN/ESLE STATEMENT.
      MESSAGE, 'ERROR: No BluSlits table found in '+ $
        'xxxx.bintabs.fits file!'
  ENDELSE

;;; find index of slit in Design table     
  desidex = WHERE(LONG(STRCOMPRESS(design.slitname, /REMOVE_ALL)) $
                                              EQ slitn, cnt) 
  
  if cnt gt 1 then begin
    message,  'Error: two or more slits with same name in Design table!'
    desidex = desidex[0]
  endif
  if cnt gt 0 then begin
;;; find index of slit in BluSlits table
    bludex =  where(bluslits.dslitid eq design[desidex].dslitid, cnt) 
;;; x and y positions are averages of corner positions
    slitx = (bluslits[bludex].slitx2 + bluslits[bludex].slitx1)/2
    slity = (bluslits[bludex].slity3 + bluslits[bludex].slity1)/2
  endif else message, 'Error: slit not found in Design table!'

  fits_close,fcb

END




