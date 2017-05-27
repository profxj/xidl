;+
; NAME:
; find_skyslits
;
; PURPOSE:
; Find numbers of sky-only slits for a mask.  
;
; CALLING SEQUENCE:
; slitn = find_skyslits, datadir
;
; INPUTS:
;
; datadir - the data directory containing the spectral data for the
;           mask. Empty string (or './') for the current directory.
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; slitn   - array containing the sky slit numbers.
;
; OPTIONAL OUTPUTS:
;
; slitx   - array containing sky slit x positions.
; slity   - guess what!
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

;----------------------------------------------------------------------


FUNCTION find_skyslits,  datadir,  slitx=slitx,  slity=slity

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
;;; DETERMINE WHICH EXTENSION IS THE Object Catalog TABLE.
  extnum = WHERE(fcb.extname EQ 'ObjectCat', cnt)
  IF cnt GT 0 THEN objcat = MRDFITS(binfile, extnum[0],/SILENT) $
  ELSE BEGIN
;;; IF NO Object Catalog TABLE IS FOUND THEN SET eflag AND FOLLOW THAT
;;; PORTION OF THE IF/THEN/ESLE STATEMENT.
      MESSAGE, 'ERROR: No ObjectCat table found in '+ $
          'xxxx.bintabs.fits file!'
  ENDELSE
;;; DETERMINE WHICH EXTENSION IS THE Slit Object Map TABLE.
  extnum = WHERE(fcb.extname EQ 'SlitObjMap', cnt)
  IF cnt GT 0 THEN slitmap = MRDFITS(binfile, extnum[0],/SILENT) $
  ELSE BEGIN
;;; IF NO Slit Object Map TABLE IS FOUND THEN SET eflag AND FOLLOW THAT
;;; PORTION OF THE IF/THEN/ESLE STATEMENT.
      MESSAGE, 'ERROR: No ObjectCat table found in '+ $
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
   
;;; FIND THE SKY SLITS. GRAB THE THIRD DIGIT OF THE OBJECT NUMBER AND
;;; USE IT TO DIFFERENTIATE SKY SLITS FROM NON-SKY SLITS.
  skydex = WHERE(STRMID(objcat.object,2,1) EQ '5', skycnt)
  IF skycnt GT 0 THEN BEGIN
    slitn = LONARR(skycnt)
    slitx = fltarr(skycnt)
    slity = fltarr(skycnt)
    FOR i=0,skycnt-1 DO BEGIN
      mapdex = WHERE(slitmap.objectid EQ objcat[skydex[i]].objectid, cnt)
      IF cnt GT 0 THEN $
        desidex = WHERE(design.dslitid EQ slitmap[mapdex].dslitid, cnt) $
      ELSE MESSAGE, 'Error: sky slit not found in slit map table!'     
      IF cnt GT 0 THEN $
        slitn[i] = LONG(STRCOMPRESS(design[desidex].slitname, /REMOVE_ALL)) $
      ELSE MESSAGE, 'Error: sky slit not found in design table!'
      if cnt gt 0 then $
        bludex =  where(bluslits.dslitid eq slitmap[mapdex].dslitid, cnt) 
      
      if cnt gt 0 then begin
        slitx[i] = (bluslits[bludex].slitx2 + bluslits[bludex].slitx1)/2
        slity[i] = (bluslits[bludex].slity3 + bluslits[bludex].slity1)/2
       endif else message, 'Error: sky slit not found in BluSlits table!'
    ENDFOR
  ENDIF ELSE MESSAGE, 'No sky slits found!'

  FITS_CLOSE, fcb

  RETURN, slitn


END
