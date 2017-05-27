

PRO add_bintab_dex
;;; FIND FILE CONTAINING BIN TABLES THAT WERE READ FROM RAW DEIMOS FRAME
;;; BY THE PROCEDURE make_bintab_file.pro. 
  binfile = FINDFILE('*.bintabs.fits*', COUNT=nfile)
;;; CHECK HOW MANY FILES WERE FOUND. IF MULTIPLE FOUND, THEN USE
;;; FIRST. IF NONE FOUND, THEN RETURN MESSAGE.
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
;;; DETERMINE HOW MANY EXTENSION ARE IN THE BINTABLE FILE.
  FITS_INFO, binfile, /SILENT, n_ext=n_ext

;;; READ-IN ALL OF THE NEEDED TABLES FROM THE BINTABLE FILE.
  FOR i=1,n_ext DO BEGIN
;;; READ-IN THE .fits HEADER
    hdr = HEADFITS(binfile, ext=i, /SILENT)
;;; PARSE THE HEADER FOR THE EXTENSION NAME AND EXTENSION TYPE.
    extname = STRCOMPRESS(fxpar(hdr, 'EXTNAME'), /REMOVE_ALL)
    type = fxpar(hdr, 'XTENSION')
;;; CHECK THAT THE PRESENT EXTENSION IS A BIN TABLE.
    yes_table = (type EQ 'BINTABLE') or (type EQ 'TABLE')
;;; NOW IF THE EXTENSION IS A BIN TABLE THEN PROCEED. OTHERWISE, GO ON
;;; TO THE NEXT EXTENSION IN THE FILE.
    IF yes_table THEN BEGIN
;;; IF EXTENSION i CONTAINS THE "DesiSlits" TABLE, THEN READ-IN THE
;;; DESIGN STRUCTURE AND SORT SO THAT THE ENTRIES IN design WILL
;;; MATCH-UP CORRECTLY WITH THE ENTRIES IN objectcat. 
      IF extname EQ 'DesiSlits' THEN design = MRDFITS(binfile, i, /SILENT)
;;; IF EXTENSION i CONTAINS THE "SlitObjMap" TABLE, THEN READ-IN THE
;;; SLIT OBJECT MAP STRUCTURE.
      IF extname EQ 'SlitObjMap' THEN slitobjmap = MRDFITS(binfile, i, /SILENT)
;;; IF EXTENSION i CONTAINS THE "MaskDesign" TABLE, THEN READ-IN THE
;;; MASK DESIGN STRUCTURE.
      IF extname EQ 'MaskDesign' THEN mbasics = MRDFITS(binfile, i, /SILENT)
;;; IF EXTENSION i CONTAINS THE "ObjectCat" TABLE, THEN READ-IN THE
;;; OBJECT CAT STRUCTURE AND SORT SO THAT THE ENTRIES IN objectcat WILL
;;; MATCH-UP CORRECTLY WITH THE ENTRIES IN design. 
      IF extname EQ 'ObjectCat' THEN objectcat = MRDFITS(binfile, i, /SILENT) 
    ENDIF
  ENDFOR
;;; NOW ITERATE THRU THE design DATA AND FIND CORRESPONDING ENTRIES IN
;;; slitobjmap AND objcat.
  ndes = N_ELEMENTS(design)
  FOR i=0,ndes-1 DO BEGIN
;;; FIND INDEX IN slitobjmap FOR THE iTH SLITLET (iTH ENTRY IN design).
    map_dex = WHERE(slitobjmap.dslitid EQ design[i].dslitid, mapcnt)
    obj_dex = INTARR(mapcnt)
    FOR j=0,mapcnt-1 DO BEGIN
        obj_dex[j] = $
          WHERE(objectcat.objectid EQ slitobjmap[map_dex[j]].objectid)
    ENDFOR
;;; CHECK FOR SNAGS.
    bad = WHERE(obj_dex EQ -1, num)
    IF num GT 0 THEN MESSAGE, 'ERROR: object not found in objcat table!'
;;; NOW CREATE THE OUTPUT STRUCTURE.
    IF i EQ 0 THEN BEGIN
      FOR j=0,mapcnt-1 DO BEGIN
        IF j EQ 0 THEN BEGIN
          slitnum = LONG(design[i].slitname)
          objnum = LONG(objectcat[obj_dex[j]].object)
          desidex = i
          mapdex = map_dex[j]
          objdex = obj_dex[j]
        ENDIF ELSE BEGIN
          slitnum = [slitnum, LONG(design[i].slitname)]
          objnum = [objnum, LONG(objectcat[obj_dex[j]].object)]
          desidex = [desidex, i]
          mapdex = [mapdex, map_dex[j]]
          objdex = [objdex, obj_dex[j]]
        ENDELSE
      ENDFOR
    ENDIF ELSE BEGIN
      FOR j=0,mapcnt-1 DO BEGIN
        slitnum = [slitnum, LONG(design[i].slitname)]
        objnum = [objnum, LONG(objectcat[obj_dex[j]].object)]
        desidex = [desidex, i]
        mapdex = [mapdex, map_dex[j]]
        objdex = [objdex, obj_dex[j]]
      ENDFOR
    ENDELSE
  ENDFOR

  output = {slitnum:slitnum, objnum:objnum, desidex:desidex, $
            mapdex:mapdex, objdex:objdex}
  FXBHMAKE, hdr, 1, 'BinTab_Index', /DATE, /INITIALIZE
  MWRFITS, output, binfile, hdr, /SILENT

END

