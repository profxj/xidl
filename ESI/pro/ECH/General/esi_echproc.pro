;+ 
; NAME:
; esi_echproc   
;     Version 1.2
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;   
;  esi_echproc, esi, indx, /IFLAT, /REDDOV
;
; INPUTS:
;   esi     -  ESI structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /SVOV    - Save OV files
;  /REDOOV  - Redo OV subtraction
;  /CLOBBER - Clobber existing image
;  /SUBSCAT - Subtract scattered light from data image
;  /NOOSCAN - Raw file has no overscan region (e.g. MIRA failure)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echproc, esi, [20L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echproc, esi, indx, CLOBBER = clobber, SUBSCAT = subscat $
                 , SVOV = svov, FLATFIL = flatfil, TWIFLATFIL = TWIFLATFIL $
                 , BIASFIL = biasfil , NOOSCAN=nooscan $
                 , CBIN = cbin, RBIN = rbin, OBJ_ID = obj_id, NOHOT = NOHOT $
                 , HOT_THRESH = HOT_THRESH

;
  if  N_params() LT 2 and n_elements(OBJ_ID) EQ 0 then begin 
      print, 'Syntax - ' + $
             'esi_echproc, esi, indx, /SVOV, /REDOOV, /CLOBBER, /SUBSCAT, OBJ_ID= [v1.2]'
      return
  endif 
  
  if n_elements(OBJ_ID) NE 0 then begin
      indx = where(esi.obj_id EQ obj_id,nindx)
      if nindx EQ 0 then begin
          print, 'esi_echproc:  No objects with ID = ', obj_id
          print, 'esi_echproc:  Returning..'
          return
      endif
  endif
  
;  Optional Keywords
  if not keyword_set( CBIN ) then cbin = esi[indx[0]].cbin
  if not keyword_set( RBIN ) then rbin = esi[indx[0]].rbin

; NIMG
  nimg = n_elements(indx)

  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      for qq=0L,nimg-1 do begin
          ;; Chk for outfil
          outfil = 'Final/f_'+esi[indx[qq]].img_root
          a = findfile(outfil+'*', count=na)
          if na NE 0 then begin
              print, 'esi_echproc: File ', outfil, ' exists, continuing...'
              esi[indx[qq]].img_final = outfil
              continue
          endif
          flg_proc = 1
      endfor
  endif else flg_proc = 1

  if flg_proc EQ 0 then begin
      print, 'esi_echproc: No files to process. Returning!'
      return
  endif

; Slit
  slit = esi[indx[0]].slit  ; Takes first value so only 1 slit width at a time

; BIAS sub if necessary
  ;if not keyword_set( REDOOV ) then bias = where(esi[indx].flg_ov EQ 0, nbias) $
  ;else begin
  ;    nbias = n_elements(indx)
  ;    bias = lindgen(nbias)
  ;endelse
  ;if nbias NE 0 then begin
  ;  endif
  
  ;; Assign filenames
  IF NOT KEYWORD_SET(BIASFIL) THEN biasfil = $
    esi_getfil('bias_fil', esi[indx[0]].mode, cbin = cbin, rbin = rbin, /name)
  IF NOT KEYWORD_SET(FLATFIL) THEN flatfil = esi[indx[0]].FLAT_FIL
  IF NOT KEYWORD_SET(TWIFLATFIL) THEN twiflatfil = esi[indx[0]].TWIFLAT_FIL
  ;; Bias subtract and save these images to OV directory
  print, 'esi_echproc: Bias subtracting'
  esi_subbias, esi, indx, /FORCE, BIASFIL = biasfil, NOHOT = NOHOT $
               , HOT_THRESH = HOT_THRESH, NOOSCAN=nooscan
  ;; Subtract scattered light (untested in new version of code)
  print, 'esi_echproc: Fixing gain'
  if keyword_set(SUBSCAT) then print, 'esi_echproc: Scatter subtracting'
  esi_echsubscat, esi, indx, GAINONLY = (not keyword_set(SUBSCAT)) $
                  , FLATFIL = flatfil
  ;; Open Flat for indx[0]
  IF NOT KEYWORD_SET(FLATFIL) OR x_chkfil(flatfil+'*') EQ 0 OR $
    strlen(strtrim(flatfil,2)) EQ 0 then begin
      print, 'esi_echproc: Flat file not found', flatfil
      print, '    *********************************'
      print, '    WE ARE NOT FLAT FIELDING IMAGES.'
      print, '    UNLESS YOU ARE PROCESSING ARCS THIS IS BAD' ;; JFH 04/08
      print, '    WE DID APPLY AN ARCHIVED GAIN CORRECTION.'
      print, '    *********************************'
  ENDIF ELSE BEGIN
      print, 'esi_echproc: Loading the Flat: ', flatfil
      img_flat = xmrdfits(flatfil, 0, fhead, /silent)
      ;;esi[indx].flat_fil = flatfil ; Kludge
      ;; Zero out bad pix
      ;;esi_echbadpix, esi[indx[0]], img_flat  ;; JFH just added this??
      gdflt = where(img_flat GT 0.)
      scatt = sxpar(fhead, 'SCATTER')
      norm = sxpar(fhead, 'NORM')
      IF scatt NE 1 OR norm NE 1 then begin
          print, 'esi_echproc: Flat not processed! Returning...'
          return
      endif
  ENDELSE
  ;; Open Twilat for indx[0]
  IF NOT KEYWORD_SET(TWIFLATFIL) OR x_chkfil(twiflatfil+'*') EQ 0 OR $
    strlen(strtrim(twiflatfil,2)) EQ 0 then begin
      print, 'esi_echproc: Twilight Flat file not found', twiflatfil
      print, '    *********************************'
      print, '    WE ARE NOT ILLUMINATION CORRECTING IMAGES '
      print, '    UNLESS YOU ARE PROCESSING ARCS THIS IS NOT OPTIMAL' ;; JFH 04/08
      print, '    *********************************'
  ENDIF ELSE BEGIN
      print, 'esi_echproc: Loading the Twilight Flat: ', twiflatfil
      illum_flat = xmrdfits(twiflatfil, 0, /silent)
      gdtwiflt = where(illum_flat GT 0.)
  ENDELSE
  


  ;IF NOT keyword_set(FLATFIL) then $
  ;  flatfil = $
  ;  esi_getfil('finflat_fil', SLIT = slit, cbin = cbin, rbin = rbin, /name)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  FOR qq=0L,nimg-1 DO BEGIN
      ;; Outfil
      outfil = 'Final/f_'+esi[indx[qq]].img_root
      ;; Read IMG
      img_fil = 'OV/ov_'+esi[indx[qq]].img_root
      IF x_chkfil(img_fil+'*') EQ 0 THEN BEGIN 
          print, 'esi_echproc: Image file does not exist!', img_fil
          stop
      ENDIF
      IF strmatch(esi[indx[qq]].FLAT_FIL, '*' + flatfil + '*') EQ 0 OR $
        strmatch(esi[indx[qq]].TWIFLAT_FIL, '*' + twiflatfil + '*') EQ 0 THEN $
        message, 'esi_ehcproc: Different flat files for different images'
      img = xmrdfits(img_fil, 0, head, /silent)
      sz_img = size(img, /dimensions)
      rnoise = esi[indx[qq]].readno
      gain   = esi[indx[qq]].gain
      ;; multiply image by gain
      img = img*gain
      ;; Create IVAR image
      ivar =  float(1.0/(abs(img - sqrt(2.0)*rnoise) + rnoise^2))
      var  =  (ivar GT 0.0)/(ivar + 3.0*(ivar EQ 0.0))
      ;;var = ((img > 1.)*esi[indx[qq]].gain + esi[indx[qq]].readno^2) 
      ;; Flat field images
      print, 'esi_echproc: Flattening'
      IF KEYWORD_SET(img_flat) THEN begin 
          img_nrm1 = x_normspec(img, img_flat, var, var_nrm, PIX = gdflt) 
          var_nrm1 = var_nrm ;; JXP bug fix?  1/16/2009
      endif ELSE BEGIN
          img_nrm1 = img
          var_nrm1 = var
      ENDELSE
      IF KEYWORD_SET(illum_flat) THEN $
        img_nrm = x_normspec(img_nrm1, illum_flat, var_nrm1, var_nrm $
                             , PIX = gdtwiflt) $
      ELSE BEGIN
          img_nrm = img_nrm1
          var_nrm = var_nrm1
      ENDELSE
      badpix = WHERE(img LE -1.0d4, nbad)
      IF nbad NE 0 THEN BEGIN 
          var_nrm[badpix] = -1.0
          img_nrm[badpix] = 0.0
      ENDIF
      delvarx, img, var
      delvarx, img_nrm1, var_nrm1
      ;; OUTPUT
      print, 'esi_echproc: Writing... ', outfil
      sxaddpar, head, 'FLAT', 'T', flatfil
      esi[indx[qq]].img_final = outfil
      mwrfits, float(img_nrm), outfil, head, /create, /silent
      mwrfits, float(var_nrm), outfil, /silent
      ;; COMPRESS
      print, 'esi_echproc: Compressing... '
      spawn, 'gzip -f '+outfil
      ;; Release
      delvarx, img_nrm, var_nrm
  endfor
  ;; DEL OV
  if not keyword_set( SVOV ) then esi_delov, esi, indx
  print, 'esi_echproc:  All done!'

  return
end


