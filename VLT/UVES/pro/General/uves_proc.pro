;+ 
; NAME:
; uves_proc   
;     Version 1.0
;
; PURPOSE:
;    Process a raw data frame.  Bias subtract and flatten.  Output is
;    written in Final/ and is a mutli-extension fits file.
;    Extension-0 is the processed image, Extension-1 is the inverse
;    variance.
;
; CALLING SEQUENCE:
;   
;  There are several ways to call uves_proc
; 
;  uves_proc, uves, indx    [processes the image with index=indx]
;  uves_proc, uves, SETUP=1, OBJ=4  [processes images with setup=1 and
;     obj_id=4]
;  rslt = uves_proc_sngl('Rawfil')    [process a single image]
;
; INPUTS:
;   uves     -  HIRES structure
;   [indx]   -  Index (or indices) of the image(s) in the HIRES
;              structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /CLOBBER - Clobber existing image
;  SETUP    - Setup number
;  EXP      - Exposure numbers
;  OBJ      - Obj_id value
;  /ARC     - Identifies the image as an Arc explicitly.  This is
;             important for bias subtraction on the red side.
;  /NOFLAT  - Do not flatten
;  FLATFIL  - Name of flat file to use.
;  /SVOV    - Save the OV subtracted image
;  /REDOOV  - Overwrite any existing OV images
;  IOUTFIL  - Name for output file (array for mult
;  /STD     - Processes standard stars corresponding to an input SETUP
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The inverse variance is given this set of values in an attempt to
;   deal with the Poisson limit
;
;   ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
;
; EXAMPLES:
;   uves_proc, uves, [20L,21L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_proc_work
;  uves_subbias_sngl
;  uves_getfil
;  uves_delov
;
; REVISION HISTORY:
;   30-Apr-2002 Written by JXP
;   19-Aug-2003 Added work routine
;   20-Feb-2005 Taken from mike_proc
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_proc_sngl, rawfil, setup, side, ARC=arc, $
                         REDOOV=redoov, CLOBBER=clobber, EXTEN=exten, $
                         SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
                         IOUTFIL=ioutfil, GAIN=gain, READNO=readno, $
                         IMG_FLAT=img_flat, FLAT_IVAR=img_flat_ivar, $
                          FRAME=frame, RBIN=rbin, CBIN=cbin, IV_CORRECT=iv_correct

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = uves_proc_sngl(rawfil, setup, side,/SVOV, /REDOOV, /CLOBBER, ' + $
        ' SETUP=, GAIN=, READNO=, SETUP=, IMG_FLAT= /NOFLAT, [v1.1]'
      return, -1
  endif 

; Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.0
  if not keyword_set( READNO ) then readno = 3.0
  if not keyword_set( SETUP ) and NOT keyword_set( FLATFIL ) and $
    not keyword_set( NOFLAT ) and NOT keyword_set( IMG_FLAT) then begin
      print, 'uves_proc_sngl: Must give info for flat file!!'
      stop
      return, -1
  endif
;  if not keyword_set( FRAME ) then begin
;      pos = strpos(rawfil, '.fits')
;      frame = long(strmid(rawfil,pos-4,4))
;  endif
;  if not keyword_set(CBIN) or not keyword_set(RBIN) then begin
;      head = xheadfits(rawfil, exten=exten)
;      cbin = round(2048./sxpar(head,'NAXIS1'))
;      rbin = round(4096./sxpar(head,'NAXIS2'))
;  endif
      

  ;; Check for existing image
  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      ;; Chk for outfil
      if keyword_set(IOUTFIL) then outfil = ioutfil else stop
;        outfil = uves_getfil('fin_fil', SUBFIL=rawfil, /name)
      if x_chkfil(outfil, /silent) NE 0 then begin
          print, 'uves_proc_sngl: File ', outfil, ' exists, returning...'
          return, -1
      endif
  endif 

  ;; Grab the info
  uves_strct, uves, FILE_LIST=[rawfil], /nomkdir, /NOFILE
      
; BIAS sub if necessary
  print, 'uves_proc_sngl: Bias subtracting'
  rslt = uves_subbias_sngl( rawfil, side, CBIN=uves.colbin, RBIN=uves.rowbin, $
                            POSCAN=uves.poscan, EXTEN=exten, /NOBIASROW)

; Flatten
     
  if NOT keyword_set ( NOFLAT )  then begin
      if NOT keyword_set( IMG_FLAT ) then begin
          if not keyword_set( FLATFIL ) then begin
              tmp = {uvesstrct}
              tmp.chip = chip
              tmp.colbin = cbin
;              flg_f = uves_getpixflat(uves[exten-1], FIL=flatfil)
              if flg_f EQ -1 then $
                img_flat = uves_getfil('nqtz_fil', setup, $
                                       CHIP=chip, FIL_NM=flatfil, /name) 
          endif
          print, 'uves_proc:  Using Flat fil ', flatfil
          img_flat = xmrdfits(flatfil, 0, fhead, /silent)
          img_flat_ivar = xmrdfits(flatfil, 1, /silent)
      endif
      
      gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                    complement=badpix, ncomplement=nbad)
  endif
      


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  ;; Outfil
  if keyword_set(IOUTFIL) then outfil = ioutfil else stop
;    outfil = uves_getfil('fin_fil', SUBFIL=rawfil, /name)
  ;; Read IMG
  img = uves_getfil('ov_fil', OBJN=uves.img_root, head=head)
  sz_img = size(img, /dimensions)
          
  ;; Create VAR image
;  ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
  ivar = 1.0/(abs(img*uves.gain  - 7.) + 7. + uves.readno^2 )   ;; JXP 2/10/06
  ;; Should consider a better method and possibly couple with
                                ; extraction
;  ivar = 1.0/((abs(img * gain)>1.) + readno^2 + 20.) 
          
  ;; Flatten
  img_nrm = img

  if NOT keyword_set(NOFLAT) then begin
      print, 'uves_proc_sngl: Flattening!'
      img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
  
      ivar_correction = 0 
      if keyword_set(img_flat_ivar)  AND keyword_set(IV_CORRECT) then begin
          ivar_correction = img_nrm^2 * ivar * $
            ((1./(img_flat_ivar + (img_flat_ivar EQ 0))) > 1.0e-6) $
            * (img_flat_ivar GT 0)
      endif
      ivar[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
      if keyword_set(IV_CORRECT) then begin
          print, 'Mean ivar_correction is ', mean(ivar_correction)
          ivar = ivar / (1 + ivar_correction)
      endif
      
  endif else print, 'uves_proc_sngl: Not flattening!'
  delvarx, img

  ;; Deal with badrows
  if keyword_set( BADROW ) then begin
      for ii=0L,n_elements(badrow)-1 do ivar[*,badrow[ii]] = 0.
  endif
          
  ;; Add in gain
  fin_img = temporary(img_nrm)*uves.gain
          
  ;; OUTPUT
  print, 'uves_proc_sngl: Writing... ', outfil
  mwrfits, fin_img, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  print, 'uves_proc_sngl: Compressing... '
  spawn, 'gzip -f '+outfil
      
  ;; DEL OV
  if not keyword_set( SVOV ) then begin
      ovfil = uves_getfil('ov_fil', OBJN=uves.img_root, /name)
      spawn, '\rm '+ovfil
  endif

  delvarx, fin_img, var_nrm, img, var, img_nrm

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_proc, uves, indx, SETUP=setup, ARC=arc, SIDE=side, $
               REDOOV=redoov, CLOBBER=clobber, OBJ=obj, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, STD=std, EXP=exp, XOFF=xoff

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'uves_proc, uves, [indx], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SIDE=, OBJ=, SETUP=, /STD, /ARC [v1.1]'
      return
  endif 
  
;  Optional Keywords
   if NOT keyword_set(redoov) then redoov = 0L

   sand = keyword_set(side) ? long(total(side)) : 0L
   mfa = uves.flg_anly
   cside = mfa AND ((keyword_set(side) ? uves.side EQ sand : 1) GT 0)

   if keyword_set(obj) then begin 
     objmask = mfa * 0
     for i=0, n_elements(obj)-1 do begin
       any = where(obj[i] EQ uves.obj_id)
       if any[0] NE -1 then objmask[any] = 1
     endfor
     objmask = objmask AND mfa
   endif else objmask = mfa AND $
              ((uves.obj_id GT 0) OR (strtrim(uves.type,2) EQ 'STD'))

   setupmask = mfa AND ((keyword_set(setup) ? $
                         uves.setup EQ setup[0] : 1) GT 0)
   stdmask = mfa AND ((keyword_set(std) ? $
                       strtrim(uves.type,2) EQ 'STD':1) GT 0)


   if keyword_set(obj) OR keyword_set(std) OR keyword_set(indx) EQ 0 then begin
     indx = where( setupmask AND stdmask AND objmask AND cside, nindx)
     print, 'uves_proc: Using obj# and overriding indx if input'
   endif 

   ;; Exposure
   if keyword_set(EXP) then begin
       indx = indx[exp]
       nindx = n_elements(indx)
   endif
   ;; Nindx

   if nindx EQ 0 then begin
     print, 'No objects found to process '+ string(fix([total(setupmask), $
                        total(stdmask), total(objmask), total(cside)])) 
     stop
     return
   endif

   fullindx = indx
   allside = uves[fullindx].side + 10*uves[fullindx].setup
   allside = allside[sort(allside)]
   uside = allside[uniq(allside)]
   flg_indx = 1
      
   ;; CHIP
  for kk=0L,n_elements(uside)-1 do begin
      cur_side = uside[kk] mod 10
      cur_setup = uside[kk] / 10


      indx = where(uves[fullindx].side EQ cur_side AND $
                   uves[fullindx].setup EQ cur_setup, nindx)
      if nindx GT 0 then indx = fullindx[indx]
      if nindx EQ 0 then continue


      ;; NIMG
      nimg = n_elements(indx)
      if nimg EQ 0 then begin
        print, 'uves_proc: No files to process. Returning!'
        continue
      endif

      wcen = uves[indx[0]].xdangl

;;;;;;;;;;;;;;;;
; Flat
      flg_fset = 0
      if NOT keyword_set( NOFLAT) then begin
          if not keyword_set( FLATFIL ) then begin
;              flg_f = uves_getpixflat(uves[indx[0]], FIL=flatfil)
;              if flg_fset EQ -1 then $
                img_flat = uves_getfil('nqtz_fil', cur_setup, $
                                        WCEN=wcen, FIL_NM=flatfil, /name) 
              flg_fset = 1
          endif
          
          print, 'uves_proc:  Using Flat fil ', flatfil
          img_flat = xmrdfits(flatfil, 0, fhead, /silent)
          img_flat_ivar = xmrdfits(flatfil, 1, /silent)
          gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                        complement=badpix, ncomplement=nbad)
          uves[indx].flat_fil = flatfil ; Kludge
      endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP

      for ii=0L,nimg-1 do begin
           ;; Chk for outfil
           if keyword_set(IOUTFIL) then outfil = ioutfil[ii] else $
             outfil =uves_getfil('fin_fil', /name, OBJN=uves[indx[ii]].obj,$
                                 FRAME=uves[indx[ii]].frame, $
                                  CHKFIL=chkf)

             if keyword_set(chkf) AND NOT keyword_set(clobber) then begin
                print, 'uves_proc: File ', outfil, ' exists, continuing...'
                uves[indx[ii]].img_final = outfil
                continue
             endif

          ;; Call proc_singl
             if cur_side EQ 2 then EXTEN=1 else EXTEN=0
             rslt = uves_proc_sngl(uves[indx[ii]].rootpth+uves[indx[ii]].img_root,$
                                   cur_setup, cur_side, $
                                   EXTEN=exten, $
                                   IMG_FLAT=img_flat, REDOOV=redoov, $
                                   IOUTFIL=outfil, GAIN=uves[indx[ii]].gain, $
                                   READNO=uves[indx[ii]].readno, /SVOV, $
                                /clobber, FLAT_IVAR = img_flat_ivar, $
                                   NOFLAT=keyword_set(noflat), $
                                   CBIN=uves[indx[ii]].colbin, $
                                   RBIN=uves[indx[ii]].rowbin $
                               )
          
          ;; Set xyoff
          if keyword_set(XOFF) then $
            uves[indx[ii]].arc_xyoff[0] = xoff
          
          ;; OUTPUT
          uves[indx[ii]].img_final = outfil
      endfor  ;; Image loop
  
      if keyword_set( FLATFIL ) and FLG_FSET then delvarx, flatfil
      ;; DEL OV
      if not keyword_set( SVOV ) then uves_delov, uves, indx
  endfor  ;; CHIP loop

  print, 'uves_proc: All done!'

  return
end


