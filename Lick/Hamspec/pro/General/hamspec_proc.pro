;+ 
; NAME:
; hamspec_proc   
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
;  There are several ways to call hamspec_proc
; 
;  hamspec_proc, hamspec, indx    [processes the image with index=indx]
;  hamspec_proc, hamspec, SETUP=1, OBJ=4  [processes images with setup=1 and
;     obj_id=4]
;  rslt = hamspec_proc_sngl('Rawfil')    [process a single image]
;
; INPUTS:
;   hamspec     -  HIRES structure
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
;             important for bias subtraction on the red chip.
;  /NOFLAT  - Do not flatten
;  FLATFIL  - Name of flat file to use.
;  IMG_FLAT - Flat field data
;  FLAT_IVAR - Flat field data inverse variance
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
;   hamspec_proc, hamspec, [20L,21L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_proc_work
;  hamspec_subbias_sngl
;  hamspec_getfil
;  hamspec_delov
;
; REVISION HISTORY:
;   30-Apr-2002 Written by JXP
;   19-Aug-2003 Added work routine
;   20-Feb-2005 Taken from mike_proc
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hamspec_proc_sngl, rawfil, ARC=arc, $
  REDOOV=redoov, CLOBBER=clobber, SETUP=setup, $
  SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
  IOUTFIL=ioutfil, GAIN=gain, READNO=readno, $
  IMG_FLAT=img_flat, FLAT_IVAR=img_flat_ivar, $
  FRAME=frame, RBIN=rbin, CBIN=cbin, $
  TRIM_CLM=trim_clm

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = hamspec_proc_sngl(rawfil, /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SETUP=, GAIN=, READNO=, SETUP=, IMG_FLAT= /NOFLAT, [v1.1]'
      return, -1
  endif 

; Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.0
  if not keyword_set( READNO ) then readno = 3.0
  if not keyword_set( SETUP ) and NOT keyword_set( FLATFIL ) and $
    not keyword_set( NOFLAT ) and NOT keyword_set( IMG_FLAT) then begin
      print, 'hamspec_proc_sngl: Must give info for flat file!!'
      stop
      return, -1
  endif
  if not keyword_set( FRAME ) then begin
      pos = strpos(rawfil, '.fits')
      frame = long(strmid(rawfil,pos-4,4))
  endif
  if not keyword_set(CBIN) or not keyword_set(RBIN) then begin
      head = xheadfits(rawfil, exten=exten)
      cbin = sxpar(head,'CBIN')
      rbin = sxpar(head,'RBIN')
  endif
      

  ;; Check for existing image
  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      ;; Chk for outfil
      if keyword_set(IOUTFIL) then outfil = ioutfil else stop
;        outfil = hamspec_getfil('fin_fil', SUBFIL=rawfil, /name)
      if x_chkfil(outfil, /silent) NE 0 then begin
          print, 'hamspec_proc_sngl: File ', outfil, ' exists, returning...'
          return, -1
      endif
  endif 

  ;; Grab the info
  hamspec_strct, hamspec, FILE_LIST=[strtrim(rawfil,2)], /nomkdir, /NOFILE
      
; BIAS sub if necessary
  print, 'hamspec_proc_sngl: Bias subtracting'
  rslt = hamspec_subbias_sngl( rawfil, NAMP=hamspec.amp, CLOBBER = REDOOV, $
                             USEBIAS=usebias, BADROW=badrow, $
                             CBIN=cbin, RBIN=rbin, FRAME=frame, GAIN_RATIO=hamspec.ratio_gain)

; Flatten
     
  if NOT keyword_set ( NOFLAT )  then begin
    if NOT keyword_set( IMG_FLAT ) then begin
        if not keyword_set( FLATFIL ) then begin
            tmp = {hamspecstrct}
            tmp.colbin = cbin
            ;flg_f = hamspec_getpixflat(hamspec[exten-1], FIL=flatfil)
            flg_f = -1
            if flg_f EQ -1 then $
              img_flat = hamspec_getfil('npixflt_fil', setup, $
                                      FIL_NM=flatfil, /name) 
        endif
        print, 'hamspec_proc:  Using Flat fil ', flatfil
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
;    outfil = hamspec_getfil('fin_fil', SUBFIL=rawfil, /name)
  ;; Read IMG
  img = hamspec_getfil('ov_fil', FRAME=frame, head=head)
  sz_img = size(img, /dimensions)
          
  ;; Create VAR image
;  ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
  ivar = 1.0/(abs(img*gain  - 7.) + 7. + readno^2 )   ;; JXP 2/10/06
;  ivar = 1.0/((img*gain  - 7.)>0. + 7. + readno^2 )   ;; JXP January 2010
  ;; Should consider a better method and possibly couple with
                                ; extraction
;  ivar = 1.0/((abs(img * gain)>1.) + readno^2 + 20.) 
          
  ;; Flatten
  img_nrm = img

  if NOT keyword_set(NOFLAT) then begin
      print, 'hamspec_proc_sngl: Flattening!'
      img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
  
      ivar_correction = 0 
      if keyword_set(img_flat_ivar) then begin
          ivar_correction = img_nrm^2 * ivar * $
            ((1./(img_flat_ivar + (img_flat_ivar EQ 0))) > 1.0e-6) $
            * (img_flat_ivar GT 0)
      endif
      ivar[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
      print, 'Mean ivar_correction is ', mean(ivar_correction)
      ivar = ivar / (1 + ivar_correction)
      ;; Need to leave the interorder gaps in!!
;      if nbad NE 0 then ivar[badpix] = 0.
      
  endif else print, 'hamspec_proc_sngl: Not flattening!'
  delvarx, img

  ;; Deal with badrows
  if keyword_set( BADROW ) then begin
      for ii=0L,n_elements(badrow)-1 do ivar[*,badrow[ii]] = 0.
  endif

  ;; Add in gain
  fin_img = temporary(img_nrm)*gain
          
  ;; OUTPUT
  print, 'hamspec_proc_sngl: Writing... ', outfil
  mwrfits, fin_img, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  print, 'hamspec_proc_sngl: Compressing... '
  spawn, 'gzip -f '+outfil
      
  ;; DEL OV
  if not keyword_set( SVOV ) then begin
      ovfil = hamspec_getfil('ov_fil', FRAME=frame, /name)
      spawn, '\rm '+ovfil
  endif

  delvarx, fin_img, var_nrm, img, var, img_nrm

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_proc, hamspec, indx, SETUP=setup, ARC=arc, $
               REDOOV=redoov, CLOBBER=clobber, OBJ=obj, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, STD=std, EXP=exp, XOFF=xoff, $
                TRIM_CLM=trim_clm

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hamspec_proc, hamspec, [indx], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' CHIP=, OBJ=, SETUP=, /STD, /ARC [v1.1]'
      return
  endif 
  
;  Optional Keywords
   if NOT keyword_set(redoov) then redoov = 0L

   mfa = hamspec.flg_anly

   if keyword_set(obj) then begin 
     objmask = mfa * 0
     for i=0, n_elements(obj)-1 do begin
       any = where(obj[i] EQ hamspec.obj_id)
       if any[0] NE -1 then objmask[any] = 1
     endfor
     objmask = objmask AND mfa
   endif else objmask = mfa AND $
              ((hamspec.obj_id GT 0) OR (strtrim(hamspec.type,2) EQ 'STD'))

   setupmask = mfa AND ((keyword_set(setup) ? $
                         hamspec.setup EQ setup[0] : 1) GT 0)
   stdmask = mfa AND ((keyword_set(std) ? $
                       strtrim(hamspec.type,2) EQ 'STD':1) GT 0)


   if keyword_set(obj) OR keyword_set(std) OR n_elements(indx) EQ 0 then begin
     indx = where( setupmask AND stdmask AND objmask , nindx)
     print, 'hamspec_proc: Using obj# and overriding indx if input'
   endif 

   ;; Exposure
   if keyword_set(EXP) then indx = indx[exp]
   ;; Nindx
   nindx = n_elements(indx)

   if nindx EQ 0 then begin
     print, 'No objects found to process '+ string(fix([total(setupmask), $
                        total(stdmask), total(objmask)])) 
     return
   endif

   fullindx = indx
   flg_indx = 1
      
   
   cur_setup = setup
   indx = where(hamspec[fullindx].setup EQ cur_setup, nindx)
   if nindx GT 0 then indx = fullindx[indx]
   if nindx EQ 0 then return
   
   ;; NIMG
   nimg = n_elements(indx)
   if nimg EQ 0 then begin
      print, 'hamspec_proc: No files to process. Returning!'
      return
   endif


;;;;;;;;;;;;;;;;
; Flat
   flg_fset = 0
   if NOT keyword_set( NOFLAT) then begin
      if not keyword_set( FLATFIL ) then begin
         flg_f = -1
                                ;flg_f = hamspec_getpixflat(hamspec[indx[0]], FIL=flatfil)
         if flg_f EQ -1 then $
            img_flat = hamspec_getfil('npixflt_fil', cur_setup, $
                                  FIL_NM=flatfil, /name) 
         flg_fset = 1
      endif
      
      print, 'hamspec_proc:  Using Flat fil ', flatfil
      img_flat = xmrdfits(flatfil, 0, fhead, /silent)
      img_flat_ivar = xmrdfits(flatfil, 1, /silent)
      gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                    complement=badpix, ncomplement=nbad)
      hamspec[indx].flat_fil = flatfil ; Kludge
   endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP

   for ii=0L,nimg-1 do begin
      ;; Chk for outfil
      if keyword_set(IOUTFIL) then outfil = ioutfil[ii] else $
         outfil =hamspec_getfil('fin_fil', /name, CHKFIL=chkf, $
                            FRAME=hamspec[indx[ii]].frame)
      
      if keyword_set(chkf) AND NOT keyword_set(clobber) then begin
         print, 'hamspec_proc: File ', outfil, ' exists, continuing...'
         hamspec[indx[ii]].img_final = outfil
         return
      endif
      
      ;; Call proc_singl
      rslt = hamspec_proc_sngl(hamspec[indx[ii]].rootpth+hamspec[indx[ii]].img_root,$
                           IMG_FLAT=img_flat, REDOOV=redoov, $
                           IOUTFIL=outfil, GAIN=hamspec[indx[ii]].gain, $
                           READNO=hamspec[indx[ii]].readno, /SVOV, $
                           /clobber, FLAT_IVAR = img_flat_ivar, $
                           NOFLAT=keyword_set(noflat), $
                           CBIN=hamspec[indx[ii]].colbin, $
                           RBIN=hamspec[indx[ii]].rowbin, $
                           FRAME=hamspec[indx[ii]].frame, $
                           TRIM_CLM=TRIM_CLM $
                          )
      
      ;; Set xyoff
      if keyword_set(XOFF) then $
         hamspec[indx[ii]].arc_xyoff[0] = xoff
      
      ;; OUTPUT
      hamspec[indx[ii]].img_final = outfil
   endfor  ;; Image loop
   
   if keyword_set( FLATFIL ) and FLG_FSET then delvarx, flatfil
   ;; DEL OV
   if not keyword_set( SVOV ) then hamspec_delov, hamspec, indx

  print, 'hamspec_proc: All done!'

  return
end


