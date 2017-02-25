;+ 
; NAME:
; mike_proc   
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
;  There are several ways to call mike_proc
; 
;  mike_proc, mike, indx    [processes the image with index=indx]
;  mike_proc, mike, SETUP=1, OBJ=4  [processes images with setup=1 and
;     obj_id=4]
;  rslt = mike_proc_sngl('Rawfil')    [process a single image]
;
; INPUTS:
;   mike     -  MIKE structure
;   [indx]   -  Index (or indices) of the image(s) in the MIKE
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
;   mike_proc, mike, [20L,21L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_proc_work
;  mike_subbias_sngl
;  mike_getfil
;  mike_delov
;
; REVISION HISTORY:
;   30-Apr-2002 Written by JXP
;   19-Aug-2003 Added work routine
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_proc_sngl, rawfil, side, ARC=arc, $
  REDOOV=redoov, CLOBBER=clobber, SETUP=setup, $
  SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
  IOUTFIL=ioutfil, GAIN=gain, READNO=readno, $
  IMG_FLAT=img_flat, FLAT_IVAR=img_flat_ivar, $
  _EXTRA=extra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = mike_proc_sngl(rawfil, [side], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SETUP=, GAIN=, READNO=, SETUP=, IMG_FLAT= /NOFLAT, [v1.1]'
      return, -1
  endif 

; Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.06
  if not keyword_set( READNO ) then readno = 3.70
  if not keyword_set( SETUP ) and NOT keyword_set( FLATFIL ) and $
    not keyword_set( NOFLAT ) and NOT keyword_set( IMG_FLAT) then begin
      print, 'mike_proc_sngl: Must give info for flat file!!'
      stop
      return, -1
  endif
  if not keyword_set( SIDE ) then begin
      head = xheadfits(rawfil)
      instr = sxpar(head,'INSTRUME')
      if strmid(instr,5,3) EQ 'Blu' then side = 1 else side = 2
  endif 
      

  ;; Check for existing image
  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      ;; Chk for outfil
      if keyword_set(IOUTFIL) then outfil = ioutfil else $
        outfil = mike_getfil('fin_fil', SUBFIL=rawfil, /name)
      if x_chkfil(outfil, /silent) NE 0 then begin
          print, 'mike_proc_sngl: File ', outfil, ' exists, returning...'
          return, -1
      endif
  endif 
      
; BIAS sub if necessary
  print, 'mike_proc_sngl: Bias subtracting'
  rslt = mike_subbias_sngl( rawfil, side, CLOBBER=REDOOV, ARC=arc, $
                            USEBIAS=usebias, BADROW=badrow, _EXTRA=extra)

; Flatten
     
  if NOT keyword_set ( NOFLAT )  then begin
    if NOT keyword_set( IMG_FLAT ) then begin
      if not keyword_set( FLATFIL ) then $
        flat_fil = mike_getfil('mflat_fil', setup, $
                               SIDE=side, FIL_NM=flatfil, HEAD=fhead, /name) 
       img_flat = xmrdfits(flatfil, 0, fhead, /silent)
       img_flat_ivar = xmrdfits(flatfil, 1, /silent)
       
    endif

    gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1 and $
                  (img_flat_ivar GT 0.), $
                complement=badpix, ncomplement=nbad)
  endif
      


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  ;; Outfil
  if keyword_set(IOUTFIL) then outfil = ioutfil else $
    outfil = mike_getfil('fin_fil', SUBFIL=rawfil, /name)
  ;; Read IMG
  img = mike_getfil('ov_fil', subfil=rawfil, head=head)
  sz_img = size(img, /dimensions)
          
  ;; Create VAR image
  ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
;  ivar = 1.0/((abs(img * gain)>1.) + readno^2 + 20.) 
          
  ;; Flatten
  img_nrm = img

  if NOT keyword_set(NOFLAT) then begin
      print, 'mike_proc_sngl: Flattening!'
    img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
  
    ivar_correction = 0 
    if keyword_set(img_flat_ivar) then begin
        ivar_correction = img_nrm^2 * ivar * ((1./(img_flat_ivar $
                                                   + (img_flat_ivar EQ 0))) $
                                              > 1.0e-6) * (img_flat_ivar GT 0)
    endif
    ivar[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
    print, 'Mean ivar_correction is ', mean(ivar_correction)
    ivar = ivar / (1 + ivar_correction)
    if nbad NE 0 then ivar[badpix] = 0.

  endif else print, 'mike_proc_sngl: Not flattening!'
  delvarx, img

  ;; Deal with badrows
  if keyword_set( BADROW ) then begin
      for ii=0L,n_elements(badrow)-1 do ivar[*,badrow[ii]] = 0.
  endif
          
  ;; Add in gain
  fin_img = temporary(img_nrm)*gain
          
  ;; OUTPUT
  print, 'mike_proc_sngl: Writing... ', outfil
  mwrfits, fin_img, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  print, 'mike_proc_sngl: Compressing... '
  spawn, 'gzip -f '+outfil
      
  ;; DEL OV
  if not keyword_set( SVOV ) then $
    spawn, '\rm '+mike_getfil('ov_fil', subfil=rawfil, /name)

  delvarx, fin_img, var_nrm, img, var, img_nrm

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_proc, mike, indx, SETUP=setup, ARC=arc, SIDE=side, $
               REDOOV=redoov, CLOBBER=clobber, OBJ=obj, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, STD=std, EXP=exp, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_proc, mike, [indx], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SIDE=, OBJ=, SETUP=, /STD, /ARC [v1.1]'
      return
  endif 
  
;  Optional Keywords
   if NOT keyword_set(redoov) then redoov = 0L

   sand = keyword_set(side) ? long(total(side)) : 0L
   mfa = mike.flg_anly
   sidemask = mfa AND ((keyword_set(side) ? mike.side AND sand : 1) GT 0)

   if keyword_set(obj) then begin 
     objmask = mfa * 0
     for i=0, n_elements(obj)-1 do begin
       any = where(obj[i] EQ mike.obj_id)
       if any[0] NE -1 then objmask[any] = 1
     endfor
     objmask = objmask AND mfa
   endif else objmask = mfa AND $
              ((mike.obj_id GT 0) OR (strtrim(mike.type,2) EQ 'STD'))

   setupmask = mfa AND ((keyword_set(setup) ? mike.setup EQ setup[0] : 1) GT 0)
   stdmask = mfa AND ((keyword_set(std) ? strtrim(mike.type,2) EQ 'STD':1) GT 0)


   if keyword_set(obj) OR keyword_set(std) OR keyword_set(indx) EQ 0 then begin
     indx = where( setupmask AND stdmask AND objmask AND sidemask, nindx)
     print, 'mike_proc: Using obj# and overriding indx if input'
   endif 

   ;; Exposure
   if keyword_set(EXP) then indx = indx[exp]
   ;; Nindx
   nindx = n_elements(indx)

   if nindx EQ 0 then begin
     print, 'No objects found to process '+ string(fix([total(setupmask), $
                        total(stdmask), total(objmask), total(sidemask)])) 
     return
   endif

   fullindx = indx
   allside = mike[fullindx].side + 10*mike[fullindx].setup
   allside = allside[sort(allside)]
   sidesetup = allside[uniq(allside)]
   flg_indx = 1
      
; SIDE
  for kk=0L,n_elements(sidesetup)-1 do begin
      cur_side = sidesetup[kk] mod 10
      cur_setup = sidesetup[kk] / 10

      indx = where(mike[fullindx].side EQ cur_side AND $
                   mike[fullindx].setup EQ cur_setup, nindx)
      if nindx GT 0 then indx = fullindx[indx]
      if nindx EQ 0 then continue

      ;; NIMG
      nimg = n_elements(indx)
      if nimg EQ 0 then begin
        print, 'mike_proc: No files to process. Returning!'
        continue
      endif


;;;;;;;;;;;;;;;;
; Flat
     if NOT keyword_set( NOFLAT) then begin
      if not keyword_set( FLATFIL ) then $
        img_flat = mike_getfil('mflat_fil', cur_setup, $
                               SIDE=cur_side, FIL_NM=flatfil) 

        img_flat = xmrdfits(flatfil, 0, fhead, /silent)
        img_flat_ivar = xmrdfits(flatfil, 1, /silent)
      gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                 complement=badpix, ncomplement=nbad)
      mike[indx].flat_fil = flatfil ; Kludge
    endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP

      for ii=0L,nimg-1 do begin
           ;; Chk for outfil
           if keyword_set(IOUTFIL) then outfil = ioutfil[ii] else $
             outfil =mike_getfil('fin_fil', subfil=mike[indx[ii]].img_root, $
                                     /name, CHKFIL=chkf)

             if keyword_set(chkf) AND NOT keyword_set(clobber) then begin
                print, 'mike_proc: File ', outfil, ' exists, continuing...'
                mike[indx[ii]].img_final = outfil
                continue
             endif

          ;; Call proc_singl
          rslt = mike_proc_sngl(mike[indx[ii]].rootpth+mike[indx[ii]].img_root,$
                                cur_side, IMG_FLAT=img_flat, REDOOV=redoov, $
                                IOUTFIL=outfil, GAIN=mike[indx[ii]].gain, $
                                READNO=mike[indx[ii]].readno, /SVOV, $
                                ARC=arc, /clobber, FLAT_IVAR = img_flat_ivar, $
                                NOFLAT=keyword_set(noflat), _EXTRA=extra)
          
          ;; Set xyoff
          arcs = where(mike.setup EQ cur_setup AND mike.side EQ cur_side AND $
                       strtrim(mike.type,2) EQ 'ARC' AND $
                       mike.flg_anly NE 0, narc)
          diff = abs( (mike[indx[ii]].date + $
                                   mike[indx[ii]].exp/2./3600./24.) - $
                                  (mike[arcs].date + mike[arcs].exp/2./3600./24.) )
          mn = min(diff, gdarc)
;;          diff = abs( (mike[indx[ii]].ut + $
;;                       mike[indx[ii]].exp/2.) - $
;;                      (mike[arcs].ut + mike[arcs].exp/2.) )
;;          b = where( diff GT 21.*3600, nb)
;;          if nb NE 0 then diff[b] = abs(diff[b] - 24.*3600.)
          mn = min(diff, gdarc)
          mike[indx[ii]].arc_xyoff = mike[arcs[gdarc]].arc_xyoff
          
          ;; OUTPUT
          mike[indx[ii]].img_final = outfil
      endfor  ;; Image loop
  
      ;; DEL OV
      if not keyword_set( SVOV ) then mike_delov, mike, indx
  endfor  ;; SIDE loop

  print, 'mike_proc: All done!'

  return
end


