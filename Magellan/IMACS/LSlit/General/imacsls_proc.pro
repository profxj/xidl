;+ 
; NAME:
; imacsls_proc   
;     Version 1.0
;
; PURPOSE:
;    Process an image (bias subtract, flatten) and create variance array.  The
;  routine also turns the DN into electrons by scaling by the gain.  This
;  routine is a mirror of mike_proc
;
; CALLING SEQUENCE:
;  imacsls_proc, imacsls, indx, SETUP=, /ARC, SIDE=, /REDOOV, /CLOBBER
;
; INPUTS:
;   imacsls -  IMACS long slit structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /REDOOV  - Redo OV subtraction
;  /CLOBBER - Clobber existing image
;  /ARC - Process arcs (no ov subtraction for now)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_proc, imacsls, [20L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function imacsls_proc_sngl, rawfil, side, ARC=arc, $
               REDOOV=redoov, CLOBBER=clobber, SETUP=setup, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, GAIN=gain, READNO=readno

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = imacsls_proc_sngl(rawfil, [side], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SIDE=, OBJ=, SETUP=, GAIN=gain, RN=rn, SETUP= [v1.0]'
      return, -1
  endif 

; Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.06
  if not keyword_set( READNO ) then readno = 3.70
  if not keyword_set( SETUP ) and NOT keyword_set( FLATFIL ) and $
    not keyword_set( NOFLAT ) then begin
      print, 'imacsls_proc_sngl: Must give info for flat file!!'
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
        outfil = imacsls_getfil('fin_fil', SUBFIL=rawfil, /name)
      if x_chkfil(outfil, /silent) NE 0 then begin
          print, 'imacsls_proc_sngl: File ', outfil, ' exists, returning...'
          return, -1
      endif
  endif 
      
; BIAS sub if necessary
  print, 'imacsls_proc_sngl: Bias subtracting'
  stop
;  rslt = imacsls_subbias_sngl( rawfil, side, CLOBBER=REDOOV, ARC=arc, /usebias)

  
; Scatter subtract (or just do gain)
;  print, 'imacsls_proc: Fixing gain'
;  if keyword_set(SUBSCAT ) then print, 'imacsls_proc: Scatter subtracting'
;  imacsls_subscat, imacsls, indx, GAINONLY=(not keyword_set(SUBSCAT)), $
;    FLATFIL=flatfil

; Flatten
  if not keyword_set( FLATFIL ) then $
    img_flat = imacsls_getfil('mflat_fil', setup, $
                           SIDE=side, FIL_NM=flatfil, HEAD=fhead) $
  else img_flat = xmrdfits(flatfil, 0, fhead, /silent)
  gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, complement=badpix, ncomplement=nbad)

  ;; Zero out bad pix
;  imacsls_badpix, img_flat
;  gdflt = where(img_flat GT 0.)
;  scatt = sxpar(fhead,'SCATTER')
;  norm = sxpar(fhead,'NORM')
;  if scatt NE 1 OR norm NE 1 then begin
;      print, 'imacsls_proc: Flat not processed! Returning...'
;      return
;  endif
          

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  ;; Outfil
  if keyword_set(IOUTFIL) then outfil = ioutfil else $
    outfil = imacsls_getfil('fin_fil', SUBFIL=rawfil, /name)
  ;; Read IMG
  img = imacsls_getfil('ov_fil', subfil=rawfil, head=head)
  sz_img = size(img, /dimensions)
          
  ;; Create VAR image
  ivar = 1.0/((abs(img * gain)>1.) + readno^2 + 20.) 
          
  ;; Flatten
  print, 'imacsls_proc_sngl: Flattening'
  img_nrm = fltarr(sz_img[0],sz_img[1])
  img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
  ivar[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
  if nbad NE 0 then ivar[badpix] = -1.
  delvarx, img
          
  ;; Add in gain
  fin_img = temporary(img_nrm)*gain
          
  ;; OUTPUT
  print, 'imacsls_proc_sngl: Writing... ', outfil
  mwrfits, fin_img, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  print, 'imacsls_proc_sngl: Compressing... '
  spawn, 'gzip -f '+outfil
      
  ;; DEL OV
  if not keyword_set( SVOV ) then $
    spawn, '\rm '+imacsls_getfil('ov_fil', subfil=rawfil, /name)

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_proc, imacsls, indx, SETUP=setup, ARC=arc, SIDE=side, $
               REDOOV=redoov, CLOBBER=clobber, OBJ=obj, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, STD=std

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'imacsls_proc, imacsls, indx, /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SIDE=, OBJ=, SETUP= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if keyword_set( OBJ ) then begin
      print, 'imacsls_proc: Using obj# and overriding indx if input'
      if not keyword_set( SETUP ) then begin
          print, 'imacsls_proc: Need to set SETUP card'
          return
      endif
      if not keyword_set( SIDE ) then side = [1L,2L]
      flg_indx = 0
  endif else begin
      if keyword_set( STD ) then  begin
          indx = where(imacsls.flg_anly NE 0 AND $
                       imacsls.setup EQ setup AND $
                       strtrim(imacsls.type,2) EQ 'STD', nindx)
          if nindx EQ 0 then stop
      endif
      ;; Set side
      allside = imacsls[indx].side
      ;;
      allside = allside[sort(allside)]
      side = allside[uniq(allside)]
      side = 1
      flg_indx = 1
  endelse
      
  if not keyword_set(setup) then setup = imacsls[indx[0]].setup

; SIDE
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]

      if flg_indx EQ 0 then begin
          indx = where(imacsls.flg_anly NE 0 AND $
                       imacsls.obj_id EQ OBJ AND imacsls.setup EQ setup AND $
                       imacsls.side EQ qq AND $
                       strtrim(imacsls.type,2) EQ 'OBJ', nindx)
          if nindx EQ 0 then continue
      endif 

      ;; NIMG
      nimg = n_elements(indx)

      flg_proc = 0
      if not keyword_set( CLOBBER ) then begin
          for ll=0L,nimg-1 do begin
              ;; Chk for outfil
              if keyword_set(IOUTFIL) then outfil = ioutfil[ll] else $
                outfil = imacsls_getfil('fin_fil', subfil=imacsls[indx[ll]].img_root, $
                                     /name, CHKFIL=chkf)
              if keyword_set(CHKF) then begin
                  print, 'imacsls_proc: File ', outfil, ' exists, continuing...'
                  imacsls[indx[ll]].img_final = outfil
                  continue
              endif
              flg_proc = 1
          endfor
      endif else flg_proc = 1
      
      if flg_proc EQ 0 then begin
          print, 'imacsls_proc: No files to process. Returning!'
          return
      endif
      
      
; BIAS sub if necessary
      if not keyword_set( REDOOV ) then $
        bias = where(imacsls[indx].flg_ov EQ 0, nbias) $
      else begin
          nbias = n_elements(indx)
          bias = lindgen(nbias)
      endelse
      if nbias NE 0 then begin
          print, 'imacsls_proc: Bias subtracting'
          imacsls_subbias, imacsls, indx[bias]
      endif
      
      ;; Flat
      if not keyword_set( STD ) then begin
          if not keyword_set( FLATFIL ) then $
            img_flat = imacsls_getfil('flat_fil', setup, $
                                      SIDE=qq, FIL_NM=flatfil, HEAD=fhead) $
          else img_flat = xmrdfits(flatfil, 0, fhead, /silent)
          print, 'imacsls_proc: Using flat ', flatfil
          gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                        complement=badpix, ncomplement=nbad)
          imacsls[indx].flat_fil = flatfil ; Kludge
      endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP

      for qq=0L,nimg-1 do begin
          ;; Outfil
          if keyword_set(IOUTFIL) then outfil = ioutfil[qq] else $
            outfil = imacsls_getfil('fin_fil', $
                                    SUBFIL=imacsls[indx[qq]].img_root, /name)
          ;; Read IMG
          img = imacsls_getfil('ov_fil', $
                               subfil=imacsls[indx[qq]].img_root, head=head)
          sz_img = size(img, /dimensions)
          
          ;; Create IVAR image
          ;;  I think this is a better unbiased estimator for low counts
          gn = imacsls[indx[qq]].gain
          rn = imacsls[indx[qq]].readno
          ivar = 1.0/(abs(img * gn - sqrt(2.0)*rn) + rn^2 ) 

          if not keyword_set( STD ) then begin
              ;; Flatten
              print, 'imacsls_proc: Flattening'
              img_nrm = fltarr(sz_img[0],sz_img[1])
              ivar_nrm = img_nrm
              img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
              ivar_nrm[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
              if nbad NE 0 then ivar[badpix] = -1.
              delvarx, img
          endif else begin
              img_nrm = temporary(img)
              ivar_nrm = temporary(ivar)
          endelse
          
          ;; Multiply in gain
          fin_img = temporary(img_nrm)*imacsls[indx[qq]].gain
          
          ;; OUTPUT
          print, 'imacsls_proc: Writing... ', outfil
          imacsls[indx[qq]].img_final = outfil
          mwrfits, fin_img, outfil, head, /create, /silent
          mwrfits, ivar_nrm, outfil, /silent
          print, 'imacsls_proc: Compressing... '
          spawn, 'gzip -f '+outfil
          ;; Release
          delvarx, fin_img, var_nrm, img, var, img_nrm
          
      endfor
  
      ;; DEL OV
      if not keyword_set( SVOV ) then imacsls_delov, imacsls, indx
  endfor

  print, 'imacsls_proc: All done!'

  return
end


