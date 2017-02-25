;+ 
; NAME:
; esi_lwdproc   
;     Version 1.0
;
; PURPOSE:
;    Process an image (bias subtract + flatten)
;      WARNING! Assumes 1 bias and 1 flat for all images
;
; CALLING SEQUENCE:
;   
;  esi_lwdproc, esi, /INTERNAL
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdproc, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdproc, esi, obj_id, FLAT=flat, BIAS=bias, CLOBBER=clobber, $
                 INDEX=index

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdproc, esi, obj_id, FLAT=, BIAS=, /CLOBBER [v1.0]'
      return
  endif 

; Index
  if not keyword_set( INDEX ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
                   esi.obj_id EQ obj_id AND esi.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_lwdproc: No images to process!', obj_id
          return
      endif
  endif else begin
      indx = obj_id
      nindx = n_elements(indx)
  endelse

  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      for qq=0L,nindx-1 do begin
          ;; Chk for outfil
          outfil = 'Final/f_'+esi[indx[qq]].img_root
          a = findfile(outfil, count=na)
          if na NE 0 then begin
              print, 'esi_lwdproc: File ', outfil, ' exists, continuing...'
          endif else flg_proc = 1
      endfor
  endif else flg_proc = 1
  
  if flg_proc EQ 0 then begin
      print, 'esi_lwdproc: No files to process. Returning!'
      return
  endif

;  Optional Keywords
  rbin = esi[indx[0]].rbin
  cbin = esi[indx[0]].cbin
  if not keyword_set( FLAT ) then flat = esi[indx[0]].flat_fil
  if not keyword_set( BIAS ) then bias = $
    esi_getfil('bias_fil', esi[indx[0]].mode, $
               cbin=cbin, rbin=rbin, /name)

;  Bias Image
  print, 'esi_lwdproc: Reading Bias ', bias
  a = findfile(bias+'*', count=na)
  if NA EQ 0 then begin
      print, 'esi_lwdproc: Bias not found! Stopping... '
      stop
  endif
  img_bias = xmrdfits(bias, range=[230,3529], /silent)
  sz_bias = size(img_bias, /dimensions)

;  Flatten
  print, 'esi_lwdproc: Reading Flat ', flat
  a = findfile(flat+'*', count=na)
  if NA EQ 0 then begin
      print, 'esi_lwdproc: Flat not found! Stopping... '
      stop
  endif
  img_flat = xmrdfits(flat, /silent)
  gdflt = where(img_flat GT 0.)
  sz_flt = size(img_flat, /dimensions)

; Loop on Images
  for q=0L,nindx-1 do begin
      ;; Chk for outfil
      outfil = 'Final/f_'+esi[indx[q]].img_root
      a = findfile(outfil+'*', count=na)
      if na NE 0 and not keyword_set( CLOBBER ) then begin
          print, 'esi_lwdproc: File ', outfil, ' exists, continuing...'
          continue
      endif
      ;; Input image
      fil = esi[indx[q]].rootpth+esi[indx[q]].img_root
      print, 'esi_lwdproc: Processing ', fil
      raw = xmrdfits(fil, 0, head, range=[230,3529], /silent, /fscale)

      ;; Bias subtract
      raw = raw[25:2072,*] - img_bias[25:2072,*]

      ;; Create VAR image
      var = ((raw>3.)*esi[indx[q]].gain + esi[indx[q]].readno^2) > 0.

      ;; Flatten
      fin_img = x_normspec(raw, img_flat, var, nrmvar, PIX=gdflt)

      ;; Add in gain
      fin_img = fin_img*esi[indx[q]].gain

      ;; Output
      sxaddpar, head, 'BIAS', 'T', bias
      sxaddpar, head, 'FLAT', 'T', flat
      esi[indx[q]].img_final = outfil
      mwrfits, fin_img, outfil, head, /create, /silent
      mwrfits, nrmvar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'esi_lwdproc: Created ', outfil
  endfor

  return
end
      
