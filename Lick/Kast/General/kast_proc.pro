;+ 
; NAME:
; kast_proc   
;     Version 1.0
;
; PURPOSE:
;    Process an image (bias subtract + flatten)
;      WARNING! Assumes 1 bias and 1 flat for all images
;
; CALLING SEQUENCE:
;   
;  kast_proc, kast, setup, obj_id
;
; INPUTS:
;   kast   -  ESI structure
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
;   kast_proc, kast 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_proc, kast, setup, side, obj_id, CLOBBER=clobber, INDEX=index

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_proc, kast, setup, side, obj_id, INDEX=, /CLOBBER [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0

; Index
  if not keyword_set( INDEX ) then begin
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.obj_id EQ obj_id AND kast.type EQ 'OBJ' AND $
                   kast.side EQ side, nindx)
      if nindx EQ 0 then begin
          print, 'kast_proc: No images to process!', obj_id
          return
      endif
  endif else begin
      indx = index
      nindx = n_elements(indx)
  endelse

  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      for qq=0L,nindx-1 do begin
          ;; Chk for outfil
          outfil = 'Final/f_'+kast[indx[qq]].img_root
          a = findfile(outfil+'*', count=na)
          if na NE 0 then begin
              print, 'kast_proc: File ', outfil, ' exists, continuing...'
          endif else flg_proc = 1
      endfor
  endif else flg_proc = 1
  
  if flg_proc EQ 0 then begin
      print, 'kast_proc: No files to process. Returning!'
      return
  endif

;  Optional Keywords
  if not keyword_set( FLAT ) then flat = kast[indx[0]].flat_fil
;  if not keyword_set( BIAS ) then bias = 'Bias/BiasS.fits'

;  Bias Image
;  print, 'kast_proc: Reading Bias ', bias
;  a = findfile(bias+'*', count=na)
;  if NA EQ 0 then begin
;      print, 'kast_proc: Bias not found! Stopping... '
;      stop
;  endif
;  img_bias = xmrdfits(bias, range=[230,3529], /silent)
;  sz_bias = size(img_bias, /dimensions)

;  Flat file
  print, 'kast_proc: Reading Flat ', flat
  a = findfile(flat+'*', count=na)
  if NA EQ 0 then begin
      print, 'kast_proc: Flat not found! Stopping... '
      stop
  endif
  img_flat = xmrdfits(flat, /silent)
  gdflt = where(img_flat GT 0.)
  sz_flt = size(img_flat, /dimensions)

; Loop on Images
  for q=0L,nindx-1 do begin
      ;; Chk for outfil
      outfil = 'Final/f_'+kast[indx[q]].img_root
      a = findfile(outfil+'*', count=na)
      if na NE 0 and not keyword_set( CLOBBER ) then begin
          print, 'kast_proc: File ', outfil, ' exists, continuing...'
          continue
      endif
      ;; Input image
      fil = kast[indx[q]].rootpth+kast[indx[q]].img_root
      print, 'kast_proc: Processing ', fil
      raw = xmrdfits(fil, 0, head, /silent, /fscale)
      ;; Bias subtract
;      raw = raw[25:2072,*] - img_bias

      ;; Create VAR image
      var = ((raw>3.)*kast[indx[q]].gain + kast[indx[q]].readno^2) > 0.

      ;; Flatten
      fin_img = x_normspec(raw, img_flat, var, nrmvar, PIX=gdflt)

      ;; Add in gain
      fin_img = fin_img*kast[indx[q]].gain

      ;; Zero out bad regions
      case kast[indx[q]].side of
          1: begin
              nrmvar[*,0:5] = -1. 
              nrmvar[*,170:*]  = -1.
              fin_img[*,0:5] = 0. 
              fin_img[*,170:*]  = 0.
          end
          2: begin
              nrmvar[*,0:20] = -1. 
              nrmvar[*,190:*]  = -1.
              fin_img[*,0:20] = 0. 
              fin_img[*,190:*]  = 0.
          end
          else: stop
      endcase
          
      ;; Output
;      sxaddpar, head, 'BIAS', 'T', bias
      sxaddpar, head, 'FLAT', 'T', flat
      kast[indx[q]].img_final = outfil
      mwrfits, fin_img, outfil, head, /create, /silent
      mwrfits, nrmvar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'kast_proc: Created ', outfil
  endfor

  return
end
      
