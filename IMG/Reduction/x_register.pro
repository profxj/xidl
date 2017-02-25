;+ 
; NAME:
; x_register   
;   Version 1.1
;
; PURPOSE:
;    Given a set of offsets between images, this routine will
;   create final images that are registered.
;
; CALLING SEQUENCE:
;   
;   x_register, img, offsets, GOODREG=, OUTPTH=
;
; INPUTS:
;   img -- String array of image names (assumes _sig extensions for sig
;          files)
;   offsets -- Integer offsets between images
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   goodreg -- Region of each image (4 int array) defining good data
;               Formatted (x0,x1,y0,y1)
;   OUTPTH  -- Output path for the images
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_regsiter, img, off, GOODREG=good
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_register, img, offsets, GOODREG=goodreg, OUTPTH=outpth

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'x_register, img, offsets, GOODREG=, OUTPTH= (v1.0) ' 
      return
  endif 
  
  nimg = n_elements(img)

;  Check offsets

  sz = size(offsets, /dimensions)
  if sz[0] NE 2 OR sz[1] NE nimg then begin
      print, 'offsets has the wrong dimension!', sz
      return
  endif

;  GOODREG


; OUTPTH

  if not keyword_set( OUTPTH ) then outpth = './'

;  Determine the size of the Final Image
      
  x1 = -10L
  x2 = 100000
  y1 = -10L
  y2 = 100000
  for i=0,nimg-1 do begin
      x1 = x1 > (goodreg[0,i] + offsets[0,i])
      x2 = x2 < (goodreg[1,i] + offsets[0,i])
      y1 = y1 > (goodreg[2,i] + offsets[1,i])
      y2 = y2 < (goodreg[3,i] + offsets[1,i])
  endfor

  ;; GOODREG
  if not keyword_set( GOODREG ) then begin
      goodreg = lonarr(4,nimg)
      flg_good = 1
  endif else flg_good = 0

; Read and copy!

  for i=0,nimg-1 do begin

      ; Read Flux img

      din = xmrdfits(img[i],0,h1, /fscale)

      ;; GOODREG
      if flg_good then begin
          for i=0,nimg-1 do begin
              sz = size(din, /dimensions)
              goodreg[1,i] = sz[0]-1
              goodreg[3,i] = sz[1]-1
          endfor
      endif

      dout = din[x1-offsets[0,i]:x2-offsets[0,i], $
                 y1-offsets[1,i]:y2-offsets[1,i]] 

      ; Write img

      fil = fileandpath(img[i])
      lfil = strlen(fil)
      Xfil = strjoin([outpth, $
                      strmid(fil,0,lfil-6),'X',$
                      strmid(fil,lfil-6,1),'.fits'])
      ; Update naxis
      sxaddpar, h1, 'NAXIS1', x2-x1+1
      sxaddpar, h1, 'NAXIS2', y2-y1+1
      writefits, Xfil, dout, h1

      if not keyword_set( NOSIG ) then begin
          limg = strlen(img[i])
          sigfil = strjoin([strmid(img[i],0,limg-5), $
                            '_sig.fits'])
          din = xmrdfits(sigfil,0,h1, /fscale)
          dout = din[x1-offsets[0,i]-1:x2-offsets[0,i]-1, $
                     y1-offsets[1,i]-1:y2-offsets[1,i]-1] 
          Sfil = strjoin([outpth, $
                      strmid(fil,0,lfil-6),'S',$
                      strmid(fil,lfil-6,1),'.fits'])
          ; Update naxis
          sxaddpar, h1, 'NAXIS1', x2-x1+1
          sxaddpar, h1, 'NAXIS2', y2-y1+1
          writefits, Sfil, dout, h1
      endif

  endfor
  
  delvarx, din, dout, h1
      
end
