;+ 
; NAME:
; mosa_regsig   
;   Version 1.1
;
; PURPOSE:
;    Given a list of images and the offsets, creates a new set of
;  registered images. 
;
; CALLING SEQUENCE:
;   mosa_regsig, img, offsets, GOODREG=, OUTPTH=, GAIN=
;
; INPUTS:
;   img -- Array of image names (assumes _sig extensions for sig
;          files)
;   offsets -- Integer offsets for the images
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   goodreg -- Region of each image (4 int array) defining good data
;               Formatted (x0,x1,y0,y1)
;   GAIN=   -- Image gain [default: 2.75]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mosa_regsig, img, off, GOODREG=good
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Jan-2004 Written by JXP (revised x_register)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mosa_regsig, img, offsets, GOODREG=goodreg, OUTPTH=outpth, GAIN=gain

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'x_register, img, offsets, GOODREG=, OUTPTH= (v1.0) ' 
      return
  endif 
  
  if not keyword_set( GAIN ) then gain = 2.75

  nimg = n_elements(img)

;  Check offsets

  sz = size(offsets, /dimensions)
  if sz[0] NE 2 OR sz[1] NE nimg then begin
      print, 'offsets has the wrong dimension!', sz
      return
  endif

;  GOODREG

  if not keyword_set( GOODREG ) then begin
      goodreg = intarr(4,nimg)
      for i=0,nimg-1 do begin
          sz = size(img[i], /dimensions)
          goodreg[1,i] = sz[0]-1
          goodreg[3,i] = sz[1]-1
      endfor
  endif

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

; Read and copy!

  for i=0,nimg-1 do begin

      ;; Read Flux img
      din = xmrdfits(img[i],0,h1, /fscale)
      dout = din[x1-offsets[0,i]:x2-offsets[0,i], $
                 y1-offsets[1,i]:y2-offsets[1,i]] 
      delvarx, din

      ;; Update gain
      sxaddpar, h1, 'GAIN', gain
      texp = sxpar(h1, 'EXPTIME')

      ;; Write img
      fil = fileandpath(img[i])
      lfil = strlen(fil)
      Xfil = strjoin([outpth, $
                      strmid(fil,0,lfil-6),'X',$
                      strmid(fil,lfil-6,1),'.fits'])
      ;; Update naxis
      sxaddpar, h1, 'NAXIS1', x2-x1+1
      sxaddpar, h1, 'NAXIS2', y2-y1+1
      sxaddpar, h1, 'TTIME', 1.

      ;; Update crpix1, crpix2
      crpix1 = sxpar(h1, 'CRPIX1') - x1 + offsets[0,i]
      crpix2 = sxpar(h1, 'CRPIX2') - y1 + offsets[1,i]
      if i EQ 0 then begin
          sv1 = crpix1
          sv2 = crpix2
      endif else begin
          if crpix1 NE sv1 then stop
          if crpix2 NE sv2 then stop
      endelse
      sxaddpar, h1, 'CRPIX1', crpix1
      sxaddpar, h1, 'CRPIX2', crpix2

      mwrfits, dout/texp, Xfil, h1, /create

      ;; Create sig image
      simg = sqrt(temporary(dout)/gain/texp)

      ;; Open exp image
      expfil = strmid(fil,0,lfil-5)+'exp.fits'
      dexp = xmrdfits(expfil, /fscale)
      dexp = dexp[x1-offsets[0,i]:x2-offsets[0,i], $
                 y1-offsets[1,i]:y2-offsets[1,i]] 

      ;; Finish sig array
      simg = simg/sqrt(temporary(dexp))
      a = where(simg LE 0.,na)
      if na NE 0 then simg[a] = 9999.

      ;; Update naxis and output
      Sfil = strjoin([outpth, $
                      strmid(fil,0,lfil-6),'S',$
                      strmid(fil,lfil-6,1),'.fits'])
      mwrfits, temporary(simg), Sfil, h1, /create

  endfor
  
end
