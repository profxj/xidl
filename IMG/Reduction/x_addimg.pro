;+ 
; NAME:
; x_addimg
;   Version 1.1
;
; PURPOSE:
;    Combines a set of images.  I think this routine is entirely
;  superseded by xcombine.
;
; CALLING SEQUENCE:
;   
;   x_addimg, img, comb, var, SCALE=, WEIGHT=, IMGINDX=, VARINDX=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   comb -- Combined image
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_addimg, all_img, comb_img
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   10-Dec-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_addimg, img, comb, var, SCALE=scale, MMEM=mmem, WEIGHT=weight, $
              IMGINDX=imgindx, VARINDX=varindx, SIGREJ=sigrej


  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'x_addimg, img, comb, [var], SCALE=, MMEM=, WEIGHT=, IMGINDX='
      print, '                            VARINDX= (v1.1)'
      return
  endif 
  
  nimg = n_elements(img)

;  Optional Keywords

  if not keyword_set( IMGINDX ) then imgindx = 0L
  if not keyword_set( VARINDX ) then varindx = 1L
  if not keyword_set( MMEM ) then mmem = 200
  if not keyword_set( WEIGHT ) then weight = replicate(1., nimg)
  if not keyword_set( SCALE ) then scale = replicate(1., nimg)
  if not keyword_set( SIGREJ ) then sigrej = 5.

;  Investigate the images


  h0 = xheadfits(img[0])
  n1_0 = sxpar(h0, 'NAXIS1')
  n2_0 = sxpar(h0, 'NAXIS2')
;  bit_0 = sxpar(h0, 'BITPIX')
  
  for i=1L,nimg-1 do begin
      h1 = xheadfits(img[i])
      n1 = sxpar(h1, 'NAXIS1')
      n2 = sxpar(h1, 'NAXIS2')

      if( n1 NE n1_0 OR n2 NE n2_0) then begin
          print, 'x_addimg: Images arent all the same dimension or type'
          return
      endif
  endfor

;  Output image
  comb = fltarr(n1_0, n2_0)
  if arg_present(var) then var = fltarr(n1_0, n2_0)

  ;  Restrict on memory
  dsiz = 4

  nmeg = (n1_0 * n2_0)*dsiz / 1024.^2   ; Number of Meg per image

  maximg =  (mmem-30) / nmeg - 2   ; Subtract 2 for output
  if maximg LE 2 then begin
      print, 'x_addimg: Get some more memory and come back later!', maximg
      return
  endif

  ; Max number of rows at a time
  nrow = fix(n2_0*float(maximg)/nimg/2.) < n2_0 ; Number of rows to read in


; Status
  print, 'Combining images: '
  for i=0,nimg-1 do print, img[i]

                                ; Loop!
  srow = 0L                     ; Starting row
  lrow = nrow-1L                ; Ending row
  while (lrow LE n2_0-1) do begin

; Status
      print, 'x_addimg: Reading rows: ', $
        strtrim(srow,2)+'-'+strtrim(lrow,2), $
	' of ', strtrim(n2_0-1,2)

      ; Create the arrays
      imgarr = fltarr(n1_0,lrow-srow+1,nimg)
      vararr = fltarr(n1_0,lrow-srow+1,nimg)

      ; Read in the arrays
      for i=0,nimg-1 do begin
          imgarr[*,*,i] = $
            xmrdfits(img[i],IMGINDX,range=[srow,lrow], /silent, /fscale)
          vararr[*,*,i] = $
            xmrdfits(img[i],VARINDX,range=[srow,lrow], /silent, /fscale)
      endfor

      ;; Scale
      if keyword_set( SCALE ) then begin
          print, 'x_addimg: Scaling'
          for i=0L,nimg-1 do imgarr[*,*,i] = imgarr[*,*,i] * scale[i]
      endif
      
      ;; CALL EXTERNAL
      print, 'x_addimg: Coadding!'
      dimvec = size(imgarr, /dimensions)
      ndim = N_elements(dimvec)
      dim = size(imgarr, /n_dim)
      ;; Allocate memory for the output array
      if (ndim GT 1) then $
        newdimvec = dimvec[ where(lindgen(ndim)+1 NE dim) ] $
      else $
        newdimvec = [1]
      newsize = N_elements(imgarr) / dimvec[dim-1]
      medarr = reform(fltarr(newsize), newdimvec)
      finvar = reform(fltarr(newsize), newdimvec)
      dum = reform(fltarr(newsize), newdimvec)

      soname = filepath('libxmath.so', $
                        root_dir=getenv('XIDL_DIR'), subdirectory='/lib')
      retval = call_external(soname, 'arrmedsclipwgt', $
                             ndim, dimvec, imgarr, vararr, long(dim), $
                             float(sigrej), medarr, finvar, float(weight), $
                             float(scale), dum, /UNLOAD)
      ;; THE PASS
      comb[*,srow:lrow] = temporary(medarr)
      if arg_present(var) then $
        var[*,srow:lrow] = temporary(finvar)
      
      ;; Release memory
      delvarx, imgarr, velarr
      
      ;; Adjust the Rows
      if(lrow EQ n2_0-1) then break
      srow = lrow+1
      lrow = lrow+nrow < (n2_0-1)
  endwhile
          

; Output
  print, 'x_addimg: All done!'

end
