;+ 
; NAME:
; mosa_fixradec   
;     Version 1.1
;
; PURPOSE:
;    Updates the header of a modified MOSA FITS file to keep the
;  WCS header in order.  The program either reads in offset values
;  in order to calculate the CRPIX values or it just reads in the
;  CRPIX values.
;
; CALLING SEQUENCE:
;  mosa_fixradec, list [offset], CRPIX=, EQUINOX=
;
; INPUTS:
;   list     -- File containing the list of images
;   [offset] -- File containing the list of offsets
;
; RETURNS:
;
; OUTPUTS:
;  Image with the updated header
;
; OPTIONAL KEYWORDS:
;  CRPIX=  -- New CRPIX1, CRPIX2 values (2-element array)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; mosa_fixradec, 'Lists/headers.list', 'Lists/cr_off.list'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Apr-2003  Written by JXP
;-
;------------------------------------------------------------------------------
pro mosa_fixradec, list, offset, CRPIX=crpix, EQUINOX=equinox

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mosa_fixradec, imglist, [offlist], CRPIX=, EQUINOX= [v1.1]'
      return
  endif 
  


  if not keyword_set( CRPIX ) then begin
      ;; Read image list
      readcol, list, oldimg, imgfil, format='a,a'
      ;; Read offsets
      readcol, offset, xoff, yoff
  endif else begin
      ;; Read image list
      readcol, list, imgfil, format='a'
  endelse


  ;; Loop
  nimg = n_elements(imgfil)

  for q=0L,nimg-1 do begin
      if not keyword_set( CRPIX ) then begin
          ;; Read old image header
          oldh = xheadfits(oldimg[q])
          
          ;; Grab old pix values
          old_cr1 = sxpar(oldh, 'CRPIX1')
          old_cr2 = sxpar(oldh, 'CRPIX2')
          
          ;; Offset
          new_cr1 = old_cr1 + xoff[q]
          new_cr2 = old_cr2 + yoff[q]
      endif else begin
          new_cr1 = crpix[0]
          new_cr2 = crpix[1]
      endelse

      ;; Open new image
      print, 'Processing Image: ', imgfil[q]
      newi = xmrdfits(imgfil[q], 0, newh,/silent)
      sxaddpar, newh, 'CRPIX1', new_cr1
      sxaddpar, newh, 'CRPIX2', new_cr2

      ;; Equinox
      if keyword_set( EQUINOX ) then sxaddpar, newh, 'EQUINOX', equinox

      ;; Write
      mwrfits, newi, imgfil[q], newh, /create

  endfor
  print, 'mosa_fixradec: All done!'
  return
end
