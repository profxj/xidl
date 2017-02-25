;+ 
; NAME:
; xdimg_darksub
;        Version 1.1
;
; PURPOSE:
;    Dark subtracts an image.
;
; CALLING SEQUENCE:
;   xdimg_darksub, struct, imgarr, /SVSTRCT
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;   imgarr   -- integer array defining the images to be process
;
; RETURNS:
;
; OUTPUTS:
;   imgarr - fits files in the dir Coadd_divided
;
; OPTIONAL KEYWORDS:
;   /svstrct  -- Save the updated structure.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_dividecoadds, nght1_strct, imgarr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-June-2007 Written by LKP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_darksub, struct, imgarr, SVSTRCT=svstrct

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_darksub, struct, imgarr, /SVSTRCT (v1.1)'
      return
  endif 

  ; find the superdark frame needed for each image.
  xdimg_finddark, struct, imgarr

  for i=0,n_elements(imgarr)-1 do begin

     img = struct[imgarr[i]].rootpth + struct[imgarr[i]].img_root
     img = mrdfits(strtrim(img,2), 0, header, /silent)
     
     if struct[imgarr[i]].flg_drk ne 1 then begin
        print, 'The necessary dark file has not been created.  Skipping.'
        continue
     endif

     darkimg = readfits(strtrim(struct[imgarr[i]].img_drk,2), /silent)
     img = img - darkimg
     statsec = struct[imgarr[i]].statsec
     med_drksub = median(img[statsec[0]:statsec[1],statsec[2]:statsec[3]])  ;take median in just a small region
     struct[imgarr[i]].flg_drksub = 1
     struct[imgarr[i]].img_drksub = 'Darks/Sub/drksub_'+strtrim(struct[imgarr[i]].img_root,2)
     struct[imgarr[i]].med_drksub = med_drksub

     ; Output
     outfil = strtrim(struct[imgarr[i]].img_drksub,2)
     sxaddpar, header, 'DRKSUB', strtrim(struct[imgarr[i]].img_drk,2)
     mwrfits, img, outfil, header

  endfor

;  Write the structure
  if keyword_set(SVSTRCT) then write_dimgstr, struct, FITS='struct.fits' 

end
