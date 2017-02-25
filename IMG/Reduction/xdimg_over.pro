;+ 
; NAME:
; xdimg_over   
;        Version 1.1
;
; PURPOSE:
;    Overscan subtracts a list of images of a structure.  Calls
;  routine xoverscan extensively.
;
; CALLING SEQUENCE:
;   
;   xdimg_over, struct, imgs, ORDR=, /ERASE, /INTER, /NOSV
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;   imgs   -- integer array defining the images to process
;
; RETURNS:
;
; OUTPUTS:
;   ovimg - fits files in the dir OV
;
; OPTIONAL KEYWORDS:
;   /erase - Erase pre-existing ov files
;   /inter - Interactive OV fitting (uses x1dfit)
;   /nosv  -- Dont save the OV image!
;   ORDR=  -- Order to fit the overscan
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_over, nght1_strct, ovimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_over, struct, ovimg, ORDR=ordr, ERASE=erase, INTER=inter, NOSV=nosv

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_over, struct, ovimg, ORDR=, /ERASE, /INTER, /NOSV (v1.1)'
      return
  endif 
  
;  Optional Keywords

;  Call xoverscan

  ovfil = struct[ovimg].rootpth + struct[ovimg].img_root
  if not keyword_set( INTER) then begin
      if keyword_set( ERASE ) then begin
          if keyword_set( ORDR ) then xoverscan, ovfil, struct[ovimg[0]].ccd, $
            medov, ORDR=ordr, /ERASE $
          else xoverscan, ovfil, struct[ovimg[0]].ccd, medov, /ERASE
      endif else begin
          if keyword_set( ORDR ) then xoverscan, ovfil, struct[ovimg[0]].ccd, $
            medov, ORDR=ordr $
          else xoverscan, ovfil, struct[ovimg[0]].ccd, medov
      endelse
  endif else begin
      if keyword_set( ERASE ) then begin
          if keyword_set( ORDR ) then xoverscan, ovfil, struct[ovimg[0]].ccd, $
            medov, ORDR=ordr, /ERASE, /INTER $
          else xoverscan, ovfil, struct[ovimg[0]].ccd, medov, /ERASE, /INTER
      endif else begin
          if keyword_set( ORDR ) then xoverscan, ovfil, struct[ovimg[0]].ccd, $
            medov, ORDR=ordr, /INTER $
          else xoverscan, ovfil, struct[ovimg[0]].ccd, medov, /INTER
      endelse
  endelse
      

;  Set ov names, flags

  struct[ovimg].flg_ov = 1
  struct[ovimg].img_ov = 'OV/ov_'+struct[ovimg].img_root
  struct[ovimg].med_ov = medov

;  Write the structure
  if not keyword_set(NOSV) then write_dimgstr, struct, FITS='struct.fits' 

end
