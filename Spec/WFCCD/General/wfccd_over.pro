;+ 
; NAME:
; wfccd_over   
;        Version 1.0
;
; PURPOSE:
;    Overscan subtracts a list of images of a structure
;
; CALLING SEQUENCE:
;   
;   wfccd_over, struct
;
; INPUTS:
;   struct -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   erase - Erase pre-existing ov files
;   inter - Interactive OV fitting
;
; OPTIONAL OUTPUTS:
;   ovimg - fits files in the dir OV
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_over, nght1_strct, ovimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_over, struct, ovimg, ORDR=ordr, ERASE=erase, INTER=inter, $
                OVSEC=ovsec, IMG=img, NOFITS=nofits

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_over, struct, ovimg, ORDR=, /ERASE, /INTER (v1.0)'
      return
  endif 
  
;  Optional Keywords

  ; xoverscan automatically takes ovsec = '[2050:2079,*]'

;  Call xoverscan

  ovfil = strtrim(struct[ovimg].rootpth,2) + $
    strtrim(struct[ovimg].img_root,2)
  xoverscan, ovfil, struct[ovimg[0]].ccd, ORDR=ordr, OVSEC=ovsec, $
    ERASE=erase, INTER=inter, OVIMG=img, NOFITS=nofits
      

;  Set ov names, flags

  if not keyword_set( NOFITS ) then begin
      struct[ovimg].flg_ov = 1
      struct[ovimg].img_ov = 'OV/ov_'+strtrim(struct[ovimg].img_root,2)
  endif


end
