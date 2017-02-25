;+ 
; NAME:
; write_dimgstr   
;    Version 1.1
;
; PURPOSE:
;    Write a dimg_struct (direct image) to an ASCII file and/or FITS file
;
; CALLING SEQUENCE:
;   
;  write_dimgstr, struct, ANONLY=anonly, OUTFIL=outfil, FITS=fits
;
; INPUTS:
;   struct        - A dimg_struct
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL  - Output file for ASCII [default is image.list]
;
; OPTIONAL KEYWORDS:
;   /ANONLY - Only print files with flg_anly NE 0   
;   FITS=   - Fits file output for the structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   write_dimgstr, nght1_strct
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-July-2001 Written by JXP
;   29-Dec-2001 Added fits option
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro write_dimgstr, struct, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'write_dimgstr, struct, ANONLY=, OUTFIL=, FITS= (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'image.list'

  nimg = n_elements(struct)
  
  close, /all
  openw, 1, 'image.list'

  for q=0,nimg-1 do begin
      
      if (not keyword_set( ANONLY ) OR struct[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a12,1x,i1,1x,a12,1x,a3,1x,f7.1,1x,a1,1x,f6.3,3i6,4i2)',$
        q, $
        struct[q].img_root, $
        struct[q].flg_anly, $
        struct[q].Obj, $
        struct[q].type, $
        struct[q].exp, $
        struct[q].filter, $
        struct[q].AM, $
        struct[q].min_raw, $
        struct[q].med_raw, $
        struct[q].max_raw, $
        struct[q].flg_ov, $
        struct[q].flg_msk, $
        struct[q].flg_skymsk, $
        struct[q].flg_final
      
  endfor

  close, 1

; FITS

  if keyword_set( FITS ) then mwrfits, struct, fits, /create


end
