;+ 
; NAME:
; write_wfccdstr   
;  Version 1.0
;
; PURPOSE:
;    Write a wfccd_struct to an ASCII file
;
; CALLING SEQUENCE:
;   
;  write_wfccdstr, struct, ANONLY=, OUTFIL=, FITS=
;
; INPUTS:
;   struct        - A wfccd_struct
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   outfil = Output file (default is image.list)
;
; COMMENTS:
;
; EXAMPLES:
;   write_wfccdstr, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro write_wfccdstr, struct, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'write_wfccdstr, struct, ANONLY=, OUTFIL=, FITS= (v1.0)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'wfccd.list'

  nimg = n_elements(struct)
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if struct[q].mask_id GE 0L then begin
          if struct[q].mask_id LT 10 then id_str = '0'+strtrim(struct[q].mask_id,2) $
          else id_str = strtrim(struct[q].mask_id,2)
      endif else id_str = ' '
      if (not keyword_set( ANONLY ) OR struct[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a6,1x,i1,1x,a12,1x,a3,a3,1x,a3,1x,i5,f4.1,2a3,i2,f5.2,2i2,i4)',$
        q, $
        struct[q].img_root, $
        struct[q].flg_anly, $
        struct[q].Obj, $
        id_str, $
        struct[q].masknm, $
        struct[q].type, $
        long(struct[q].exp), $
        struct[q].gain, $
        struct[q].filter, $
        struct[q].grism, $
        struct[q].aperpos, $
        struct[q].AM, $
        struct[q].flg_ov, $
        struct[q].flg_final, $
        struct[q].nslits
        
      
  endfor

  close, 1

  if keyword_set( FITS ) then mwrfits, struct, fits, /create


end
