;+ 
; NAME:
; wfccd_bias   Version 1.0
;
; PURPOSE:
;    Creates Bias frame given structure
;
; CALLING SEQUENCE:
;   
;   wfccd_bias, struct, SVOV=svov
;
; INPUTS:
;   struct -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   bias - fits file in the dir BIAS named 'Bias.fits'
;
; OPTIONAL KEYWORDS:
;   SVOV - save ov files
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_bias, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_bias, struct, SVOV=svov

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'wfccd_bias, struct, SVOV= (v1.0)'
      return
  endif 
  
;  Optional Keywords
  
;  Find the Bias frames

  bias = where(struct.type EQ 'ZRO' AND struct.flg_anly NE 0, nbias)

;  Overscan

  ovflt = where(struct[bias].flg_ov EQ 0, nov)
  if nov NE 0 then wfccd_over, struct, bias[ovflt], ORDR=4

; Create directory if necessary

  a = findfile('Bias/..', count=count)
  if count EQ 0 then file_mkdir, 'Bias'

  outfil = 'Bias/Bias.fits'

; Status
  print, 'Combining images: '
  for i=0,nbias-1 do print, struct[bias[i]].img_ov
  print, '             into ', outfil

; Combine Images
  xcombine, struct[bias].img_ov, comb, head

; Output
  writefits, outfil, comb, head

; Delete the ov files

  if not keyword_set( SVOV ) then wfccd_delov, struct, bias

  print, 'All done with Bias frame!'

end
