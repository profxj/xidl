;+ 
; NAME:
; xdimg_bias   
;    Version 1.1
;
; PURPOSE:
;    Creates Bias frame given structure
;
; CALLING SEQUENCE:
;   
;   xdimg_bias, struct, /SVOV, OUTFIL=
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;    A Fits file in the dir BIAS named 'Bias.fits'
;
; OPTIONAL KEYWORDS:
;   SVOV - Save ov files (default is to delete them)
;   OUTFIL - Name of Bias file (default is 'Bias/Bias.fits')
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_bias, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XDIMG_OVER
;  XCOMBINE
;  MWRFITS
;  XDIMG_DELOV
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_bias, struct, SVOV=svov, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_bias, struct, /SVOV, OUTFIL= (v1.1)'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OUTFIL ) then begin
      ; Create directory if necessary
      a = findfile('Bias/..', count=count)
      if count EQ 0 then file_mkdir, 'Bias'
      outfil = 'Bias/Bias.fits'
  endif
  
;  Find the Bias frames

  bias = where(struct.type EQ 'ZRO' AND struct.flg_anly NE 0, nbias)

;  Overscan

  ovflt = where(struct[bias].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, bias[ovflt], ORDR=4


; Status
  print, 'Combining images: '
  for i=0,nbias-1 do print, struct[bias[i]].img_ov
  print, '             into ', outfil

; Combine Images
  xcombine, struct[bias].img_ov, comb, head

; Output
  mwrfits, comb, outfil, head, /create

; Delete the ov files

  if not keyword_set( SVOV ) then xdimg_delov, struct, bias

  print, 'All done with Bias frame!'

end
