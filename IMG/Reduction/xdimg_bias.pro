;+ 
; NAME:
; xdimg_bias   
;    Version 1.1
;
; PURPOSE:
;    Creates Bias frame given direct image structure (which must
;   contain a number of bias files!)
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
;   OUTFIL - Name of Bias file [default is 'Bias/Bias.fits']
;
; OPTIONAL KEYWORDS:
;   SVOV - Save ov files (default is to delete them)
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
; Note: if there's only one bias frame, i think this code will crash!!
; nbias is never checked before xcombine.

;  Overscan

  ovflt = where(struct[bias].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, bias[ovflt], ORDR=4

; Print a warning message if there's only one bias frame.
  if nbias eq 1 then print, 'WARNING: You dont want to run xcombine on only one image.'

; Status
  print, 'Combining images: '
  for i=0,nbias-1 do print, struct[bias[i]].img_ov
  print, '             into ', outfil


; Combine images only if there's more than one of them.
  if nbias gt 1 then begin
    ; Be Careful!
    ; The file name may be padded with spaces because everytime struct is created it does this to strings! Yikes!
     for hh=0,nbias-1 do begin
        struct[bias[hh]].img_ov=strtrim(struct[bias[hh]].img_ov,2)
     endfor
     xcombine, struct[bias].img_ov, comb, head
endif

; Output
  if nbias gt 1 then begin
     mwrfits, comb, outfil, head, /create
  endif else begin
     comb=mrdfits(strtrim(struct[bias].img_ov,2),0,head)  ;not sure if i need to get a header as well.
     mwrfits, comb, outfil, head, /create
  endelse

; Delete the ov files
  if not keyword_set( SVOV ) then xdimg_delov, struct, bias

; Resave the updated structure so that you can pick up where you left off easily.
  mwrfits, struct, 'struct.fits', /create

  print, 'All done with Bias frame!'

end
