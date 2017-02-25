;+ 
; NAME:
; mslitstrct__define
;   Version 1.1
;
; PURPOSE:
;  This structure is used to save info about a slit in a mask.  In
;  particular, one saves the edges of the slit.
;
; CALLING SEQUENCE:
;   tmp = {mslitstrct}
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;  Written by JXP
;-
;----------------------------------------------------------------------
pro thrustrct__define

; Version 1.0

  tmp = {thrustrct, $
         observatory: ' ', $
         instrument: ' ',        $      ; ID number
         date: ' ',   $         ; DDMMYYYY (e.g. 20JAN2004)
         grating: ' ',   $      ; Field name
         gratepos: 0,   $      ; Field name
         central_wave: 0., $    ; Central wavelength of grating
         blocking: '', $        ; Blocking filter
         airmass: 0., $
         spec_bin: 0, $
         slit_width: 0., $
         std_name: '', $
         ra: ' ',    $              ; RA (J2000)
         dec: ' ',    $             ; DEC (J2000)
         conditions: ' ',   $       ; x offset from the QSO (arcsec)
         filename: '', $
         frame: 0L $
         }

end
  
         
