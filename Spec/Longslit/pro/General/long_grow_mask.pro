;+
; NAME:
;  long_combspec 
;  Version 1.1
;
; PURPOSE:
;  Combine the spectra from multiple exposures taken through the same
;  mask or longslit.
;
; CALLING SEQUENCE:
;  newmask = LONG_GROW_MASK(mask, ngrow)
;
; INPUTS:
;   mask       -- input mask
;   ngrow      -- amount to grow mask by
;
; COMMENTS:

; EXAMPLES:
;
; PROCEDURES CALLED:
;  
;
; REVISION HISTORY:
;   15-October-2007  Created by JFH
;-
;------------------------------------------------------------------------------


FUNCTION LONG_GROW_MASK, mask, ngrow
;; convention here is 1=good, 0 =bad
nmask = n_elements(mask)
newmask = mask
bad = WHERE(mask EQ 0, nbad)
IF nbad GT 0 THEN BEGIN
    FOR kk = 0L, nbad-1L DO BEGIN
        ind = bad[kk]
        IF mask[ind] EQ 0 THEN $
          newmask[(ind-ngrow) > 0:(ind+ngrow) < (nmask-1L)] = 0
    ENDFOR
ENDIF

RETURN, newmask
END

