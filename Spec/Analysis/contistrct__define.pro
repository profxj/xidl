;+ 
; NAME:
; xsp1dstrct__define
;    Version 1.0
;
; PURPOSE:
;  This routine defines the structure for individual object spectra
;
; CALLING SEQUENCE:
;   tmp = {xsp1dstrct}
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
;   02-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------
pro contistrct__define

;  This routine defines the structure for individual object spectra

  tmp = { contistrct, $
          npts: 0L, $
          xval: fltarr(100), $
          yval: fltarr(100), $
          msk: lonarr(100) }

   return
end
