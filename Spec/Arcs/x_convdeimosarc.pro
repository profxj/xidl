;+
; NAME:
;   x_convdeimosarc
;    Version 1.0
;
; PURPOSE:
;  Simple program to convert the DEIMOS line list to my own format
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------
pro x_convdeimosarc, fil, outfil

   ; Need 3 parameters
   if (N_params() LT 2) then begin
      print, 'x_convdeimosarc, fil, outfil, [v1.0]'
      return
   endif

   ;; Read in fil
   readcol, fil, wv, cnt, qual, elm, FORMAT='D,I,A,A'
   nwv = n_elements(wv)

   ;; Parse?
   writecol, outfil, wv, replicate(1,nwv), elm, $
     FMT='(f10.4,1x,i1,1x,a2)'

   return
   
end
;------------------------------------------------------------------------------
