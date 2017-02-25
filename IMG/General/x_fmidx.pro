;+ 
; NAME:
; x_fmidx
;    Version 1.1
;
; PURPOSE:
;    Converts frame number to index -- Useful for ccd###
;
; CALLING SEQUENCE:
;   
; indx = x_fmidx(struct, frameno)
;
; INPUTS:
;   frameno - Frameno [array] of the CCD image
;
; RETURNS:
;   indx -- Index of the frame
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
;   idx = x_fmidx(n1, 203)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_fmidx, struct, frame

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'idx = x_fmidx, struct, frame, [v1.1]'
      return, -1
  endif 

  ;
  nfrm = n_elements(frame)
  rval = lonarr(nfrm)

  for i=0L,nfrm-1 do begin
      idx = where(struct.frame EQ frame[i], count)
      case count of
          1: rval[i] = idx[0]
          else: rval[i] = -1
      endcase
  endfor
  ; Return
  if nfrm EQ 1 then return, rval[0] else return, rval
end
  
      
