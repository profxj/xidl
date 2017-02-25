;+ 
; NAME:
; x_slittoobj   
;    Version 1.0
;
; PURPOSE:
;    Simple program to move some tags from the slit structure into the
;    object structure
;
; CALLING SEQUENCE:
;   x_slittoobj, objstr, nobj, slitstr, clm
;
; INPUTS:
;   objstr      - Object sturcture
;   nobj        - Index of the object
;   slitstr     - Slit structure
;   clm         - Column where the object was identified
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_slittoobj, objstr, nobj, slitstr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_slittoobj, objstr, nobj, slitstr, clm


;  Error catching
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_slittoobj, objstr, nobj, slitstr, clm [v1.1]'
    return
  endif 


;  Optional Keywords

;  Simple copy

  objstr[nobj].field = slitstr.field
  objstr[nobj].slit_id = slitstr.id
  objstr[nobj].xcen = clm

;
  return
end
