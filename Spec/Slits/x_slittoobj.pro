;+ 
; NAME:
; x_slittoobj   
;    Version 1.0
;
; PURPOSE:
;    Given the slitstr and the map, find slit positions in the
;    original image
;
; CALLING SEQUENCE:
;   
;   x_slittoobj, slitstr, map
;
; INPUTS:
;   img         - Flux image
;   slitstr     - Slit structure
;   map         - y-distortion map (fits is ok)
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
             'x_slittoobj, objstr, nobj, slitstr, clm [v1.0]'
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
