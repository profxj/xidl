;+ 
; NAME:
; x_objnumid   
;    Version 1.0
;
; PURPOSE:
;    Given the slitstr and the map, find slit positions in the
;    original image
;
; CALLING SEQUENCE:
;   
;   x_objnumid, slitstr, map
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
;   x_objnumid, num
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_objnumid, num


;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'string = x_objnumid(num) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

;  Simple copy

  case num of 
      0: return, 'a'
      1: return, 'b'
      2: return, 'c'
      3: return, 'd'
      4: return, 'e'
      5: return, 'f'
      6: return, 'g'
      7: return, 'h'
      8: return, 'i'
      9: return, 'i'
     10: return, 'i'
     11: return, 'j'
     12: return, 'k'
      else: stop
  endcase

  return, -1
end
