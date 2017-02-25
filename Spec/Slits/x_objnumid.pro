;+ 
; NAME:
; x_objnumid   
;    Version 1.1
;
; PURPOSE:
;    Turns an integer into a letter for Obj identification
;
; CALLING SEQUENCE:
;   
;   idval = x_objnumid( num )
;
; INPUTS:
;  num - An integer
;
; RETURNS:
;  idval -- A letter (0='a', 1='b', etc.)
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
;   idval = x_objnumid( 1 )
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
  val = strarr(n_elements(num))
  for qq=0L, n_elements(num)-1 do begin
      case num[qq] of 
          0: val[qq]= 'a'
          1: val[qq]= 'b'
          2: val[qq]= 'c'
          3: val[qq]= 'd'
          4: val[qq]= 'e'
          5: val[qq]= 'f'
          6: val[qq]= 'g'
          7: val[qq]= 'h'
          8: val[qq]= 'i'
          9: val[qq]= 'j'
          10: val[qq]= 'k'
          11: val[qq]= 'l'
          12: val[qq]= 'm'
          13: val[qq]= 'n'
          14: val[qq]= 'o'
          15: val[qq]= 'p'
          16: val[qq]= 'q'
          17: val[qq]= 'r'
          18: val[qq]= 's'
          19: val[qq]= 't'
          20: val[qq]= 'u'
          else: stop
      endcase
  endfor
  if n_elements(num) EQ 1 then val = val[0]

  return, val
end
