;+ 
; NAME:
; x_chkfil   
;    Version 1.0
;
; PURPOSE:
;    Given an array of strings, return an array of unique values +
;    number
;
; CALLING SEQUENCE:
;   
; uniq = x_chkfil(strings, COUNT=count)
;
; INPUTS:
;   strings - Array of strings
;
; RETURNS:
;   uniq  - Array of unique members
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   sort - Sort the output
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   flg = x_chkfil( lbls, count=count)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_chkfil, fil, COUNT=count, SILENT=silent

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'flg = x_chkfil(fil, COUNT=) (v1.0)'
      return, -1
  endif 

  a = findfile(fil, COUNT=count)
  case count of 
      0: begin
          if not keyword_set( SILENT ) then $
            print, 'x_chkfil: '+strtrim(fil,2)+' not found!'
          return, 0
      end
      1: return, 1
      else: return, 2
  endcase

end
