;+ 
; NAME:
; wfccd_getobjnm
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;   indx = wfccd_getobjnm(wffspec, obj_nm)
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
;   indx = wfccd_getobjnm(wffspec, obj_nm)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function wfccd_getobjnm, wffspec, obj_nm, OBJSTR=OBJSTR

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_getobjnm, fspec_fil, /OBJSTR [v1.0]'
    return, -1
  endif 

; Make list

  list = strarr(n_elements(wffspec))
  for q=0L,n_elements(wffspec)-1 do $
    list[q] = strtrim(wffspec[q].slit_id,2)+strtrim(wffspec[q].obj_id,2)

  ; Search

  if not keyword_set( OBJ_NM ) then obj_nm = x_guilist(list, INDX=indx) $
  else indx = (where(list EQ obj_nm))[0]

  return, indx
end
