;+ 
; NAME:
; x_getobjnm
;    Version 1.1
;
; PURPOSE:
;  Get the index of the object name OR return the list of obj_nm.
;  This is a pretty specific routine.  If obj_nm is not specified the
;  routine launches a GUI to have the user choose.
;
; CALLING SEQUENCE:
;   indx = x_getobjnm(objstr, [obj_nm], LST=)
;
; INPUTS:
;  objstr  -- Object structure
;  [obj_nm] -- Name of object for which the index is desired
;
; RETURNS:
;   indx -  Index of the object name OR list of obj_nm
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /LST - Return the list instead of the index
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  indx = x_getobjnm(objstr)
;  list = x_getobjnm(objstr, /LST)
;
;
; PROCEDURES/FUNCTIONS CALLED:
; x_guilist
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_getobjnm, objstr, obj_nm, LST=lst

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'indx = x_getobjnm(fspec_fil, obj_nm, /list) [v1.0]'
    return, -1
  endif 

; Make list

  list = strarr(n_elements(objstr))
  for q=0L,n_elements(objstr)-1 do $
    list[q] = strtrim(objstr[q].slit_id,2)+strtrim(objstr[q].obj_id,2)

  if not keyword_set(LST) then begin
      ;; Search
      if not keyword_set( OBJ_NM ) then obj_nm = x_guilist(list, INDX=indx) $
      else indx = (where(list EQ obj_nm))[0]
      return, indx
  endif else return, list

end
