;+
; NAME:
;   verbset
;
; PURPOSE:
;   set the verbosity level used by vprint
;
; CALLING SEQUENCE:
;   verbset, verblevel   ; to set verbosity level
;   verbset              ; to query verbosity level
; 
; INPUTS:
;   verblevel - verbosity level (usually an integer, 0..5 or so)
;
; COMMON BLOCKS:
;   verbosity_level
;
; EXAMPLES:
;   verbset, 1
;
; COMMENTS:
;
; REVISION HISTORY:
;
;       Sat Feb 23 17:46:38 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
pro verbset, verblevel_in

  common verbosity_level, verblevel
  
  if n_params() EQ 0 then begin 
     if n_elements(verblevel) EQ 0 then begin 
        print, 'VERBSET:  verbosity is not set; set with verbset, vlevel'
     endif else begin 
          
        print, 'VERBSET:  verbosity is set to:', verblevel, $
          format='(A,I3)'
     endelse 
     
     return
  endif 

  verblevel_in = fix(verblevel_in)

  if verblevel_in LT 0 or verblevel_in GT 100 then message, 'be reasonable'

  verblevel = verblevel_in
  vprint, 1, 'Setting verbosity level to:', verblevel

  return
end
