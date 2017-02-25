;+ 
; NAME:
; x_chkwindow
;   Version 1.1
;
; PURPOSE:
;  Check to see if a window is open
;
; CALLING SEQUENCE:
;   
;   x_chkwindow, win_id, get_all_active=
;
; INPUTS:
;   win_id
;
; RETURNS:
;   dat       - Byt array describing open windows
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_chkwindow, win_id, get_all_active = gaa


;
;  if  N_params() LT 1  then begin 
;    print,'Syntax - ' + $
;             'x_psopen, psfile, /MAXS, _EXTRA= [V1.1]'
;    return
;  endif 

  ;; 
  device, window_state = ws
  if arg_present(gaa) then $
    gaa = where(ws gt 0)
  
  if n_elements(win_id) gt 0 then begin
      out = bytarr(n_elements(ws))
      out[win_id] = 1B
      return, 1B - (out - ws)[win_id]
  endif else return, -1B
  
end
