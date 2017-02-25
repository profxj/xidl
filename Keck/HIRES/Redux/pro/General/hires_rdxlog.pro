;+ 
; NAME:
; hires_rdxlog   
;     Version 1.1
;
; PURPOSE:
;    Simple code to help build a log file of the pipeline
;
; CALLING SEQUENCE:
;  hires_rdxlog, logfil, mssg
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_rdxlog, logfil, mssg, OPHDR=ophdr, CLHDR=clhdr

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_rdxlog, logfil, mssg, [v1.1]'
      return
  endif 
  
  ;; File
  fil = hires_getfil(logfil)
  if keyword_set(OPHDR) then splog, '------  '+systime()+'  ---------------', $
    filename=fil, /append, /close
  splog, mssg, filename=fil, /append, /close
  if keyword_set(CLHDR) then splog, '------  '+systime()+'  ---------------', $
    filename=fil, /append, /close

  return

end
