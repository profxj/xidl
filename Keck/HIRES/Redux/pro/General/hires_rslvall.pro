;+ 
; NAME:
; hires_rslvall  
;     Version 1.1
;
; PURPOSE:
;   Resolves a number of key codes for HIRES redux
;
; CALLING SEQUENCE:
;   
;  hires_rslvall
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
;   hires_rslvall
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_rslvall

  compile_opt strictarr
  ;; Subbias
  resolve_routine, 'hires_subbias', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_proc', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Arcs
  resolve_routine, 'hires_procarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_fitarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_fit2darc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_tracearc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_fittrcarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'hires_mkaimg', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'hires_arcxyoff', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'hires_allarc', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'hires_arcalign', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Flats
  resolve_routine, 'x_trcflat', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_fntobj', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_ordermask', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'slitflat_qa', /no_recompile, /either, /COMPILE_FULL_FILE

  ;; Other
  resolve_routine, 'hires_slitlen', /no_recompile, /either, /COMPILE_FULL_FILE

end
