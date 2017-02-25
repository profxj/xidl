;+ 
; NAME:
; uves_rslvall  
;     Version 1.1
;
; PURPOSE:
;   Resolves a number of key codes for MIKE redux
;
; CALLING SEQUENCE:
;   
;  uves_rslvall
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
;   uves_rslvall
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_rslvall

  compile_opt strictarr
  ;; Subbias
  resolve_routine, 'uves_subbias', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_proc', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Arcs
  resolve_routine, 'uves_procarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_fitarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_fit2darc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_tracearc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_fittrcarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_mkaimg', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'uves_arcxyoff', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'uves_allarc', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'uves_arcalign', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Flats
  resolve_routine, 'x_trcflat', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_fntobj', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_ordermask', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'slitflat_qa', /no_recompile, /either, /COMPILE_FULL_FILE

  ;; Other
;  resolve_routine, 'uves_slitlen', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'uves_fitscatt', /no_recompile, /either, /COMPILE_FULL_FILE

end
