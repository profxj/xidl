;+ 
; NAME:
; mike_rslvall  
;     Version 2.0
;
; PURPOSE:
;   Resolves a number of key codes for MIKE redux.  The main problem
;   is that IDL will not find functions 'hidden' inside of main
;   routines.  This routine compiles the majority of these functions.
;
; CALLING SEQUENCE:
;  mike_rslvall
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
;   mike_rslvall
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_rslvall

  compile_opt strictarr
  ;; Subbias
  resolve_routine, 'mike_subbias', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_proc', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Arcs
  resolve_routine, 'mike_procarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_fitarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_fit2darc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_tracearc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_fittrcarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_mkaimg', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_arcxyoff', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_allarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'mike_arcalign', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Slit flat
  resolve_routine, 'slitflat_qa', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'frac_order', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_trcflat', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_fntobj', /no_recompile, /either, /COMPILE_FULL_FILE

  ;; Correct double-precision routine for bsplines
  ;cd, getenv('MIKE_DIR')+'/pro/Object/', current=back
  ;resolve_routine, 'bsplvn', /either
  ;cd,  back


end
