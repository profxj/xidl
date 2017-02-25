;+ 
; NAME:
; hires_tstrdx   
;     Version 1.1
;
; PURPOSE:
;  Routine used to check the HIRES pipeline.  For internal use only.
;
; CALLING SEQUENCE:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Feb-2006 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_tstrdx 

;
;  if  N_params() LT 1  then begin 
;      print,'Syntax - ' + $
;        'hires_tstrdx, afil, ffil, ofil, CHIP=, SETUP= [v1.1]'
;      return
;  endif 

  if not keyword_set(CHIP) then chip = [1L,2L,3L]

  ;; Create structure
  hires_strct, hires
  nhires = n_elements(hires)
  

  ;; Setup
  hires_setup, hires, XTOLER=0.0015

  ;; Setup
  setup = 1

  ;; Flats
  hires_allflat, hires, setup, CHIP

  ;; Arcs
  hires_allarc, hires, setup, CHIP, _EXTRA=extra

  ;; Slitflat
  hires_slitflat, hires, setup, CHIP, _EXTRA=extra

  ;; Object
  obj_id = 1L
  for qq=0L,n_elements(CHIP)-1 do begin
      ;; Process
      hires_proc, hires, SETUP=setup, OBJ=obj_id, CHIP=chip[qq]
      ;; Trace
      hires_fntobj, hires, setup, obj_id, chip[qq], CHK=chk
      ;; Skysub
      hires_skysub, hires, setup, obj_id, chip[qq], CHK=chk
      ;; Extract
      hires_extract, hires, setup, obj_id, chip[qq], CHK=chk, _EXTRA=EXTRA
      ;; Combine
      hires_combspec, hires, setup, obj_id, chip[qq], /noflux
  endfor


  return
end

