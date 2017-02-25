;+
; NAME: hires_allobj   
;     Version 1.1
;
; PURPOSE:
;    Run all of the processing steps on a single object.  This
;    includes processing, CR rejection, tracing, sky subtraction and
;    extraction.  Future versions will allow multiple (or all
;    objects).  
;
;   Note -- This method is not recommended presently.
;
;
; CALLING SEQUENCE:
;    hires_allobj, hires, setup, obj_id, chip, [exp], /PROCALL
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1) or Green (2) and/or Red (3) side
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /ALLPROC - Perform all steps
;  /SKYCHK  - Check Final sky image
;  /RESCHK  - Check Final extracted image
;  /DOPROC  - Do image processing
;  /DOFNT   - Do object tracing
;  /DOSKY   - Do sky sub
;  /DOBOX   - Do extraction
;  /NOCR    - Turn off CR rejection routine (hires_objcr)
;  /CLOBBER - Overwrite previous images (mainly the processed image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_allobj, hires, 1L, 1L
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-Apr-2004 Written by JXP 
;
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_allobj, hires, setup, obj_id, chip, exp, PROCALL=PROCALL, $
                 RESCHK=reschk, SKYCHK=skychk, TRCCHK=trcchk, CHKALL=chkall, $
                 DOSKY=dosky, DOBOX=dobox, NOCR=nocr, CLOBBER=clobber, $
                 DOCOMB=docomb, DOFLUX=doflux, DO1D=do1d, _EXTRA=EXTRA, $
                 DOFNT=dofnt, DOPROC=doproc, BOXCHK=boxchk, STD=std

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_allobj, hires, setup, obj_id, chip, [exp], /PROCALL, /RESCHK'
      print, '   /SKYCHK, /TRCCHK, /DOSKY, /DOCOMB, /NOCR, /STD [v1.1]'
      return
  endif 

;  Optional Keywords

  if keyword_set( PROCALL ) then begin
      DOPROC = 1
      DOCR = 1
      DOFNT = 1
      DOSKY = 1
      DOEXT = 1
  endif

  if keyword_set( NOCR ) then docr = 0

  if keyword_set( CHKALL ) then begin
      CRCHK=1
      RESCHK=1
      SKYCHK=1
      TRCCHK=1
      COMBCHK=1
      BOXCHK=1
  endif

  if NOT keyword_set(chip) then chip = [1L,2L,3L]
  trimtype = strtrim(hires.type,2)

  if keyword_set(obj_id) then this_id = obj_id $
  else begin 
      check_id = where(hires.flg_anly NE 0 AND hires.setup EQ setup AND $
                       (trimtype EQ 'OBJ' OR trimtype EQ 'STD'))
      if check_id[0] EQ -1 then begin
          print, 'No obj_id found'
          return
      endif
      sort_id = hires[check_id[sort(hires[check_id].obj_id)]].obj_id
      this_id = sort_id[uniq(sort_id)]
  endelse

  for ii=0, n_elements(this_id)-1 do begin

    for ss=0L,n_elements(chip)-1 do begin
      qq = chip[ss]

      ;; Index
      indx = where(hires.flg_anly NE 0 AND hires.setup EQ setup AND $
                   hires.chip EQ qq  AND $
                   hires.obj_id EQ this_id[ii] AND $
                   (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)
  
      ;; Exposures
      if not keyword_set(exp) then this_exp = lindgen(nindx) $
      else this_exp = exp
      nexp = n_elements(this_exp)

      ;; Process
      if keyword_set( DOPROC ) then $
        hires_proc, hires, SETUP=setup, obj=this_id[ii], CHIP=qq, EXP=this_exp, $
        CLOBBER=clobber

      ;; CR
      if keyword_set( DOCR ) and nexp GT 1 then $
        hires_objcr, hires, setup, this_id[ii], qq, CHK=crchk

      ;; Find and Trace
      if keyword_set( DOFNT ) then $
        hires_fntobj, hires, setup, this_id[ii], qq, this_exp, CHK=TRCCHK, $
        STD=std

      ;; Skysub
      if keyword_set( DOSKY ) and not keyword_set( STD ) then $
        hires_skysub, hires, setup, this_id[ii], qq, this_exp, $
          FCHK=SKYCHK, chk=SKYCHK

      ;; Extract
      if keyword_set( DOEXT ) then $
        hires_extract, hires, setup, this_id[ii], qq, this_exp, $
          CHK=BOXCHK, STD=std, _EXTRA=EXTRA

      ;; Flux
;      if keyword_set( DOFLUX ) then $
;        hires_flux, hires, setup, this_id[ii], qq, _EXTRA=EXTRA

      ;; Combine
;      if keyword_set( DOCOMB ) then $
;        hires_combspec, hires, setup, this_id[ii], qq, this_exp, CHK=COMBCHK, $
;        _EXTRA=EXTRA

      ;; 1D
;      if keyword_set( DO1D ) then $
;        hires_1dspec, hires, setup, this_id[ii], qq, this_exp, CHK=COMBCHK, $
;        _EXTRA=EXTRA

    endfor
  endfor
;  DONE
  print, 'hires_allobj: All done! '
  return
end


