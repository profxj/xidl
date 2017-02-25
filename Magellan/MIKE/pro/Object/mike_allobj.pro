;+
; NAME: mike_allobj   
;     Version 1.1
;
; PURPOSE:
;    Run all of the processing steps on a single object.  This
;    includes processing, CR rejection, tracing, sky subtraction and
;    extraction.  Future version will allow multiple (or all objects)
;    and also fluxing.
;
;
; CALLING SEQUENCE:
;   
;    mike_allobj, mike, setup, obj_id, side, [exp], /PROCALL
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) or Red (2) side
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
;  /NOCR    - Turn off CR rejection routine (mike_objcr)
;  /CLOBBER - Overwrite previous images (mainly the processed image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_allobj, mike, 1L, 1L
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

pro mike_allobj, mike, setup, obj_id, side, exp, PROCALL=PROCALL, $
                 RESCHK=reschk, SKYCHK=skychk, TRCCHK=trcchk, CHKALL=chkall, $
                 DOSKY=dosky, DOBOX=dobox, NOCR=nocr, CLOBBER=clobber, $
                 DOCOMB=docomb, DOFLUX=doflux, DO1D=do1d, _EXTRA=EXTRA, $
                 DOFNT=dofnt, DOPROC=doproc, BOXCHK=boxchk

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_allobj, mike, setup, obj_id, side, [exp], /PROCALL, /RESCHK' + $
        '/SKYCHK, /BOXCHK, /TRCCHK, /DOSKY, /DOBOX, /DOCOMB, /NOCR [v1.1]' 
      return
  endif 

;  Optional Keywords

  if keyword_set( PROCALL ) then begin
      DOPROC = 1
      DOCR = 1
      DOFNT = 1
      DOSKY = 1
      DOBOX = 1
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

  if NOT keyword_set(side) then side = [1L,2L]
  trimtype = strtrim(mike.type,2)

  if keyword_set(obj_id) then this_id = obj_id $
  else begin 
    check_id = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   (trimtype EQ 'OBJ' OR trimtype EQ 'STD'))
    if check_id[0] EQ -1 then begin
       print, 'No obj_id found'
       return
    endif
    sort_id = mike[check_id[sort(mike[check_id].obj_id)]].obj_id
    this_id = sort_id[uniq(sort_id)]
 endelse

  for ii=0, n_elements(this_id)-1 do begin

    for ss=0L,n_elements(side)-1 do begin
      qq = side[ss]

      ;; Index
      indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   mike.side EQ qq  AND $
                   mike.obj_id EQ this_id[ii] AND $
                   (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)

      ;; Exposures
      if not keyword_set(exp) then this_exp = lindgen(nindx) $
      else this_exp = exp
      nexp = n_elements(this_exp)

      ;; Process
      if keyword_set( DOPROC ) then $
        mike_proc, mike, indx, CLOBBER=clobber
;        mike_proc, mike, SETUP=setup, obj=this_id[ii], SIDE=qq, EXP=this_exp, $
;        CLOBBER=clobber

      ;; CR
      if keyword_set( DOCR ) and nexp GT 1 then $
        mike_objcr, mike, setup, this_id[ii], qq, CHK=crchk

      ;; Find and Trace
      if keyword_set( DOFNT ) then $
        mike_fntobj, mike, setup, this_id[ii], qq, this_exp, CHK=TRCCHK

      ;; Skysub
      if keyword_set( DOSKY ) then $
        mike_skysub, mike, setup, this_id[ii], qq, this_exp, $
          FCHK=SKYCHK, chk=SKYCHK

      ;; Extract
      if keyword_set( DOBOX ) then $
        mike_box, mike, setup, this_id[ii], qq, this_exp, RESCHK=RESCHK, $
          CHK=BOXCHK, _EXTRA=EXTRA

      ;; Flux
      if keyword_set( DOFLUX ) then $
        mike_flux, mike, setup, this_id[ii], qq, _EXTRA=EXTRA

      ;; Combine
      if keyword_set( DOCOMB ) then $
        mike_combspec, mike, setup, this_id[ii], qq, this_exp, CHK=COMBCHK, $
        _EXTRA=EXTRA

      ;; 1D
      if keyword_set( DO1D ) then $
        mike_1dspec, mike, setup, this_id[ii], qq, this_exp, CHK=COMBCHK, $
        _EXTRA=EXTRA

    endfor
  endfor
;  DONE
  print, 'mike_allobj: All done! '
  return
end


