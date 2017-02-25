;+ 
; NAME:
; mike_tracearc   
;     Version 1.1
;
; PURPOSE:
;    To trace the arc lines in each order (individually) and fit a
;    straight line to each one.  The following steps are taken:
;    1.  Scattered light is removed from the image
;    2.  All significant arc lines are identified (5 sigma)
;    3.  trace_crude is used to trace the lines 
;    4.  trace_crude is reapplied to only those lines which are
;    entirely in the order
;    5.  xy2traceset is used to fit a straight line to each arc line
;    6.  Only the good lines are saved for 2D fitting in a structure
;    which is written to disk
;
; CALLING SEQUENCE:
;   
;  mike_tracearc, mike, setup, obj_id, [side], /CLOBBER, INIO=
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;  all_xset -- A structure describing the fits to the arc lines
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   INIO      - Initial order (for debugging)
;   /CLOBBER  - Overwrite previous fits
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  mike_tracearc, mike, setup, obj_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  trace_crude
;  xy2traceset
;  mike_tracearc_work
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mike_tracearc_work, arc_fil, setup, side, $
                             slopecoeff=slopecoeff, ximage=ximage, $
                             INIO=inio, ALL_XSET=all_xset, CLOBBER=clobber, $
                             SHFTPRM=shftprm
;                             XOFF=xoff
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_tracearc_sngl(arc_fil, setup, side, ' + $
        'SLOPECOEFF=, INIO=, /CLOBBER, ALL_XSET=, XIMAGE= ) [v1.1]'
      return, -1
  endif 

;  Optional Keywords
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.
  if not keyword_set( SHFTPRM ) then shftprm = [0., 0.]
;  if not keyword_set( XOFF ) then xoff =  0.

; Loop on side

  ;; OUTFIL
  out_fil = mike_getfil('arc_trc', subfil=arc_fil, /name, CHKFIL=chkf)
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
      print, 'mike_tracearc: Arc fit file exists. ' + $
        'Continuing..'
      return, out_fil
  endif

  ;; Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  nordr = n_elements(ordr_str)
          
  ;; Arc file
  qafil = mike_getfil('qa_tracearc', setup, SUBFIL=arc_fil) 

  rslt = x_echtrcarc(arc_fil, ordr_str, out_fil, $
                     CLOBBER=clobber, INIO=inio, QAFIL=qafil)

  return, out_fil
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_tracearc, mike, setup, obj_id, side, INIO=inio, CLOBBER=clobber
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_tracearc, mike, setup, obj_id, [side], INIO=, /CLOBBER'
      print, '     [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.


; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]
      ;; SIDE
      if qq EQ 1 then print, 'mike_tracearc: Tracing BLUE arc' $
      else print, 'mike_tracearc: Tracing RED arc'

      ;; Grab all obj indices
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                   strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_tracearc: No Obj found!  Returning' 
          continue
      endif

      arcfil_all= mike[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

      ;; Loop on exposures
      for mm=0L,n_elements(auniq)-1 do begin

          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          stop  ;;  Need to add shftprm !!
          rslt = mike_tracearc_work(arc_fil, setup, qq, $
                  CLOBBER=clobber, INIO=inio)
;                  CLOBBER=clobber, XOFF=mike[idx].arc_xyoff[0], INIO=inio)
          if size(rslt,/tname) NE 'STRING' then stop
      endfor

  endfor

  print, 'mike_tracearc:  All done!'
  return
end

