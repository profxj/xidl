;+ 
; NAME:
; hires_tracearc   
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
;  hires_tracearc, hires, setup, [obj_id, chip], /CLOBBER, INIO=
;
; INPUTS:
;   hires     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  Output structure describing the fits to the arc lines
;
; OPTIONAL KEYWORDS:
;   INIO      - Initial order (for debugging)
;   /CLOBBER  - Overwrite previous fits
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  hires_tracearc, hires, setup, obj_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  trace_crude
;  xy2traceset
;  x_echtracearc
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_tracearc, hires, setup, obj_id, chip, INIO=inio, CLOBBER=clobber, $
                    ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_tracearc, hires, setup, obj_id, [chip], INIO=, /CLOBBER'
      print, '     ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
  if keyword_set(ARCFIL) then begin
      ;; Grab Chip from name
      case strmid(arcfil,9,1) of
          'S': chip = -1L
          'B': chip = 1L
          'G': chip = 2L
          'R': chip = 3L
          else: stop
      endcase
  endif
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.


; Loop on chip
  for kk=0L,n_elements(chip)-1 do begin
      qq = chip[kk]
      ;; CHIP
      case qq of
          -1: print, 'hires_tracearc: Tracing Single chip'
          1: print, 'hires_tracearc: Tracing BLUE arc' 
          2: print, 'hires_tracearc: Tracing GREEN arc' 
          3: print, 'hires_tracearc: Tracing RED arc' 
      endcase

      if not keyword_set(ARCFIL) then begin
          ;; Grab all obj indices
          indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                       hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                       (strtrim(hires.type,2) EQ 'OBJ' OR $
                        strtrim(hires.type,2) EQ 'STD'), nindx)
          if nindx EQ 0 then begin
              print, 'hires_tracearc: No Obj found!  Returning' 
              continue
          endif

          arcfil_all= hires[indx].arc_fil
          asort = sort(arcfil_all)
          auniq = uniq(arcfil_all[asort])
      endif else begin
          arcfil_all = [arcfil]
          auniq = [0L]
          asort = [0L]
          pos = strpos(arcfil, '.fits')
          frame = long(strmid(arcfil,pos-4,4))
          indx = where(hires.frame EQ frame AND hires.chip EQ chip)
      endelse

      ;; ORD_STR
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip)

      ;; Loop on exposures
      for mm=0L,n_elements(auniq)-1 do begin

          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; Outfil
          pos = strpos(arc_fil, '.fits')
          frame = long(strmid(arc_fil,pos-4,4))

          ;;  Check for outfil
          out_fil = hires_getfil('arc_trc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          ;;  Check for outfil
          if chkf NE 0 and not keyword_set( CLOBBER ) then begin
              print, 'hires_tracearc: Arc fit file exists. ' + $
                'Continuing..'
              continue
          endif
;          stop  ;;  Need to add shftprm !!
          ;; QA
          qafil = hires_getfil('qa_tracearc', setup, CHIP=chip, $
                                 /name, CHKFIL=chkf,  FRAME=frame)
          rslt = x_echtrcarc(arc_fil, ordr_str, out_fil, $
                             CLOBBER=clobber, INIO=inio, QAFIL=qafil)
          if size(rslt,/tname) NE 'STRING' then stop
      endfor

  endfor

  print, 'hires_tracearc:  All done!'
  return
end

