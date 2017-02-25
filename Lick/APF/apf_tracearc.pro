;+ 
; NAME:
; apf_tracearc   
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
;  apf_tracearc, apf, setup, [obj_id], /CLOBBER, INIO=
;
; INPUTS:
;   apf     -  MIKE structure
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
;  apf_tracearc, apf, setup, obj_id
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

pro apf_tracearc, apf, setup, obj_id, INIO=inio, CLOBBER=clobber, $
                    ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_tracearc, apf, setup, obj_id, INIO=, /CLOBBER'
      print, '     ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.


; Loop on chip
  if not keyword_set(ARCFIL) then begin
     ;; Grab all obj indices
     indx = where(apf.flg_anly NE 0 AND $
                  apf.obj_id EQ obj_id AND apf.setup EQ setup AND $
                  (strtrim(apf.type,2) EQ 'OBJ' OR $
                   strtrim(apf.type,2) EQ 'STD'), nindx)
     if nindx EQ 0 then begin
        print, 'apf_tracearc: No Obj found!  Returning' 
        return
     endif
     
     arcfil_all= apf[indx].arc_fil
     asort = sort(arcfil_all)
     auniq = uniq(arcfil_all[asort])
  endif else begin
     arcfil_all = [arcfil]
     auniq = [0L]
     asort = [0L]
     pos = strpos(arcfil, '.fits')
     frame = long(strmid(arcfil,pos-5,5))
     indx = where(apf.frame EQ frame)
  endelse

  ;; ORD_STR
  ordr_str = apf_getfil('ordr_str', setup)
  
  ;; Loop on exposures
  for mm=0L,n_elements(auniq)-1 do begin
     
     arc_fil = arcfil_all[asort[auniq[mm]]]
     idx = indx[asort[auniq[mm]]]
     
     ;; Outfil
     pos = strpos(arc_fil, '.fits')
     frame = long(strmid(arc_fil,pos-5,5))
     
     ;;  Check for outfil
     out_fil = apf_getfil('arc_trc', setup, $
                          /name, CHKFIL=chkf,  FRAME=frame)
     ;;  Check for outfil
     if chkf NE 0 and not keyword_set( CLOBBER ) then begin
        print, 'apf_tracearc: Arc fit file exists. ' + $
               'Continuing..'
        return
     endif
;          stop  ;;  Need to add shftprm !!
     ;; QA
     qafil = apf_getfil('qa_tracearc', setup, $
                        /name, CHKFIL=chkf,  FRAME=frame)
     rslt = x_echtrcarc(arc_fil, ordr_str, out_fil, $
                        CLOBBER=clobber, INIO=inio, QAFIL=qafil)
     if size(rslt,/tname) NE 'STRING' then stop
  endfor


  print, 'apf_tracearc:  All done!'
  return
end

