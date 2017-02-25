;+ 
; NAME:
; uves_fittrcarc
;     Version 1.1
;
; PURPOSE:
;   To fit the slope of the arc lines as a function of order number
;   and y position on the CCD.  This information is then used to
;   construct a 2D wavelength image.  The fitting routine is the usual
;   least-squares algorithm with two rounds of rejection.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   uves     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  Fits file with the coefficients of the 2D fit.  Filename like
;  'Arcs/TRC/Arc_mb0539_F.fits' 
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Plots residuals
;  /CLOBBER -- Overwrite previous solution
;  /ORDRCLOB -- Overwrite arc_m in the order structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_fittrcarc_work -- Main driver
;  uves_getfil
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_fittrcarc_work, arc_fil, setup, side, CHK=chk, CLOBBER=clobber, $
                      ORDRCLOB=ordrclob, NYCOEFF=nycoeff, NOCOEFF=nocoeff
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = uves_fittrcarc_work( arc_fil, setup, side, /CHK, /CLOBBER' + $
        ' /ORDRCLOB, NYCOEFF=, NOCOEFF=) [v2.0]'
      return, -1
  endif 

;  Optional Keywords

  if keyword_set( CLOBBER ) then ordrclob = 1

  if side EQ 1 then begin
      if NOT keyword_set(nycoeff) then nycoeff = 2
      if NOT keyword_set(nocoeff) then nocoeff = 3
  endif else begin
;      if NOT keyword_set(nycoeff) then nycoeff = 3
;      if NOT keyword_set(nocoeff) then nocoeff = 6
      if NOT keyword_set(nycoeff) then nycoeff = 2
      if NOT keyword_set(nocoeff) then nocoeff = 3
  endelse

  ;; wcen
  wcen = uves_getwcen(arc_fil,FRAME=frame, /AFIL)
  
; Loop on side
  out_fil = uves_getfil('arc_fittrc', setup, FRAME=frame, WCEN=wcen, /name, CHKFIL=chkf)
  if CHKF NE 0 then begin
      print, 'uves_fittrcarc: File exists: ', out_fil
      if not keyword_set( CLOBBER ) then begin
          print, 'uves_fittrcarc_work: Use /CLOBBER to overwrite'
          return, out_fil
      endif else print, 'uves_fittrcarc_work: Overwriting', out_fil
  endif
          
  ;; Read in order structure
  ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen, FIL_NM=ordr_fil)
  nordr = n_elements(ordr_str)

  ;; Grab Arc Trace
  trc_fil = uves_getfil('arc_trc', setup, WCEN=wcen, FRAME=frame, /name, CHKFIL=chkf)

  qafil = uves_getfil('qa_fittrcarc', setup, WCEN=wcen,  FRAME=frame, /name)

  rslt = x_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, qafil, $
                     CHK=chk, CLOBBER=clobber, $
                     ORDR_FIL=ordr_fil, _EXTRA=extra) 
          
  return, out_fil
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_fittrcarc, uves, setup, obj_id, side, CHK=chk, CLOBBER=clobber
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_fittrcarc, uves, setup, obj_id, [side], /CHK [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]

; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]

      ;; SIDE
      if qq EQ 1 then print, 'uves_fittrcarc: Fitting BLUE trcarc' $
      else print, 'uves_fittrcarc: Fitting RED trcarc'

      ;; Grab all obj indices
      indx = where(uves.flg_anly NE 0 AND uves.side EQ qq AND $
                   uves.obj_id EQ obj_id AND uves.setup EQ setup AND $
                   (strtrim(uves.type,2) EQ 'OBJ' OR $ 
                   strtrim(uves.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'uves_fittrcarc: No Obj found!  Returning' 
          continue
      endif

      arcfil_all= uves[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

      ;; LOOP on exposures
      for mm=0L,n_elements(auniq)-1 do begin
          arc_fil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; Main driver
          rslt = uves_fittrcarc_work(arc_fil, setup, qq, $
                                     CHK=chk, CLOBBER=clobber) 
          if size(rslt,/tname) NE 'STRING' then stop
      endfor
  endfor

  print, 'uves_fittrcarc: All done!'

  return
end
