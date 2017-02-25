;+ 
; NAME:
; uves_setarcfil
;     Version 1.1
;
; PURPOSE:
;  Set the arc file name of a given setup
;
; CALLING SEQUENCE:
;   uves_setarcfil, uves, setup
;
; INPUTS:
;   uves  -  HIRES structure
;   setup  -  Integer defining setup
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
;   uves_setarcfil, uves, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_getfil
;
; REVISION HISTORY:
;   28-Sep-2005 Written by JXP
;
;  Usage:
;-
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uves_setarcfil, uves, setup

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_setarcfil, uves, setup [v1.1]'
      return
  endif 

  obj = where(uves.type EQ 'OBJ' AND $
              uves.setup EQ setup AND $
              uves.flg_anly NE 0,nobj)
  arcs = where(uves.type EQ 'ARC' AND $
               uves.setup EQ setup AND $
               uves.flg_anly NE 0,narc)

  if narc EQ 0 then begin
      print, 'uves_setarcfil: No arcs with this setup! Returning..'
      stop
      return
  endif

  for qq=0L,nobj-1 do begin
      ;; Find arcfil
      arc_fil = uves_getarcfil(uves, obj[qq], arcs, ARC_img=arc_img)
      uves[obj[qq]].arc_fil = arc_fil
      uves[obj[qq]].arc_img = arc_img
  endfor

  print, 'uves_setarcfil: All done!'

  return
end

