;+ 
; NAME:
; hires_setarcfil
;     Version 1.1
;
; PURPOSE:
;  Set the arc file name of a given setup
;
; CALLING SEQUENCE:
;   hires_setarcfil, hires, setup
;
; INPUTS:
;   hires  -  HIRES structure
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
;   hires_setarcfil, hires, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_getfil
;
; REVISION HISTORY:
;   28-Sep-2005 Written by JXP
;
;  Usage:
;-
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_setarcfil, hires, setup

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_setarcfil, hires, setup [v1.1]'
      return
  endif 

  obj = where(hires.type EQ 'OBJ' AND $
              hires.setup EQ setup AND $
              hires.flg_anly NE 0,nobj)
  arcs = where(hires.type EQ 'ARC' AND $
               hires.setup EQ setup AND $
               hires.flg_anly NE 0,narc)

  if narc EQ 0 then begin
      print, 'hires_setarcfil: No arcs with this setup! Returning..'
      stop
      return
  endif

  for qq=0L,nobj-1 do begin
      ;; Find arcfil
      arc_fil = hires_getarcfil(hires, obj[qq], arcs, ARC_img=arc_img)
      hires[obj[qq]].arc_fil = arc_fil
      hires[obj[qq]].arc_img = arc_img
  endfor

  print, 'hires_setarcfil: All done!'

  return
end

