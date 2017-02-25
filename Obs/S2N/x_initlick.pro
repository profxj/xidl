; NAME:
; x_initkeck
;    Version 1.0
;
; PURPOSE:
;    Initialize a telescope structure for Keck
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21-Mar-2008 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_initlick, licktel

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initlick, lickstr  [v1.0]'
      return
  endif 

  licktel = {telestruct}
  licktel.AREA  = 63617. ; cm^2 -- 3m telescope with 10% central obscuration
  licktel.PLATE_SCALE = 1.379 ; "/mm  [Should confirm; secondary dependent]
  licktel.name = 'Lick-3m'

  return
end
