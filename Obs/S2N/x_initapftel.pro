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
;   20-Oct-2005 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_initapftel, apftel

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initapftel, apfstr  [v1.0]'
      return
  endif 

  apftel = {telestruct}
  apftel.AREA  = 42053. ; cm^2 -- 10m telescope with 7.9% central obscuration
  apftel.PLATE_SCALE = 5.8241 ; "/mm

  return
end
