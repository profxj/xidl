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
pro x_initkeck, kecktel

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initkeck, keckstr  [v1.0]'
      return
  endif 

  kecktel = {telestruct}
  kecktel.AREA  = 723674. ; cm^2 -- 10m telescope with 7.9% central obscuration
  kecktel.PLATE_SCALE = 1.379 ; "/mm

  return
end
