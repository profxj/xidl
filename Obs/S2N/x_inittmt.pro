; NAME:
; x_inittmt
;    Version 1.0
;
; PURPOSE:
;    Initialize a telescope structure for TMT
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Jan-2006 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_inittmt, tmttel

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_initkeck, keckstr  [v1.0]'
      return
  endif 

  tmttel = {telestruct}
  tmttel.AREA  = 7068580.  ; cm^2 -- 30m telescope with no obscuration
  tmttel.PLATE_SCALE = 0.4583 ; "/mm

  return
end
