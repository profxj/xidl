;+ 
; NAME:
; x_ccdinf   
;   Version 1.0
;
; PURPOSE:
;    Returns arcpix and orientation of ccd+tel combination.  This
;   routine saves all of that key info in one spot, although it 
;   is primarily used only for direct imaging.
;
; CALLING SEQUENCE:
;   x_ccdinf, ccd, tel, arcpix, [orient], SAT=
;
; INPUTS:
;   ccd   -  Name of CCD ('SITe1', 'LRISR', 'WFTek5')
;   tel   -  Name of Telescope ('LCO-40', 'Keck', 'LCO-100')
;
; RETURNS:
;   arcpix --  Arcsec per pixel
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   orient - Orientation of the CCD:
;;                   -2 = SITe1 on LCO-40  (E up, N left)
;  SAT=   -- Saturation limit of the CCD
;
; COMMENTS:
;
; EXAMPLES:
;   x_ccdinf, 'SITe1', 'LCO-40', arcpix, SAT=sat
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_ccdinf, ccd, tel, arcpix, orient, SAT=sat

  ; 
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_ccdinf, ccd, tel, arcpix, [orient], SAT= (v1.1)'
      return
  endif 

  
; Set Orientation

  case tel of 
      'LCO-40' : begin
          case ccd of
              'SITe1' : begin
                  orient = -2     ; E up, N left
                  arcpix = 0.6964
                  sat = 26000.
              end
              'SITe3' : begin
                  orient = -3     ; S up, E left
                  arcpix = 0.4349    ; WCS solution 4/29/02
                  sat = 26000.
              end
              else : begin
                  print, 'This ccd, tel combo not set up: ', ccd, tel
                  return
              end
          endcase
      end
      'LCO-100' : begin
          case ccd of
              'Tek5' : begin     ; S up, E left
                  orient = -3
                  arcpix = 0.259
                  sat = 25000.
              end
              'WFTek5' : begin     ; Update this!
                  orient = -3
                  arcpix = 0.753
                  sat = 25000.
              end
              else : begin
                  print, 'This ccd, tel combo not set up: ', ccd, tel
                  return
              end
          endcase
      end
      'Keck' : begin
          case ccd of
              'LRISR' : begin
                  orient = 3   ; S up, E right
                  arcpix = 0.213
                  sat = 50000.
              end
              else : begin
                  print, 'This ccd, tel combo not set up: ', ccd, tel
                  return
              end
          endcase
      end
      'KPNO-4m' : begin
          case ccd of
              'MOSA' : begin
                  orient = 1   ; N up, E left
                  arcpix = 0.258
                  sat = 50000.
              end
              else : begin
                  print, 'This ccd, tel combo not set up: ', ccd, tel
                  return
              end
          endcase
      end
      else : begin
          print, 'This tel not supported: ', tel
          return
      end
  endcase


return
end
  
      
