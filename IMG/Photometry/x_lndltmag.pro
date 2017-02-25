;+ 
; NAME:
; x_lndltmag   
;    Version 1.1
;
; PURPOSE:
;    Returns the Landolt magnitude for a given filter
;
; CALLING SEQUENCE:
;   
; mag = x_lndltmag(filter, lndltstr, SIGMAG=)
;
; INPUTS:
;   filter - String name of the filter
;   lndlstr - Structure of the landolt star
;
; RETURNS:
;   mag - magnitude
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  SIGMAG = Error in the magnitude
;
; COMMENTS:
;
; EXAMPLES:
;   mag = x_lndltmag('R', landolt, SIGMAG=sig)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------
function x_lndltmag, filter, landolt, SIGMAG=sigmag

;
  if  N_params() LT 2  then begin 
      print, 'Syntax - ' +$
        'mag = x_lndltmag(filter, landolt, SIGMAG=)'
      return, -1
  endif 

; Parse on filter

  case filter of 
      'U' : begin
          mag = landolt.V + landolt.BV + landolt.UB
          sigmag = landolt.sig_UB
      end
      'B' : begin
          mag = landolt.V + landolt.BV 
          sigmag = landolt.sig_BV
      end
      'V' : begin
          mag = landolt.V
          sigmag = landolt.sig_V
      end
      'R' : begin
          mag = landolt.V - landolt.VR
          sigmag = landolt.sig_VR
      end
      'I' : begin
          mag = landolt.V - landolt.VI
          sigmag = landolt.sig_VI
      end
      else : begin
          print, 'x_lndltmag: ', filter, ' is not a standard Landolt fitler!'
          return, -1
      end
  endcase
  
  return, mag

end

