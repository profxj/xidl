;+ 
; NAME:
; fit_convxptox
;    Version 1.1
;
; PURPOSE:
;    Convert the fitted coefficients of a POLY fit from the normalized
;    values to 'sensible' ones.
;
; CALLING SEQUENCE:
;   newcoeff = fit_convxptox( calib, FFIT=, NRM= )
;
; INPUTS:
;  calib -- FIT structure
;
; RETURNS:
;  newcoeff -- Unnormalized coefficients
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FFIT= -- Coefficients [default: *calib.ffit]
;  NRM=  -- Normalization values [default: calib.nrm]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------
function fit_convxptox, calib, FFIT=ffit, NRM=nrm

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'fit_convxptox, claib, FFIT=, NRM=, [v1.1]'
    return, -1
  endif 

  if not keyword_set(FFIT) then ffit = *calib.ffit
  if not keyword_set(NRM) then nrm = calib.nrm

  dum1 = 2.d/nrm[1]
  dum2 = -2.d*nrm[0]/nrm[1]

  coeff = dblarr(n_elements(ffit))

  tot = 0.d
  for i=1,n_elements(ffit)-1 do begin
      case i of  
          1 : coeff[1] = ffit[i]*2/nrm[1]
          2 : begin
              coeff[1] = coeff[1] - ffit[i]*8.*nrm[0]/(nrm[1]^2) 
              coeff[2] = ffit[i]*4/(nrm[1]^2) 
          end
          3 : begin
              coeff[1] = coeff[1] + ffit[i]*24*(nrm[0]^2)/(nrm[1]^3)
              coeff[2] = coeff[2] - ffit[i]*24*nrm[0]/(nrm[1]^3)
              coeff[3] = ffit[i]*8/(nrm[1]^3)
          end
          4 : begin
              coeff[1] = coeff[1] - ffit[i]*64*(nrm[0]^3)/(nrm[1]^4)
              coeff[2] = coeff[2] + ffit[i]*96*(nrm[0]^2)/(nrm[1]^4)
              coeff[3] = coeff[3] - ffit[i]*64*nrm[0]/(nrm[1]^4)
              coeff[4] = ffit[i]*16/(nrm[1]^4)
          end
      endcase
      tot = tot + ffit[i] * (dum2^i)
  endfor
  coeff[0] = ffit[0] + tot

  return, coeff

end



