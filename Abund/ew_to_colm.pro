;+ 
; NAME:
; ew_to_colm   
;   Version 1.1
;
; PURPOSE:
;    Calculates the column density given EW and rest wave (or
;    vice-versa).  Assumes the transition lies on the linear COG.
;
; CALLING SEQUENCE:
;   colm = ew_to_colm([lambda], [EW], /RVRS)
;
; INPUTS:
;   lambda  - Rest Wavelength (Ang)
;   EW       - EW (mA) [or column density (linear)]
;
; RETURNS:
;   colm   - Column density (linear)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SILENT -- Suppress written output
;  /RVRS   -- Take a column density input and output the EW assuming
;             the linear COG
;  TOLER=  -- Tolerance on wavelength [default: 0.1]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   colm = ew_to_colm([1215.6701],[50.])
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function ew_to_colm, lambda, EW, RVRS=rvrs, SILENT=silent, TOLER=toler

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'ew = ew_to_colm([labmda], [EW], /RVRS, /SILENT, TOLER=) [v1.1]'
    return, -1
  endif 

 ;; Optional Keywords

 nlmb = n_elements(lambda)
 if n_elements(EW) NE nlmb then return, -1
 Nval = dblarr(nlmb)

 if not keyword_set(TOLER) then close = 1

 ;; Loop
 for i=0L,nlmb-1 do begin
     ;; f-value
     getfnam, lambda[i], fval, CLOSE=close, NEWWV=newwav, TOLER=toler

     ;; Assume weak limit
     if not keyword_set( RVRS ) then $
       Nval[i] = (EW[i]*1e-3 / lambda[i]) / 8.85249e-13 / fval / $
       (lambda[i] * 1e-8) $
     else Nval[i] = EW[i]/((1e-3 / lambda[i]) / 8.85249e-13 / fval / $
                           (lambda[i] * 1e-8)) 
       
     ;; Print
     if not keyword_set(SILENT) then begin
         if not keyword_set(RVRS) then $
           print, 'ew_to_colm: lambda = ', newwav, ' N = ', alog10(Nval[i])$
         else $
           print, 'ew_to_colm: lambda = ', newwav, ' EW(mA) = ', Nval[i]
     endif
 endfor

 return, Nval

end
     
