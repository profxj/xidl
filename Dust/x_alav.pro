;+ 
; NAME:
; x_alav
;   Version 1.0
;
; PURPOSE:
;    Report A_lambda/A(V) given lambda for a user specific extinction
;    law. Default is the Cardelli 1989 paramaterization of the MW.
;    The 'laws' are contained in .dat files in XIDL_DIR/Dust
;
; CALLING SEQUENCE:
;   
;   dat = x_alav(lambda, RV=, /SMC, /CALZ)
;
; INPUTS:
;   lambda  -- Wavelengths to evaluate A at (angstroms)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  RV= -- Value of R_V
;  /SMC -- Use the SMC law
;  /CALZ -- Use the Calzetti law
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
;   2006 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_alav, lambda, RV=rv, SMC=smc, CALZ=calz

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'alav = x_alav(lambda, RV=, /SMC, /CALZ) [V1.0]'
    return, -1
  endif 

	;  Optional Keywords

  ;; Milky Way
  if not (keyword_set(SMC) or keyword_set(CALZ)) then $
    return, x_getalav(lambda,RV=rv)

  if keyword_set(SMC) then $
    readcol, getenv('XIDL_DIR')+'/Dust/smcN.dat', wv, AlAV, /sile $
  else readcol, getenv('XIDL_DIR')+'/Dust/calzetti.dat', wv, AlAV, /sile
  return, interpol(alav, wv, lambda)
end

