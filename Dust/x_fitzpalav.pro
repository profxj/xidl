;+ 
; NAME:
; x_fitzpalav
;   Version 1.0
;
; PURPOSE:
;    Report A_lambda/A(V) given lambda for a Fitzpatrick & Massa
;    parameterized extinction law.
;
; CALLING SEQUENCE:
;   
;   AlAV = x_alav(lambda, RV=, c3=, c4=)
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
;  Written   W. Landsman        Raytheon  STX   October, 1998
;  2008 Revised by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fitzpalav, lambda, R_V=r_v, c3=c3, c4=c4

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'alav = x_alav(lambda, R_V=, c3=, c4= ) [V1.0]'
    return, -1
  endif 

  if not keyword_set(R_V) then R_V = 3.1
  x = 10000./ lambda
  curve = x*0.

  if N_elements(x0) EQ 0 then x0    =  4.596  
  if N_elements(gamma) EQ 0 then gamma =  0.99	
  if N_elements(c3) EQ 0 then c3    =  3.23	
  if N_elements(c4) EQ 0 then c4   =  0.41    
  if N_elements(c2) EQ 0 then c2    = -0.824 + 4.717/R_V
  if N_elements(c1) EQ 0 then c1    =  2.030 - 3.007*c2


; Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
; R-dependent coefficients
 
  xcutuv = 10000.0/2700.0
  xspluv = 10000.0/[2700.0,2600.0]
  iuv = where(x ge xcutuv, N_UV, complement = iopir, Ncomp = Nopir)
  IF (N_UV GT 0) THEN xuv = [xspluv,x[iuv]] ELSE  xuv = xspluv

  yuv = c1  + c2*xuv
  yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma)^2)
  yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)^2+0.05644*((xuv>5.9)-5.9)^3)
  yuv = yuv + R_V
  yspluv  = yuv[0:1]            ; save spline points

  IF (N_UV GT 0) THEN curve[iuv] = yuv[2:*] ; remove spline points
 
; Compute optical portion of A(lambda)/E(B-V) curve
; using cubic spline anchored in UV, optical, and IR

  xsplopir = [0,10000.0/[26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]]
  ysplir   = [0.0,0.26469,0.82925]*R_V/3.1 
  ysplop   = [poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ), $
              poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ), $
              poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ), $
              poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, $ 
                          -4.45636e-05] ) ]
  
  ysplopir = [ysplir,ysplop]
  
  if (Nopir GT 0) then $
    curve[iopir] = CSPLINE([xsplopir,xspluv],[ysplopir,yspluv],x[iopir])

 ; Now apply extinction correction to input flux vector

;   Alambda = ebv*curve ;; Alambda 
   Extcurve = curve/R_V ;; Alambda/AV
;   if N_params() EQ 3 then flux = flux * 10.^(0.4*curve) else $
;        funred = flux * 10.^(0.4*curve)       ;Derive unreddened flux

;   ExtCurve = Curve - R_V  ;; E(l-V)/E(B-V)

   return, Extcurve

end

