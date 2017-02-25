;+ 
; NAME:
; x_getalav
;   Version 1.0
;
; PURPOSE:
;    Report A_lambda/A(V) given lambda.  
;      Adopts the Cardelli 1989 paramaterization
;
; CALLING SEQUENCE:
;   
;   A_lambda = x_getalav(lambda, RV=)
;
; INPUTS:
;   lambda -- Wavelength to calculate extiction at
;
; RETURNS:
;   A_lambda -- Extinction (relative to A_V)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  R_V=  --  The R_V value for dust
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_getalav, lambda, RV=rv

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'alav = x_getalav(lambda, RV=) [V1.0]'
    return, -1
  endif 

	;  Optional Keywords
  if not keyword_set( RV ) then RV = 3.1

  nx = n_elements(lambda)
  alav = dblarr(nx)
  ax = dblarr(nx)
  bx = dblarr(nx)

  ;; 
  x = 1.d4 / lambda
  y = (x - 1.82)

  ;; IR
  b = where(x LT 1.1, nb)
  if nb NE 0 then begin
      ax[b] = 0.574*(x[b]^1.61)
      bx[b] = -0.527*(x[b]^1.61)
  endif 

  b = where(x GE 1.1 AND x LE 3.3, nb)
  ;; Optical
  if nb NE 0 then begin 
      ax[b] = 1. + 0.17699*y[b] - 0.50447*(y[b]^2) - 0.02427*(y[b]^3) $
        + 0.72085*(y[b]^4) $
        + 0.01979*(y[b]^5) - 0.7753*(y[b]^6) + 0.33*(y[b]^7)
      bx[b] = 1.41338*y[b] + 2.28305*(y[b]^2) + 1.07233*(y[b]^3) $
        - 5.38434*(y[b]^4) $
        -0.62251*(y[b]^5) + 5.30260*(y[b]^6) - 2.09002*(y[b]^7)
  endif

  b = where(x GE 3.3 AND x LE 8, nb)
  ;; UV
  if nb NE 0 then begin 
      ;; Fa, Fb
      Fa = fltarr(nb)
      Fb = fltarr(nb)
      gg = where(x[b] GT 5.9,ngg,complement=hh,ncomplement=nhh)
      if ngg NE 0 then begin ;; x>5.9
          Fa[gg] = -0.04473*(x[b[gg]]-5.9)^2 - 0.009779*$
            (x[b[gg]]-5.9)^3
          Fb[gg] = 0.2130*(x[b[gg]]-5.9)^2 + 0.1207*$
            (x[b[gg]]-5.9)^3
      endif
      
      ;; Fb
      ax[b] = 1.752 - 0.316*x[b] - 0.104/( (x[b]-4.67)^2 + $
                                           0.341) + Fa
      bx[b] = -3.090 + 1.825*x[b] + 1.206/ $
        ( (x[b]-4.62)^2 + 0.263) + Fb
  endif

  ;; Far-UV
  b = where(x GT 8 AND x LE 10, nb)
  if nb NE 0 then begin 
      ax[b] = 1.073 - 0.628*(x[b]-8) + 0.137*(x[b]-8)^2 - 0.070*(x[b]-8)^3 
      bx[b] = 13.67 + 4.257*(x[b]-8) - 0.420*(x[b]-8)^2 + 0.374*(x[b]-8)^3
  endif

  alav = ax + bx/RV

  return, alav

end
