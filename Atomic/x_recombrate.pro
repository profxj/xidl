;+ 
; NAME:
; x_recombrate
;   Version 1.0
;
; PURPOSE:
;    Calculates the recombination rate for an ion given 
;  values specifying the ion and the gas temperature.
;
; CALLING SEQUENCE:
;   
; recomb = x_recombrate( Z, ion, T )
;
; INPUTS:
;   Z      - Atomic number
;   ion    - Number of free electrons
;   T      - Energy of photon (K)
;
; RETURNS: 
;  rate -- s^-1 cm^3
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;    CII->CI
;    MgII->MgI
;    FeII->FeI
;
; EXAMPLES:
;   
;      print, x_recombrate(12, 1, 1000.)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   5-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_recombrate, Z, ion, T

  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'x_recombrate( Z, ion [lower], T) (v1.1)'
      return, -1
  endif 
;
  ;; Optional keywords
  

  ;; Constants
  c = x_constants()

  ;; Element
  case Z of
      6: begin
          case ion of 
              1: begin ; CII -> CI 
                  ;; Radiative (from Cloudy; ion_recomb.c)
                  k = [7.651e-09,0.8027,1.193e-03,9.334e12]
                  tt = sqrt(T / k[2])
                  alpha_r = k[0] / $
                         (tt * (tt+1)^(1-k[1]) * (1+sqrt(t/k[3]))^(1+k[1]) )
                  ;; Dielectronic
                  alpha_d = 0.  ; Could be wrong...
                  ;; Total
                  alpha = alpha_r + alpha_d
              end
              else: stop
          endcase
      end
      12: begin
          case ion of 
              1: begin  ; MgII -> MgI
                  ;; Radiative
                  c_r = 1.4d-13
                  e_r = -0.855
                  alpha_r = c_r * (T / 1d4)^e_r ; Radiative
                  ;; Dielectronic
                  c_d1= 4.49d-4
                  c_d2 = 0.0021
                  c_T1 = -5.01d4
                  c_T2 = -2.61d4
                  alpha_d = c_d1 * T^(-1.5) * exp(c_T1/T) *(1 + c_d2*exp(c_T2/T))
                  ;; Total
                  alpha = alpha_r + alpha_d
              end
		else: stop
	endcase
	end
      26: begin
          case ion of 
              1: begin ; FeII -> FeI  (from Cloudy)
                  ;; Radiative
                  k = [8.945e-9, 0.2156, 4.184e-2, 5.353e13]
                  tt = sqrt(T / k[2])
                  alpha_r = k[0] / $
                         (tt * (tt+1)^(1-k[1]) * (1+sqrt(t/k[3]))^(1+k[1]) )
                  ;; Dielectronic
                  te = T * c.k / c.eV
                  d = 2.2e-4 * exp(-5.12/te) + 1e-4 * exp(-12.9/te)
                  alpha_d = d * T^(-1.5)
                  ;; Total
                  alpha = alpha_r + alpha_d
              end
              else: stop
          endcase
      end
      else: stop
  endcase

  return, alpha

end
