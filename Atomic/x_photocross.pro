;+ 
; NAME:
; x_photocross
;   Version 1.1
;
; PURPOSE:
;    Calculates the photoionization cross-section given by Verner &
;    Ferland (1996)
;
; CALLING SEQUENCE:
;   
; cross = x_photocross( Z, ion, E )
;
; INPUTS:
;   Z      - Atomic number
;   ion    - Number of bound electrons
;   E      - Energy of photon (eV)
;
; RETURNS: 
;  cross -- Cross-section in cm^2
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Allowed ions:
;     CI, CIV, MgI, FeI
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   5-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_photocross, Z, ion, E, ETH=eth

  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'x_photocross( Z, ion, E [eV]) (v1.1)'
      return, -1
  endif 
;
  ;; Optional keywords
  

  ;; Constants
  c = x_constants()

  ;; Element
  case Z of
      1: begin
          case ion of 
              1: begin
                  Eth = 13.60
                  E0 = 4.298e-1
                  s0 = 5.475e4
                  ya = 3.288e1
                  P = 2.963
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      2: begin
          case ion of 
              1: begin
                  Eth = 24.59
                  E0 = 13.61
                  s0 = 9.492e2
                  ya = 1.469
                  P = 3.188
                  yw = 2.039
                  y0 = 4.434e-1
                  y1 = 2.136
              end
              2: begin
                  Eth = 54.42
                  E0 = 1.72
                  s0 = 1.369e4
                  ya = 3.288e1
                  P = 2.963
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      6: begin
          case ion of 
              1: begin
                  Eth = 11.26
                  E0 = 2.144
                  s0 = 5.027e2
                  ya = 6.216e1
                  P = 5.101
                  yw = 9.157e-2
                  y0 = 1.133
                  y1 = 1.607
              end
              2: begin
                  Eth = 24.38
                  E0 = 4.058e-1
                  s0 = 8.709
                  ya = 1.261e2
                  P = 8.578
                  yw = 2.093
                  y0 = 4.929e1
                  y1 = 3.234
              end
              3: begin
                  Eth = 47.89
                  E0 = 4.614
                  s0 = 1.539e4
                  ya = 1.737
                  P = 1.593e1
                  yw = 5.922
                  y0 = 4.378e-3
                  y1 = 2.528e-2
              end
              4: begin
                  Eth = 64.49
                  E0 = 3.506
                  s0 = 1.068e2
                  ya = 1.436e1
                  P = 7.457
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      7: begin
          case ion of 
              1: begin
                  Eth = 14.53
                  E0 = 4.034
                  s0 = 8.235e2
                  ya = 80.33
                  P = 3.928
                  yw = 9.097e-2
                  y0 = 8.598e-1
                  y1 = 2.325
              end
              2: begin
                  Eth = 29.60
                  E0 = 6.128e-2
                  s0 = 1.944
                  ya = 8.163e2
                  P = 8.773
                  yw = 10.43
                  y0 = 4.28e2
                  y1 = 2.03e1
              end
              3: begin
                  Eth = 47.45
                  E0 = 2.42e-1
                  s0 = 9.375e-1
                  ya = 2.788e2
                  P = 9.156
                  yw = 1.85
                  y0 = 1.877e2
                  y1 = 3.999
              end
              4: begin
                  Eth = 77.47
                  E0 = 5.494
                  s0 = 1.69e4
                  ya = 1.714
                  P = 17.06
                  yw = 7.904
                  y0 = 6.415e-3
                  y1 = 1.937e-2
              end
              5: begin
                  Eth = 97.89
                  E0 = 4.471
                  s0 = 8.376e1
                  ya = 3.297e1
                  P = 6.003
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      8: begin
          case ion of 
              1: begin
                  Eth = 13.62
                  E0 = 1.24
                  s0 = 1.745e3
                  ya = 3.784
                  P =  1.764e1
                  yw = 7.589e-2
                  y0 = 8.698
                  y1 = 1.271e-1
              end
              7: begin
                  Eth = 7.393e2
                  E0 = 8.709e1
                  s0 = 1.329e2
                  ya = 2.535e1
                  P =  2.336
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              8: begin
                  Eth = 8.714e2
                  E0 = 2.754e1
                  s0 = 8.554e2
                  ya = 3.288e1
                  P =  2.963
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      12: begin
          case ion of 
              1: begin
                  Eth = 7.646
                  E0 = 1.197e1
                  s0 = 1.372e8
                  ya = 2.228e-1
                  P = 1.574d1
                  yw = 2.805e-1
                  y0 = 0.
                  y1 = 0.
              end
              2: begin
                  Eth = 15.04
                  E0 = 8.139
                  s0 = 3.278
                  ya = 4.341e7
                  P = 3.610
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      13: begin
          case ion of 
              2: begin
                  Eth = 18.83
                  E0 = 2.048e-1
                  s0 = 6.948e-2
                  ya = 5.675e2
                  P = 9.049
                  yw = 4.615e-1
                  y0 = 9.149e1
                  y1 = 6.565e-1
              end
              3: begin
                  Eth = 28.45
                  E0 = 1.027e1
                  s0 = 4.915
                  ya = 1.99e6
                  P = 3.477
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      14: begin
          case ion of 
              2: begin
                  Eth = 16.35
                  E0 = 2.556
                  s0 = 4.140
                  ya = 1.337e1
                  P = 11.91
                  yw = 1.57
                  y0 = 6.634
                  y1 = 1.272e-1
              end
              3: begin
                  Eth = 33.49
                  E0 = 1.659e-1
                  s0 = 5.79e-4
                  ya = 1.474e2
                  P = 1.336e1
                  yw = 8.626e-1
                  y0 = 9.613e1
                  y1 = 6.442e-1
              end
              4: begin
                  Eth = 45.14
                  E0 = 12.88
                  s0 = 6.083
                  ya = 1.356e6
                  P = 3.353
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      26: begin
          case ion of 
              1: begin
                  Eth = 7.902
                  E0 = 5.461e-2
                  s0 = 3.062d-1
                  ya = 2.671e7
                  P = 7.923
                  yw = 20.69
                  y0 = 138.2
                  y1 = 2.481e-1
              end
              26: begin
                  Eth = 9.278e3
                  E0 = 2.932e2
                  s0 = 8.099e1
                  ya = 3.288e1
                  P = 2.963
                  yw = 0.
                  y0 = 0.
                  y1 = 0.
              end
              else: stop
          endcase
      end
      else: stop
  endcase

  x = E/E0 - y0
  y = sqrt(x^2 + y1^2)

  F = ((x-1.)^2 + yw^2) * y^(0.5*P - 5.5) * $
    (1 + sqrt(y/ya) )^(-1.*P)

  sigma = s0 * F * 1e-18 ; cm^2

  ;; Energy threshold
  low = where(E LT Eth,nlow)
  if nlow NE 0 then sigma[low] = 0.

  return, sigma
end
