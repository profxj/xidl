;+ 
; NAME:
; ism_uvradfield
;   Version 1.0
;
; PURPOSE:
;    Calculates the far-UV radiation field at a specified wavelength
;    using the formalism from GPW80
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   lambda -- Wavelength for the evaluation
;   [model] -- 1: GPW80
;   /FNU -- Return Fnu instead of Flambda
;
; RETURNS:
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2006 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function ism_uvradfield, lambda, model, FNU=fnu

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'flux = ism_uvradfield(lambda (Ang), [model]) [v1.0]'
      return, -1
  endif 

  ;; Optional keywords
  if not keyword_set(MODEL) then model = 1

  ;; 
  c = x_constants()

  case model of
      1: begin ; Gondhalekar et al. 1980 (GPW80)
          obs = [ $
                [930., 8.31], $
                [975., 14.82], $
                [1025., 14.57], $
                [1075., 15.93], $
                [1125., 18.16], $
                [1175., 16.67], $
                [1220.,  9.81], $
                [1270.,  17.54], $
                [1325.,  17.78], $
                [1380.,  15.69], $
                [1400.,  14.94], $
                [1420.,  14.94], $
                [1440.,  15.44], $
                [1460.,  15.69], $
                [1480.,  15.31], $
                [1500.,  15.07], $
                [1520.,  14.33], $
                [1540.,  13.46], $
                [1560.,  13.09], $
                [1580.,  13.46], $
                [1600.,  12.84], $
                [1620.,  12.47], $
                [1640.,  12.84], $
                [1660.,  13.34], $
                [1680.,  13.59], $
                [1700.,  13.34], $
                [1720.,  12.60], $
                [1780.,  11.58], $
                [1800.,  11.82], $
                [1820.,  11.41], $
                [1840.,  11.04], $
                [1860.,  10.51], $
                [1880.,  10.55], $
                [1900.,  10.10], $
                [1920.,   9.76], $
                [1940.,   9.72], $
                [1960.,   9.72], $
                [1980.,   9.66], $
                [2000.,   9.31], $
                [2020.,   8.88], $
                [2040.,   8.69], $
                [2060.,   8.41], $
                [2080.,   8.23], $
                [2100.,   8.14], $
                [2120.,   7.89], $
                [2180.,   7.51], $
                [2200.,   7.36], $
                [2220.,   7.36], $
                [2240.,   7.14], $
                [2260.,   7.27], $
                [2280.,   7.27], $
                [2300.,   7.14], $
                [2320.,   6.92], $
                [2340.,   6.72], $
                [2360.,   6.56], $
                [2380.,   6.36], $
                [2400.,   6.37], $
                [2420.,   6.43], $
                [2440.,   6.56], $
                [2460.,   6.52], $
                [2480.,   6.48], $
                [2500.,   6.50], $
                [2520.,   6.22], $
                [2740.,   5.38] ]
          ;; Interpolate
          flux = interpol(obs[1,*], obs[0,*], lambda) * 1e-7
      end
      else: stop
  endcase

  ;; Fnu?

  if keyword_set(FNU) then flux = flux * lambda^2 * 1e-8 / c.c

  return, flux
end

