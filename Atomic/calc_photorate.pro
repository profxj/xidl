;+ 
; NAME:
; calc_photorate
;   Version 1.1
;
; PURPOSE:
;    Calculates the photoionization rate given a radiation field and
;      energy range to integrate across
;
; CALLING SEQUENCE:
;   
; rate = calc_photorate( Z, ion, Erange, uvfield )
;
; INPUTS:
;   Z       - Atomic number
;   ion     - Number of free electrons
;   Erange  - Energy range of photons (eV)
;   uvfield - Flag for prescription to calcualte the UV field
;                1: GPW80
;
; RETURNS: 
;  rate -- Photoionization rate (s^-1)
;   
; OUTPUTS:
;   
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
;   5-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function calc_photorate, Z, ion, Erang, uvfield, NTAB=ntab, PLOT=plot

  if  N_params() LT 4  then begin 
      print, 'Syntax - ' +$
        'rate = calc_photorate( Z, ion, Erng (eV), uvfield) (v1.1)'
      return, -1
  endif 
;
  ;; Optional keywords
  if not keyword_set(NTAB) then ntab = 10000L

  ;; Constants
  c = x_constants()

  ;; Create E range 
  if keyword_set(LOG) then $
    E_array = Erang[0] * exp( alog(Erang[1]/Erang[0]) * dindgen(ntab+1)/ntab)$
  else $
    E_array = Erang[0] + (Erang[1]-Erang[0])*dindgen(ntab+1)/ntab

  nu = E_array * c.eV / c.h

  ;; UV Flux
  case UVFIELD of 
      1: begin ; GPW80
          lambda = c.c/nu * 1e8 ; Ang
          uvfield = ism_uvradfield(lambda, 1, /fnu)
      end
      else: stop
  endcase
          
  ;; Photoionization cross-section
  phcross = x_photocross( Z, ion, E_array )

  ;; 
  integr = phcross * uvfield / (E_array*c.eV)
  if keyword_set(PLOT) then begin
      yplt = integr
      gd = where(integr GT 0., complement=bad, ncomplement=nbad)
      yplt[gd] = alog10(integr[gd])
      if nbad NE 0 then yplt[bad] = -50.
      x_splot, E_array, yplt, /bloc
  endif
  rate = int_tabulated(nu, integr, /double)

  return, rate
end
