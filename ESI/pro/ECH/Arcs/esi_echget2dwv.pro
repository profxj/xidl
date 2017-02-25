;+ 
; NAME:
; esi_echget2dwv
;     Version 1.1
;
; PURPOSE:
;    Grab a wavelength from the 2D solution
;
; CALLING SEQUENCE:
;   
;  wav = esi_echget2dwv(arc_str, yval, ordr)
;
; INPUTS:
;   arc_str -  The 2D fit structure
;   yval    -  Array of values along the order
;   ordr    -  Array of Order numbers
;
; RETURNS:
;
; OUTPUTS:
;  Array of wavelength values (logarithmic)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_fit2darc, esi, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_echget2dwv, arc_str, yval, ordr
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'wv = esi_echget2dwv(arc_str, yval, ordr) [v1.1]'
      return, -1
  endif

  ;; Create the arrays
  work2d = dblarr(n_elements(yval),arc_str.ny*arc_str.no)
  pix_nrm = 2. * (double(yval) - arc_str.nrm[0])/arc_str.nrm[1]
  t_nrm = 2. * (replicate(ordr,n_elements(yval)) $
                - arc_str.nrmt[0])/arc_str.nrmt[1]
  worky = fchebyshev(pix_nrm[*], arc_str.ny)
  workt = fchebyshev(t_nrm[*], arc_str.no)
  
  for i=0,arc_str.no-1 do begin
      for j=0,arc_str.ny-1 do begin
          work2d[*,j*arc_str.no+i] = worky[*, j] * workt[*,i]
      endfor
  endfor
  wav = pix_nrm*0.d
  wav[*] = work2d # arc_str.res

  return, wav

end

