;+ 
; NAME:
; x_slitprofile_return
;     Version 1.1
;
; PURPOSE:
;    Return a slit profile for a MIKE or HIRES order as a function of
;      slit position and row
;
; CALLING SEQUENCE:
;   
;    profile = x_slitprofile_return(slit_frac, y, ordr_p)  
;
; INPUTS:
;
;    slit_frac   : array corresponding to fractional position along slit,
;                  where -1,1 are left and right edges respectively.
;    y           : An array with corresponding row positions
;    ordr_p      : Order structure with profile information, p0, p1
;
; RETURNS:
; 
;    profile     : Slit profile for each pixel in the arrays slit_frac & y
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-2003     Written by SMB
;-
;------------------------------------------------------------------------------
function x_slitprofile_return, slit_frac, y, ordr_p

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'profile = x_slitprofile_return(slit_frac, y, ordr_p)  [v1.1]'
      return, -1
  endif 
  nrow = n_elements(ordr_p.lhedg)
  ynorm      = (2.0* y- nrow)/nrow
  
  
  p0 = ordr_p.profile0
  p1 = ordr_p.profile1
  
  if total(finite(p1) EQ 0) GT 0 then p1[*] = 0. 
  
  npoints = n_elements(p0)
  xp      = (npoints-1)/2 + 100.0 * slit_frac
  inter0 = interpolate(p0, xp)
  inter1 = interpolate(p1, xp)
  
  profile = (inter0 + inter1 * ynorm) * (xp GE 0 AND xp LE npoints-1)
  return, profile
end

