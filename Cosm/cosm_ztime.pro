;+ 
; NAME:
; cosm_ztime
;
; PURPOSE:
;   Calculates the redshift given the age (z=0 corresponds to t=0yr)
;
; CALLING SEQUENCE:
;  z  = cosm_ztime(t, /init)
;
; INPUTS:
;   t  -- Time
;
; RETURNS:
;   z  -- Redshift
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /INIT  -- Initializes the cosmology to the default values
;  H0=    -- Hubbles constant (km/s/Mpc)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   z = cosm_ztime(1e10)
;   
;
; PROCEDURES CALLED:
;   cosm_common
;
; REVISION HISTORY:
;   14-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function cosm_funcztime, z
common cosmolgy_cmmn, cosm_dm, cosm_k, cosm_h, cosm_Ob, cosm_L, cosm_r
common cosmolgy_fztime, age

  return, cosm_time(z) - age
end


function cosm_ztime, t, H0=h0, INIT=init, _EXTRA=extra

common cosmolgy_cmmn, cosm_dm, cosm_k, cosm_h, cosm_Ob, cosm_L, cosm_r
common cosmolgy_fztime, age

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'z = cosm_ztime(t, H0=, /INIT) [v1.1]'
    return, -1
  endif 

  if keyword_set( INIT ) or keyword_set(EXTRA) then $
    cosm_common, H0=h0, _EXTRA=extra
  if keyword_set( H0 ) then cosm_h = h0
  ;; Constants
  c = x_constants()

  nt = n_elements(t)
  z = fltarr(nt)
  for jj=0L,nt-1 do begin
      ;; Iterate
      age = t[jj];/c.kpc/1e3*1e5*c.yr
      z[jj] = fx_root([1.,50.,0.],'cosm_funcztime',/double)
  endfor

  if nt EQ 1 then return, z[0] else return, z

end

