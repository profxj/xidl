;+ 
; NAME:
; cosm_time
;
; PURPOSE:
;   Calculates the age of the universe at an arbitrary redshift where
;  the interval is from z=0 to z_i.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   z_i  -- Redshift
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /INIT  -- Initializes the cosmology to default values (0.7, 0.3, 75)
;  /W06MAP  -- Initializes the cosmology to WMAP06 (0.72, 0.28, 73)
;  H0=    -- Hubbles constant (km/s/Mpc)
;  /START  -- Calculate the age from z=Infinity to z_i
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   age = cosm_time(1.)
;   
;
; PROCEDURES CALLED:
;   cosm_common
;
; REVISION HISTORY:
;   02-Jul-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function cosm_inttime, z
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

  distintg = 1./ cosm_hubble(z) / (1.+z)
  return, distintg
end

function cosm_b_inttime, v
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

  t = exp(v) - 1
  distintg = 1./ cosm_hubble(t)
  return, distintg
end

function cosm_time, z, H0=h0, INIT=init, _EXTRA=extra, START=start

  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'age = cosm_time(z, H0=, /INIT, /W06MAP, /START) [v1.1]'
    return, -1
  endif 

  if keyword_set( INIT ) or keyword_set(EXTRA) then $
    cosm_common, H0=h0, _EXTRA=extra
  if keyword_set( H0 ) then cosm_h = h0

  nz = n_elements(z)
  time = fltarr(nz)
  if not keyword_set(START) then begin
      for jj=0L,nz-1 do time[jj] = qromb('cosm_inttime', 0., z[jj]>0., /double) 
  endif else begin
      time = qromo('cosm_b_inttime', alog(z+1), /midexp, /double)
  endelse

      ;; Pass back
  if nz EQ 1 then time = time[0]

  c = x_constants()
  return, time*c.kpc*1e3/1e5/c.yr

end

