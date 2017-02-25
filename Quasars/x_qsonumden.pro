;+ 
; NAME:
; x_qsonumden
;
; PURPOSE:
;    Gives #qsos per comoving Mpc-3 at a given redshift for 
;    a given magnitude limit (at 1450A).  Uses the LF from Hopkins et
;    al. 2007.  Assumes the input Magnitude is at 1450
;
; CALLING SEQUENCE:
;   
; INPUTS:
;   z    -- Redshift for evaluation
;   Mlim -- Limiting magnitude for the calculation
;
; RETURNS:
;  NUM_DEN -- Number density of quasars per Mpc^3 at z to Mlim
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   13-Feb-2008 Written by JXP
;-
;------------------------------------------------------------------------------
function x_qsonumden, z, Mlim, FLG=flg

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'numqso = x_qsonumden(z, Mlim, FLG=) [v1.0]'
    return, -1
  endif 

  if not keyword_set(FLG) then FLG = 0  ;; HRH07  (1=Willott)

  c = x_constants()
  num_den = dblarr(n_elements(Mlim))

  qso_lf = x_qsolumfunc(z, FLG, M1450_eval=M1450)

  ;; Convert M_B -> M_1450

  ;; Sum 
  dM = abs(median(M1450-shift(M1450,1)))  ;; I hope Mlim is regular!
  num_M = dM * total( qso_lf, /cumul)

  ;; Interpolate
  num_den = interpol(num_M, M1450, Mlim)
  ;x_psclose
  ;!p.multi=[0,1,1]
  ;stop

  return, num_den
end
