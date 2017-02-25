;+ 
; NAME: 
; eval_powerfn   
;    Version 1.1
;
; PURPOSE:
;    Evaluate f(N,X) from a power-law model
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;  log f(N,X)
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   May-2011 Written by JXP
;   Oct-2013 Added hook to Spline
;-
;------------------------------------------------------------------------------
function eval_powerfn, powerfn_strct, lgNval, z, VEL_ARRAY=vel_array, $
                       FNSTR=fnstr

  c = x_constants()

  ;; Grab the write member of the structure
  if not keyword_set(FNSTR) then begin
     fidx = where(powerfn_strct.zmnx[0] LT z and powerfn_strct.zmnx[1] GT z, nidx)
     if nidx EQ 0 then stop
     fnstr = powerfn_strct[fidx]
  endif

  ;; Spline?
  szs = size(fnstr, /str)
  case szs.structure_name of
     'SPLINEFNSTRCT': begin
        logfNX = eval_splinefn(fnstr, lgNval, z, vel_array=vel_array)
        return, logfNX
     end
     'POWERFNSTRCT': 
     else: stop
  endcase

  ;; Evaluate f(N,X)
  n_nval = n_elements(lgNval)

  ;; Find the right power laws
  sz = n_elements(fnstr.pivots)
  Nval_array = lgNval # replicate(1., sz) 
  Npiv_array = replicate(1., n_nval) # fnstr.pivots

  fpiv1 = Nval_array - Npiv_array
  fpiv2 = Nval_array - shift(Npiv_array,0,-1)

  flg = (fpiv1 GT 0.d) AND (fpiv2 LE 0.d)  ;; Could be round off troubles here
  idx = where(flg, nidx) 
  if nidx NE n_nval then stop

  idx2 = idx / n_nval  ;; This works, honest


  ;; Create a grid if desired
  if keyword_set(VEL_ARRAY) then begin
     z_val = z + (1+z) * vel_array/(c.c/1e5)
     lgNHI_grid = (fnstr.fn_pivot[idx2] + (lgNval-fnstr.pivots[idx2])*fnstr.beta[idx2]) # $
                replicate(1., n_elements(z_val))
     z_grid1 = 10.^( fnstr.gamma[idx2] # alog10(1+z_val))  ;; (1+z)^gamma
     z_grid2 = ((1./(1+fnstr.zpivot[idx2]))^fnstr.gamma[idx2]) # $
               replicate(1.,n_elements(z_val)) 
     logfNX = lgNHI_grid + alog10(z_grid1*z_grid2) 
  endif else begin
     ;; One redshift for evaluation
     logfNX = fnstr.fn_pivot[idx2] + (lgNval-fnstr.pivots[idx2])*fnstr.beta[idx2] + $
              fnstr.gamma[idx2] * alog10((1+z)/(1+fnstr.zpivot[idx2]))
  endelse
              

  return, logfNX
end
