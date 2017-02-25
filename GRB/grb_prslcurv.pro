;+ 
; NAME:
;  grb_prslcurv
;   Version 1.1
;
; PURPOSE:
;    Parses a data file of GRB light curves and returns a structure
;    which parameterizes it.
;
; CALLING SEQUENCE:
;
; INPUTS:
;     fil -- Data file of light curves
;
; RETURNS:
;   struct -- Structure describing afterglow light curves
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NAM=  -- Name of GRB of interest.  Leave blank to get an array of
;           structures for the various light curves.
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
;   17-Feb-2006 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function grb_prslcurv, fil, NAM=nam

;
;  if  N_params() LT 2  then begin 
;    print,'Syntax - ' + $
;             'grb = grb_prslcurv()'
;    return, -1
;  endif 

  ;; Optional keywords
  if not keyword_set(FIL) then fil = getenv('XIDL_DIR')+'/GRB/grb_lightcrv.dat'

  readcol, fil, nam, z, alpha, beta, t0, mag, magwv, refnum, $
           FORMAT='A,F,F,F,F,F,F,A', /silen
  ngrb = n_elements(z)
  grb = replicate({grbafterglow}, ngrb)

  grb.nam = nam
  grb.z = z
  grb.alpha = alpha
  grb.beta = beta
  grb.t0 = t0
  grb.mag = mag
  grb.mag_wv = magwv
  
  if not keyword_set(REFFIL) then reffil = $
    getenv('XIDL_DIR')+'/GRB/grb_lightcrv.ref'
  readcol, reffil, ref_idx, ref_str, format='L,A'

  ;; References
  for qq=0L,ngrb-1 do begin
      ;; Parse
      refn = long(strsplit(refnum[qq], ':', /extract))
      nref = n_elements(refn)
      grb[qq].refnum[0:nref-1] = refn
      ;; Loop
      for ii=0L,nref-1 do begin
          mt = where(ref_idx EQ refn[ii],nmt)
          if nmt NE 1 then stop
          grb[qq].ref[ii] = ref_str[mt[0]]
      endfor
  endfor


  return, grb
end
