;+ 
; NAME:
; x_maxsngpow   
;    Version 1.1
;
; PURPOSE:
;  Run a Maximum Likelihood fit for a single power-law N^alpha
;
; CALLING SEQUENCE:
;  alpha = x_maxsngpow( arr, NMIN=nmin )
;
; INPUTS:
;   Array
;
; RETURNS:
;
; OUTPUTS:
;  val == Best fit values
;
; OPTIONAL KEYWORDS:
;  NOISE -- Helps with the KS test for data with discrete values
;           (logarithmic value only!)
;  NMIN= -- Minimum value for alpha
;  NMAX= -- Maximum value for alpha
;  CL=   -- Confidence limit for the value
;  KSPROB= -- KS propbability that this model is acceptible
;
; OPTIONAL OUTPUTS:
; ERR= -- Error in the best fitted value (corrsponds to CL)
;
; COMMENTS:
;
; EXAMPLES:
;   x_maxsngpow
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function x_maxsngpow_kscumf, x, EXTRA=extra

common x_maxsngpow_cmmn, sngpow_k, sngpow_maxb, sngpow_nmin

    return, sngpow_k / (sngpow_maxb+1) * (x^(sngpow_maxb+1.) - $
                                          sngpow_nmin^(sngpow_maxb+1.) ) 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_maxsngpow, arr, Nmin=nmin, SMM=smm, ERR=err, KSPROB=ksprob, $
                      PLOT=plot, NOISE=noise, NMAX=nmax, ARNG=arng, CL=cl

common x_maxsngpow_cmmn

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'alpha = x_maxsngpow( arr, NMIN=) [v1.1]'
      return, 0.
  endif 
  if not keyword_set(NMIN) then nmin = 1.

  narr = n_elements(arr)
  if not keyword_set( NMAX ) then begin
      smm =  total( alog( arr / nmin ) )
      maxb = -1.*(1 + narr / smm)
  endif else begin
      if not keyword_set(ARNG) then stop
      if not keyword_set(NSTEP) then nstep = 10000L
      aval = arng[0] + dindgen(nstep)*(arng[1]-arng[0])/float(nstep)
      lgL = aval*total(alog(arr)) - $
            narr*alog( (Nmax^(aval+1) - Nmin^(aval+1))/(aval+1))
      ;; Spline
      maxL = x_maxspln(aval, lgL, SPMX=maxb)
  endelse

  ;; Error Analysis
  if arg_present(ERR) then begin
      if not keyword_set(CL) then cl = 0.683
      if not keyword_set(NMAX) then begin
          if not keyword_set(ARNG) then stop
          if not keyword_set(NSTEP) then nstep = 10000L
          aval = arng[0] + dindgen(nstep)*(arng[1]-arng[0])/float(nstep)
          lgL = aval*total(alog(arr)) + narr*alog(-1.*aval-1.) - $
            narr*alog( Nmin^aval )
          maxL = maxb*total(alog(arr)) + narr*alog(-1.*maxb-1.) - $
            narr*alog( Nmin^maxb )
      endif
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd EQ 0 then stop
      ;; Integrate
      err = x_maxcl(aval[gd], lgl[gd], maxb, CL=cl)
  endif

  ;; KS Test
  if arg_present(KSPROB) and not keyword_set(NMAX) then begin
      sngpow_k = -1.* (maxb+1.) / Nmin^(maxb+1.)  ; Normalized to unity
      sngpow_maxb = maxb
      sngpow_nmin = nmin
      if keyword_set( NOISE ) then darr = arr * 10^(noise*randomu(-2442,narr)) $
      else darr = arr
      ;; 
      ksone, darr, 'x_maxsngpow_kscumf', d, ksprob;, PLOT=plot 
  endif

  return, maxb

end
