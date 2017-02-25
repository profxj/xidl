;+ 
; NAME:
; x_poisscl  
;    Version 1.1
;
; PURPOSE:
;  Given the number of objects obsrved and confidence limits, this
;  program returns the values corresponding to the confidence limits.
;
; CALLING SEQUENCE:
;  val = x_poisscl( x, cl, SIGMA= )
;
; INPUTS:
;   x = Number of objects detected
;   [cl] =  Confidence limits
;
; RETURNS:
;  val == Best fit values
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   SIGMA= Number of sigma (in lieu of cl)
;  /VERBOSE 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   print, x_poisscl(40, sigma=2)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Dec-2004 Written by JXP
;   12-Apr-2012 More properly defining CL
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_poisscl, x, icl, SIGMA=sigma, NSPL=nspl, SILENT=silent, VERBOSE=verbose

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'prob = x_poisscl( x, [cl], SIGMA=, /VERBOSE {/SILENT} ) [v1.1]'
      return, 0.
  endif 
  if not keyword_set( NSPL ) then nspl = 1000L
  if keyword_set( SIGMA ) then begin
     cl = gauss_pdf(sigma)
  endif else begin
     CL = iCL + (1-iCL)/2.
  endelse
  if keyword_set(VERBOSE) then print, 'x_poisscl: Using a CL value of ', cl

  if not keyword_set(CL) then begin
     print, 'x_poisscl:  You need to set CL or SIGMA'
     return, 0.
  endif

  n = round(x) 
  val = dblarr(2)

  ;; Large value?  If so, then use Gaussian
  if n GT 120 then begin
      S = abs(gauss_cvf(cl))
      val[0] = float(n+1)*(1 - 1./(9*(n+1)) + S/(3*sqrt(n+1.)))^3
      val[1] = float(n) - S*sqrt(double(n)) + (S^2 - 1)/3.
      return, val
  endif

  ;; Positive limit
  ;; Normalize
;  norm = 1. - igamma(x+1,x)  ;; Note x! divided out

  ;; I think there is a bug in here for zero events!!  Check against Gehrels 1986 before proceeding
  if not keyword_set( MAXX ) then maxx = (n+10.*sqrt(n)) > 10.
  lup = n + (maxx-n)*dindgen(nspl)/float(nspl)
  
  luparr = lup # replicate(1.d, n+1)
  narr = replicate(1., nspl) # dindgen(n+1)

  summarr = luparr^narr * exp(-1.d*luparr) / gamma(narr+1.d)
  if n GE 1 then totup = total(summarr,2) else totup = summarr

  if min(totup) GT (1.-cl) then stop

  ;; Spline
  splin = spl_init(lup, totup, /double)
  val[0] = x_fndspln(lup,totup,1.-cl,splin,TOLER=1d-3*(1-cl))

  if n LE 0 then return, val

  ;; Negative limit
  if not keyword_set( MINX ) then minx = (x - 10*sqrt(x)) > 1e-5
  llow = x - (x-minx)*findgen(nspl)/float(nspl)
  llowarr = llow # replicate(1.d, n)
  narr = replicate(1., nspl) # dindgen(n)

  summarr = llowarr^narr * exp(-1.d*llowarr) / gamma(narr+1)
  if n GT 1 then totlow = total(summarr,2) else totlow = summarr

  if max(totlow) LT cl then stop

  ;; Spline
  splin = spl_init(llow, totlow, /double)
  val[1] = x_fndspln(llow,totlow,cl,splin,TOLER=1d-5*cl)

  return, val

end
