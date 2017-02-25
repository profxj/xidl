;+ 
; NAME:
; lls_calcabnd
;  V1.1
;
; PURPOSE:
;    Calculates the abundance relative to the Sun, [X/Y]
;
; CALLING SEQUENCE:
;   lls_calcabnd, lls, nn, sys, X, Y, ans, sig
;
; INPUTS:
;   lls -- DLA structure array
;   nn  -- Index of the structure
;   X   -- Atomic number of first element
;   Y   -- Atomic number of second element
;
; RETURNS:
;
; OUTPUTS:
;  ans --  [X/Y]
;  sig --  error in [X/Y]
;
; OPTIONAL KEYWORDS:
;  /NOSIGY -- Do not include error in Y in the calculation
;     This is only useful when dealing with [X/H]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lls_allabd, lls
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   21-Jul-2008 Written by JXP
;- 
;------------------------------------------------------------------------------
pro lls_calcabnd, lls, nn, sys, X, Y, ans, sig, NOSIGY=nosigy

  if (N_params() LT 5) then begin 
    print,'Syntax - ' + $
             'lls_calcabnd, lls, nn, sys, X, Y, ans, sig, /NOSIGY [v1.1]'
    return
  endif 


  ;; Get abundances
  getabnd, nm, X, Xabnd, flag=1
  getabnd, nm, Y, Yabnd, flag=1

  ;; Calculate [X/Y]

  if X NE 1 then begin
      logX = alog10(lls[nn].elm[X].clm)
      logXsig = sqrt(((1./(alog(10.0)*lls[nn].elm[X].clm))^2)* $
                     lls[nn].elm[X].sigclm^2)
  endif else begin
      logX = lls[nn].elm[X].clm
      logXsig = lls[nn].elm[X].sigclm
  endelse
  if Y NE 1 then begin
      logY = alog10(lls[nn].elm[Y].clm)
      logYsig = sqrt(((1./(alog(10.0)*lls[nn].elm[Y].clm))^2)* $
                     lls[nn].elm[Y].sigclm^2)
  endif else begin
      logY = lls[nn].elm[Y].clm
      logYsig = lls[nn].elm[Y].sigclm
  endelse

  ans = logX - logY - Xabnd + Yabnd
  
  ;;     Error
  
  if keyword_set(NOSIGY) then logYsig = 0.
  sig = sqrt(logXsig^2 + logYsig^2)


  return
end

