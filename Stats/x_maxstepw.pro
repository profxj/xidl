;+ 
; NAME:
; x_maxstepw
;    Version 1.1
;
; PURPOSE:
;  Run a Stepwise Maximum Likelihood analysis.  
;  NOT OPERATIONAL
;
; CALLING SEQUENCE:
;  x_maxstepw, arr, bins, guess, ans, ERR=err
;
; INPUTS:
;   Array
;   Bins -- Bin width (uniform value)
;   Guess
;
; RETURNS:
;
; OUTPUTS:
;  ans -- Answer
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_maxstepw
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function stepw_func


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_maxstepw, arr, bin, guess, ans, ERR=err, NITER=niter, BMNX=bmnx, $
                MINV=minv

  ; 
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_maxstepw, arr, bins, guess, ans, ERR=, [v1.1]'
      return
  endif 
  stop  ;; Work in progress!

  ;; Keywords
  if not keyword_set(NITER) then niter = 100L
  if not keyword_set(toelr) then toler = 1e-5
  if not keyword_set(MINV) then stop

  if not keyword_set(BMNX) then begin
      bmnx = dblarr(2)
      bmnx[0] = min(arr, max=mx)
      nbin = fix((mx-bmnx[0])/bin) + 1L
      bmnx[1] = bmnx[0] + nb*bin
  endif else nbin = round((bmnx[1]-bmnx[0])/bin) 
  if n_elements(guess) NE nbin then stop

  narr = n_elements(arr)

  ;; Setup
  val = bmnx[0] + findgen(nbin)*bin + bin/2.

  phi = guess
  for qq=0L,niter-1 do begin
      ;; Iterate on phik
      for jj=0L,nbin-1 do begin
          widx = where(abs(arr-val[jj]) LT bin/2.,nw)
          if nw EQ 0 then phi[jj] = 0. else begin
              num = float(nw)
              ;; The next is only good for DLA
              ;hidx = where(val[jj] GT MINV, nh)
              num2 = total( phi*bin )
              ;; Denom (only good for DLA)
              denom = float(narr)
              phi[jj] = num*num2/denom/bin
          endif
      endfor
  endfor
              
  return

end
