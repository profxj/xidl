;+ 
; NAME:
; x_poissdraw   
;    Version 1.1
;
; PURPOSE:
;   Pass back the number of objects given a random number
;  between 0 and 1 and an average number.
;
; CALLING SEQUENCE:
;  nobj = x_poissdraw( random, avg )
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_poisson
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-May-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_poissdraw, rand, avg

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'prob = x_poissdraw(rand, avg, TOL= ) [v1.1]'
      return, 0L
  endif 

  if not keyword_set(NMAX) then nmax = round(10.*avg) > 10

  nval = n_elements(rand)
  values = lonarr(nval)
  cum_tot = total(x_poisson(lindgen(nmax), avg), /cumul)

  pos = where(rand GT cum_tot[0], npos)
  for ss=0,npos-1 do begin
      idx = where(cum_tot GT rand[pos[ss]], nidx) 
      values[pos[ss]] = idx[0]
  endfor

  ;; Scalar?
  if nval EQ 1 then values = values[0]

  ;; Return
  return, values

end
