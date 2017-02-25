;+ 
; NAME:
; x_mknewlls
;    Version 1.1
;
; PURPOSE:
;   Creates the absorption line list for an LLS
;
; CALLING SEQUENCE:
;   
;   lines = x_mknewlls(zabs, NABS=, BVAL=, SET=)
;
; INPUTS:
;  zabs -- Redshift of the LLS
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NABS -- N(HI) value (cm^-2) [default: 18.]
;  BVAL -- b value (km/s) [default: 30]
;  SET  -- Index of the LLS lines (groups them together)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lines = x_mknewlls(2., NABS=17.3) 
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_mknewlls, zabs, NABS=nabs, BVAL=bval, SET=set

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'lines = x_mknewlls(zabs, NAbs=, b=, set=) [v1.0]'
    return, -1
  endif 

; Optional keywords
  if not keyword_set(BVAL) then bval = 30.
  if not keyword_set(Nabs) then nabs = 18.0
  if not keyword_set(SET) then set = 0

  tmp1 = { newabslinstrct }

  ;; Open Line list
  if not keyword_set( tlist ) then $
    tlist = getenv('XIDL_DIR')+'/LLS/Lines/LLS_std.lst'
  readcol, tlist, wr, skipline=1, /sile

  lines = replicate(tmp1, n_elements(wr))

  ;; Get wrest, f, gamma
  for i=0L,n_elements(wr)-1 do begin
      tmp = x_setline(wr[i])
      lines[i] = tmp
  endfor

  ;; All values
  lines.N = Nabs
  lines.set = set
  lines.b = bval
  lines.zabs = zabs
  lines.ion = 'H I'

  return, lines

end
