;+ 
; NAME:
; x_mkslitqa
;     Version 1.1
;
; PURPOSE:
;   Creates the QA structure used to save info from the slitflat
;   routine for final printing via qa_profile
;    
; CALLING SEQUENCE:
;  qa = x_mkslitqa( x, z, mask, ordr, lower=, upper= )
;
; INPUTS:
;  x -- slit_frac
;  z -- Profile
;  mask -- Profile mask
;  ordr -- Order structure (parsed to include only one order)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; LOWER=
; UPPER=
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;     temp_qa = x_mkslitqa(slit_frac[xsort], profile[xsort], profile_mask, $
;                    ordr_str[q])
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2004 Written by SB
;   Feb-2005  Ported to XIDL by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_mkslitqa, x, z, mask, ordr, lower=lower, upper=upper, $
                     per05=per05,  per95=per95
      
  if not keyword_set(PER05) then per05 = 0.05
  if not keyword_set(PER95) then per95 = 0.95
   n = n_elements(ordr[0].profile0)
   if NOT keyword_set(lower) then begin
     lower = findgen(n)/100. - 1.255
     upper = findgen(n)/100. - 1.245
   endif
   cen = 0.5*(lower+upper)

   mid   = ordr[0].profile0
   edgel = ordr[0].profile0-ordr[0].profile1
   edger = ordr[0].profile0+ordr[0].profile1
  
   tt = { order : ordr[0].order, $
          cen   : cen, $
          mid   : mid, $
          edgel : edgel, $
          edger : edger, $
          median: lower*0.0, $
          p05   : lower*0.0, $
          p95   : lower*0.0, $
          ab    : [1.0, 0.0], $
          Npix  : 0L,  $
          Nrej  : 0L, $
          chi2  : 0.0  } 
     
   for i=0, n-1 do begin
     any = where(x GE lower[i] AND x LT upper[i], nany)
     if nany GT 10 then begin
       md = median(z[any]) 
       s = sort(z[any])
       tt.median[i] = md
       tt.p95[i] = (z[any[s]])[long(nany*per95)]
       tt.p05[i] = (z[any[s]])[long(nany*per05)]
     endif
   endfor

return, tt
end


