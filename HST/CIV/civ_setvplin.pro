;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_setvplin.pro             
; Author: Kathy Cooksey                      Date: 23 Mar 2008
; Project: MLSS: HST survey of CIV with Xavier Prochaska
; Description: From civcandstrct, create array of structures
;              to be used with x_voigt (i.e. x_setline)
; Input: 
;   strct_fil - randomly generated civcandstrct (from 
;               civ_genprof)
; Optional:
;   dvabs - apply civcandstrct.zsig (which is components)
; Output: 
;   vplin - array of x_voigt-compatible structures (either
;           returned as structure or as FITS file)
; Example:
;   civ_setvplin,'randciv.fits',vplin,/dvabs
; History:
;   23 Mar 2008 - created by KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_setvplin,strct_fil,vplin,dvabs=dvabs
  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_setvplin: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent,structyp='civcandstrct') 
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  ;; Instantiate VPLIN struct and expand to remaining ions
  for ii=0,nstrct-1 do begin
     gd = where(strct[ii].wrest gt 0.,ngd)
     if ngd eq 0 then stop,'civ_setvplin: no elements'

     tmp = x_setline(strct[ii].wrest[gd[0]]) ;includes fval
     lin = replicate(tmp,ngd)
     for jj=1,ngd-1 do lin[jj] = x_setline(strct[ii].wrest[gd[jj]])
     lin.zabs = strct[ii].zabs[gd]
     if keyword_set(dvabs) then lin.zabs = lin.zabs + strct[ii].zsig[gd]

     lin.n = strct[ii].ncolm[gd]
     lin.b = strct[ii].b[gd]
     if keyword_set(allvp) then allvp = [allvp,lin] $
     else allvp = lin
  endfor                        ; nstrct loop

  ;; Return info
  if size(vplin,/type) eq 7 then begin
     mwrfits,allvp,vplin,/create,/silent
     print,'civ_setvplin: created ',vplin
  endif else vplin = allvp

end
