;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_synthnoise.pro               
; Author: Kathy Cooksey                      Date: 15 Apr 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Randomly insert Voigt profiles into spectrum with
;              appropriate noise
; Input: 
;   inspec -- structure of original {wave,flux,error}
; Optional:
;   conti -- if array, then add noise to it, 
;            if /conti, then assume normalized (continuum = 1)
;            else, use original spectrum with absorption masked out
;   seed -- random number seed
;   /debug -- flag to stop at key moments
; Output: 
;   outspec -- new structure with Voigt profiles {wave,flux,error}
; Example:
;   civ_synthnoise,inspec,vpfx,outspec,/conti
; History:
;    6 Apr 2008 -- created by KLC, taken from spec_addlin
;   21 Apr 2008 -- mask out absorption features to preserve error
;   26 Aug 2008 -- streamline to be just adding noise
;   13 Oct 2008 -- Better check for flooring (compare to neighbors)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_synthnoise,inspec,varflr,vpfx,outspec,conti=conti,$
                   seed=seed,debug=debug
if (N_params() lt 4) then begin
    print,'Syntax - '+ $
      'civ_synthnoise,inspec,varflr,vpfx,outspec,[conti=conti,'+$
      '  seed=seed,debug=debug]'
    return
endif 
nwspec = inspec

npix = n_elements(inspec.wave)
if n_elements(vpfx) ne npix then $
  stop,'civ_synthnoise: Voigt profile and input spectrum different size'
wvlim = [min(inspec.wave), max(inspec.wave)]

;; Generate new error array (accounting for floor)
nwerr = sqrt(nwspec.error^2 - (nwspec.error^2-varflr)*(1-vpfx))

;; Prevent error being higher than neighbors
gap = where(vpfx lt 1.,ngap,complement=unabs,ncomplement=nunabs)
lo = where(gap-1 ne shift(gap,1),nlo)
hi = where(gap+1 ne shift(gap,-1),nhi)
if nlo ne nhi then stop,'civ_synthnoise: mismatch bounds'
bnd = lonarr(nlo,2)
bnd[*,0] = gap[lo]
bnd[*,1] = gap[hi]
for ii=0,nlo-1 do begin
   ;; Compare to neighboring error array pixels
   cent = round(mean(bnd[ii,*]))
   srt = sort(abs(cent-unabs))
   erunabs = median(nwerr[unabs[srt[0:10]]]) ; might consider using replaced error
;   erunabs = median(inspec.error[bnd[ii,0]:bnd[ii,1]])
   erabs = median(nwerr[bnd[ii,0]:bnd[ii,1]])
   ;; Return to original
   if erabs gt erunabs then $
      nwerr[bnd[ii,0]:bnd[ii,1]] = inspec.error[bnd[ii,0]:bnd[ii,1]]
endfor                          ; loop features

;; Generate new flux array (adding scatter as necessary)
if keyword_set(conti) then begin
    ;; Conti is an array
    if n_elements(conti) eq npix then $
      nwfx = vpfx * conti + randomu(seed,npix,/normal)*abs(nwerr) $
      ;; Assume normalized
    else nwfx = vpfx + randomu(seed,npix,/normal)*abs(nwerr) 
endif else nwfx = vpfx * nwspec.flux ;more scatter

;; Account for gaps and/or bad pixels
bd = where(nwspec.error eq -1.,nbd)
if nbd ne 0 then begin 
    nwerr[bd] = -1.
    nwfx[bd] = 0.
endif 
if keyword_set(debug) then stop,'civ_synthnoise debug: created new flux'

;; Create output structures
outspec = {wave:nwspec.wave,flux:nwfx,error:nwerr}

if keyword_set(debug) then stop,'civ_synthnoise debug: about to exit'
end
