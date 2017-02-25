;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_synthew.pro               
; Author: Kathy Cooksey                      Date: 29 Aug 2008
; Project: 
; Description: 
; Input: 
; Optional:
; Output: 
; Example:
; History:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_synthew,vp_fil,strct_fil,wave=wave,totew=totew

if size(strct_fil,/type) eq 7 then civcand = xmrdfits(strct_fil,1,/silent) $
else civcand = strct_fil
nciv = n_elements(civcand)
gd = where(civcand[0].wrest gt 0.,ngd)
unq = uniq(civcand[0].wrest[gd],sort(civcand[0].wrest[gd]))
nion = n_elements(unq)
ncomp = ngd/nion

if size(vp_fil,/type) eq 7 then begin
   if keyword_set(wave) then vpfx = xmrdfits(vp_fil,1,/silent) $
   else vpfx = x_readspec(vp_fil,wav=wave) 
endif else begin
   vpfx = vp_fil
   if not keyword_set(wave) then stop,'civ_synthew: must set wavelenght array'
endelse 

dwv = wave - shift(wave,1)
dwv[0] = dwv[1]

for ii=0,nciv-1 do begin
   if keyword_set(totew) then begin
      for jj=0,nion-1 do begin
         if ncomp ne nion then begin
            ;; Must determine wavelength limits
            indx = lindgen(ncomp)*nion+jj
            wv_lim = [min(civcand[ii].wv_lim[indx,0]),$
                      max(civcand[ii].wv_lim[indx,1])]
         endif else wv_lim = [civcand[ii].wv_lim[jj,0],civcand[ii].wv_lim[jj,1]]

         dum = min(wv_lim[0]-wave,pmin,/absolute) 
         dum = min(wv_lim[1]-wave,pmax,/absolute) 
         civcand[ii].ew[jj] = total((1-vpfx[pmin:pmax])*dwv[pmin:pmax])*$
                              1000./(1+civcand[ii].zabs[jj])
      endfor                    ; loop nion
   endif else begin
      for kk=0,ncomp-1 do begin
         dum = min(civcand[ii].wv_lim[kk,0]-wave,pmin,/absolute) 
         dum = min(civcand[ii].wv_lim[kk,1]-wave,pmax,/absolute) 
         civcand[ii].ew[kk] = total((1-vpfx[pmin:pmax])*dwv[pmin:pmax])*$
                              1000./(1+civcand[ii].zabs[kk])
      endfor                    ; loop ncomp
   endelse 
endfor                          ; loop candidates

if size(strct_fil,/type) eq 7 then mwrfits,civcand,strct_fil,/create,/silent $
else strct_fil = civcand
end
