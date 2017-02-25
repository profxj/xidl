;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; search_overlap.pro        
; Author: Kathy Cooksey                      Date: 14 Nov 2005
; Project: HST Archive Metal-Line System Survey with 
;          Jason X. Prochaska
; Description: Accepts structure of the form as those created
;              in searchspec, and attempts to sort through
;              features that may have been detected and fit
;              several times and to divide regions with
;              several features in to distinct ones
; Input: 
;   strct -- searchspec formatted structure
; Optional:
;   /divid -- try to divide regions
;   /debug -- stops at select points to aid in debugging
; Output: 
;   returns: array of indices of sorted and separated features
; Optional output:
;   elim= -- return indices of those features eliminated 
; Example:
;   rslts = search_overlap(strct,/divid,elim=elim)
; History:
;   17-Nov-05  created by KLC
;   15-Dec-05  cleaned up search_best logic order and flags;
;              added search_chkbnd, which does a lot
;   15-Oct-10  Renamed to civ_search_...
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function civ_search_best, substrct, flag=flag, ngood=ngood, good=good
;;flag: (1) chi^2 ~ 1; (2) if substrct.flag all 2, take largest;
;;(4) chi^2 != -9.99; (8) whichever substrct.flag = 4 (Gaussian fit)
num = n_elements(substrct)
mask = replicate(1b,num)
best = -1


;; flag == 4 better than flag == 2 
bttr = where(substrct.flg eq 4 and mask eq 1b,nbttr,$
             complement=bd,ncomplement=nbd)
if nbttr eq 1 then begin
    best = bttr[0]
    flag = 8
    goto, done
endif
if nbd ne 0 and nbd ne num then mask[bd] = 0b

gd = where(mask eq 1b,ngd)
if ngd eq 0 then $
  stop,'search_overlap: problematic result--out of possibilities'

;; chi^2 != -9.99, means only one Gaussian actually fitted
bst = where((-9.99-substrct.chisq) lt 0.0001 and mask eq 1b,nbst,$
            complement=bd,ncomplement=nbd)
if nbst eq 1 then begin
    best = bst[0]
    flag = 4
    goto, done
endif 
if nbd ne 0 and nbd ne ngd then mask[bd] = 0b

gd = where(mask eq 1b,ngd)
if ngd eq 0 then $
  stop,'search_overlap: problematic result--out of possibilities'


;; if some completely inside one, exclude smaller (only if boxcar)
dum = where(substrct.flg eq 2 and mask eq 1b,ndum)
if ndum eq ngd then begin
    wvmx = max(substrct[dum].wv_lim[1],iwvmx)
    wvmn = min(substrct[dum].wv_lim[0],iwvmn)
    if iwvmx eq iwvmn then begin
        best = dum[iwvmx]
        flag = 2
        goto, done
    endif 
endif 
    

;; chi^2 ~ 1. 
mn = min(abs(substrct[gd].chisq-1.),imn)
best = gd[imn]
flag = 1                     ;viable ones

done:

;;Only if resorted to chi^2 are there other possibilities
;;Equivalent to rechecking mask
if flag eq 1 then good = gd else good = best
ngood = n_elements(good)

return,best
end                             ;function search_best


function search_group, substrct, nsets,debug=debug
sets = replicate(0,n_elements(substrct))
wvmin = min(substrct.wv_lim[0],iwvmn)
wvmax = max(substrct.wv_lim[1],iwvmx)
diff = substrct.wv_lim[1]-substrct.wv_lim[0]
;gd = max(diff,ngd)
;sets[ngd] = -1                  ; overarching group --> straddle group

wvlomx = max(substrct.wv_lim[0],iwvlomx)
wvhimn = min(substrct.wv_lim[1],iwvhimn)
if keyword_set(debug) then stop,'search_group debug: set range'
if wvlomx ge wvhimn then begin ;>= 2 distinct groups
    gd = where(substrct.wv_lim[0] ge wvmin and substrct.wv_lim[1] le wvhimn,$
               ngd)             ; lower group
    sets[gd] = 1    
    mx1 = max(diff[gd],imx1)
    mx1 = gd[imx1]


    gd = where(substrct.wv_lim[0] ge wvlomx and substrct.wv_lim[1] le wvmax,$
               ngd)             ; upper group
    sets[gd] = 2               
    mx2 = max(diff[gd],imx2)
    mx2 = gd[imx2]


    gd = where(substrct.wv_lim[1] ge wvlomx and $
               substrct.wv_lim[0] le wvhimn,ngd) ; straddle group
    if ngd ne 0 then sets[gd] = 3


    ;; leftovers
    gd = where(substrct.wv_lim[0] ge substrct[mx1].wrest and sets eq 0,ngd)
    if ngd ne 0 then sets[gd] = 2 ; more 'upper' like
    gd = where(substrct.wv_lim[1] le substrct[mx2].wrest and sets eq 0,ngd)
    if ngd ne 0 then sets[gd] = 1 ; more 'lower' like
    
    if keyword_set(debug) then stop,'search_group debug: grouped'
endif 

nsets = max(sets)
if keyword_set(debug) then stop,'search_group debug: exiting'
return,sets
end                             ;function search_group


function search_chkbnd,strct,ngd,nbd=nbd,bd=bd
;;Require that the wavelength bounds be at least one sigma (b/sqrt(2))
;;from the centroid
dwv1 = strct.wv_lim[1]-strct.wrest
dwv0 = strct.wrest-strct.wv_lim[0]
oneb = strct.doppler/(2.998e5*sqrt(2)) * strct.wrest

gd = where(dwv1 ge oneb and dwv0 ge oneb,ngd,complement=bd,ncomplement=nbd)

return,gd

end                             ;function search_chkbnd



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function civ_search_overlap, strct, divid=divid, elim=elim, debug=debug

srt = sort(strct.wrest)
auto = strct[srt]

nauto = n_elements(auto)
mask = replicate(0b,nauto)
det = -1
ndet = 0

gd = search_chkbnd(auto,ngd,nbd=nbd,bd=bd)
if ngd eq 0 then stop,'search_overlap: no good lines'
mask[gd] = 1b

done = 0
ii = gd[0]

while not done do begin
    ;;Centroid within limits of another
    gd = where(auto.wrest ge auto[ii].wv_lim[0] and $
               auto.wrest le auto[ii].wv_lim[1] and mask eq 1b,ngd)
    if keyword_set(debug) then stop,'search_overlap: begin loop ',ii

    case ngd of 
        0: stop,'search_overlap: problematic result--no match!'
        1: begin
            det = [det,gd[0]]
            ndet++
            mask[gd[0]] = 2b
        end
        else: begin
            ;;Ambiguous, explore further
            best = civ_search_best(auto[gd],ngood=ngood,good=good)
            if keyword_set(divid) and ngood gt 1 then begin
                ;;Attempt to divide further
                grp = search_group(auto[gd[good]],ngrp)
                if ngrp ne 0 then begin
                    for jj=1,max(grp) do begin
                        gg = where(grp eq jj,ngg)
                        if ngg gt 1 then begin
                            tmp = min(abs(auto[gd[good[gg]]].chisq-1.),igg)
                            igg = gd[good[gg[igg]]]
                        endif else igg = gd[good[gg[0]]]
                        det = [det,igg]
                        ndet++
                        mask[igg] = 2b
                    endfor 
                    ;; Take best of both groups
                    ;g1 = where(grp eq 1,ng1)
                    ;if ng1 gt 1 then begin
                    ;    tmp = min(abs(auto[gd[good[g1]]].chisq-1.),ig1) 
                    ;    ig1 = gd[good[g1[ig1]]]
                    ;endif else ig1 = gd[good[g1[0]]]
                    ;    
                    ;g2 = where(grp eq 2,ng2)
                    ;if ng2 gt 1 then begin
                    ;    tmp = min(abs(auto[gd[good[g2]]].chisq-1.),ig2) 
                    ;    ig2 = gd[good[g2[ig2]]]
                    ;endif else ig2 = gd[good[g2[0]]]
                    ;
                    ;; Straddle group (need check)
                    ;g3 = where(grp eq 3, ng3)
                    ;if ng3 gt 1 then begin
                    ;    tmp = min(abs(auto[gd[good[g3]]].chisq-1.),ig3)
                    ;    ig3 = gd[good[g3[ig3]]]
                    ;endif else ig3 = gd[good[g3[0]]]
                    ;
                    ;det = [det,ig1,ig2]
                    ;ndet = ndet + 2
                    ;mask[[ig1,ig2]] = 2b
                endif else begin ; end possible to group
                    ;;Just take the indistinguishable best (for the
                    ;;default, should make gd[good] -> gd[best])
                    det = [det,gd[good]]
                    ndet = ndet + ngood
                    mask[gd[good]] = 2b
                endelse         ;end not possible to group           
                if keyword_set(debug) then stop,'search_overlap debug: '+$
                  'several in overlap region with good chi^2'
            endif else begin    ;end /divid and necessary to divide
                ;;Take best
                det = [det,gd[good]]
                ndet = ndet + ngood
                mask[gd[good]] = 2b
            endelse 
        end
    endcase
    if keyword_set(debug) and ngd gt 1 then $
      stop,'search_overlap debug: end loop ',ii

    gd = where(mask eq 1b,ngd)
    if ngd eq 0 then done = 1 else ii = gd[0]
endwhile ;done

;;clean up
det = det[1:ndet]
elim = where(mask ne 2b,nelim)
if nelim ne 0 then elim = srt[elim]

;;Would need to sort
;test = where(mask eq 2b)
;diff = det-test
;test = where(diff ne 0,ntest)
;if ntest ne 0 then stop,'search_overlap: discrepency'

if keyword_set(debug) then stop,'search_overlap debug: last stop'

return,srt[det]
end
