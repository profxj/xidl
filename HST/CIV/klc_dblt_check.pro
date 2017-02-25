;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; dblt_check.pro                
; Author: Kathy Cooksey                      Date: 29 Sep 2005
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Blind search (e.g. no Lya necessary) for given
;              doublet in a searchspec structure
; Input: 
;   auto - searchspec structure 
;   pair - doublet structure, as given by dblt_retrieve
; Optional:
;   debug - stop occasionally
;   zqso - redshift of quasar, constrains search
;   srchewlim = instrument file name; place upper limits on pair
;               and append to auto
; Output: 
; Example:
; History:
;   20 Mar 2006 - working version, KLC
;   22 Jan 2007 - added /span
;   17 May 2007 - Match to same instrument first, search with
;                 little leniency around expected range
;   20 Dec 2007 - Only mask out one per loop; 
;                 match only on observed wavelength
;    7 Feb 2008 - Add srchewlim
;   22 Feb 2008 - search for pair I inside pair II as well
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function klc_dblt_check,auto,pair,debug=debug,zqso=zqso,zgal=zgal,$
                    srchewlim=srchewlim

;;Parameters
if not keyword_set(zqso) then zqso = 100. ;ridiculous limit
if not keyword_set(zgal) then zgal = 0.
if keyword_set(srchewlim) then begin 
    readcol,srchewlim,spec,dw,w0,format='a,f,f',/silent
    ninstr = n_elements(spec)
    wvlim = fltarr(ninstr,2)
    snr = fltarr(ninstr)
    for ii=0,ninstr-1 do begin
        test = file_search(spec[ii],count=ntest)
        if ntest eq 0 then begin
            wvlim[ii,*] = [-1,-1]
            continue
        endif 

        if ii le 6 then $       ;FUSE
          fx = x_readspec(spec[ii],wav=wv,inflg=3) $
        else fx = x_readspec(spec[ii],wav=wv) ;stis
        gd = where(fx ne 0,ngd)
        wv = wv + dw[ii]*w0[ii]/wv

        wvlim[ii,*] = [min(wv[gd],max=mx),mx]
        if ii lt 7 then gd = measureSN(spec[ii],inflg=3,/silent) $
        else begin              ;HST
            espec = strmid(spec[ii],0,strpos(spec[ii],'f.fits'))+'e.fits'
            gd = measureSN(spec[ii],err=espec,/silent)
        endelse 
        snr[ii] = gd[2]
    endfor
    autoxtr = 0
endif 

;;Constants
sol = 2.998e5                   ;km/s
rt2 = sqrt(2)
rtpi = sqrt(!pi)
nauto = n_elements(auto.wrest)

;;Variables
doublet = {ion:pair.ion,$
           flag:0, $
           zabs:0d, $ 
           ew:0., $             ;EW ratio (ewI/ewII)
           b:0., $             ;doppler parameter ratio (dopbI/dopbII)
           dv:0. $          ;velocity separation ratio (dv/dvexpected)
          }
dblts = replicate(doublet,nauto)

;;Define region where possible to detect given pair (and must give
;;buffer region for poor wavelength centroids
;;Search for 1000 km/s less than z = 0
gdI = where(auto.wrest ge pair.wvI*(1.+zgal) and $
            auto.wrest le pair.wvI*(1.+zqso),ngdI)

if ngdI eq 0 then begin
    print,'dblt_check: '+pair.ion+' not detectable'
    if keyword_set(debug) then stop
    return,dblts
endif 

;; assumed redshift
dblts[gdI].zabs = auto[gdI].wrest/pair.wvI - 1. 
;; expected wavelength of pair II
wvIIexp = fltarr(nauto)
wvIIexp[gdI] = pair.wvII*(1.+dblts[gdI].zabs) 
;; delta(wavelength) limits for pair I (maybe asymmetric)
;dwvI = fltarr(nauto,2)
;dwvI[gdI,*] = [auto[gdI].wrest-auto[gdI].wv_lim[0],$
;               auto[gdI].wv_lim[1]-auto[gdI].wrest]
;; wavelength bounds for pair II
wv_limII = fltarr(nauto,2)
;wv_limII[gdI,*] = [wvIIexp[gdI]-dwvI[gdI,0],$
;                   wvIIexp[gdI]+dwvI[gdI,1]]
wv_limII[gdI,0] = pair.wvII*auto[gdI].wv_lim[0]/pair.wvI
wv_limII[gdI,1] = pair.wvII*auto[gdI].wv_lim[1]/pair.wvI

mask = replicate(0,nauto)

;;Loop over potential pair I
for ii=0,ngdI-1 do begin
    jj = gdI[ii]
    ;;Separation within tolerance and in same instrument
    ;; Pair II wobs inside pair I limits 
    ;; OR pair II limits enclosing pair I redshift/wobs
    gdII = where(((auto.wrest ge wv_limII[jj,0] and $
                   auto.wrest le wv_limII[jj,1]) or $
                  (pair.wvI*auto.wv_lim[0]/pair.wvII le auto[jj].wrest and $
                   pair.wvI*auto.wv_lim[1]/pair.wvII ge auto[jj].wrest)) and $
                 auto.instr eq auto[jj].instr and mask eq 0, ngdII)
    case ngdII of
        0: begin                
            ;; Matches perhaps in diff channels
            gdII = where(((auto.wrest ge wv_limII[jj,0] and $
                           auto.wrest le wv_limII[jj,1]) or $
                          (pair.wvI*auto.wv_lim[0]/pair.wvII le auto[jj].wrest and $
                           pair.wvI*auto.wv_lim[1]/pair.wvII ge auto[jj].wrest)) and $
                         mask eq 0, ngdII)

            ;; One one match in another instrument
            if ngdII eq 1 then begin
                mask[gdII] = mask[gdII] + 1
                dblts[jj].flag = gdII[0] ;1b 
            endif 

            if ngdII gt 1 then begin
                ;; Multiple matches in other instruments
                dum = min(abs(wvIIexp[jj]-auto[gdii].wrest),igdII)
                gdII = gdII[igdII]
                dblts[jj].flag =  gdII[0] ;1b ;yes
                mask[gdII] = mask[gdII] + 1
            endif 
        end                     ;no match in instrument
        1: begin
            mask[gdII] = mask[gdII] + 1
            dblts[jj].flag = gdII[0] ;1b ;; One match in current instrument
        end 
        else: begin             
            ;; Multiple matches in current instrument, perhaps
            ;; puzzling
            dum = min(abs(wvIIexp[jj]-auto[gdII].wrest),igdII)
            gdII = gdII[igdII]
            dblts[jj].flag = gdII[0] ;1b ;yes
            mask[gdII] = mask[gdII] + 1
        end                     ;more than one match in instrument
    endcase 
    
    if dblts[jj].flag ne 0 then begin
        dblts[jj].ew = auto[jj].ew/auto[gdII].ew
        dblts[jj].b = auto[jj].doppler/auto[gdII].doppler
        dblts[jj].dv = sol*(auto[gdII].wrest/wvIIexp[jj] - 1.)
    endif 

    ;;;;;;;;;;;;;;;;;;;;;;;
    ;; Include upper limits
    ;;;;;;;;;;;;;;;;;;;;;;;
    if dblts[jj].flag eq 0 and keyword_set(srchewlim) then begin
        mtch = where(wvlim[*,0] le wv_limII[jj,0] and $
                     wvlim[*,1] ge wv_limII[jj,1],nmtch)
        if nmtch eq 0 then goto,skip
       
        if nmtch gt 1 then begin
            mx = max(snr[mtch],imx) ;use highest S/N spectra
            mtch = mtch[imx] 
        endif else mtch = mtch[0]

        if mtch lt 7 then $     ;FUSE
          fx=x_readspec(spec[mtch],inflg=3,wav=wv,sig=er) $
        else begin              ;STIS
            espec = strmid(spec[mtch],0,strpos(spec[mtch],'f.fits'))+'e.fits'
            fx=x_readspec(spec[mtch],wav=wv,sig=er,fil_sig=espec)
        endelse 
        dwv = wv - shift(wv,1)
        dwv[0] = dwv[1]

        ;; Identify expected pair II region
        gd = min(wv-wv_limII[jj,0],imn,/absolute)
        gd = min(wv-wv_limII[jj,1],imx,/absolute)

        ;; Modify 1-flux when flux < -|error| and flux > 1+|error|
        rng = lindgen(imx-imn+1)+imn
        nfx = where(fx[rng] lt -abs(er[rng]),nnfx)
        pfx = where(fx[rng] gt 1.+abs(er[rng]),npfx)
        modfx = 1.-fx[rng]
        if nnfx ne 0 then modfx[nfx] = 1+2*abs(er[rng[nfx]])
        if npfx ne 0 then modfx[pfx] = -2*abs(er[rng[pfx]])
        ew = total(modfx*dwv[rng])
        sigew = sqrt(total((er[rng]*dwv[rng])^2))
;        if abs(ew) lt sigew then goto,skip ;non-detection

        ;; Test if ratio within reason
        rtoexp = (pair.wvI*pair.fI)/(pair.wvII*pair.fII)
        rto = auto[jj].ew/ew
        sigrto = abs(rto)*sqrt((sigew/ew)^2+(auto[jj].sigew/auto[jj].ew)^2)
        if rto lt 1.-sigrto or rto gt rtoexp+sigrto then goto,skip

        ;; Load information to append to auto
        tmp = auto[jj]          ;preserves redshift, misc. searchspec stuff
        tmp.wrest = wvIIexp[jj]
        tmp.wv_lim = wv_limII[jj,*]
        tmp.ew = ew
        tmp.sigew = sigew
        tmp.instr = 2^mtch
        tmp.flg = 5

        if keyword_set(autoxtr) then autoxtr = [autoxtr,tmp] $
        else autoxtr = tmp


        ;; Load information to pass with dblts structure
        dblts[jj].flag = nauto + n_elements(autoxtr) - 1
        dblts[jj].ew = auto[jj].ew/ew
        dblts[jj].b = 0.
        dblts[jj].dv = 0.

        skip:
    endif                       ;keyword set srchewlim


    if keyword_set(debug) then begin
        print,string('pair:',dblts[jj].zabs,auto[jj].wrest,$
                     auto[gdII].wrest,format='(a5,2x,f8.6,2(2x,f9.4))')
        if size(debug,/type) eq 7 then $
          x_specplot,debug,zin=dblts[jj].zabs,/lls,/block
    endif                       ;debug
    
endfor                          ;loop over gdI

;; Modify input auto structure
if keyword_set(srchewlim) and keyword_set(autoxtr) then begin
   print,'dblt_check: returning appended auto.-search structure'
   auto = [auto,autoxtr]
endif 

bd = where(mask gt 1,nbd)
if nbd ne 0 then stop,'dblt_check: mask element > 1'

if keyword_set(debug) then stop,'dblt_check debug: finished'

return, dblts
end
