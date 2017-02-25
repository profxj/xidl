;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_cleanspec.pro               
; Author: Kathy Cooksey                      Date: 26 Aug 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Remove strong absorption lines and compute flux
;              floor
; Input: 
;   inspec -- structure of original {wave,flux,error}
; Optional:
;   search_fil -- clean out features from searchspec
;   instr_flg -- instrument flag number 
;   seed -- random number seed
;   /debug -- flag to stop at key moments
; Output: 
;   outspec -- new structure with Voigt profiles {wave,flux,error}
; Example:
;   civ_cleanspec,inspec,outspec
; History:
;   26 Aug 2008 -- created by KLC, taken from civ_synthnoise
;    1 Sep 2008 -- added search_fil
;   16 Sep 2008 -- added instr_flg
;   13 Oct 2008 -- MC neighbors into search_fil results
;   12 Jan 2010 -- measure variance floor with error > 0.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@measureSN
pro civ_cleanspec,inspec,outspec,varflr=varflr,search_fil=search_fil,$
                  instr_flg=instr_flg,seed=seed,lls=lls,debug=debug
  if (N_params() lt 2) then begin
     print,'Syntax - '+ $
           'civ_cleanspec,inspec,outspec,[varflr=varflr,'+$
           '  seed=seed,debug=debug]'
     return
  endif 

  nwspec = inspec
  npix = n_elements(inspec.wave)
  wvlim = [min(inspec.wave), max(inspec.wave)]

  ;; Prepare for adding scatter to line-center saturated Voigt profiles
  ;; Isolate median flux
  snr = measureSN(inspec.flux,err=inspec.error,wave=inspec.wave,$
                  dev=dev,inflg=4,/silent)

  ;; Determine floor in higher S/N portion of spectrum (basically
  ;; exclude edges)
  lo = where(inspec.wave le wvlim[0]+0.25*(wvlim[1]-wvlim[0]))
  gd = where(inspec.flux[lo]/inspec.error[lo] ge snr[2])
  if gd[0] ne -1 then lo = median(inspec.wave[lo[gd]]) $
  else lo = min(inspec.wave)
  hi = where(inspec.wave ge wvlim[1]-0.25*(wvlim[1]-wvlim[0]))
  gd = where(inspec.flux[hi]/inspec.error[hi] ge snr[2])
  if gd[0] ne -1 then hi = median(inspec.wave[hi[gd]]) $
  else hi = max(inspec.wave)

  ;; Find bottom of strong absorption features in high S/N region
  trough = where(inspec.wave ge lo and inspec.wave le hi and $
                 (inspec.flux lt snr[1] or inspec.flux lt 0.05) and $
;                 inspec.error ne -1.,ntrough)
                 inspec.error gt 0.,ntrough)
  if ntrough ne 0 then $
     varflr = mean(inspec.error[trough]^2) > median(inspec.error[trough]^2) $
  else varflr = 0.05            ; assumption

  if keyword_set(search_fil) then begin
     ;; Mask out automatically detected features
     if size(search_fil,/type) eq 7 then $
        srch_strct = xmrdfits(search_fil,1,/silent) $
     else srch_strct = search_fil

     ;; Make sure features found
     if size(srch_strct,/type) ne 8 then goto,QUIT
     
     ;; Apply instrument flag
     if keyword_set(instr_flg) then begin
        gd = where((srch_strct.instr and instr_flg) ge instr_flg, ngd) 
        if ngd gt 0 then srch_strct = srch_strct[gd] $
        else goto,QUIT
     endif 

     ;; Cull and sort
     srt = sort(srch_strct.wv_lim[0])
     srch_strct = srch_strct[srt]
     
     ;; Find regions outside of features
     nsrch = n_elements(srch_strct)
     msk = replicate(0,npix)
     for ii=0,nsrch-1 do begin
        dum = min(nwspec.wave-srch_strct[ii].wv_lim[0],pmin,/absolute)
        dum = min(nwspec.wave-srch_strct[ii].wv_lim[1],pmax,/absolute)
        msk[pmin:pmax] = 1
     endfor 

     ;; Fix LLS and prevent nearly infinite loops
     bd = where(msk eq 1,nbd)
     if keyword_set(lls) and nbd ne 0 then begin
        bnd1 = where(bd-shift(bd,1) ne 1,nbnd1)
        bnd2 = where(bd-shift(bd,-1) ne -1,nbnd2)
        dum = max(nwspec.wave[bnd2]-nwspec.wave[bnd1],ills)
        ;; New "continuum" fit (perfectly straight line at median flux
        ;; in LLS)
        nwspec.flux[bnd1[ills]:bnd2[ills]] = $
           1. + randomn(seed,bnd2[ills]-bnd1[ills]+1)*$
           nwspec.error[bnd1[ills]:bnd2[ills]]
        ;; Match S/N properties to make tests realistic
        ;; (looks horrible)
        mdsnr = median(inspec.flux[bnd1[ills]:bnd2[ills]]/$
                       inspec.error[bnd1[ills]:bnd2[ills]])
        nwspec.error[bnd1[ills]:bnd2[ills]] = $
           nwspec.error[bnd1[ills]:bnd2[ills]]/mdsnr
        ;; Set non-LLS region
        msk[bnd1[ills]:bnd2[ills]] = -1
        unlls = where(msk ge 0)
        ;; Set new S/N properties
        snr = measureSN(nwspec.flux[unlls],err=inspec.error[unlls],$
                        wave=inspec.wave[unlls],dev=dev,inflg=4,/silent)
        ;; Indicate that region already handled
        msk[bnd1[ills]:bnd2[ills]] = 0
     endif                      ; /lls
     
     ;; Random in only regions within acceptable flux limits
     unabs = where(msk eq 0,nunabs,complement=gap,ncomplement=ngap)
     nrng = 10L
     cnt = 0L
     while ngap ne 0 do begin
        num = long((randomu(seed,ngap)-0.5)*nrng) ; draw around pixel
        rng = gap + num
        
        ;; If outside limits, restore original for next round
        bd = where(rng gt npix-1 or rng lt 0,nbd)
        if nbd ne 0 then rng[bd] = gap[bd]
        nwspec.flux[gap] = inspec.flux[rng]
        nwspec.error[gap] = inspec.error[rng]

        ;; New check
        bd = where(msk[rng] eq 1,ngap)  ; used in-line value
        if ngap ne 0 then gap = gap[bd] ; apply from rng to gap

        cnt = cnt + 1
        if (cnt mod 3) eq 0 then nrng = nrng*2 ; may need to go farther
        if cnt gt 1e7 then begin
           stop,'civ_cleanspec: masking going poorly, about to abort'
           ngap = 0
        endif 
     endwhile                                  ; loop gaps
  endif                                        ; search_fil masking
QUIT: 

     ;; Mask out strong absorption regions
  if keyword_set(unlls) then begin
     unabs = where(nwspec.flux[unlls] ge snr[0]-2*dev[0] and $
                   nwspec.error[unlls] ne -1.,nunabs,complement=gap,$
                   ncomplement=ngap) 
     if ngap ne 0 then gap = unlls[gap]
     if nunabs ne 0 then unabs = unlls[unabs]
  endif else $
     unabs = where(nwspec.flux ge snr[0]-2*dev[0] and $
                   nwspec.error ne -1.,nunabs,complement=gap,ncomplement=ngap)


  ;; LLS must be fixed (wasn't before)
  if keyword_set(lls) and ngap ne 0 and $
     not keyword_set(search_fil) then begin
     stop,'civ_cleanspec: fixing LLS w/o search_fil is imperfect'

     bnd1 = where(gap-shift(gap,1) ne 1,nbnd1)
     bnd2 = where(gap-shift(gap,-1) ne -1,nbnd2)
     dum = max(nwspec.wave[bnd2]-nwspec.wave[bnd1],ills)
     ;; New "continuum" fit (perfectly straight line at median flux
     ;; in LLS)
     nwspec.flux[bnd1[ills]:bnd2[ills]] = $
        1. + randomn(seed,bnd2[ills]-bnd1[ills]+1)*$
        nwspec.error[bnd1[ills]:bnd2[ills]]
     ;; Match S/N properties to make tests realistic
     ;; (looks horrible)
     mdsnr = median(inspec.flux[bnd1[ills]:bnd2[ills]]/$
                    inspec.error[bnd1[ills]:bnd2[ills]])
     nwspec.error[bnd1[ills]:bnd2[ills]] = $
        nwspec.error[bnd1[ills]:bnd2[ills]]/mdsnr
     ;; Set non-LLS region
     msk = replicate(0,npix)
     msk[bnd1[ills]:bnd2[ills]] = -1
     unlls = where(msk ge 0)
     ;; Set new S/N properties
     snr = measureSN(nwspec.flux[unlls],err=inspec.error[unlls],$
                     wave=inspec.wave[unlls],dev=dev,inflg=4,/silent)
     ;; Set next step
     unabs = where(nwspec.flux[unlls] ge snr[0]-2*dev[0] and $
                   nwspec.error[unlls] ne -1.,nunabs,complement=gap,$
                   ncomplement=ngap)
     if ngap ne 0 then gap = unlls[gap]
     if nunabs ne 0 then unabs = unlls[unabs]
  endif                         ;/lls

  nrng = 10L                    ; area from which to draw
  cnt = 0L
  while ngap ne 0 do begin
     num = long((randomu(seed,ngap)-0.5)*nrng) ; draw around pixel
     rng = gap + num

     ;; If outside limits, restore original for next round
     bd = where(rng gt npix-1 or rng lt 0,nbd)
     if nbd ne 0 then rng[bd] = gap[bd] 
     nwspec.flux[gap] = inspec.flux[rng]
     nwspec.error[gap] = inspec.error[rng]

     ;; New S/N
     if keyword_set(unlls) then begin
        snr = measureSN(nwspec.flux[unlls],err=nwspec.error[unlls],$
                        wave=nwspec.wave[unlls],dev=dev,inflg=4,/silent)
                                ;: Test
        gap = where(nwspec.flux[unlls] lt snr[0]-2*dev[0] or $
                    nwspec.error[unlls] eq -1.,ngap)
        if ngap ne 0 then gap = unlls[gap]
     endif else begin
        snr = measureSN(nwspec.flux,err=nwspec.error,wave=nwspec.wave,$
                        dev=dev,inflg=4,/silent)
                                ;: Test
        gap = where(nwspec.flux lt snr[0]-2*dev[0] or $
                    nwspec.error eq -1.,ngap)
     endelse 
     cnt = cnt + 1
     if (cnt mod 3) eq 0 then nrng = nrng * 2 ;may need to go farther
  endwhile                                    ; ngap ne 0

  if keyword_set(debug) then begin
     print,'civ_cleanspec debug: displaying region included for noise'
     x_splot,nwspec.wave,nwspec.flux,ytw=nwspec.error,$
             xthr=nwspec.wave[unabs],ythr=nwspec.flux[unabs],/block
  endif 
  
  ;; Create output structures
  outspec = nwspec

  if keyword_set(debug) then stop,'civ_cleanspec debug: about to exit'
end
