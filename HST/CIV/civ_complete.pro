;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_complete.pro               
; Author: Kathy Cooksey                      Date: 26 Aug 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Driver program for completeness tests
; Input: 
;   true_list -- string array or ASCII file listing input profiles
;                (civcandstrct)
;   rec_list -- string array or ASCII file listing output profiles
;               (civcandstrct)
; Optional:
;   ostrct_fil -- name of output completeness comparison file
;   /silent -- don't print messages
;   extra -- pass to civ_complete_hist
; Output: 
;   ostrct_fil -- name of output completeness comparison file
;                 ext=1 = completeness comparison structure
;                 ext=2 = completeness arrays
; External Calls:
; Example:
; History:
;   26 Aug 2008  created by KLC
;   13 Oct 2008  fix case where everything above limit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_genprof                    ; need civ_genprof_config()

function civ_complete_reorg,trufits
;; Total components into first structure arrays
nmx = n_elements(trufits[0].wrest)

nabs = n_elements(trufits)
gd = where(trufits[0].wrest gt 0.)
unq = uniq(trufits[0].wrest[gd],sort(trufits[0].wrest[gd]))
nion = n_elements(unq)

nwfits = trufits            ;civ_instantstrct(nabs)
for ii=0,nabs-1 do begin
    gd = where(trufits[ii].wrest gt 0.,ngd)
    ncomp = ngd/nion

    for jj=0,ncomp-1 do begin
        indx = lindgen(ncomp)*nion+jj
        nwfits[ii].ncolm[jj] = alog10(total(10^trufits[ii].ncolm[indx]))
        nwfits[ii].ew[jj] = total(trufits[ii].ew[indx])
        nwfits[ii].wv_lim[jj,0] = min(trufits[ii].wv_lim[indx,0])
        nwfits[ii].wv_lim[jj,1] = max(trufits[ii].wv_lim[indx,1])
    endfor                      ; loop components

    ;; Eliminate confusion
    nwfits[ii].wrest[nion:nmx-1] = 0.
    nwfits[ii].ion[nion:nmx-1] = 0.
    nwfits[ii].zabs[nion:nmx-1] = 0.
    nwfits[ii].zsig[nion:nmx-1] = 0.
    nwfits[ii].ew[nion:nmx-1] = 0.
    nwfits[ii].sigew[nion:nmx-1] = 0.
    nwfits[ii].ncolm[nion:nmx-1] = 0.
    nwfits[ii].signcolm[nion:nmx-1] = 0.
    nwfits[ii].b[nion:nmx-1] = 0.
    nwfits[ii].sigb[nion:nmx-1] = 0.
    nwfits[ii].flg_colm[nion:nmx-1] = 0.
    nwfits[ii].wv_lim[nion:nmx-1,*] = 0.
    
endfor                          ; loop lines
return,nwfits

end                             ; civ_complete_reorg()



pro civ_complete_loadstrct,arr,nwarr,cindx,dindx
;; Transfer data (only common to civcandstrct)
;; cindx = array indice for true or recovered
;; dindx = array indice for doublet
nwarr.zabs[cindx,dindx] = arr.zabs[dindx]
nwarr.ew[cindx,dindx] = arr.ew[dindx]
nwarr.sigew[cindx,dindx] = arr.sigew[dindx]
nwarr.ncolm[cindx,dindx] = arr.ncolm[dindx]
nwarr.signcolm[cindx,dindx] = arr.signcolm[dindx]
nwarr.b[cindx,dindx] = arr.b[dindx]
nwarr.sigb[cindx,dindx] = arr.sigb[dindx]
nwarr.wv_lim[cindx,dindx,0] = arr.wv_lim[dindx,0]
nwarr.wv_lim[cindx,dindx,1] = arr.wv_lim[dindx,1]
nwarr.instr[cindx,dindx] = arr.instr[dindx]
end                             ;civ_complete_loadstrct



function civ_complete_hist,compltstrct_fil,binncolm=binncolm,binz=binz,$
                           binew=binew,list=list,nlim=nlim,ewlim=ewlim,$
                           zlim=zlim,view=view,savfits=savfits,$
                           combine_strct=combine_strct,$
                           debug=debug,silent=silent
  ;; Histograming
  ;;   binncolm -- resolution for searching column density (default: 0.1)
  ;;   binz -- resolution for binning redshift (default: 0.005 or
  ;;           1500km/s)
  ;;   zlim[] -- forces histogram to same redshift bins
  if (N_params() lt 1) then begin
     print,'Syntax - '+$
           'rslt = civ_complete_hist(compltstrct_fil,[binncolm=,binz=,binew=,/list,zlim=,/view,/savfits,/debug,/silent])'
     return,-1
  endif

  maxno = 9999.d                ; some large number

  ;; Parameters
  ;; Safety values
  config_fil = getenv('MLSS_DIR')+'/pro/civ_genprof.config'
  test = file_search(config_fil,count=ntest)
  if ntest ne 0 then config_strct = civ_genprof_config(config_fil)

  ;; Compile results
  if not keyword_set(binncolm) then begin
     if ntest ne 0 then binncolm = config_strct.binncolm $
     else binncolm = 0.1
  endif 
  if not keyword_set(binz) then begin 
     if ntest ne 0 then binz = config_strct.binz $
     else binz = 0.005          ; 1500./3e5
  endif 
  if not keyword_set(binew) then binew = 10. ; mA

  if not keyword_set(list) then begin
     ;; Read in single file
     if size(compltstrct_fil,/type) eq 7 then $
        fil = xmrdfits(compltstrct_fil,1,/silent) $
     else fil = compltstrct_fil
     nfil = 1

     if not keyword_set(zlim) then $
        zlim = [min(fil.zabs[0,0],max=mx),mx]
     if not keyword_set(nlim) then $
        nlim = [min(fil.ncolm[0,0],max=mx),mx]
     if not keyword_set(ewlim) then $
        ewlim = [min(fil.ew[0,0],max=mx),mx]
  endif else begin
     ;; Multiple files to combine
     if n_elements(compltstrct_fil) eq 1 then $
        readcol,compltstrct_fil,fil,format='a',comment=';' $
     else fil = compltstrct_fil
     nfil = n_elements(fil)

     if not keyword_set(zlim) then zlim = [-0.01,1.01] ; more than enough
     if not keyword_set(nlim) then begin
        if keyword_set(config_strct) then nlim = config_strct.nlim $
        else nlim = [12.7,14.3]
     endif 
     if not keyword_set(ewlim) then ewlim = [10.,2400.] ; mA
  endelse 

  ;; Determine bins and bounds
  nbinz = ceil((zlim[1] - zlim[0])/binz)
  locz = zlim[0] + binz*dindgen(nbinz)
  dz = replicate(binz,nbinz)
  dz[nbinz-1] = zlim[1]- locz[nbinz-1]
  nbinn =  ceil((nlim[1]-nlim[0])/binncolm)
  locn = nlim[0] + binncolm*dindgen(nbinn)
  nbinew = ceil((ewlim[1]-ewlim[0])/binew)
  locw = ewlim[0] + binew*dindgen(nbinew)
  if keyword_set(debug) then $
     stop,'civ_complete_hist() debug: measured bins and bounds'


  ;; Instantiate structures
  cmplt95 = replicate(maxno,nbinz,2)
  cmplt90 = replicate(maxno,nbinz,2)
  histn = replicate(0d,nbinn,2) 
  histw = replicate(0d,nbinew,2) 


  for zz=0,nbinz-1 do begin

     for ff=0,nfil-1 do begin
        ;; Read in
        if size(fil,/type) eq 7 then $
           strct = xmrdfits(fil[ff],1,/silent) $
        else strct = fil

        sub = where(strct.zabs[0,0] ge locz[zz] and $
                    strct.zabs[0,0] lt locz[zz]+binz,nsub)
        if nsub ne 0 then begin
           if size(compltstrct,/type) eq 8 then $ ; Concatenate
              compltstrct = [compltstrct,strct[sub]] $
           else compltstrct = strct[sub] ; Instantiate       
        endif else begin
           if keyword_set(debug) then $
              stop,'civ_complete_hist() debug: no features in redshift bin ',$
                   locz[zz]
           continue
        endelse 

        ;; Warning
        if min(strct.zabs[0,*],max=mx) lt zlim[0] or mx gt zlim[1] then $
           stop,'civ_complete_hist(): some redshifts out of zlim ',$
                fil[ff]
     endfor                     ; loop nfil

     if size(compltstrct,/type) ne 8 then begin
        ;; These bins must be masked out b/c now redshift coverage
        if not keyword_set(silent) then $
           print,'civ_complete_hist(): no features or coverage in locz =',$
                 locz[zz]
        cmplt95[zz,*] = -maxno
        cmplt90[zz,*] = -maxno
        compltstrct = 0         ; reset for next loop
        continue                ; nothing
     endif                      ; no items in locz
     nsub = n_elements(compltstrct)
     
     ;; Histogram all column density
     hist = histogram(compltstrct.ncolm[0,0],binsize=binncolm,$
                      min=nlim[0],max=nlim[1])
     histn[*,0] = histn[*,0] + hist
     ;; Histogram all EW
     hist = histogram(compltstrct.ew[0,0],binsize=binew,$
                      min=ewlim[0],max=ewlim[1])
     histw[*,0] = histw[*,0] + hist

     ;; Want CIV 1548 to be 3sig feature
     mtch = where(compltstrct.flg[0] ne 5 and $
                  compltstrct.flg[1] ne 5,nmtch)
     if nmtch ne 0 then begin
        ;; Histogram recovered line column density
        hist = histogram([compltstrct[mtch].ncolm[0,0]],binsize=binncolm,$
                         min=nlim[0],max=nlim[1])
        histn[*,1] = histn[*,1] + hist
        ;; Histogram recovered line EW
        hist = histogram([compltstrct[mtch].ew[0,0]],binsize=binew,$
                         min=ewlim[0],max=ewlim[1])
        histw[*,1] = histw[*,1] + hist
     endif                      ; nmtch ne 0

     case nmtch of
        0: begin                ; do nothing
           if keyword_set(debug) then $
              stop,'civ_complete_hist() debug: no recovered lines ' + $
                   'for z-bin ',locz[zz]
        end 
        1: begin
           prcnt = nmtch/double(nsub) 
           if prcnt ge 0.9 then begin
              cmplt90[zz,0] = compltstrct[mtch[0]].ncolm[0,0]
              cmplt90[zz,1] = compltstrct[mtch[0]].ew[0,0]
           endif 
           if prcnt ge 0.95 then begin
              cmplt95[zz,0] = compltstrct[mtch[0]].ncolm[0,0]            
              cmplt95[zz,1] = compltstrct[mtch[0]].ew[0,0]            
           endif 
           if keyword_set(debug) then $
              stop,'civ_complete_hist() debug: one recovered line '+$
                   'for z-bin ',locz[zz]
        end
        else: begin
           histncolm = histogram(compltstrct.ncolm[0,0],$
                                 binsize=binncolm,loc=locncolm)
           nnbins = n_elements(locncolm)
           histncolmrec = histogram(compltstrct[mtch].ncolm[0,0],$
                                    binsize=binncolm,$
                                    min=locncolm[0],max=locncolm[nnbins-1])
           histew = histogram(compltstrct.ew[0,0],$
                              binsize=binew,loc=locew)
           newbins = n_elements(locew)
           histewrec = histogram(compltstrct[mtch].ew[0,0],$
                                 binsize=binew,$
                                 min=locew[0],max=locew[newbins-1])

           ;; Interpolate completeness limit
           prcntncolm = replicate(0.d,nnbins)
           bdn = where(histncolm le 0,nbdn,complement=gdn,ncomplement=ngdn)
           if nbdn ne 0 then prcntncolm[bdn] = 0.
           if ngdn ne 0 then prcntncolm[gdn] = $
              histncolmrec[gdn]/double(histncolm[gdn])
           prcntew = replicate(0.d,nnbins)
           bdew = where(histew le 0,nbdew,complement=gdew,ncomplement=ngdew)
           if nbdew ne 0 then prcntew[bdew] = 0.
           if ngdew ne 0 then prcntew[gdew] = $
              histewrec[gdew]/double(histew[gdew])

           ;; Column Density Limit
           test = where(prcntncolm gt 0.9,ntest)
           if ntest gt 0 then begin
              gdn = where(prcntncolm gt 0.,ngdn)
              case ngdn of
                 0: if keyword_set(debug) then stop,$ ; do nothing
                    'civ_complete_hist() debug: nothing recovered for z-bin ',$
                    locz[zz]
                 1: begin
                    ;; Can't extrapolate; just test
                    if prcntncolm[gdn[0]] ge 0.9 then $
                       cmplt90[zz,0] = locncolm[gdn[0]] 
                    if prcntncolm[gdn[0]] ge 0.95 then $
                       cmplt95[zz,0] = locncolm[gdn[0]]
                    if keyword_set(debug) then stop,$
                       'civ_complete_hist() debug: one recovered line, '+$
                       'strict test z-bin ',locz[zz]
                 end 
                 else: begin
                    ;; 4th-order polynomial w/Poisson error
                    rslt = svdfit(locncolm[gdn],prcntncolm[gdn],5,yfit=nfit,$
                                  measure_errors=sqrt(histncolmrec[gdn]))
                    ;; Solve for poly(N) = 0.9
                    tmp = rslt
                    tmp[0] = rslt[0]-0.9
                    soln = fz_roots(tmp,/double)
                    rln = real_part(soln)
                    imn = imaginary(soln)
                    gd = where(imn eq 0. and rln ge min(locncolm[gdn],max=mx) $
                               and rln le mx,ngd)
                    if ngd eq 0 then begin
                       if keyword_set(debug) then stop,$
                          'civ_complete_hist(): no 90% column density limit ',$
                          locz[zz] $
                       else print,$
                          'civ_complete_hist(): no 90% column density limit ',$
                          locz[zz] 
                    endif else cmplt90[zz,0] = min(rln[gd])
                    ;; Solve for poly(N) = 0.95
                    tmp = rslt
                    tmp[0] = rslt[0]-0.95
                    soln = fz_roots(tmp,/double)
                    rln = real_part(soln)
                    imn = imaginary(soln)
                    gd = where(imn eq 0. and rln ge min(locncolm[gdn],max=mx) $
                               and rln le mx,ngd)
                    if ngd eq 0 then begin
                       if keyword_set(debug) then stop,$
                          'civ_complete_hist(): no 95% column density limit ',$
                          locz[zz] $
                       else print,$
                          'civ_complete_hist(): no 95% column density limit ',$
                          locz[zz] 
                    endif else cmplt95[zz,0] = min(rln[gd])
                 end            ; else (4th order polynomial)
              endcase           ; ngdn (# recovered) 
           endif                ; test

           ;; EW limit
           test = where(prcntew gt 0.9,ntest)
           if ntest gt 0 then begin
              gdew = where(prcntew gt 0.,ngdew)
              case ngdew of
                 0: if keyword_set(debug) then stop,$ ; do nothing
                    'civ_complete_hist() debug: nothing recovered for z-bin ',$
                    locz[zz]
                 1: begin
                    ;; Can't extrapolate; just test
                    if prcntew[gdew[0]] ge 0.9 then $
                       cmplt90[zz,1] = locew[gdew[0]] 
                    if prcntew[gdew[0]] ge 0.95 then $
                       cmplt95[zz,1] = locew[gdew[0]]
                    if keyword_set(debug) then stop,$
                       'civ_complete_hist() debug: one recovered line, '+$
                       'strict test z-bin ',locz[zz]
                 end 
                 else: begin
                    ;; 4th-order polynomial with Poisson errors
                    rslt = svdfit(locew[gdew],prcntew[gdew],5,yfit=ewfit,$
                                  measure_errors=sqrt(histewrec[gdew]))
                    ;; Solve for poly(EW) = 0.9
                    tmp = rslt
                    tmp[0] = rslt[0]-0.9
                    soln = fz_roots(tmp,/double)
                    rln = real_part(soln)
                    imn = imaginary(soln)
                    gd = where(imn eq 0. and rln ge min(locew[gdew],max=mx) $
                               and rln le mx,ngd)
                    if ngd eq 0 then begin
                       if keyword_set(debug) then $
                          stop,'civ_complete_hist(): no 90% EW limit ',$
                               locz[zz] $
                       else $
                          print,'civ_complete_hist(): no 90% EW limit ',$
                                locz[zz] 
                    endif else cmplt90[zz,1] = min(rln[gd])
                    ;; Solve for poly(EW) = 0.95
                    tmp = rslt
                    tmp[0] = rslt[0]-0.95
                    soln = fz_roots(tmp,/double)
                    rln = real_part(soln)
                    imn = imaginary(soln)
                    gd = where(imn eq 0. and rln ge min(locew[gdew],max=mx) $
                               and rln le mx,ngd)
                    if ngd eq 0 then begin
                       if keyword_set(debug) then $
                          stop,'civ_complete_hist(): no 95% EW limit ',$
                               locz[zz] $
                       else $
                          print,'civ_complete_hist(): no 95% EW limit ',$
                                locz[zz] 
                    endif else cmplt95[zz,1] = min(rln[gd])
                 end            ; else
              endcase           ; ngdew (# recovered)
           endif                ; test 

           if keyword_set(view) then $
              x_splot,locncolm,histncolm,ytwo=histncolmrec,psym1=10,psym2=10,$
                      title=string(locz[zz],format='(f5.3)'),/block
        end
     endcase                    ; nmtch (# recovered) 
     
     ;; Process for NaN and +/-Inf
     if not finite(cmplt90[zz,0]) then begin
        print,'civ_complete: non-finite column density 90% limit'
        cmplt90[zz,0] = maxno
     endif 
     if not finite(cmplt90[zz,1]) then begin
        print,'civ_complete: non-finite EW 90% limit'
        cmplt90[zz,1] = maxno
     endif 

     if not finite(cmplt95[zz,0]) then begin
        print,'civ_complete: non-finite column density 95% limit'
        cmplt95[zz,0] = maxno
     endif 
     if not finite(cmplt95[zz,1]) then begin
        print,'civ_complete: non-finite EW 95% limit'
        cmplt95[zz,1] = maxno
     endif 

     ;; Reset
     compltstrct = 0
  endfor                        ; loop nbinz

  ;; Post-process and blanked out bins
  for ii=0,4-1 do begin
     case ii of
        0: tmp = cmplt90[*,0]
        1: tmp = cmplt90[*,1]
        2: tmp = cmplt95[*,0]
        3: tmp = cmplt95[*,1]
        else: stop,'civ_complete_hist(): this is weird'
     endcase

     bd = where(tmp ge maxno,nbd)
     gd = where(abs(tmp) lt maxno,ngd)
     if ngd eq 0 then begin
        print,'civ_complete_hist(): no good data'
        continue
     endif 
     for jj=0,nbd-1 do begin
        ;; Find neighboring good values
        done = 0
        ilo = bd[jj] - 1
        while (done ne 1) and (ilo ge 0) do begin
           mtch = where(ilo eq gd,nmtch)
           if nmtch eq 1 then done = 1 $
           else ilo = ilo - 1
        endwhile
        done = 0
        ihi = bd[jj] + 1
        while (done ne 1) and (ihi le nbinz-1) do begin
           mtch = where(ihi eq gd,nmtch)
           if nmtch eq 1 then done = 1 $
           else ihi = ihi + 1
        endwhile

        ;; Assess value
        sub = [ilo,ihi]
        if ilo ge 0 and ihi le nbinz-1 then begin
           rslt = interpol(tmp[sub],locz[sub],locz[bd[jj]])
           tmp[bd[jj]] = rslt
        endif else begin        ; ilo and ihi in bounds
           if ilo lt 0 then tmp[bd[jj]] = tmp[ihi]
           if ihi gt nbinz-1 then tmp[bd[jj]] = tmp[ilo]
        endelse                 ; ilo or ihi out of bounds
     endfor

     ;; Max out gaps in redshift coverage
     bd = where(tmp le -maxno,nbd)
     if nbd ne 0 then $
        tmp[bd] = maxno

     ;; Store
     case ii of 
        0: cmplt90[*,0] = tmp
        1: cmplt90[*,1] = tmp
        2: cmplt95[*,0] = tmp
        3: cmplt95[*,1] = tmp
        else: stop,'civ_complete_hist(): this is weird'
     endcase

     ;; Write out (safety)
     if keyword_set(savfits) then begin
        ofil = savfits+'_z'+strtrim(zz,2)+'.fits'
        mwrfits,compltstrct,ofil,/create,/silent
        strct = {locz:locz[zz], dz:dz[zz], $
                 cmplt90:cmplt90[zz,*],cmplt95:cmplt95[zz,*]}
        mwrfits,strct,ofil,/silent
        strct = {locn:locncolm,histn:histn[*,0],histnrec:histn[*,1],$
                 locew:locw,histew:histw[*,0],histewrec:histw[*,1]}
        mwrfits,strct,ofil,/silent ; continually updating    
        if not keyword_set(silent) then $
           print,'civ_complete_hist(): created ',ofil
     endif
  endfor                        ; loop blanked out bins

  strct = {locz:locz, dz:dz, cmplt90:cmplt90, cmplt95:cmplt95}
  combine_strct = {locn:locncolm,histn:histn[*,0],histnrec:histn[*,1],$
                   locew:locw, histew:histw[*,0],histewrec:histw[*,1]}
  if keyword_set(debug) then stop,'civ_complete debug: about to return'
  return,strct
end                             ; civ_complete_hist()


pro civ_complete_plot,strct_fil,psfil,option,siiv=siiv,$
                      csize=csize,lthick=lthick,$
                      xrng=xrng,yrng=yrng,binhist=binhist,$
                      rebin=rebin,pltnhist=pltnhist,pltnum=pltnum,_extra=extra
  ;; Plot 90% and 95% completeness limits vs redshift
  if not keyword_set(option) then option = 0 ; Column density

  ;; Read in structures
  compltstrct = xmrdfits(strct_fil,1,/silent)
  if keyword_set(rebin) then $
     histstrct = civ_complete_hist(compltstrct,_extra=extra) $ ;binncolm,binz
  else histstrct = xmrdfits(strct_fil,2,/silent)

  ;; Recover parameters
  binz = histstrct.locz[1]-histstrct.locz[0]
  if not keyword_set(binhist) then begin
     if option eq 0 then binhist = 0.1 $
     else binhist = 10.
  endif 
  nbinz = n_elements(histstrct.locz)
  zlim = [histstrct.locz[0],histstrct.locz[nbinz-1]]
  gd = where(histstrct.cmplt90[*,option] gt 0 and $
             histstrct.cmplt95[*,option] gt 0,ngd)
  if ngd ne 0 then $
     nlim = [min([histstrct.cmplt90[gd,option],$
                  histstrct.cmplt95[gd,option]],max=mx),mx] $
  else begin
     if option eq 0 then nlim = [12.7,14.3] $
     else nlim = [10.,150.] ; mA
  endelse

  ;; Plot (and params)
  if keyword_set(psfil) then begin
     if stregex(psfil,'ps',/boolean) eq 0 then psfil = psfil+'.ps'
     x_psopen,psfil,/maxs
  endif 
  if keyword_set(pltnhist) then !p.multi = [2,1,2] $
  else !p.multi = [1,1,1]
  if keyword_set(pltnum) then !x.margin = [8,6] else !x.margin = [8,3]
  !y.margin = [3,1] 
  if keyword_set(pltnhist) then !y.omargin = [1,0]
  clr = getcolor(/load)
  if not keyword_set(csize) then csize = 2
  if not keyword_set(lthick) then lthick = 2.
  if not keyword_set(xrng) then xrng = [zlim[0]-0.5*binz,zlim[1]+0.5*binz]
  if not keyword_set(yrng) then yrng = [nlim[0]-0.05,nlim[1]+0.05]

  if keyword_set(siiv) then begin
     if option eq 0 then ytitle = 'log !8N!X(Si!E+3!N) Limit' $
     else ytitle = '!8W!X!D!8r!X,1393!N Limit (m'+STRING("305B) +')'
  endif else begin
     if option eq 0 then ytitle = 'log !8N!X(C!E+3!N) Limit' $
     else ytitle = '!8W!X!D!8r!X,1548!N Limit (m'+STRING("305B) +')'
  endelse 
  xtitle = '!8z!X'

  ;; Frame
  if keyword_set(pltnum) then ysty = 9 else ysty = 1
  plot,xrng,yrng,/nodata,ystyle=ysty,/xstyle,background=clr.white,$
       color=clr.black,ytitle=ytitle,xtitle=xtitle,charsize=csize,thick=lthick
  ;; 90%
  oplot,histstrct.locz,histstrct.cmplt90[*,option],psym=10,color=clr.black,$
        thick=lthick
  ;; 95%
  oplot,histstrct.locz,histstrct.cmplt95[*,option],psym=10,color=clr.red,$
        thick=lthick

  if keyword_set(pltnum) then begin
     histz = histogram(compltstrct.zabs[0,option],binsize=binz,$
                       min=zlim[0],max=zlim[1])

     ;; # per z-bin
     oplot,histstrct.locz,histz/double(total(histz))*100.*$
           (yrng[1]-yrng[0])+yrng[0],psym=10,color=clr.blue,thick=lthick
     ;; New axis is inverted equation to scale histogram numbers
     axis,/yaxis,yrange=[0,total(histz)/100.],charsize=csize*0.75,/ystyle,$
          ytitle='Number per bin',color=clr.black
  endif 

  if keyword_set(pltnhist) then begin
     if keyword_set(siiv) then begin
        if option eq 0 then xtitle = 'log !8N!X(Si!E+3!N)' $
        else xtitle = '!8W!X!D!8r!X,1393!N (m'+STRING("305B) +')'
     endif else begin
        if option eq 0 then xtitle = 'log !8N!X(C!E+3!N)' $
        else xtitle = '!8W!X!D!8r!X,1548!N (m'+STRING("305B) +')'
     endelse 
     ytitle = 'Number'
     
     ;; Frame
     if option eq 0 then $
        histn = histogram(compltstrct.ncolm[0,0],binsize=binhist,min=yrng[0],$
                          max=yrng[1],locations=locn) $
     else $
        histn = histogram(compltstrct.ew[0,0],binsize=binhist,min=yrng[0],$
                          max=yrng[1],locations=locn)      
     plot,yrng,[0,max(histn)+5],/nodata,/ystyle,/xstyle,background=clr.white,$
          color=clr.black,ytitle=ytitle,xtitle=xtitle,$
          charsize=csize,thick=lthick

     ;; All
     oplot,locn,histn,psym=10,color=clr.black,thick=lthick

     ;; Recovered
     gd = where(compltstrct.flg[0] ne 5 and compltstrct.flg[1] ne 5,ngd)
     if ngd ne 0 then begin
        if option eq 0 then $
           histn = histogram(compltstrct[gd].ncolm[0,0],binsize=binhist,$
                             min=yrng[0],max=yrng[1],locations=locn) $
        else $
           histn = histogram(compltstrct[gd].ew[0,0],binsize=binhist,$
                             min=yrng[0],max=yrng[1],locations=locn) 
        oplot,locn,histn,psym=10,color=clr.red,thick=lthick
     endif                      ; recovered 
     
  endif                         ; /pltnhist

  if keyword_set(psfil) then begin
     x_psclose
     print,'civ_complete_plot: created ',psfil
  endif 

end                             ; civ_complete_plot


pro civ_complete_strongabs,cmplt_fil,search_fil,spec_fil,siiv=siiv,$
                           instr_flg=instr_flg
  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  ;; Account for strong Lya absorbers and the pathlength blocked
  maxno = 9999.d                  ; some large number
  npix_buff = 5L                  ; grouping of saturated pixels

  ;; Read-in completeness file
  if size(cmplt_fil,/type) eq 7 then $
     cmplt_strct = xmrdfits(cmplt_fil,2,/silent) $
  else cmplt_strct = cmplt_fil
  nbinz = n_elements(cmplt_strct.locz)
  nwlocz = cmplt_strct.locz
  nwdz = cmplt_strct.dz
  nwcmplt90 = cmplt_strct.cmplt90
  nwcmplt95 = cmplt_strct.cmplt95

  ;; This already accounts for the pertinent CIV coverage
  zlim = dblarr(2)
  zlim[0] = min(cmplt_strct.locz-cmplt_strct.dz)
  zlim[1] = max(cmplt_strct.locz+cmplt_strct.dz)
  
  ;; Automatically-detected features
  if size(search_fil,/type) eq 7 then $
     srch_strct = xmrdfits(search_fil,1,/silent) $
  else srch_strct = search_fil
  if keyword_set(instr_flg) then begin
     gd = where((srch_strct.instr and instr_flg) ge instr_flg,nstrct)
     if nstrct eq 0 then begin
        ;; Exit and change nothing
        print,'civ_complete_strongabs: no absorbers in ',search_fil
        return
     endif 
     srch_strct = srch_strct[gd]
  endif 

  ;; Forcing the issue around CIV
  civ = dblt_retrieve(dblt_name)

  ;; Assume everything HST format
  fil_sig = strmid(spec_fil,0,strpos(spec_fil,'_f.fits'))+$
            '_e.fits'
  fx = x_readspec(spec_fil,wav=wv,sig=er,fil_sig=fil_sig,npix=npix)

  ;; Maximize range: zlo <= zobs(if 1550) and zobs(if 1548) <= zhi
  gd = where(srch_strct.wrest/civ.wvII-1 ge zlim[0] and $
             srch_strct.wrest/civ.wvI-1. le zlim[1],nstrct)

  if nstrct eq 0 then begin
     ;; Exit and change nothing
     print,'civ_complete_strongabs: no absorbers blocking '+dblt_name+$
           ' in ',search_fil
     return
  endif 
  ;; Subset of structure
  srch_strct = srch_strct[gd]

  for ii=0,nstrct-1 do begin
     ;; Measure percent of feature that is "strong"
     dum = min(srch_strct[ii].wv_lim[0]-wv,pmin,/absolute)
     dum = min(srch_strct[ii].wv_lim[1]-wv,pmax,/absolute)
     nrng = pmax - pmin + 1
     rng = lindgen(nrng) + pmin
     sat = where(fx[rng] le er[rng],nsat) ; within one sigma of zero
     
     while nsat ne 0 do begin
        ;; Group by three pixels
        gd = sat[0]
        if gd[0] gt 0 then gd = [(gd-npix_buff > 0),gd] ; buffer

        if nsat gt 1 then done = 0 else begin
           sat = -1
           nsat = 0
           done = 1
        endelse 
        imx = 1
        while (done eq 0) do begin
           if rng[sat[imx]]-max(rng[gd]) lt 2*npix_buff then begin
              gd = [gd,sat[imx]] 
              if imx lt (nsat-1) then begin
                 imx = imx + 1 
              endif else begin
                 sat = -1
                 nsat = 0
                 done = 1
              endelse 
           endif else begin
              sat = sat[imx:nsat-1]
              nsat = n_elements(sat)
              done = 1
           endelse 
        endwhile 
        ngd = n_elements(gd)
        if gd[ngd-1] lt nrng-1 then begin
           gd = [gd,(gd[ngd-1]+npix_buff < nrng-1)] ; buffer
           ngd = ngd + 1
        endif 

        for jj=0,1 do begin
           ;; Mask out for both 1548 and 1550
           if jj eq 0 then zciv = wv/civ.wvI-1. $
           else zciv = wv/civ.wvII-1.

           ;; Re-sample strongest part (make contiguous)
           pmin = rng[gd[0]]
           pmax = rng[gd[ngd-1]]

           if zciv[pmax] lt nwlocz[0] or $
              zciv[pmin] gt nwlocz[nbinz-1]+nwdz[nbinz-1] then continue
           
           ;; Bins locz defined on left-hand side (preserve orientation)
           dum = min(nwlocz-zciv[pmin],imn,/absolute)
           if nwlocz[imn] gt zciv[pmin] and imn ne 0 then imn = imn-1
           dum = min(nwlocz-zciv[pmax],imx,/absolute)
           if nwlocz[imx] gt zciv[pmax] and imx ne 0 then imx = imx-1
           
           ;; 90% completeness
           tmp0 = [nwcmplt90[0:imn,0],maxno,nwcmplt90[imx:nbinz-1,0]]
           tmp1 = [nwcmplt90[0:imn,1],maxno,nwcmplt90[imx:nbinz-1,1]]
           nnwbinz = n_elements(tmp0)
           nwcmplt90 = dblarr(nnwbinz,2)
           nwcmplt90[*,0] = tmp0
           nwcmplt90[*,1] = tmp1
           ;; 95% completeness
           tmp = dblarr(nbinz+3,2)
           tmp0 = [nwcmplt95[0:imn,0],maxno,nwcmplt95[imx:nbinz-1,0]]
           tmp1 = [nwcmplt95[0:imn,1],maxno,nwcmplt95[imx:nbinz-1,1]]
           nwcmplt95 = dblarr(nnwbinz,2)
           nwcmplt95[*,0] = tmp0
           nwcmplt95[*,1] = tmp1
           ;; dz (to the right)
           if imn gt 0 then begin
              if imx lt nbinz-1 then $
                 nwdz = [nwdz[0:imn-1],zciv[pmin]-nwlocz[imn],$
                         zciv[pmax]-zciv[pmin],$
                         nwlocz[imx+1]-zciv[pmax],nwdz[imx+1:nbinz-1]] $
              else $
                 nwdz = [nwdz[0:imn-1],zciv[pmin]-nwlocz[imn],$
                         zciv[pmax]-zciv[pmin],$
                         nwlocz[imx]+nwdz[nbinz-1]-zciv[pmax]] 
           endif else begin
              if imx lt nbinz-1 then $
                 nwdz = [zciv[pmin]-nwlocz[imn],zciv[pmax]-zciv[pmin],$
                         nwlocz[imx+1]-zciv[pmax],nwdz[imx+1:nbinz-1]] $
              else $
                 nwdz = [zciv[pmin]-nwlocz[imn],zciv[pmax]-zciv[pmin],$
                         nwlocz[imx]+nwdz[nbinz-1]-zciv[pmax]] 
           endelse 
           ;; locz (defined on left)
           if imx lt nbinz-1 then $
              nwlocz = [nwlocz[0:imn],zciv[pmin],zciv[pmax],$
                        nwlocz[imx+1:nbinz-1]] $
           else $
              nwlocz = [nwlocz[0:imn],zciv[pmin],zciv[pmax],$
                        nwlocz[nbinz-1]]
           nbinz = nnwbinz
        endfor                  ; loop CIV doublet
     endwhile                   ; while nsat > 0
  endfor                        ; loop nstrct
  
  ;; Store results
  nwcmplt_strct = {locz:nwlocz, dz:nwdz, cmplt90:nwcmplt90, cmplt95:nwcmplt95}
  if size(cmplt_fil,/type) eq 7 then begin
     tmp = xmrdfits(cmplt_fil,1,/silent)
     mwrfits,tmp,cmplt_fil,/create,/silent
     mwrfits,nwcmplt_strct,cmplt_fil,/silent
     print,'civ_complete_strongabs: overwrote ',cmplt_fil
  endif else cmplt_fil = nwcmplt_strct

end                             ; civ_complete_strongabs


pro civ_complete_combine, compltstrct_fil, ostrct_fil, list=list,$
                          debug=debug
  ;; Combine overlapping spectra by evaluating which has better
  ;; (lower) completeness limts
  maxno = 9999.d

  if keyword_set(list) then $
     readcol,compltstrct_fil,fil,format='a' $
  else fil = compltstrct_fil
  nfil = n_elements(fil)

  ;; Determine coverage
  zlim = dblarr(nfil,2)
  for ii=0,nfil-1 do begin
     strct = xmrdfits(fil[ii],2,/silent)
     nbin = n_elements(strct.locz)
     zlim[ii,0] = strct.locz[0]
     zlim[ii,1] = strct.locz[nbin-1]+strct.dz[nbin-1]
  endfor 

  ;; Sort, to move left to right (low-z to high-z)
  srt = sort(zlim[*,0])
  fil = fil[srt]
  zlim[*,0] = zlim[srt,0]
  zlim[*,1] = zlim[srt,1]

  for ii=0,nfil-1 do begin
     strct = xmrdfits(fil[ii],2,/silent)
     nbin = n_elements(strct.locz)
     if ii eq 0 then begin
        ;; Set up first file
        locz = strct.locz
        dz = strct.dz
        cmplt90_n = strct.cmplt90[*,0]
        cmplt90_ew = strct.cmplt90[*,1]
        cmplt95_n = strct.cmplt95[*,0]
        cmplt95_ew = strct.cmplt95[*,1]
        nbin_tot = nbin
        continue
     endif 

     ;; Constrain maximum overlap region
     lap = where(strct.locz+strct.dz ge locz[0] and $
                 strct.locz le $
                 locz[nbin_tot-1]+dz[nbin_tot-1],nlap,$
                 complement=gd,ncomplement=ngd)

     ;; Non-overlapping portion plus seam
     if ngd ne 0 then begin
        if keyword_set(debug) then $
           stop,'civ_complete_combine: add non-overlapping portion ',fil[ii]
        locz = [locz,locz[nbin_tot-1]+dz[nbin_tot-1],strct.locz[gd]]
        dz = [dz,strct.locz[gd[0]]-(locz[nbin_tot-1]+dz[nbin_tot-1]),strct.dz[gd]]
           nbin_tot = nbin_tot + ngd + 1
        if nlap eq 0 then begin
           ;; Gap between spectra
           cmplt90_n = [cmplt90_n,maxno,strct.cmplt90[gd,0]]
           cmplt90_ew = [cmplt90_ew,maxno,strct.cmplt90[gd,1]]
           cmplt95_n = [cmplt95_n,maxno,strct.cmplt95[gd,0]]
           cmplt95_ew = [cmplt95_ew,maxno,strct.cmplt95[gd,1]]
        endif else begin
           ;; Partial bin
           cmplt90_n = [cmplt90_n,strct.cmplt90[gd[0],0],strct.cmplt90[gd,0]]
           cmplt90_ew = [cmplt90_ew,strct.cmplt90[gd[0],1],strct.cmplt90[gd,1]]
           cmplt95_n = [cmplt95_n,strct.cmplt95[gd[0],0],strct.cmplt95[gd,0]]
           cmplt95_ew = [cmplt95_ew,strct.cmplt95[gd[0],1],strct.cmplt95[gd,1]]

           ;; Adjust lapped bin
           strct.dz[lap[nlap-1]] = locz[nbin_tot-ngd-1] - strct.locz[lap[nlap-1]]
        endelse
     endif                      ; handle non-overlaping portion
     
     ;; Overlapping portion
     if nlap ne 0 then begin
        ;; Determine borders
        dum = min(strct.locz[lap[0]]-locz,ilo,/absolute)
        if locz[ilo] gt strct.locz[lap[0]] $
           and ilo ne 0 then $
           ilo = ilo - 1
        dum = min(strct.locz[lap[nlap-1]]+strct.dz[lap[nlap-1]]-locz,ihi,/absolute)

        rng = ilo + lindgen(ihi-ilo+1)
        if keyword_set(debug) then $
           stop,'civ_complete_combine: overlapping portion range ',fil[ii]

        ;; Test solely on 95% completeness in logN
        replace = 0
        gd_nw = where(strct.cmplt95[lap,0] lt maxno,ngd_nw)
        gd_old = where(cmplt95_n[rng] lt maxno,ngd_old)

        if ngd_old eq 0 then begin
           ;; Replace stored region with new 'lap' region
           replace = 1
        endif 
        if ngd_nw ne 0 and ngd_old ne 0 then begin
           if median(strct.cmplt95[lap[gd_nw],0]) lt $
              median(cmplt95_n[rng[gd_old]]) then begin
              ;; Replace stored region with new 'lap' region
              replace = 1
           endif 
           ;; if ngd_nw eq 0 then keep original; so do nothing
        endif

        if replace eq 1 then begin
           ;; Actually replace b/c new completeness is better
           if keyword_set(debug) then $
              stop,'civ_complete_combine: actually replace section ',fil[ii]

           ;; Account for partial bin
           locz = [locz[0:ilo],strct.locz[lap[0:nlap-1]],locz[ihi:nbin_tot-1]]
           if ilo ne 0 then $
              dz = [dz[0:ilo-1],strct.locz[lap[0]]-locz[ilo],strct.dz[lap[0:nlap-1]],$
                    dz[ihi:nbin_tot-1]] $
           else dz = [strct.locz[lap[0]]-locz[ilo],strct.dz[lap[0:nlap-1]],$
                      dz[ihi:nbin_tot-1]]
           cmplt90_n = [cmplt90_n[0:ilo],strct.cmplt90[lap[0:nlap-1],0],$
                        cmplt90_n[ihi:nbin_tot-1]]
           cmplt90_ew = [cmplt90_ew[0:ilo],strct.cmplt90[lap[0:nlap-1],1],$
                         cmplt90_ew[ihi:nbin_tot-1]]
           cmplt95_n = [cmplt95_n[0:ilo],strct.cmplt95[lap[0:nlap-1],0],$
                        cmplt95_n[ihi:nbin_tot-1]]
           cmplt95_ew = [cmplt95_ew[0:ilo],strct.cmplt95[lap[0:nlap-1],1],$
                         cmplt95_ew[ihi:nbin_tot-1]]
           nbin_tot = nbin_tot + nlap - (ihi-ilo) + 1
        endif                   ; repl
     endif                      ; handle overlapping portion

;     x_splot,locz,cmplt95_n,psym1=10,xtwo=strct.locz,ytwo=strct.cmplt95[*,0],psym2=10,/block
  endfor                        ; loop nfil
  
  cmplt90 = dblarr(nbin_tot,2)
  cmplt90[*,0] = cmplt90_n
  cmplt90[*,1] = cmplt90_ew
  cmplt95 = dblarr(nbin_tot,2)
  cmplt95[*,0] = cmplt95_n
  cmplt95[*,1] = cmplt95_ew
  strct = {locz:locz, dz:dz, cmplt90:cmplt90, cmplt95:cmplt95}
  if size(ostrct_fil,/type) eq 7 then begin
     ;; Write out (the histogram must be the 2nd extension)
     mwrfits,strct,ostrct_fil,/create,/silent
     mwrfits,strct,ostrct_fil,/silent
     print,'civ_complete_combine: created ',ostrct_fil
  endif else ostrct_fil = strct
  
end                             ; civ_complete_combine


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_complete,true_list,rec_list,search_fil=search_fil,spec_fil=spec_fil,$
                 ostrct_fil=ostrct_fil,silent=silent,$
                 _extra=extra
  if (N_params() lt 2) then begin
     print,'Syntax - '+$
           'civ_complete,true_list,rec_list,[' + $
           '           ostrct_fil=,_extra=,/silent]'
     return
  endif

  ;; Parameters
  if not keyword_set(ostrct_fil) then ostrct_fil = 'civ_complete.fits'

  ;; Read in
  if size(true_list,/dimension) eq 0 then $
     readcol,true_list,trufits,format='a',/silent $
  else trufits = true_list
  if size(rec_list,/dimension) eq 0 then $
     readcol,rec_list,recfits,format='a',/silent $
  else recfits = rec_list

  nreal = n_elements(trufits)
  if nreal ne n_elements(recfits) then $
     stop,'civ_complete: mismatch array sizes for comparison'


  ;; Compile information
  tmpstrct = { wrest:dblarr(2), $
               zabs:dblarr(2,2), $
               ew:dblarr(2,2), $
               sigew:dblarr(2,2), $
               b:dblarr(2,2), $
               sigb:dblarr(2,2), $ 
               ncolm:dblarr(2,2), $
               signcolm:dblarr(2,2), $
               wv_lim:dblarr(2,2,2), $
               flg:replicate(0,2), $
               flg_sys:replicate(0,2), $ 
               instr:lonarr(2,2), $ 
               fil:strarr(2) $
             }

  for ii=0,nreal-1 do begin
     tru = xmrdfits(trufits[ii],1,/silent)
     nutru = civ_complete_reorg(tru)
     gd = where(tru[0].wrest gt 0.,ngd)
     if ngd eq 0 then stop,'civ_complete: no values in array'
     unq = uniq(tru[0].wrest[gd],sort(tru[0].wrest[gd]))
     nion = n_elements(unq)
     ncomp = ngd/nion

     nsub = n_elements(tru)
     tmp = replicate(tmpstrct,nsub)
     civ_complete_loadstrct,nutru,tmp,0L,0L
     civ_complete_loadstrct,nutru,tmp,0L,1L
     tmp.wrest[0] = nutru.wrest[0]
     tmp.fil[0] = trufits[ii]

     test = file_search(recfits[ii],count=ntest)
     if ntest eq 0 then begin
        ;; Assume nothing recovered
        for kk=0,nion-1 do tmp.flg[kk] = 5
        
     endif else begin
        rec = xmrdfits(recfits[ii],1,/silent)
        
        tmp.wrest[1] = nutru.wrest[1]
        tmp.fil[1] = recfits[ii]
        
        if size(rec,/type) eq 8 then begin
           for jj=0,nsub-1 do begin
              for kk=0,nion-1 do begin
                 ;; if embedded line within bounds of auto-line
                 mtch = where(nutru[jj].wrest[kk]*(1+nutru[jj].zabs[kk]) ge $
                              rec.wv_lim[kk,0] and $
                              nutru[jj].wrest[kk]*(1+nutru[jj].zabs[kk]) le $
                              rec.wv_lim[kk,1],nmtch)
                 case nmtch of
                    0: begin    ;not detected
                       tmp[jj].flg[kk] = 5
                       if nutru[jj].ncolm[0] ge 14. then print,$
                          'civ_complete: logN >= 14 feature not recovered, b =', $
                          nutru[jj].b[kk]
                    end         ; nmtch = 0
                    1: begin
                       dum = tmp[jj]
                       civ_complete_loadstrct,rec[mtch],dum,1L,kk
                       tmp[jj] = dum
                    end 
                    else: begin
                       if not keyword_set(silent) then begin
                          print,'civ_complete: more than one match'
                          print,jj,nutru[jj].wrest[kk]*(1+nutru[jj].zabs[kk]),$
                                nutru[jj].wv_lim[kk,0],$
                                nutru[jj].wv_lim[kk,1],nutru[jj].ew[kk]
                          printcol,mtch,rec[mtch].wrest[kk]*(1+rec[mtch].zabs[kk]),$
                                   rec[mtch].wv_lim[kk,0],rec[mtch].wv_lim[kk,1],$
                                   rec[mtch].ew[kk]
                       endif    ; print
                       
                       ;; Accept the closest
                       tmp[jj].flg[kk] = 3 ;lower limit (because split)
                       dum = min(nutru[jj].wrest[kk]*(1+nutru[jj].zabs[kk])-$
                                 rec[mtch].wrest[kk],imtch,/absolute)
                       mtch = mtch[imtch]
                       dum = tmp[jj]
                       civ_complete_loadstrct,rec[mtch],dum,1L,kk
                       dum.flg_sys = rec[mtch].flg_sys
                       tmp[jj] = dum
                    end
                 endcase        ; nmtch
                 
                 if keyword_set(debug) then $
                    stop,'civ_complete debug: finished added line ',jj
              endfor                                  ; loop nion
           endfor                                     ; loop nsub
        endif else for kk=0,nion-1 do tmp.flg[kk] = 5 ; nothing recovered
     endelse                                          ; lines recovered

     ;; Save the results
     if keyword_set(compltstrct) then compltstrct = [compltstrct,tmp] $
     else compltstrct = tmp

     if keyword_set(debug) then $
        stop,'civ_complete debug: finished file ',fitfil[ii]
  endfor                        ; loop nreal
  
  mwrfits,compltstrct,ostrct_fil,/create,/silent      ; ext=1
  strct = civ_complete_hist(compltstrct,_extra=extra) ; binncolm,binz
  mwrfits,strct,ostrct_fil,/silent                    ; ext=2
  if not keyword_set(silent) then print,'civ_complete: created ',ostrct_fil

  ;; Account for strong liens
  if keyword_set(search_fil) and keyword_set(spec_fil) then $
     civ_complete_strongabs,ostrct_fil,search_fil,spec_fil,$
                            _extra=extra ; /siiv, instr_flg

end
