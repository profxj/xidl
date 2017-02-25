;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_genprof.pro             
; Author: Kathy Cooksey                      Date: 23 Mar 2008
; Project: MLSS--HST survey of CIV with Xavier Prochaska
; Description: Generate civcandstrct of randomly generated
;              values of z, [N, b, dv, ncomp]
; Input: 
;   nabs - number of absorption lines to make for each wrest
;   zlim - bounds of randomly generated redshifts
;   nlim - bounds of randomly generated column densities (logN)
;   blim - bounds of randomly generated Doppler parameters (km/s)
;   nlimobs - bounds of observed column densities that define f(N)
;   nabsobs - number of absorbers found in nlimobs (normalize f(N))
;   aobs - power-law exponent f(N) = b*N^(-aobs)
; Optional:
;   seed - IDL randomu random number generation seed
;   dvabs - bounds of randomly generated velocity offsets for 
;           components (km/s)
;   ndv - number of components, used for each input wrest, and
;         cannot exceed number civcandstrct holds
;   debug - stop at opportune times for debugging purposes
; Output: 
;   strct_fil - name of output civcandstrct, with nabs elements,
;               with n_elements(wrest) times ndv each
; Example:
;   civ_genprof,'randciv.fits',50,[0,1],[12,15],[10,40],
;      [12.7,14.5],21,1.83,,dvabs=[-200,200],ndv=5
;   [Danfort & Shull 2007ph: 12.7 < ncolm < 14.5, alpha = 1.83+/-0.32,
;                            delta(z) = 2.22, nabsobs = 21 for CIV]
; History:
;   23 Mar 2008 - created by KLC
;   21 Apr 2008 - draw ncolm from frequency distribution 
;                 (power-law)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@spec_ewlin
pro civ_genprof_plot,strct_fil,fwhmkms,savfil=savfil,fwhm=fwhm,$
                     _extra=extra
;; Display profiles
  if (N_params() lt 2) then begin
     print,'Syntax - '+ $
           'civ_genprof_plot,strct_fil,fwhmkms,[savfil=savfil,'+$
           '   fwhm=fwhm,_extra=extra]'
     return
  endif

  if not keyword_set(fwhm) then fwhm = 2 ; assume STIS
  done = 0

  ;; Read structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_setvplin: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent) 
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  ii = 0
  while not done do begin
     civ_setvplin,strct[ii],vplin,_extra=extra ;e.g., /dvabs

     ;; Make (logarithmic) wavelength array with buffer
     bmx = max(vplin.b)
     nmx = max(vplin.n)
     zlo = min(vplin.zabs,max=zhi)
;     dv = (zhi-zlo)*2.998e5
     dv = 2.998e5 * ((1+zhi)/(1+zlo) - 1)
     wvlo = min(vplin.wrest,max=wvhi)*(1.+zlo-10*bmx/2.998e5)
     wvhi = wvhi*(1+zhi+10*bmx/2.998e5)
     cdelt1 = alog10(1.+fwhmkms/2.998e5)
     npix = ceil(alog10(wvhi/wvlo)/cdelt1)+1
     wave = 10^(alog10(wvlo)+dindgen(npix)*cdelt1)

     ;; Generate Voigt profile and plot
     vpfx = x_voigt(wave,vplin,fwhm=fwhm,_extra=extra) ;e.g., nosmooth
     unq = uniq(vplin.wrest,sort(vplin.wrest))
     ncomp = n_elements(vplin)/n_elements(unq)
     ttl = 'logN = '+string(nmx,format='(f5.2)')+'; b = '+$
           strtrim(round(bmx),2)+' km/s; ncomp = '+strtrim(ncomp,2) + $
           '; dv = '+strtrim(round(dv),2)+' km/s'
     x_splot,wave,vpfx,title=ttl,ymnx=[-0.1,1.1],/block

     if ii eq nstrct-1 then done = 1 else ii = ii+1
  endwhile

end                             ;civ_genprof_plot


function civ_genprof_sub,strct_fil,nabs,zlim,indx=indx,savfil=savfil,$
                         ncomp=ncomp,noblend=noblend,seed=seed,$
                         nlim=nlim,blim=blim
  ;; DO NOT pass config_fil values

  if (N_params() lt 1) then begin
     print,'Syntax - '+ $
           'civ_genprof_sub,strct_fil,nabs,zlim,[indx=indx,'+$
           '   savfil=savfil,ncomp=ncomp,noblend=noblend,seed=seed,'+$
           '   nlim=nlim,blim=blim]'
     return,-1
  endif
  
  if keyword_set(seed) then rseed = seed

  ;; Select the desired set of (multi-component) lines
  ;; ncomp = mean of Poisson distribution of components
  ;; Open structure
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_setvplin: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent) 
  endif else strct = strct_fil
  nstrct = n_elements(strct)
  msk = replicate(1,nstrct)

  ;; Extract based on redshift
  if keyword_set(zlim) then begin
     if n_elements(zlim) eq 2 then $
        gd = where(strct.zabs[0] ge zlim[0] and $
                   strct.zabs[0] le zlim[1],ngd) $ 
     else gd = where(strct.zabs[0] le zlim[0],ngd)
     if ngd eq 0 then $
        stop,'civ_genprof_sub: cannot span the redshift limits'
     strct = strct[gd]
     nstrct = ngd

     ;; Update mask
     msk = msk[gd]

     if ngd lt nabs then begin
        ;; RETURN
        print,'civ_genprof_sub: redshift cut yields <= number requested'
        return,strct
     endif 
  endif 

  ;; Assume each element has same number of ions
  gd = where(strct[0].wrest gt 0.,ngd)
  unq = uniq(strct[0].wrest[gd],sort(strct[0].wrest[gd]))
  wrest = strct[0].wrest[gd[unq]]
  nion = n_elements(wrest)
  ndv = ngd/nion                ;number of repetitions 

  ;; Eliminate blends
  if keyword_set(noblend) then begin
     ;; Store bounds of features
     wobs = fltarr(nion*nstrct)
     wv_lim = fltarr(nion*nstrct,2)
     for ii=0,nstrct-1 do begin
        for jj=0,nion-1 do begin
           wobs[jj*nstrct+ii] = wrest[jj]*(1+strct[ii].zabs[jj])
           wv_lim[jj*nstrct+ii,0] = $
              min(strct[ii].wv_lim[nion*lindgen(ndv)+jj,0])
           wv_lim[jj*nstrct+ii,1] = $
              max(strct[ii].wv_lim[nion*lindgen(ndv)+jj,1])
        endfor                  ;loop nion
     endfor                     ;loop nstrct
     
     ii = 0
     while ii le nion*nstrct-1 do begin
        ;; Centroids not lost in other bounds (less stringent)
        if noblend gt 1 then $
           bd = where(wobs[ii] gt wv_lim[*,0] and $
                      wobs[ii] lt wv_lim[*,1] and $
                      msk[(ii mod nion)] eq 1,nbd) $
                  ;; Wings do not overlap (more stringent)
        else bd = where(wv_lim[ii,1] gt wv_lim[*,0] and $
                        wv_lim[ii,0] lt wv_lim[*,1] and $
                        msk[(ii mod nion)] eq 1,nbd)
        ;; Always includes self
        if nbd gt 1 then begin
           msk[(ii mod nion)] = 0 ; remove system
           ii = (ii mod nion) + 1 ; next system
        endif else ii = ii + 1    ; next line in list
     endwhile                     ;loop nion*nstrct

     gd = where(msk eq 1,ngd)
     if ngd eq 0 then $
        stop,'civ_genprof_sub: cannot unblend lines'
     strct = strct[gd]
     nstrct = ngd

     ;; Update mask
     msk = msk[gd]

     if ngd le nabs then begin
        ;; RETURN
        print,'civ_genprof_sub: deblending yields <= number requested'
        return,strct
     endif 
  endif                         ;/noblend

  ;; Column density cut
  if keyword_set(nlim) then begin
     if n_elements(nlim) eq 2 then $
        gd = where(strct.ncolm[0] ge nlim[0] and $
                   strct.ncolm[1] le nlim[1],ngd) $
     else gd = where(strct.ncolm[0] le nlim[0],ngd)
     if ngd eq 0 then $
        stop,'civ_genprof_sub: cannot span the column density limits'
     strct = strct[gd]
     nstrct = ngd

     ;; Update mask
     msk = msk[gd]

     if ngd le nabs then begin
        ;; RETURN
        print,'civ_genprof_sub: column density cut yields <= number requested'
        return,strct
     endif 
  endif                         ;nlim

  ;; Doppler parameter cut
  if keyword_set(blim) then begin
     brng = fltarr(nstrct,2)
     for ii=0,nstrct-1 do begin
        brng[ii,0] = min(strct[ii].b[0:nion*ndv-1],max=mx)
        brng[ii,1] = mx
     endfor 
     if n_elements(blim) eq 2 then $
        gd = where(brng[*,0] ge blim[0] and brng[*,1] le blim[1],ngd) $ 
     else gd = where(brng[*,1] le blim[0],ngd)
     if ngd eq 0 then $
        stop,'civ_genprof_sub: cannot span the Doppler parameter limits'
     strct = strct[gd]
     nstrct = ngd

     ;; Update mask
     msk = msk[gd]

     if ngd le nabs then begin
        ;; RETURN
        print,'civ_genprof_sub: Doppler parameter cut yields <= ' + $
          'number requested'
        return,strct
     endif 
  endif                         ;blim

  ;; Randomly extract from larger array (assume checked nstrct >= nabs
  ;; earlier) 
  num = randomu(rseed,nstrct,/double)
  sub = sort(num)
  indx = sub[0:nabs-1]
  substrct = strct[indx]

  ;; Consider number of components (ncomp)
  mx = n_elements(substrct[0].zabs)
  if keyword_set(ncomp) then begin
     if ncomp gt ndv then $
        stop,'civ_genprof_sub: desired no. components larger than sample'
     ;; Randomly define number of components 
     num = randomu(rseed,nabs,poisson=ncomp) ;mean = ncomp
     num = long(num) < ndv                  ;cap 
     gd = where(num eq 0,ngd,complement=bd,ncomplement=nbd)
     ;; Want no no-components (num = 0)
     while (ngd ne 0) do begin
        ;; Replace any zeros
        tmp = randomu(rseed,ngd,poisson=ncomp) ; mean = ncomp
        num[gd] = long(tmp) < ndv              ; cap
        gd = where(num eq 0,ngd)
     endwhile 
  endif else begin              ;/ncomp
     gd = lindgen(nabs)         ; all single component
     num = replicate(1,nabs)
  endelse 
  
  ;; Trim components as necessary
  for ii=0,nabs-1 do begin
     if num[ii] eq ndv then continue
     ilo = num[ii]*nion 
     ihi = ndv*nion
     ;; Preserve freq. distr. f(N) and sum ncolm back in
     for jj=0,nion-1 do begin
         rng = lindgen((ihi-ilo)/nion)*nion+ilo+jj
         substrct[ii].ncolm[jj] = 10^substrct[ii].ncolm[jj] + $
           total(10^substrct[ii].ncolm[rng])
         substrct[ii].ncolm[jj] = alog10(substrct[ii].ncolm[jj])
     endfor                     ; loop nion

     ;; Erase
     substrct[ii].ion[ilo:mx-1] = ''
     substrct[ii].wrest[ilo:mx-1] = 0.
     substrct[ii].wv_lim[ilo:mx-1,*] = 0.
     substrct[ii].instr[ilo:mx-1] = 0
     substrct[ii].ew[ilo:mx-1] = 0.
     substrct[ii].sigew[ilo:mx-1] = 0.
     substrct[ii].ncolm[ilo:mx-1] = 0.
     substrct[ii].b[ilo:mx-1] = 0.
     substrct[ii].zabs[ilo:mx-1] = 0.
     substrct[ii].zsig[ilo:mx-1] = 0.
  endfor                        ;loop nabs

  return,substrct
end ; civ_genprof_sub()


function civ_genprof_config,config_fil
cfg = {targlist:'',$
       nlim:fltarr(2), $
       blim:fltarr(2), $
       dvabs:fltarr(2), $
       ndv:0L, $
       nlimobs:fltarr(2), $
       nabsobs:0L, $
       aobs:0., $
       dzunblck:0.,$
       binncolm:0.,$
       binz:0. }

if not keyword_set(config_fil) then return,cfg

readcol,config_fil,val,descript,format='a,a',delimiter='|'
descript = strtrim(descript,2)
nval = n_elements(val)

for ii=0,nval-1 do begin
   case descript[ii] of
      'targlist': cfg.targlist = val[ii]
      'nlim': begin
         prs = strsplit(val[ii],/extract)
         cfg.nlim = [float(prs[0]),float(prs[1])]
      end 
      'blim': begin
         prs = strsplit(val[ii],/extract)
         cfg.blim = [float(prs[0]),float(prs[1])]
      end 
      'dvabs': begin
         prs = strsplit(val[ii],/extract,count=nprs)
         if nprs eq 2 then $
            cfg.dvabs = [float(prs[0]),float(prs[1])] $
         else cfg.dvabs = [-float(prs[0]),float(prs[0])]
      end 
      'ndv': cfg.ndv = long(val[ii])
      'nlimobs': begin
         prs = strsplit(val[ii],/extract)
         cfg.nlimobs = [float(prs[0]),float(prs[1])]
      end 
      'nabsobs': cfg.nabsobs = long(val[ii])
      'aobs': cfg.aobs = float(val[ii])
      'dzunblck': cfg.dzunblck = float(val[ii])
      'binncolm': cfg.binncolm = float(val[ii])
      'binz': cfg.binz = float(val[ii])
   endcase
endfor                          ;loop nval

return,cfg
end ; civ_genprof_config()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_genprof,strct_fil,nabs,zlim,nlim,blim,nlimobs,nabsobs,aobs,$ 
                siiv=siiv,$
                seed=seed,dvabs=dvabs,ndv=ndv,lya=lya,calcew=calcew,$
                config_fil=config_fil,debug=debug
if (N_params() lt 8) and not keyword_set(config_fil) then begin
    print,'Syntax - '+ $
      'civ_genprof,strct_fil,nabs,zlim,nlim,blim,'+$
      '   nlimobs,nabsobs,aobs,[seed=seed,dvabs=dvabs,'+$
      '   ndv=ndv,lya=lya,debug=debug]'
    return
endif                           ; param prompt

if (n_params() lt 3) then stop,'civ_genprof: must set first 3 params'

if keyword_set(siiv) then dblt_name = 'SiIV' $
else dblt_name = 'CIV'

if not keyword_set(config_fil) then begin
    cfg = civ_genprof_config(config_fil)
    cfg.nlim = nlim
    cfg.blim = blim
    cfg.nlimobs = nlimobs
    cfg.nabsobs = nabsobs
    cfg.aobs = aobs
    if keyword_set(dvabs) then begin
       if n_elements(dvabs) eq 2 and dvabs[0] ge 0 then $
          stop,'civ_genprof: cfg.dvabs[0] must be < 0'
       if n_elements(dvabs) eq 1 then cfg.dvabs = [-dvabs,dvabs] $
       else cfg.dvabs = dvabs
    endif                       ; dvabs=
    if keyword_set(ndv) then cfg.ndv = ndv
endif else begin               ; load structure
   if size(config_fil,/type) eq 8 then cfg = config_fil $
   else cfg = civ_genprof_config(config_fil)
endelse 

;; Default Values
if keyword_set(seed) then rseed = seed

if n_elements(zlim) ne 2 then stop,'civ_genprof: zlim must be range'

;; Over-ride
if keyword_set(nlim) then cfg.nlim = nlim
if keyword_set(blim) then cfg.blim = blim

;; Frequency distribution coefficient: 
;; f(N) = m/(dN*dX) where m = # absorbers in bin dN 
;; let f(N) = b*N^(-a)
;; m = int(b*N^(-a),Nmin,Nmax) where limits set by observing range
;; m = b/(1-a)*(Nmax^(1-a) - Nmin^(1-a))
;; b = m*(1-a)*(Nmax^(1-a) - Nmin^(1-a))^(-1)
if abs(cfg.aobs) le 1. then $
  stop,'civ_genprof: power-law exponent cannot <= 1'
if cfg.nlimobs[0] lt 100 or cfg.nlimobs[1] lt 100 then $
  nrng = 10.^cfg.nlimobs else nrng = cfg.nlimobs
eobs = 1.-abs(cfg.aobs)
bobs = cfg.nabsobs*eobs*(nrng[1]^eobs - nrng[0]^eobs)^(-1.) ;not in log
;; Normalize probability distribution
;; P(N)dN = f(N)dN / int(f(N)dN,Nmin,Nmax)
;; int(f(N)dN,Nmin,Nmax) = b/(1-a)*[Nmax^(1-a) - Nmin^(1-a)] = C
nrng = 10.^cfg.nlim
cobs = bobs/eobs*(nrng[1]^eobs-nrng[0]^eobs) ;not in log

;; Rest wavelength (here is where /lya could be added)
civ = dblt_retrieve(dblt_name)
wrest = [civ.wvI,civ.wvII]
strct = civ_instantstrct(nabs)
nion = n_elements(wrest)
getfnam,wrest,fval,nam

;; Store name
strct[*].ion[0:nion-1] = nam
strct[*].wrest[0:nion-1] = wrest

;; Generate and store redshift 
num = randomu(rseed,nabs,/double)
strct[*].zabs[0] = num*(zlim[1]-zlim[0])+zlim[0]
for ii=1,nion-1 do strct[*].zabs[ii] = strct[*].zabs[0]

for ii=0,nion-1 do begin
    print,'civ_genprof: generating profiles for ',nam[ii],'...'
    ;; Random values: (wrest must be in logical order, w/CIV 1548
    ;; first and 1550 2nd and anything else next)
    ;; But for now, only CIV

    ;; Column Density
    if ii ge 1 then begin
        ;; Compare to previous
        prs = strsplit(nam[ii-1],/extract)
        if stregex(nam[ii],strtrim(prs[0],2),/boolean) then test = 1 $
        else test = 0
    endif else test = 0

    if test eq 1 then begin
        ;; Only on partner line (1550)
        strct[*].ncolm[ii] = strct[*].ncolm[ii-1] 
        strct[*].b[ii] = strct[*].b[ii-1]
    endif else begin
        ;; Generate new values
        ;; Column density
        if n_elements(cfg.nlim) eq 2 then begin
            num = randomu(rseed,nabs,/double) ;CDF
            ;; Solve cumulative probability distribution to get Nx
            ;; int(P(N)dN,Nmin,Nx) = R (random number)
            ;; Nx = (R * C*(1-a)/b + Nmin^(1-a))^(1/(1-a))
            ncolm = (num * cobs*eobs/bobs + nrng[0]^eobs)^(1./eobs) ;not log
            strct[*].ncolm[ii] = alog10(ncolm)
           
            ;; Sanity check
            if ((cobs*eobs/bobs + nrng[0]^eobs)^(1./eobs)-nrng[1]) gt 1e10 or $
            ((nrng[0]^eobs)^(1./eobs)-nrng[0]) gt 1e10 then $
               stop,'civ_genprof: CDF not normalized properly'
         endif else strct[*].ncolm[ii] = cfg.nlim[0]
        
        ;; Doppler Parameter
        if n_elements(cfg.blim) eq 2 then begin
            num = randomu(rseed,nabs,/double)
            strct[*].b[ii] = num*(cfg.blim[1]-cfg.blim[0])+cfg.blim[0]
        endif else strct[*].b[ii] = cfg.blim[0]
    endelse                     ; new column density for new ion

    ;; Offset in redshift/velocity
    if keyword_set(cfg.ndv) then begin
        print,'... adding components'
        ;; Need to sample +/- limits for all instances
        ;; Random select in velocity space and transform to redshift
        for jj=1,cfg.ndv-1 do begin
            strct[*].zabs[jj*nion+ii] = strct[*].zabs[ii] ; fill in
            if test eq 1 then begin
                strct[*].zsig[jj*nion+ii] = strct[*].zsig[jj*nion+ii-1]
                strct[*].ncolm[jj*nion+ii] = strct[*].ncolm[jj*nion+ii-1]
                strct[*].b[jj*nion+ii] = strct[*].b[jj*nion+ii-1]
            endif else begin
                num = randomu(rseed,nabs,/double)
                ;; Assume cfg.dvabs[0] < 0 and cfg.dvabs[1] > 0
                if n_elements(cfg.dvabs) eq 2 then zcmp = $
                  ((num*(cfg.dvabs[1]-cfg.dvabs[0])+$
                    cfg.dvabs[0])/2.998e5 + 1)* $
                  (1+strct[*].zabs[ii]) - 1 $
                else zcmp = $
                  (num*cfg.dvabs[0]-0.5*cfg.dvabs[0])/2.998e5*$
                   (1+strct[*].zabs[ii])-1
                strct[*].zsig[jj*nion+ii] = zcmp - strct[*].zabs[ii]


                ;; Column Density
                if n_elements(cfg.nlim) eq 2 then begin
                    num = randomu(rseed,nabs,/double)
                    ;; Component must have column density less than
                    ;; first (main) component in order to make z
                    ;; sensible.  And must reduce column density
                    ;; to preserve freq. distr. (total column density)
                    nwncolm = num*(strct[*].ncolm[jj*nion+ii-2]-cfg.nlim[0])+$
                      cfg.nlim[0]
                    prevncolm = alog10(10^(strct[*].ncolm[jj*nion+ii-2]) - $
                                       10^nwncolm)
                    
                    ;; Check that primary component (0th element) is
                    ;; >= cfg.nlim[0]
                    ;; Otherwise, remaining components just have to be
                    ;; consistent 
                    if jj eq 1 then bd = where(prevncolm lt cfg.nlim[0] or $
                                               finite(prevncolm) eq 0,nbd) $
                    else bd = where(prevncolm lt cfg.nlim[0]*frac or $
                                    finite(prevncolm) eq 0,nbd)
                    frac = 1.d
                    ibd = 0
                    while nbd ne 0 do begin
                       ;; Keep all component colm. density > minimum
                       ;; Unless never converging, then decrease limits
                       num = randomu(rseed,nbd,/double)
                       nwncolm[bd] = num*(strct[bd].ncolm[jj*nion+ii-2] - $
                                          cfg.nlim[0]*frac)+cfg.nlim[0]*frac
                       prevncolm[bd] = $
                          alog10(10^(strct[bd].ncolm[jj*nion+ii-2]) - $
                                 10^nwncolm[bd])
                       if jj eq 1 then bd = where(prevncolm lt cfg.nlim[0] or $
                                                  finite(prevncolm) eq 0,nbd) $
                       else bd = where(prevncolm lt cfg.nlim[0]*frac or $
                                       finite(prevncolm) eq 0,nbd)
                       ibd = ibd + 1
                       if (ibd mod 10) eq 0 then frac = frac - 0.02
                       if frac le 0. then stop,'civ_genprof: impossible'
                    endwhile    ; testing component ncolm
                    strct[*].ncolm[jj*nion+ii-2] = prevncolm
                    strct[*].ncolm[jj*nion+ii] = nwncolm
                endif else strct[*].ncolm[jj*nion+ii] = cfg.nlim[0]

                ;; Doppler Parameter
                if n_elements(cfg.blim) eq 2 then begin
                    num = randomu(rseed,nabs,/double)
                    strct[*].b[jj*nion+ii] = num*(cfg.blim[1]-$
                                                  cfg.blim[0])+cfg.blim[0]
                endif else strct[*].b[jj*nion+ii] = cfg.blim[0]
            endelse             ; generate new component
        endfor                  ; loop num components
    endif                       ; cfg.dvabs set

    print,'... finished ',nam[ii]
    print,''
    if keyword_set(debug) then $
      stop,'civ_genprof debug: finished random z, N, b'
endfor                          ;nion loop


;; Wavelength bounds (depend on z and b)
;; Treating Voigt profile like Gaussian and including +/-3*dopb
if keyword_set(calcew) then print,'civ_genprof: computing limits and EW' $
else print,'civ_genprof: computing limits'

if keyword_set(cfg.ndv) then begin
    ;; Repeat common values throughout civcandstrct
    for jj=1,cfg.ndv-1 do begin     ;skip centroid (already covered)
        ilo = jj*nion
        ihi = (jj+1)*nion - 1
        strct[*].ion[ilo:ihi] = nam
        strct[*].wrest[ilo:ihi] = wrest
        ;; ncolm and b randomized above
    endfor 

    ;; Calculate and store EW
    if keyword_set(calcew) then $
       for ii=0,nion*cfg.ndv-1 do begin
       elem = x_setline(wrest[(ii mod nion)])
       for jj=0,nabs-1 do begin
          elem.n = strct[jj].ncolm[ii]
          elem.b = strct[jj].b[ii]
          strct[jj].ew[ii] = spec_ewlin(elem,/silent)
       endfor                   ; loop nabs
    endfor                      ; loop nion    
    
    strct[*].wv_lim[0:nion*cfg.ndv-1,0] = strct[*].wrest[0:nion*cfg.ndv-1]*$
      (1+strct[*].zabs[0:nion*cfg.ndv-1])*$
      (1-3*strct[*].b[0:nion*cfg.ndv-1]/2.998e5)
    strct[*].wv_lim[0:nion*cfg.ndv-1,1] = strct[*].wrest[0:nion*cfg.ndv-1]*$
      (1+strct[*].zabs[0:nion*cfg.ndv-1])*$
      (1+3*strct[*].b[0:nion*cfg.ndv-1]/2.998e5)
endif else begin
    ;; No components; just copy
    ;; Calculate and store EW
    if keyword_set(calcew) then $
       for ii=0,nion-1 do begin
       elem = x_setline(wrest[ii])
       for jj=0,nabs-1 do begin
          elem.n = strct[jj].ncolm[ii]
          elem.b = strct[jj].b[ii]
          strct[jj].ew[ii] = spec_ewlin(elem,/silent)
       endfor                   ; loop nabs
    endfor                      ; loop nion

    strct[*].wv_lim[0:nion-1,0] = strct[*].wrest[0:nion-1]*$
      (1+strct[*].zabs[0:nion-1])*$
      (1-3*strct[*].b[0:nion-1]/2.998e5)
    strct[*].wv_lim[0:nion-1,1] = strct[*].wrest[0:nion-1]*$
      (1+strct[*].zabs[0:nion-1])*$
      (1+3*strct[*].b[0:nion-1]/2.998e5)
endelse                         ; cfg.ndv not set

if size(strct_fil,/type) eq 7 then begin
    mwrfits,strct,strct_fil,/create,/silent
    print,'civ_genprof: created ',strct_fil
endif else strct_fil = strct
end
