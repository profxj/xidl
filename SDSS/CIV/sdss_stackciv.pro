;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_stackciv.pro               
; Author: Kathy Cooksey                      Date: 24 Apr 2012
; Project: SDSS Metal-line survey 
; Description: Stack spectra of absorbers in given structure.
; Input: 
;   civstrct_fil -- sdsscivstrct structure of absorbers with
;                   doublets in first two array elements
;
; Optional Input:
;   /debug -- print some info and plot at end.
;   /clobber -- overwrite outfil if exists.
;   gwave= -- global, log-linear, rest wavelength array for stack.
;   wvmnx= -- 2-element array of rest wavelength bounds for 
;             global wavelength solution (default: 900--9000 Ang)
;             with SDSS pixel scale.
;   /final -- use *_FINAL tags in civstrct.
;   wvmsk= -- [N,2] array of observed wavelength regions to
;             mask out (e.g., sky lines) but will not exclude
;             doublet(s) wavelength bounds.
;   wvnrm= -- 2-element array of quasar emitted wavelengths to 
;             normalize to 1; recommended: 1450, 1470 Ang. 
;             Strange to combine with /conti= (warning).
;   /conti -- normalize each spectrum by its continuum fit. 
;   cflg= -- overrides the sdsscontistrct value to select continuum
;   cmplt_fil= -- completeness structure(s) used to weight spectra
;                 by completeness fraction; does not work with
;                 /ivarwgt or /litwgt. Must be one file per ndblt.
;   /civobs_corr -- use civstrct to exclude pathlength blocked
;                   by absorbers. Only works with cmplt_fil=.
;   /litwgt -- "light-weight" by normalizing median
;              inverse variance to unity.
;   /ivarwgt -- inverse-variance weighted mean/median.
;   /median -- median flux per pixel instead of weighted mean;
;              re-binned input spectra can be weighted by completeness,
;              /ivarwgt, or /litwgt.
;   percentile= -- 2-element array with percentiles to store for
;                  median stack [default: 15.9%, 84.1% (like +/-1
;                  Gaussian sigma]
;   /reflux -- re-collapse the stack with sdss_stackciv_stack() and
;              does what /refit does
;   /refit -- re-normalize and find lines; calls sdss_stackciv_fitconti()
;   /reerr -- re-analyze the error from Monte-Carlo bootstrapping with
;             function sdss_stackciv_errmc() (unless /qerr also set)
;   /qerr -- quick error analysis (not reasonable for median)
;   ndblt= -- number of doublets that were selected for the sample 
;             and so have to explicitly loop for wvmsk= and 
;             cmplt_fil=. Doublets assumped to be in elements with
;             even index (e.g., 2*lindgen(n_elements(ndblt))).
;             [default: 1]
;   _extra= -- other options sent to other routines
;
; Output: 
;   outfil -- name of output file for SDSS-formatted stacked
;             spectrum in 0th extension, sdsscontistrct in 1st
;   
; Optional Output:
;
; Example:
;    IDL>civstr=sdss_getcivstrct(/default,zlim=[2.9,3.1],ewlim=[0.6,!values.f_infinity])
;    IDL>sdss_stackciv,civstr,'civ2p9z3p1.fit',/debug,/clobber,cmplt_fil='inputs/SNRge4/gz/cmpltrec_CIV_1z6cum.fit',wvmsk=transpose([[5575.,5583],[5888,5895],[6296,6308],[6862,6871]]),wvnrm=[1450.,1470.]
;
; History:
;   25 Apr 2012  created by KLC
;   21 Jul 2014  change /nowgt to /litwgt (flip usage), KLC
;   27 Jul 2014  weighted median implemented, KLC
;                also ndblt= functionality, KLC
;   29 Jul 2014  use MAD instead of std err on median, KLC
;   30 Jul 2014  MAD doesn't work; hack something else, KLC
;   31 Jul 2014  Create functions to stack and MC error, KLC
;   11 Jul 2017  Rebin variance; change e.g., sigma ne 0. to gt 0.,
;                enable ivarwgt, change median /qerr, KLC
;   21 Jul 2017  Revamp to enable /sigew in *errmc(), KLC
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@sdss_fndlin                    ; resolve sdss_fndlin_fitspl()

function sdss_stackciv_fitconti, spec_fil, wave=wave, nlmax=nlmax,$
                                 lin_fil=lin_fil, dvlin=dvlin, _extra=extra
  ;; Use as many other functions as possible
  if n_params() ne 1 then begin
     print,'sdss_stackciv_fitconti( spec_fil, [wave=, nlmax=, lin_fil=, dvlin=, _extra=])'
     return,-1
  endif 

  ;; Defaults
  if not keyword_set(lin_fil) then $
     lin_fil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls_stack.lst'
  linstr = x_setllst(lin_fil,0) ; lls.lst for
  nlin = (size(linstr,/dim))[0] > 1
  if not keyword_set(dvlin) then dvlin = 600. ; km/s like sdss_dblfitconti
  c = 299792.458                              ; km/s
  if not keyword_set(nlmax) then nlmax = 500

  if size(spec_fil,/type) eq 7 then begin
     parse_sdss,spec_fil,flux,wave,sig=sigma,npix=npix
  endif else begin
     flux = spec_fil[*,0]
     sigma = spec_fil[*,2]
     if not keyword_set(wave) then $
        stop,'sdss_stackciv_fitconti stop: must set wavelength array.'
     npix = (size(flux,/dim))[0] > 1
  endelse

  ;; May need a larger continuum structure
  cstrct = sdss_expandcontistrct(npix+1,nlin=nlmax) ; stores lsnr
  cstrct.cflg = sdss_getcflg(/spl)
  cindx = sdss_getcflg(/spl,/index)

  gdpix = where(sigma gt 0.)
  cstrct.snr_conv[gdpix,cindx] = flux[gdpix]/sigma[gdpix] ; sleight of hand

  ;; Mask out lines and record "centroids"
  premask = replicate(1,npix)
  dlim = dvlin / c
  for ll=0,nlin-1 do begin
     if abs(linstr[ll].wave-1215.6701) lt 1e-4 then $ ; wider for Lya
        sub = where(abs(wave-linstr[ll].wave) lt 2*dlim*linstr[ll].wave) $
     else sub = where(abs(wave-linstr[ll].wave) lt dlim*linstr[ll].wave)
     if sub[0] ne -1 then premask[sub] = 0
  endfor

  ;; _extra= includes everyn=, sset=, maxrej=, lower=, upper=, nord=,
  ;; /groupbadpix, bsplmask=, /debug, /silent
  ;; Force bsline_iterfit() to keep premask b/c we *know* where the
  ;; lines are to avoid.
  everyn_lcl = 75 ; wider breakpoint spacing; "stiffer"
  if keyword_set(_extra) then begin
     extra0 = extra             ; preserve originl
     
     ;; Unless the user explicitly passed in values, pass in values
     ;; better suited to stack normalization than
     ;; sdss_fndlin_fitspline() defaults
     tags = tag_names(extra)
     if not stregex(tags,'EVERYN',/boolean) then $
        extra = create_struct(extra0,'EVERYN',everyn_lcl)  
  endif else extra = {EVERYN:everyn_lcl}
  cstrct = sdss_fndlin_fitspline(wave, flux, sigma, 0.0, $
                                 premask=premask, /nopca, /sticky, $
                                 cstrct_fil=cstrct, _extra=extra)
  if keyword_set(extra0) then extra = extra0 $ ; restore
  else undefine, extra                         ; remove from downstream

  ;; Find absorption lines
  ;; _extra= includes lsnr=, /debug
  mask = sigma * 0.
  mask[gdpix] = 1. 
  fconti = cstrct.conti[0:cstrct.npix-1,cindx]
  ivconti = 1./fconti
  fx = flux * ivconti
  sig = sdss_calcnormerr(flux, sigma, cstrct, cflg=cstrct.cflg, $
                         baderrval=9.e9)
  ;; _extra includes lsnr=, /debug, nfind=
  cstrct = sdss_fndlin_srch(cstrct, cstrct.cflg, wave=wave, flux=fx, $
                            sigma=sig, mask=mask, _extra=extra)

  ;; Check that close lines are separated
  srt = sort(linstr.wave)
  dvlls_lo = abs(c*(shift(linstr[srt].wave,1)-linstr[srt].wave)/$
                 linstr[srt].wave)
  dvlls_hi = abs(c*(shift(linstr[srt].wave,-1)-linstr[srt].wave)/$
                 linstr[srt].wave)
  chk = where(dvlls_lo lt 800. or dvlls_hi lt 800.,nchk) ; gets MgII
  ;; nchk will be even; and order will be e.g.,
  ;;   dvlls_lo  dvlls_hi      wave  ion
  ;; -1275.9026 370.27962 1036.3367  CII
  ;; -369.82284 13398.364 1037.6167  OVI
  ;; -26802.108 723.75852 1190.4158  SiII
  ;; -722.01543 1864.1660 1193.2897  SiII*
  ;; -9611.0882 506.88761 1302.1685  OI
  ;; -506.03201 6932.3651 1304.3702  SiII
  ;; -4161.0136 498.62298 1548.1950  CIV
  ;; -497.79504 11150.821 1550.7700  CIV
  ;; -21032.051 769.64922 2796.3520  MgII
  ;; -767.67838 5286.0843 2803.5310  MgII
  for ii=0,nchk*0.5-1 do begin
     ilo = srt[chk[2*ii]]
     ihi = srt[chk[2*ii]+1]
     mnlo = min(cstrct.centroid[*,cindx]-linstr[ilo].wave,imnlo,/abs)
     mnhi = min(cstrct.centroid[*,cindx]-linstr[ihi].wave,imnhi,/abs)
     ;; Must check whether found the right line; and check whether
     ;; this stil divides
     if imnlo eq imnhi and (abs(mnlo/linstr[ilo].wave*c) lt 800. and abs(mnhi/linstr[ihi].wave*c) lt 800.) then begin
        ;; Definitely was detected based on lsnr= cut but it's
        ;; the same line and don't want that; divide it up
        cstrct.centroid[imnlo,cindx] = linstr[ilo].wave
        if cstrct.ncent[cindx] eq nlmax then begin
           ;; replace from the back, which might blow away ones I
           ;; want... 
           print,'sdss_stackciv_fitconti(): need more space in abslin strct. Expanding.'
           tmp_cstrct = sdss_expandcontistrct(npix+1,nlin=nlmax+nchk*0.5) ; big enough
           cstrct = sdss_cpstrct(cstrct, tmp_cstrct)
        endif 
        ;; add until full
        cstrct.centroid[cstrct.ncent[cindx],cindx] = linstr[ihi].wave
        cstrct.ncent[cindx]++
     endif                      ; imnlo == imhi
  endfor                        ; loop ii=nchk*0.5
  
  ;; Another check might be that all lines in linstr have a
  ;; corresponding centroid... 
  if cstrct.ncent[cindx] eq 0 then $
     print,'sdss_stackciv_fitconti(): WARNING!!!: should find lines.'

  ;; Calculate EW
  ;; _extra= includes /keepwvlim, /debug, /plot
  cstrct = sdss_fndlin_calcew(cstrct, wave=wave, flux=flux, sigma=sigma,$
                              _extra=extra)
  
  return, cstrct
end                             ; sdss_stackciv_fitconti()


function sdss_stackciv_stack, gstrct, median=median, percentile=percentile, $
                              fonly=fonly

  if n_params() ne 1 then begin
     print,'Syntax - sdss_stackciv_stack( gstrct, [/median, percentile=, '
     print,'                             /fonly] )'
     return, -1
  endif

  ngpix = (size(gstrct.gwave,/dim))[0] > 1

  ;; Collapse
  fdat = fltarr(ngpix,7,/nozero)
  fdat[*,1] = total(gstrct.gnspec,2)
  ;; Going to follow SDSS (pre-SDSS-III) spSpec format, differences
  ;; noted in square brackets (ext=0, what parse_sdss can handle):
  ;; fdat[*,0] = flux
  ;; fdat[*,1] = number of spectra per pixel [not continuum-subtracted
  ;;             spectrum] 
  ;; fdat[*,2] = 1-sigma uncertainty; weighted error
  ;;             propogation (very small) if mean stack or
  ;;             median absolute deviation (MAD; badness!) if median stack
  ;; fdat[*,3] = weight (or median weight) [not mask array]
  ;; fdat[*,4] = continuum [not in SDSS spSpec]
  ;; fdat[*,5] = lower percentile [not in SDSS spSpec]
  ;; fdat[*,6] = upper percentile [not in SDSS spSpec]
  if not keyword_set(percentile) then begin
     ;; Consider these like "errors" but just alternate percentiles
     ;; to save (where median is 50%)
     percentile = replicate(gauss_pdf(1.),2) ; ~84.1%; +1 Gaussian sigma
     percentile[0] = 1. - percentile[0]      ; ~15.9%; -1 Gaussian sigma
;     percentile = [0.25,0.75]
  endif

  if keyword_set(median) then begin

     ;; Estimate error straight from the statistics
     for pp=0L,ngpix-1 do begin
        ;; Use only CIV that contribute meaninful to this pixel
        gdciv = where(finite(gstrct.gflux[pp,*]),ngdciv) 

        if ngdciv eq 0 then begin
           fdat[pp,[0,2,3,5,6]] = 0. ; must instantiate flux
           ;; [4] is continuum
        endif else begin

           if ngdciv eq 1 then begin
              fdat[pp,0] = gstrct.gflux[pp,gdciv]
              ;; [1] is number of spec set above
              if keyword_set(fonly) then continue; save time

              fdat[pp,2] = sqrt(gstrct.gvariance[pp,gdciv]) ; lack of better to do
              fdat[pp,3] = gstrct.gweight[pp,gdciv]
              ;; [4] is continuum
              fdat[pp,[5,6]] = 0 ; percentiles not calculable
           endif else begin

              fdat[pp,0] = sdss_medianw(reform(gstrct.gflux[pp,gdciv]),$
                                        reform(gstrct.gweight[pp,gdciv]),/even,$
                                        /silent)

              if keyword_set(fonly) then continue ; save time


              ;; median weight per pix... hope to be something related to
              ;; completeness or 1 
              fdat[pp,3] = median(reform(gstrct.gweight[pp,gdciv]),/even) 

              ;; KLC changed 11 Jul 2017 to just be median absolute
              ;; deviation (MAD)---weighted as needed; now works with
              ;; /conti,/extra
              fdat[pp,2] = sdss_medianw($
                           reform(abs(gstrct.gflux[pp,gdciv]-fdat[pp,0])), $
                           reform(gstrct.gweight[pp,gdciv]),/even,/silent) ;$

              ;; Percentiles for error
              fdat[pp,5] = sdss_medianw(reform(gstrct.gflux[pp,gdciv]),$
                                        reform(gstrct.gweight[pp,gdciv]),$
                                        /even,percent=percentile[0],$
                                        /silent)
              fdat[pp,6] = sdss_medianw(reform(gstrct.gflux[pp,gdciv]),$
                                        reform(gstrct.gweight[pp,gdciv]),$
                                        /even,percent=percentile[1],$
                                        /silent)

           endelse              ; more than one value

        endelse                 ; no good values 
     endfor                     ; loop=ngpix

  endif else begin

     ;; This is the definition of a weighted average, where the
     ;; weight could be 1's and hence just an average 
     ;; fbar = sum(fi*wi)/sum(wi)
     ;; Propogate errors and var(fbar) = sum(var(fi)*wi^2)/(sum(wi))^2
     ;; If /ivarwgt and /litwgt weren't set in sdss_stackciv, this 
     ;; reduces to the regular mean.
     twgt = 1. / total(gstrct.gweight,2)
     fdat[*,0] = total(gstrct.gflux*gstrct.gweight,2,/nan) * twgt
     ;; [1] is number of spec set above
     if not keyword_set(fonly) then begin ; save time
        fdat[*,2] = total(gstrct.gvariance*gstrct.gweight^2,2,/nan) * twgt^2
        fdat[*,2] = sqrt(fdat[*,2])
        fdat[*,3] = total(gstrct.gweight,2) ; weight
        ;; [4] is continuum
        
        for pp=0L,ngpix-1 do begin
           ;; Use only CIV that contribute meaninful to this pixel
           gdciv = where(finite(gstrct.gflux[pp,*]),ngdciv) 

           if ngdciv le 1 then fdat[pp,[5,6]] = 0. $ ; like /median
           else begin
              ;; Going to be real stupid for low ngdciv... basically
              ;; assuming almost line or something (and this is not
              ;; quite even with median, which uses the weights)
              srt = gdciv[sort(gstrct.gflux[pp,gdciv])]
              cumhist = lindgen(ngdciv)/float(ngdciv)
              fprob = interpol(gstrct.gflux[pp,srt], cumhist, percentile)
              fdat[pp,5] = fprob[0]
              fdat[pp,6] = fprob[1]
           endelse              ; ngdciv > 0
        endfor                  ; pp=ngpix
        
     endif

  endelse 

  if not keyword_set(fonly) then $
     gstrct.percentile = percentile ; posterity
  
  ;; fdat[*,4] = conti in sdss_stackciv (main)

  return, fdat
end                             ; sdss_stackciv_stack()



function sdss_stackciv_errmc, fdat, gstrct0, fexcl=fexcl, sigew=sigew, $
                              cstrct_resmpl=cstrct_resmpl, $
                              niter=niter, seed=seed, oseed=oseed, $
                              minmax=minmax, _extra=extra
  if n_params() ne 2 then begin
     print,'Syntax -- sdss_stackciv_errmc(fdat, gstrct0, [fexcl=, /sigew, '
     print,'                              niter=, seed=, oseed=, _extra=])'
     return, -1
  endif

  ;; fexcl= fraction to exclude (if < 1) else total number

  if not keyword_set(niter) then niter = 1000
  ;; /minmax stores the resampled flux array when top-fexcl of the
  ;; input sample (gstrct0.ewabs) are excluded in the second-to-last
  ;; spectrum in gstrct0.gflux_resample and similar for the
  ;; bottom-fexcl in the last spectrum 
  if keyword_set(minmax) then niter_xtr = 2 else niter_xtr = 0
  if not keyword_set(fexcl) then fexcl = 0.25 ; 25%
  if keyword_set(seed) then oseed = seed
  
  ngpix = n_elements(fdat[*,0])
  if ngpix ne n_elements(gstrct0.gwave) then $
     stop,'sdss_stackciv_errmc() stop: input have different number of pixels'
  if keyword_set(sigew) then begin
     if size(sigew,/type) eq 8 then cstrct0 = sigew $ ; assume
     else cstrct0 = $
        sdss_stackciv_fitconti(fdat, wave=gstrct0.gwave, _extra=extra)
     ;; _extra includes lin_fil=, dvlin=, and lots of other stuff
     ;; (see function)
     
     ;; Force lsnr=0.1 because we want to find lines maximally (even
     ;; with hokey error estimate, as happens with median-stacks w/o
     ;; this function)
     lsnr_lcl = 0.1
     if keyword_set(_extra) then begin
        extra0 = extra          ;  preserve original
        tags = tag_names(extra)
        test = where(stregex(extra,'LSNR',/boolean),ntest) ; might be LSNR2=
        if ntest ne 0 then extra.lsnr = lsnr_lcl $
        else extra = create_struct(extra0,'LSNR',lsnr_lcl) ; append
     endif else extra = {LSNR:lsnr_lcl}
  endif
  
  nspec = max(fdat[*,1])         ; largest number of spectra per pixel
  if fexcl lt 1 then nexcl = round(fexcl*nspec) $
  else nexcl = fexcl

  ;; Need to save the values of the flux from the stack
  fx_resmpl = fltarr(ngpix,niter+niter_xtr,/nozero) ; last two are extremum

  ;; Bootstrapping
  for ii=0L,niter+niter_xtr-1 do begin

     gstrct = gstrct0           ; reset

     if ii ge niter then begin
        if (size(gstrct0.ewabs,/dim))[0] gt 1 then $
           stop,'sdss_stackciv_errmc() stop: Monte-Carlo bootstrap does not work for more than one doublet-based stacks'
        srt = sort(gstrct0.ewabs[0,*]) ; don't know how to handle doublets yet
        if ii eq niter then begin
           ;; exclude *top* quartile of input absorbers and resample
           ;; for remaining 75% to probe *lower* extreme of equivalent
           ;; widths
           srt = srt[0:nspec-nexcl-1]
        endif else begin
           ;; ii == niter+1 --> exclude *bottom* quartile
           srt = srt[nexcl:nspec-1]
        endelse
        iexcl = srt[sort(randomu(oseed,nspec-nexcl))] ; indexes gstrct0
        iresmpl = srt[sort(randomu(oseed,nspec-nexcl))]
     endif else begin
        ;; Exclude and re-sample with replacement
        iexcl = sort(randomu(oseed,nspec))
        iresmpl = sort(randomu(oseed,nspec))
     endelse

     
     ;; Replace (have to loop; don't know why)
     for rr=0L,nexcl-1 do begin
        gstrct.ewabs[iexcl[rr]] = gstrct0.ewabs[iresmpl[rr]]
        gstrct.zabs[iexcl[rr]] = gstrct0.zabs[iresmpl[rr]]
        gstrct.gflux[*,iexcl[rr]] = gstrct0.gflux[*,iresmpl[rr]]
        gstrct.gvariance[*,iexcl[rr]] = gstrct0.gvariance[*,iresmpl[rr]]
        gstrct.gweight[*,iexcl[rr]] = gstrct0.gweight[*,iresmpl[rr]]
        gstrct.gnspec[*,iexcl[rr]] = gstrct0.gnspec[*,iresmpl[rr]]
        gstrct.medsnr_spec[iexcl[rr],*] = gstrct0.medsnr_spec[iresmpl[rr],*]
     endfor                     ; loop rr=nloop

     new_fdat = sdss_stackciv_stack(gstrct, /fonly, median=gstrct.median)

     if keyword_set(sigew) then begin
        ;; need error array for fitting continuum and finding lines
        ;; so have to get it
        new_fdat[*,2] = sdss_stackciv_errmc(new_fdat, gstrct, $
                                            sigew=0, minmax=0, $
                                            fexcl=fexcl, niter=niter, $
                                            seed=oseed, oseed=oseed)
        

        ;; _extra includes lin_fil=, dvlin=, and lots of other stuff
        ;; (see function)... but definitely includes LSNR = lsnr_lcl
        cstrct = sdss_stackciv_fitconti(new_fdat, wave=gstrct0.gwave, $w
                                        _extra=extra)
        if ii eq 0 then cstrct_resmpl = cstrct $
        else begin
           cstrct_resmpl = [cstrct_resmpl,cstrct] ; wonder how RAM intensive this is
           if (ii mod 50) eq 0 then begin
              print,'sdss_stackciv_errmc(/sigew): iter = ',ii+1
              save,/all,filename='sdss_stackciv_errmc_sigew'+$
                   string(ii+1,format='(i05)')+'.sav'
           endif
        endelse

     endif                                          ; /sigew
        
     fx_resmpl[*,ii] = new_fdat[*,0] ; save just the flux
     
  endfor                        ; loop ii=niter+2

  ;; Analyze just the distribution of the flux
  if keyword_set(gstrct0.median) then begin
     ;; now use the MAD
     tmp = reform(fdat[*,0])
     flux = rebin(tmp,ngpix,niter)                                     ; dupliclate to get right dimensionality
     error = median( abs(fx_resmpl[*,0:niter-1] - flux), dim=2, /even) ; excl extrema
  endif else $
     error = stddev( fx_resmpl[*,0:niter-1], dim=2 )

  ;; Should we do the percentiles? 
  
  oseed = oseed[0]              ; for next iteration

  ;; Modify output (cannot save structures within structure, so
  ;; cstrct_resmpl returned separately)
  gstrct0 = create_struct(gstrct0,'GFLUX_RESAMPLE',fx_resmpl)
  
;  x_splot,gstrct.gwave,fdat[*,0],ytwo=fdat[*,2],psym1=10,psym2=10,ythr=error,psym3=10,/block
;
;  idx = lindgen(10)
;  x_splot,gstrct.gwave,fdat[*,0],psym1=10,ytwo=fx_resmpl[*,idx[0]],psym2=10,ythr=fx_resmpl[*,idx[1]],psym3=10,yfou=fx_resmpl[*,idx[2]],psym4=10,yfiv=fx_resmpl[*,idx[3]],psym5=10,ysix=fx_resmpl[*,idx[4]],psym6=10,ysev=fx_resmpl[*,idx[5]],psym7=10,yeig=fx_resmpl[*,idx[6]],psym8=10

  if keyword_set(sigew) then stop
  return, error

end                             ; sdss_stackciv_errmc()


pro sdss_stackciv_jackknife, stack_fil, oroot, fjk=fjk, clobber=clobber, _extra=extra
  ;; "Deleted-d jackknife" or "group jackknife" which is more robust
  ;; for non-smooth estimators like median
  ;; (http://www.stat.berkeley.edu/~hhuang/STAT152/Jackknife-Bootstrap.pdf
  ;; and
  ;; https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/9-6.pdf).
  ;; But really this isn't a true "deleted-d jackknife" because
  ;; it does not excise *every* subsample of size f_jk*n_obj.
  ;; This instead systematically removes f_jk*n_obj in increasing rest
  ;; equivalent widths.
  if n_params() ne 2 then begin
     print,'Syntax -- sdss_stackciv_jackknife, stack_fil, oroot, [fjk=, /clobber, _extra=]'
     return
  endif

  ;; Timing
  tstart = systime(/seconds)
  
  ;; Read in (shortcut on the processing)
  fdat0 = xmrdfits(stack_fil,0,hdr0,/silent)
  cstrct0 = xmrdfits(stack_fil,1,/silent)
  gstrct0 = xmrdfits(stack_fil,2,/silent)
  srt0 = sort(gstrct0.ewabs)

  nabs = fxpar(hdr0, 'NABS')    ; = (size(gstrct0.ewabs,/dim))[0]
  if not keyword_set(fjk) then begin
     ;; Optimal choice is sqrt(nabs) < nexcl < nabs
     nexcl = ceil(sqrt(nabs))   ; ceil() to ensure above lower limit
     fjk = nexcl/float(nabs)    ; posterity, should be roughly 1/sqrt(nabs) 
  endif else $
     nexcl = round(fjk*nabs)    ; round() fairest down the line
  nloop = ceil(1/fjk)           ; definitely make it to the end
  istart = 0L

  for ii=0L,nloop-1 do begin
     istop = (istart + nexcl - 1)
     if ii eq nloop-1 then istop = nabs - 1 ; all the way to end

     ;; Extract information
     if istart eq 0 then sub = srt0[istop+1:*] $ ; exclude istart:istop
     else if istop eq nabs-1 then sub = srt0[0:istart-1] $
     else sub = [srt0[0:istart-1],srt0[istop+1:*]]
     nsub = nabs - (istop-istart+1) 

     ewmin = min(gstrct0.ewabs[srt0[istart:istop]],max=ewmax)
     ofil = oroot+string(ewmin,ewmax,istop-istart+1,$
                         format="('_',f4.2,'w',f4.2,'_n',i04)")+'.fit'
     test = file_search(ofil+'*',count=ntest)
     if ntest ne 0 and not keyword_set(clobber) then begin
        print,'sdss_stackciv_jackknife: will not clobber ',ofil
        istart = istop + 1      ; next loop 
        continue                ; skip write
     endif
     
     gstrct = {$                        ; this has to match what's in gstrct
              median:gstrct0.median,$    ; whatever sdss_stackciv_errmc() needs
              percentile:gstrct0.percentile,$
              ewabs:gstrct0.ewabs[sub],$
              zabs:gstrct0.zabs[sub],$
              gwave:gstrct0.gwave,$
              gflux:gstrct0.gflux[*,sub],$
              medsnr_spec:gstrct0.medsnr_spec[sub,*],$ ; [<f>,<sig>,<f/sig>]
              gvariance:gstrct0.gvariance[*,sub],$
              gweight:gstrct0.gweight[*,sub],$
              gnspec:gstrct0.gnspec[*,sub] $ ; counter and may divide
              }                             ; will end up gflux_resample too
     
     
     ;; _extra= includes percentile=, /fonly
     fdat = sdss_stackciv_stack(gstrct, median=gstrct.median, $
                                _extra=extra)
     ;; _extra= includes fexcl=, niter=, seed=, oseed= and others
     ;; passed through (typically when /sigew)
     fdat[*,2] = sdss_stackciv_errmc(fdat, gstrct, sigew=0, minmax=0, $
                                     _extra=extra)
     
     ;; _extra includes lin_fil=, dvlin=, and lots of other stuff
     ;; (see function)
     cstrct = sdss_stackciv_fitconti(fdat, wave=gstrct.gwave, $
                                     _extra=extra)
          
     ;; Update header
     hdr = hdr0
     sxaddpar,hdr,'NABS',nsub   ; not whole set
     sxaddpar,hdr,'ZMED',median(gstrct.zabs,/even)
     sxaddpar,hdr,'ZMEAN',mean(gstrct.zabs)
     sxaddpar,hdr,'ZMIN',min(gstrct.zabs,max=mx)
     sxaddpar,hdr,'ZMAX',mx
     sxaddpar,hdr,'EWMED',median(gstrct.ewabs,/even)
     sxaddpar,hdr,'EWMEAN',mean(gstrct.ewabs)
     sxaddpar,hdr,'EWMIN',min(gstrct.ewabs,max=mx) ; incl. in stack
     sxaddpar,hdr,'EWMAX',mx
     sxaddpar,hdr,'EWAVE_JK',mean(gstrct0.ewabs[srt0[istart:istop]]),$
              'Mean EW excluded in jackknife' ; new keywords; "EWMEAN_JK" too long
     sxaddpar,hdr,'EWMED_JK',median(gstrct0.ewabs[srt0[istart:istop]],/even),$
              'Median EW excluded in jackknife'
     sxaddpar,hdr,'EWMIN_JK',ewmin,'Min EW excluded in jackknife'
     sxaddpar,hdr,'EWMAX_JK',ewmax,'Max EW excluded in jackknife'

     ;; Write file (must match sdss_stackciv output)
     mwrfits,fdat,ofil,hdr,/create,/silent ; ext = 0
     mwrfits,cstrct,ofil,/silent           ; ext = 1
     mwrfits,gstrct,ofil,/silent           ; ext = 2
     spawn,'gzip -f '+ofil
     print,'sdss_stackciv_jackknife: created ',ofil

     ;; Setup for next loop
     istart = istop + 1         ; non-overlapping
  endfor                        ; loop ii=nloop

  ;; Final messages
  if not keyword_set(silent) then $
     print, 'sdss_stackciv_jackknife: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_stackciv: Elapsed clock time for '+strtrim(nloop,2)+$
        ' iterations (m) = ',dt/60.
  
end                             ; sdss_stackciv_jackknife


function sdss_stackciv_jackknife_stats, stack_list, refstack_fil, $
                                        lin_fil=lin_fil, wrt_ref=wrt_ref,$
                                        _extra=extra
  if n_params() ne 2 then begin
     print,'Syntax -- sdss_stackciv_jackknife_stats(stack_list, refstack_fil, [lin_fil=, /wrt_ref, _extra=])'
     return,-1
  endif

  if not keyword_set(lin_fil) then $
     lin_fil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls_stack.lst'
  linstr = x_setllst(lin_fil,0) ; lls.lst for
  nlin = (size(linstr,/dim))[0] > 1

  if n_elements(stack_fil) eq 1 then $
     readcol,stack_list,stack_fil,format='(a)',/silent $
  else stack_fil = stack_list
  nfil = (size(stack_fil,/dim))[0]

  ;; Setup output (dynamically sized)
  rslt = {stack_fil:stack_fil, $
          nabs:lonarr(nfil), $    ; will use average
          ewave:fltarr(nfil,2), $ ; [mean,median] of excluded in jackknife
          ewlim:fltarr(nfil,2), $ ; [min,max] of excluded
          median:-1, $            ; 0: mean, 1: median
          lin_fil:lin_fil, $
          ion:linstr.name, $
          wrest:linstr.wave, $
          ewref:fltarr(nlin,2), $ ; [EW, variance] --> sigma at end
          nref:0L, $
          percentile:fltarr(2), $ ; be consistent with sdss_stackciv
          wrt_ref:keyword_set(wrt_ref), $ 
          ewcdf:fltarr(nlin,3), $ ; [ewref or mean/median(EW), low-, high-sigma]
          ewion:fltarr(nlin,nfil), $
          sigewion:fltarr(nlin,2,nfil), $
          ewion_excl:fltarr(nlin,nfil,2), $ ; excluding current; mean & variance --> sigma at end
          ewion_est:fltarr(nlin,2), $       ; estimate of mean & variance -- sigma at end
          ewion_bias:fltarr(nlin,2), $      ; bias estimate of above
          ewion_estcorr:fltarr(nlin,2) $    ; bias-corrected estimators
         }

  ;; _extra= includes dwvtol=
  stackstr_ref = sdss_mkstacksumm(refstack_fil, lin_fil=lin_fil, _extra=extra)
  rslt.ewref[*,0] = stackstr_ref.ew
  rslt.ewref[*,1] = stackstr_ref.sigew^2
  rslt.nref = stackstr_ref.nabs
  if keyword_set(wrt_ref) then $
     rslt.ewcdf[*,0] = rslt.ewref[*,0]
  rslt.percentile = stackstr_ref.percentile 
  stackstr = sdss_mkstacksumm(stack_fil, lin_fil=lin_fil, _extra=extra)
  rslt.nabs = stackstr.nabs
  rslt.ewave = transpose(stackstr.ewave[2:3]) ; from EWAVE_JK and EWMED_JK header keywords
  rslt.ewlim = transpose(stackstr.ewlim[2:3]) ; from EWMIN_JK and EWMAX_JK

  ;; Sanity check on uniformity (stackstr uniformity handled by
  ;; sdss_mkstacksumm())
  if stackstr_ref.median ne stackstr[0].median then begin
     print,'sdss_stackciv_jackknife_stats(): ERROR!!! reference stack and jackknife stacks not all mean or median; exiting.'
     return,-1
  endif
  rslt.median = stackstr_ref.median
  if rslt.median then iave = 1 else iave = 0

  if stackstr_ref.percentile[0] ne stackstr[0].percentile[0] or $
     stackstr_ref.percentile[1] ne stackstr[0].percentile[1] then begin
     print,'sdss_stackciv_jackknife_stats(): ERROR!!! reference stack and jackknife stacks do not use same percentile; exiting.'
     return,-1
  endif
  rslt.percentile = stackstr_ref.percentile
        

  ;; Aggregate statistics per line
  for ll=0,nlin-1 do begin
     ;; Using median of objects in stack (zabs[1]) and extracting
     ;; _extra= includes zrng=, dztol=, dwvtol=
     ;; Use skip_null=0 to (1) ensure all lines accounted for in order
     ;; and (2) zeros *are* information for the estimators
     iondat = sdss_getstackdat(stackstr, stackstr.zave[iave], linstr[ll].name, $
                               skip_null=0, /nosrt, $
                               dztol=max(stackstr.zabs[iave],min=mn)-mn,$
                               _extra=extra)

     rslt.ewion[ll,*] = iondat.ydat
     if not keyword_set(wrt_ref) then begin
        ;; weighted by number bin
        wgt = rslt.nref-rslt.nabs
        val = transpose(rslt.ewion[ll,*])
        if rslt.median then $
           rslt.ewcdf[ll,0] = sdss_medianw(val,wgt,/even) $
        else $
           rslt.ewcdf[ll,0] = total(val*wgt)/total(wgt)
     endif 
     rslt.sigewion[ll,0,*] = iondat.sigydat[*,0] ; likely symmetric
     rslt.sigewion[ll,1,*] = iondat.sigydat[*,1]

     ;; Look at statistics excluding one stack each time
     ;; Following Wikipedia entry but also see _The Jackknife, the
     ;; Bootstrap, and Other Resampling Plans_ by Efron (1982):
     for ff=0,nfil-1 do begin
        ;; find subsample less each sample
        if ff eq 0 then rng = 1 + lindgen(nfil-1) $ ; [1:*]
        else begin
           if ff eq nfil-1 then rng = lindgen(nfil-1) $     ; [0:nfil-1]
           else rng = [lindgen(ff),ff+1+lindgen(nfil-ff-1)] ; gap
        endelse

        ;; Estimator (e.g., mean) <x_(i)> of excluding i-th subsample
        ;; and variance estimation (e.g., MAD) (this latter may not be
        ;; statistically sound)
        if rslt.median then begin
           ;; <x_(i)>
           rslt.ewion_excl[ll,ff,0] = median(rslt.ewion[ll,rng],/even)
           ;; <varx_(i)>
           rslt.ewion_excl[ll,ff,1] = median(abs(rslt.ewion[ll,rng] - $
                                                 rslt.ewion_excl[ll,ff,0]),/even)
        endif else begin
           ;; <x_(i)> (mean() will divide by (nfil-1))
           rslt.ewion_excl[ll,ff,0] = mean(rslt.ewion[ll,rng]) 
           ;; <varx_(i)> (nfil-2 to be variance of *sample* population)
           rslt.ewion_excl[ll,ff,1] = total((rslt.ewion[ll,rng] - $
                                             rslt.ewion_excl[ll,ff,0])^2)/(nfil-2.)
        endelse

     endfor                     ; loop ff=nfil

     ;; Jackknife statistics (most comments in mean-estimator section)
     if rslt.median then begin
        ;; Estimator of expectation value is median of whole sample: <x_(.)>
        rslt.ewion_est[ll,0] = median(rslt.ewion_excl[ll,*,0],/even)

;        ;; Estimator of variance (need some kind of scaling wrt nfil?
;        ;; or should this be the same equation as the mean-estimator
;        ;; section? ... YES): <varx_jk>
;        rslt.ewion_est[ll,1] = median(abs(rslt.ewion_excl[ll,*,0] - $
;                                          rslt.ewion_est[ll,0]),/even)
        
     endif else begin
        ;; Jackknife estimator of expectation value (e.g., mean:
        ;; <x_(.)> = 1/n SUM( <x_(i)>, i=1, n ) )
        ;; Mean or median of deleted averages; should equal estimator
        ;; of whole sample.
        rslt.ewion_est[ll,0] = mean(rslt.ewion_excl[ll,*,0])

        ;; If had median-specific <varx_jk>, would need one here too
     endelse
;     print,rslt.ion[ll],median(rslt.ewion_excl[ll,*,0],/even),mean(rslt.ewion_excl[ll,*,0])

     ;; Jackknife estimator of variance
     ;; <varx_jk> = (n-d)/d SUM( (<x_(i)> - <x_(.)>)^2, i=1, n )
     ;; where d is number excluded (_Essential Statistical
     ;; Interference: Theory and Methods_ by Boos & Stefanksi (2013)).
     ;; If this is used here, maybe all other (n-1) should be (n-d)? 
;     if rslt.median then d = median(rslt.nref-rslt.median,/even) $
;     else d = mean(rslt.nref-rslt.median) 
     ;; or coefficient is (n-1)/n (and this is hard to reconcile
     ;; with (n-d)/d in limit d goes to 1)
     rslt.ewion_est[ll,1] = (nfil-1)*mean((rslt.ewion_excl[ll,*,0] - $
                                           rslt.ewion_est[ll,0])^2)
     
     ;; Bias estimator for mean and variance
     ;; Bias = (n - 1) ( <x_(.)> - x_ref )
     ;; where <x_j> is the verage of the "leave-one-out" estimates and
     ;; x_ref 
     rslt.ewion_bias[ll,0] = (nfil-1)*( rslt.ewion_est[ll,0] - $
                                        rslt.ewref[ll,0] )
     rslt.ewion_bias[ll,1] = (nfil-1)*( rslt.ewion_est[ll,1] - $
                                        rslt.ewref[ll,1] )
     
     ;; Jackknife bias-corrected estiamtes
     ;; <x> = n x_ref - (n - 1) <x_(.)>
     rslt.ewion_estcorr[ll,0] = nfil*rslt.ewref[ll,0] - $
                                (nfil-1)*rslt.ewion_est[ll,0]
     rslt.ewion_estcorr[ll,1] = nfil*rslt.ewref[ll,1] - $
                                (nfil-1)*rslt.ewion_est[ll,1]


     ;; CDF manually constructed around references (has failure modes)
     ewrng = [0.99*min(rslt.ewion[ll,*]-sqrt(rslt.sigewion[ll,0,*])),$
              1.01*max(rslt.ewion[ll,*]+sqrt(rslt.sigewion[ll,1,*]))]     
     dloc = (ewrng[1]-ewrng[0])/(nfil-1.) 
     loc = dloc*findgen(nfil) + ewrng[0] ; matches what histogram(loc=) would return
     hist = lonarr(nfil)                    ; zeros
     for ff=0,nfil-1 do begin
        sub = where(rslt.ewion[ll,*] ge loc[ff] and $
                    rslt.ewion[ll,*] lt loc[ff]+dloc,nsub)
        if nsub eq 0 then continue
        hist[ff] += total(rslt.nref-rslt.nabs[sub]) ; number contributing
     endfor                                         ; loop ff=nfil
     cdf_ion = total(hist,/cum)/float(rslt.nref)
     ;; OR 
;     srt_ewion = sort(rslt.ewion[ll,*])
;     cdf_ion = total(rslt.nref-rslt.nabs[srt_ewion],/cum)/float(rslt.nref) ; 1/nabs to 1.
;     loc = rslt.ewion[ll,srt_ewion]
     ;; if /wrt_ref, then ewcdf finds error estimate with respect to
     ;; reference, otherwise, with respect to mean/median of jackknife
     ;; sample
     fprob = interpol(loc, cdf_ion, rslt.percentile)
     rslt.ewcdf[ll,1] = rslt.ewcdf[ll,0] - fprob[0]
     rslt.ewcdf[ll,2] = fprob[1] - rslt.ewcdf[ll,0]

  endfor                        ; loop ll=nlin

  ;; Turn all variances to sigma (does not apply to CDF)
  rslt0 = rslt                  ; preserve for debugging
  rslt.ewref[*,1] = sqrt(rslt.ewref[*,1])
  if not rslt.median then $ ; don't sqrt MAD
     rslt.ewion_excl[*,*,1] = sqrt(rslt.ewion_excl[*,*,1])
  rslt.ewion_est[*,1] = sqrt(rslt.ewion_est[*,1])
  rslt.ewion_bias[*,1] = sqrt(rslt.ewion_bias[*,1]) ; is this fair?
  rslt.ewion_estcorr[*,1] = sqrt(rslt.ewion_estcorr[*,1])
  ;; rslt.ewcdf[*,1:2] already sigma
  
  return, rslt
end                             ; sdss_stackciv_jackknife_stats


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_stackciv, civstrct_fil, outfil, debug=debug, clobber=clobber, $
                   gwave=gwave, wvmnx=wvmnx, wvnrm=wvnrm, final=final,$
                   wvmsk=wvmsk, cmplt_fil=cmplt_fil, litwgt=litwgt, $
                   ivarwgt=ivarwgt, ndblt=ndblt, $ 
                   civobs_corr=civobs_corr, conti=conti, cflg=cflg, $
                   median=median, percentile=percentile, reflux=reflux, $
                   refit=refit, reerr=reerr, qerr=qerr, _extra=extra

  if N_params() LT 2 then begin 
     print,'Syntax - sdss_stackciv, civstrct_fil, outfil, [/debug, /clobber, '
     print,'                   gwave=, wvmnx=, wvnrm=, /final, wvmsk=, '
     print,'                   cmplt_fil=, /civobs_corr, /conti, cflg=, /litwgt, '
     print,'                   /ivarwgt, /median, percentile=, /reflux, /refit,'
     print,'                   /reerr, ndblt=, /qerr, _extra=]' 
     return
  endif

  ;; Timing
  tstart = systime(/seconds)

  ;; Check file and clobber option
  test = file_search(outfil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     if keyword_set(refit) or keyword_set(reflux) or keyword_set(reerr) then $
        goto, begin_fit $       ; SKIP
     else begin
        print,'sdss_stackciv: file exists; will not clobber ',outfil
        return                  ; EXIT
     endelse 
  endif 

  sdssdir = sdss_getsdssdir()
  if not keyword_set(ndblt) then ndblt = 1

  if not keyword_set(gwave) and not keyword_set(wvmnx) then $
     wvmnx = [900.,9000.d]      ; Lyman limit to basically highest SDSS lambda
  
  ;; _extra= includes /default, rating=, zlim=, dvqso=, dvgal=,
  ;; ewlim=, /noBAL, dblt_name=
  civstr = sdss_getcivstrct(civstrct_fil,_extra=extra)
  nciv = (size(civstr,/dim))[0]

  tags = tag_names(civstr)
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
     wvlimtag = (where(tags eq 'WVLIM_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
     wvlimtag = (where(tags eq 'WVLIM_ORIG'))[0]
  endelse

  if keyword_set(cmplt_fil) then begin
     ;; Sanity check
     if keyword_set(ivarwgt) or keyword_set(litwgt) then begin
        print,'sdss_stackciv: WARNING!!! does not seem sensible to completeness-weight *and* /ivarwgt or /litwgt; stopping'
        stop
     endif
     
     ;; Completeness correct each absorber... probably only works when
     ;; *not* normalizing the spectra
     ;; _extra= includes /grid, /ewmax
     if n_elements(cmplt_fil) ne ndblt then $
        stop,'sdss_stackciv stop: number of completeness files DNE ndblt'

     czw = fltarr(nciv,ndblt,/nozero)
     for dd=0,ndblt-1 do begin
        if keyword_set(civobs_corr) then begin
           if ndblt eq 1 then civcorr = civstr $
           else $  
              ;; Must rotate location of "primary" doublet
              civcorr = sdss_cpstrct(civstr,{sdsscivstrct},nshift=-2*dd)
           
           czw[*,dd] = sdss_getdxw(cmplt_fil[dd], civcorr.(ztag)[0], $
                                   civcorr.(ewtag)[0], $
                                   sigewabs=civcorr.(sigewtag)[0], /czw, $
                                   final=final, civobs_corr=civcorr, $
                                   _extra=extra)
        endif else begin
           ;; KLC: JS found this issue so temporarily use civstr
           ;; instead of civcorr until further thought, 30 May 2015
           czw[*,dd] = sdss_getdxw(cmplt_fil[dd], civstr.(ztag)[0], $
                                   civstr.(ewtag)[0], $
                                   sigewabs=civstr.(sigewtag)[0], /czw, $
                                   final=final, civobs_corr=0, $
                                   _extra=extra)
           civcorr = civstr     ; for test statement printing below
        endelse 
        

        ;; Sanity check
        test = where(czw[*,dd] eq 0. or finite(czw[*,dd]) eq 0,ntest)
        if ntest ne 0 then begin
           ;; Seriously, I don't understand how this gets
           ;; triggered but that's b/c I don't understand
           ;; sdss_getdxw()
           print,''
           print,'sdss_stackciv: WARNING!!! 0 or non-finite C(z,W):'
           printcol,civcorr[test].qso_name,civcorr[test].wrest[0],$
                    civcorr[test].(ztag)[0],civcorr[test].(ewtag)[0],$
                    czw[test,dd],format='(5x,a14,1x,f9.4,1x,f7.5,1x,f5.2,f5.3)'
           print,''
           czw[test,dd] = !values.f_nan
        endif
     endfor                     ; loop dd=ndblt

  endif 
  
  
  if keyword_set(gwave) then begin
     ;; Global, log-linear, rest wavelength array for stack
     ngpix = (size(gwave,/dim))[0] 
     wvmnx = [gwave[0],gwave[ngpix-1]]
     pixscale = (gwave[1]-gwave[0])/gwave[0]
  endif else begin
     ;; Set up gwave
     pixscale = sdss_getspecpixscale(/loglam)/alog(10) ; this is resolution of 
     ngpix = round((alog10(wvmnx[1]/wvmnx[0]))/pixscale) + 1L
     gwave = 10.^(alog10(wvmnx[0]) + dindgen(ngpix)*pixscale)
  endelse 
  
  ;; This is memory intensive; try floats for now (update changes here
  ;; to trimmed gstrct = tmp below)
  gstrct = {$
           median:keyword_set(median),$ ; whatever sdss_stackciv_errmc() needs
           percentile:fltarr(2),$       ; stored in fdat[*,[5,6]]
           ewabs:civstr.(ewtag)[2*indgen(ndblt)],$ ; [ndblt, nabs]
           zabs:civstr.(ztag)[2*indgen(ndblt)],$
           gwave:gwave,$
           gflux:fltarr(ngpix,nciv,/nozero),$
           medsnr_spec:fltarr(nciv,3,/nozero),$ ; [<f>,<sig>,<f/sig>]
           gvariance:fltarr(ngpix,nciv,/nozero),$
           gweight:fltarr(ngpix,nciv,/nozero),$
           gnspec:fltarr(ngpix,nciv,/nozero) $ ; counter and may divide
           }

  if keyword_set(wvmsk) then nwvmsk = (size(wvmsk[*,0],/dim))[0] > 1

  
  ;; Get spectra names
  spec_fil = sdss_getname(civstr,/strct,dir=specdir)
  if keyword_set(conti) then $  ; get continuum names for normalized-then-stack
     ;; _extra may include /extrap to get *abslinx.fit
     conti_fil = sdss_getname(civstr,/strct,dir=cdir,/abslin,_extra=extra)


  ;; Should sort QSOs to save time in reading in
  for ff=0L,nciv-1 do begin
     if ff eq 0 then begin
        parse_sdss,sdssdir+specdir[ff]+spec_fil[ff],flux,wave,sig=sigma,npix=npix
        if keyword_set(conti) then $ ; read in continuum fits
           cstrct = xmrdfits(sdssdir+cdir[ff]+conti_fil[ff],1,/silent)
     endif else begin
        ;; Save number of reading in by checking if current file same
        ;; as previous
        if spec_fil[ff] ne spec_fil[ff-1] then begin
           parse_sdss,sdssdir+specdir[ff]+spec_fil[ff],flux,wave,sig=sigma,npix=npix
           if keyword_set(conti) then $ 
              cstrct = xmrdfits(sdssdir+cdir[ff]+conti_fil[ff],1,/silent)
        endif
     endelse

     if keyword_set(wvmsk) then begin
        ;; Considering masking out sky lines as in McDonald et al. (2006)
        ;; and Ivashchenko et al. (2011) so long as absorber not in them
        ;; 5575--5583 Ang, 5888-5895 Ang, 6296--6308 Ang, 6862--6871 Ang.
        ;; sdss_getskylinwave() only has 5579 and 6302.
        ;; But make sure full doublet is included always.
        for ii=0,nwvmsk-1 do begin
           for dd=0,ndblt-1 do begin
              sub = where(wave ge wvmsk[ii,0] and wave le wvmsk[ii,1] and $
                          not (wave ge civstr[ff].(wvlimtag)[dd,0] and $
                               wave le civstr[ff].(wvlimtag)[dd,1]) and $ ; wvI
                          not (wave ge civstr[ff].(wvlimtag)[dd+1,0] and $
                               wave le civstr[ff].(wvlimtag)[dd+1,1])) ; wvII
              if sub[0] ne -1 then sigma[sub] = 0.                     ; exclude
           endfor                                                      ; loop dd=ndblt
        endfor                                                         ; loop ii=nwvmsk
     endif                                                             ; wvmsk=


     ;; Always shift to rest wavelength of absorber 
     rwave = wave / (1. + civstr[ff].(ztag)[0]) 

     if keyword_set(conti) then begin
        ;; Normalize
        if cstrct.npix ne npix then $
           stop,'sdss_stackciv stop: npix != cstrct.npix'
        if keyword_set(cflg) then cstrct.cflg = cflg ; could normalize by other continuum model 

        cindx = fix(alog(cstrct.cflg)/alog(2))
        gdpix = where(sigma gt 0. and $ ; avoid floating point errors
                      cstrct.conti[0:cstrct.npix-1,cindx] gt 0.,ngdpix)
        ;; Note: should never have conti be negative (but
        ;; /eig,/extrap could) but going to call this non-physical
        nwsig = sdss_calcnormerr(flux,sigma,cstrct)
        nwfx = flux * 0.        ; set array
        nwfx[gdpix] = flux[gdpix]/cstrct.conti[gdpix,cindx]

        ;; Copy over
        flux = nwfx
        sigma = nwsig
     endif else gdpix = where(sigma gt 0.,ngdpix) ; avoid floating point errors


     if keyword_set(wvnrm) then begin
        ;; Consider normalizing all spectra with arithmetic mean flux in
        ;; all pixels within emitted wavelengths 1450--1470 Ang due to
        ;; similarity of quasar spectra (Press et al. 1993, Zheng et
        ;; al. 1997), as done in Ivashchenko et al. (2014, MNRAS, 437,
        ;; 4, 3343)
        rwv_qso = wave / (1. + civstr[ff].z_qso)
        gdnrm = where(rwv_qso[gdpix] ge wvnrm[0] and $
                      rwv_qso[gdpix] le wvnrm[1] and $
                      sigma gt 0., $ ; excl. wvmsk and bad conti
                      ngdnrm)
        count = 0
        while ngdnrm lt 10 do begin
           case count of
              0: begin
                 print,'sdss_stackciv: does not span wvnrm= region: ',$
                       spec_fil[ff],civstr[ff].z_qso
                 gdnrm = where(rwv_qso[gdpix] ge 1250.,ngdnrm) ; matches sdss_fndciv
              end
              1: begin
                 print,'sdss_stackciv: does not span >= 1250 Ang: ',$
                       spec_fil[ff],civstr[ff].z_qso
                 gdnrm = where(rwv_qso[gdpix] ge 1215.,ngdnrm) 
              end
              2: begin
                 print,'sdss_stackciv: does not span >= 1215 Ang: ',$
                       spec_fil[ff],civstr[ff].z_qso
                 gdnrm = gdpix
                 ngdnrm = ngdpix > 10 ; have to get out of loop
              end
           endcase
           count++
        endwhile                ; ngdnrm < 10

        if keyword_set(median) then $
           norm = 1./median(flux[gdpix[gdnrm]],/even) $ ; robust to absorption
        else norm = 1./mean(flux[gdpix[gdnrm]],/nan)    ; paper rec.

        if keyword_set(conti) then $
           print,'sdss_stackciv: NOTICE!!! continuum-normalized and '+$
                 'normalized in specified wavelength range. Overkill?'
     endif else norm = 1. 
     flux = flux * norm
     sigma = sigma * norm

     ;; Inverse variance for weighting
     ivar = sigma * 0. 
     ivar[gdpix] = 1./sigma[gdpix]^2     

     ;; Rebin (faster than x_specrebin())
     fxnew = rebin_spectrum(flux[gdpix], rwave[gdpix], gstrct.gwave)
     varnew = rebin_spectrum(sigma[gdpix]^2, rwave[gdpix], gstrct.gwave)
     gdpixnew = where(varnew gt 0.,ngdpixnew)
     ivarnew = varnew*0.      ; make right size
     ivarnew[gdpixnew] = 1./varnew[gdpixnew]

     ;; Instantiate weights (for mean and median stacks)
     if keyword_set(ivarwgt) then weight = ivarnew $   
     else weight = replicate(1.,ngdpixnew)
     
     ;; "Light-Weighting": normalize inverse variance to median 1.
     ;; Like Weiner et al. (2009)
     if keyword_set(litwgt) then $
        medivar = median(ivar[gdpix],/even) $ ; go back to original
     else medivar = 1.
     weight = weight / medivar ; median of this won't necessarily be 1...

     if keyword_set(cmplt_fil) then begin
        ;; Median or weighted mean can get here
        ;; Multiple completeness assumes independent 
        ;; Should fold in error... 
        gd = where(finite(czw[ff,*]) eq 1,ngd)
        if ngd eq 0 then weight = 0. * weight $ ; really wrong...
        else begin
           cwgt = czw[ff,gd[0]]
           for dd=1,ngd-1 do cwgt = cwgt*czw[ff,gd[dd]]
           weight = weight / cwgt ; completeness-weighted
        endelse
     endif 

     ;; Add to global arrays
     gdvar = where(ivarnew gt 0.,complement=bdvar) ; avoid infinite values

     gstrct.gweight[gdvar,ff] = weight[gdvar]
     gstrct.gflux[gdvar,ff] = fxnew[gdvar]
     gstrct.gvariance[gdvar,ff] = 1./ivarnew[gdvar]
     gstrct.gnspec[gdvar,ff] = 1
     
     ;; Store for later; sub should never fail so long as working
     ;; with metal-line systems outside Lya forest
     sub = where(ivarnew gt 0. and $
                 gstrct.gwave*(1.+civstr[ff].(ztag)[0])/(1.+civstr[ff].z_qso) $
                 ge 1250.)            ; match sdss_fndlin
     if sub[0] eq -1 then sub = gdvar ; fall back
     sig_tmp = sqrt(gstrct.gvariance[sub,ff])
     gstrct.medsnr_spec[ff,*] = [median(fxnew[sub],/even),$
                                 median(sig_tmp,/even),$
                                 median(fxnew[sub]/sig_tmp,/even)]
     
     if bdvar[0] ne -1 then begin
        ;; Must instantiate; and exclude from median() by making
        ;; NaN 
        gstrct.gweight[bdvar,ff] = 0.
        gstrct.gflux[bdvar,ff] = !values.f_nan
        gstrct.gvariance[bdvar,ff] = !values.f_nan
        gstrct.gnspec[bdvar,ff] = 0
     endif 
     
     if keyword_set(debug) and $
        (ff mod 100) eq 1 then print,'sdss_stackciv: count =',ff

  endfor                        ; loop ff=nciv


  ;; Collapse and do faux-error assessment because want percentiles
  fdat = sdss_stackciv_stack(gstrct, median=gstrct.median, $
                             percentile=percentile)

  ;; Trim leading and trailing (fdat[*,1] is number of spectra per pixel)
  gd = where(fdat[*,1] ne 0. and finite(fdat[*,2])) ; spectra added to pixels
  if gd[0] ne -1 then begin
     gapstrt = where(gd ne shift(gd,1)+1, ngap)
     gapstop = where(gd ne shift(gd,-1)-1)
     istrt = gd[gapstrt[0]]
     istop = gd[gapstop[ngap-1]]
     ngpix = istop - istrt + 1
     fdat = fdat[istrt:istop,*]
     ;; save space by trimming (keep synched with gstrct = {} def above)
     tmp = {$
           median:gstrct.median,$ ; whatever sdss_stackciv_errmc() needs
           percentile:gstrct.percentile,$
           ewabs:gstrct.ewabs,$
           zabs:gstrct.zabs,$
           gwave:gstrct.gwave[istrt:istop],$
           gflux:gstrct.gflux[istrt:istop,*],$
           medsnr_spec:gstrct.medsnr_spec,$ ; [<f>,<sig>,<f/sig>]
           gvariance:gstrct.gvariance[istrt:istop,*],$
           gweight:gstrct.gweight[istrt:istop,*],$
           gnspec:gstrct.gnspec[istrt:istop,*] $ ; counter and may divide
           }

     gstrct = tmp 
  endif 
  
  ;; Make basic header and reproduce SDSS-like info
  ;; Information
  fxhmake, header, fdat
  sxaddpar,header,'NABS',nciv,'Number of absorbers in stack'
  sxaddpar,header,'WVION',civstr[0].wrest[0],'Ion rest wavelength'
  sxaddpar,header,'ZMED',median(civstr.(ztag)[0],/even),'Median ion redshift'
  sxaddpar,header,'ZMEAN',mean(civstr.(ztag)[0]),'Mean ion redshift'
  sxaddpar,header,'ZMIN',min(civstr.(ztag)[0],max=mx),'Min ion redshift'
  sxaddpar,header,'ZMAX',mx,'Max ion redshift'
  sxaddpar,header,'EWMED',median(civstr.(ewtag)[0],/even),'Median ion rest EW'
  sxaddpar,header,'EWMEAN',mean(civstr.(ewtag)[0]),'Mean ion rest EW'
  sxaddpar,header,'EWMIN',min(civstr.(ewtag)[0],max=mx),'Min ion rest EW'
  sxaddpar,header,'EWMAX',mx,'Max ion rest EW'

  ;; Options
  sxaddpar,header,'MEDIAN',keyword_set(median),'0: mean; 1: median flux'
  ;; Normalized?
  sxaddpar,header,'WVNRM',keyword_set(wvnrm),'Normalize spectra at wavelength region'
  if keyword_set(wvnrm) then begin
     sxaddpar,header,'WVNRM0',wvnrm[0],'Lower norm wave bound'
     sxaddpar,header,'WVNRM1',wvnrm[1],'Upper norm wave bound'
  endif
  ;; Inverse-variance weighted
  sxaddpar,header,'IVARWGT',keyword_set(ivarwgt),'0: not inverse-var weighted; 1: is'
  ;; "Light-weighted"
  sxaddpar,header,'LITWGT',keyword_set(litwgt),'Norm med inverse-var to 1'
  if size(cmplt_fil,/type) eq 7 then begin
     for dd=0,ndblt-1 do begin
        prs = strsplit(cmplt_fil[dd],'/',/extract,count=nprs)
        numstr = strtrim(dd+1,2)
        sxaddpar,header,'CMPLTFIL'+numstr,prs[nprs-1],'File '+$
                 numstr+' for completeness-corr weight' 
     endfor                     ; loop dd=ndblt
  endif else sxaddpar,header,'CMPLTFIL',keyword_set(cmplt_fil),$
                      'Completeness-corr weight'
  sxaddpar,header,'NDBLT',ndblt,'Number of doublets for wvmsk and cmpltfil'

  ;; Mask sky lines
  root = 'WVMSK'
  sxaddpar,header,root,keyword_set(wvmsk),'Mask out regions'
  if keyword_set(wvmsk) then begin
     for ii=0,nwvmsk-1 do begin
        iistr = strtrim(ii,2)
        sxaddpar,header,root+iistr+'0',wvmsk[ii,0],'Lower mask wave bound '+iistr
        sxaddpar,header,root+iistr+'1',wvmsk[ii,1],'Upper mask wave bound '+iistr
     endfor                     ; loop ii=nwvmsk
  endif                         ; wvmsk=
  
  ;; Wavelength solution
  sxaddpar,header,'COEFF0',alog10(gstrct.gwave[0]),'Center wavelength (log10) of first pixel'
  sxaddpar,header,'COEFF1',pixscale,'Log10 dispersion per pixel'
  sxaddpar,header,'CRVAL1',alog10(gstrct.gwave[0]),'Iraf zero point'
  sxaddpar,header,'CD1_1',pixscale,'Iraf dispersion'
  sxaddhist,'Zeroth dimen is flux',header,/comment
  sxaddhist,'First dimen is Nspec per pix',header,/comment
  sxaddhist,'Second dimen is sigma',header,/comment
  sxaddhist,'Third dimen is weight per pix',header,/comment     
  sxaddhist,'Fourth dimen is continuum',header,/comment
  sxaddpar,header,'PERLOW',percentile[0],'Lower percentile in fifth dimen'
  sxaddpar,header,'PERHIGH',percentile[1],'Upper percentile in sixth dimen'
  sxaddhist,'Fifth dimen is lower percentile',header,/comment
  sxaddhist,'Sixth dimen is upper percentile',header,/comment


  begin_fit:
  if keyword_set(reflux) then begin
     header = xheadfits(outfil,/silent)
     gstrct = xmrdfits(outfil,2,/silent) 
     print,'sdss_stackciv: re-collapsing the stack and finding in ',outfil
     fdat = sdss_stackciv_stack(gstrct, median=gstrct.median, $
                                percentile=gstrct.percentile) 
  endif else if keyword_set(refit) or keyword_set(reerr) then begin
     ;; Read in
     fdat = xmrdfits(outfil,0,header,/silent)
;     civstr = xmrdfits(outfil,2,/silent)
     gstrct = xmrdfits(outfil,2,/silent) 
     print,'sdss_stackciv: re-fitting and finding in ',outfil
  endif

  if keyword_set(qerr) then begin
     if keyword_set(reerr) then begin
        print,"sdss_stackciv: /reerr with /qerr doesn't do anything" 
        sxdelpar,header,'NITER' ; no longer applicable, if exists
        sxdelpar,header,'FEXCL'
     endif
  endif else begin
     print,'sdss_stackciv: Monte-Carlo bootstrap error (re-)estimate'
     if not keyword_set(niter) then niter = 1000
     if not keyword_set(fexcl) then fexcl = 0.25 ; 25%
     ;; _extra= includes seed=, oseed=, sigew=, and stuff for
     ;; sdss_stackciv_fitconti()
     fdat[*,2] = sdss_stackciv_errmc(fdat, gstrct, niter=niter, $
                                     fexcl=fexcl, cstrct_resmpl=cstrct_resmpl, $
                                     _extra=extra)
     ;; cstrct_resmpl only returned if _extra= includes sigew=
     ;; can be passed to sdss_mkstacksumm() for further analysis
     
     sxaddpar,header,'NITER',niter,'Number of MC iterations'
     sxaddpar,header,'FEXCL',fexcl,'Frac or num exclude each iter'
  endelse


  ;; Conti and line-finding always done with /reflux, /refit, /reerr
  ;; _extra= includes lin_fil=, dvlin=, and lots of other stuff (see
  ;; function) 
  cstrct = sdss_stackciv_fitconti( fdat, wave=gstrct.gwave, debug=debug, $
                                   _extra=extra)
  fdat[*,4] = cstrct.conti[0:cstrct.npix-1,fix(alog(cstrct.cflg)/alog(2))]

  ;; Write
  mwrfits,fdat,outfil,header,/create,/silent ; ext = 0
  mwrfits,cstrct,outfil,/silent              ; ext = 1
;  mwrfits,civstr,outfil,/silent              ; ext = ...
  mwrfits,gstrct,outfil,/silent ; ext = 2
;  mwrfits,cstrct_resmpl,outfil,/silent ; ext = 3  ; when writable to FITS
  spawn,'gzip -f '+outfil
  print,'sdss_stackciv: created ',outfil

  if keyword_set(cstrct_resmpl) then begin
     ;; Can't save array of structures to FITS
     tmpfil = strmid(outfil,0,strpos(outfil,'.',/reverse_search))+'_cstrct_resmpl.sav'
     save,cstrct_resmpl,filename=tmpfil
     print,'sdss_stackciv: created ',cstrct_resmpl
  endif

  ;; Final messages
  if not keyword_set(silent) then $
     print, 'sdss_stackciv: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_stackciv: Elapsed clock time (m) = ',dt/60.

  ;; Plot
  if keyword_set(debug) then begin
     x_specplot, outfil, ytwo=fdat[*,4], inflg=5, /lls, zin=1.e-6,/block
     stop
  endif 


end
