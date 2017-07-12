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
;                  median stack [default: 25th and 75th]
;   /refit -- re-normalize and find lines; calls sdss_stackciv_stack() 
;   /qerr -- quick error analysis (not correct)
;   /reerr -- re-analyze the error from Monte Carlo bootstrapping
;             function sdss_stackciv_errmc()
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
  cstrct.snr_conv[gdpix,cindx] = flux[gdpix]/sigma[gdpix] 

  ;; Mask out lines and record "centroids"
  premask = replicate(1,npix)
  dlim = dvlin / c
  for ll=0,nlin-1 do begin
     sub = where(abs(wave-linstr[ll].wave) lt dlim*linstr[ll].wave)
     if sub[0] ne -1 then premask[sub] = 0
  endfor

  ;; _extra= includes everyn=, sset=, maxrej=, lower=, upper=, nord=,
  ;; /groupbadpix, /sticky, bsplmask=, /debug, /silent
  cstrct = sdss_fndlin_fitspline(wave, flux, sigma, 0.0, $
                                 premask=premask, /nopca, $
                                 cstrct_fil=cstrct, _extra=extra)

  ;; Find absorption lines
  ;; _extra= includes lsnr=, /debug
  mask = sigma * 0.
  mask[gdpix] = 1. 
  fconti = cstrct.conti[0:cstrct.npix-1,cindx]
  ivconti = 1./fconti
  fx = flux * ivconti
  sig = sdss_calcnormerr(flux, sigma, cstrct, cflg=cstrct.cflg, $
                         baderrval=9.e9)
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
                              wvnrm=wvnrm, fonly=fonly

  if n_params() ne 1 then begin
     print,'Syntax - sdss_stackciv_stack( gstrct, [/median, percentile=, '
     print,'                              /wvnrm, /fonly] )'
     return, -1
  endif

  ;; Collapse
  ;; Going to follow SDSS (pre-SDSS-III) spSpec format (ext=0, what
  ;; parse_sdss can handle):
  ;; fdat[*,0] = flux
  ;; fdat[*,1] = number of spectra per pixe [not continuum-subtracted
  ;;             spectrum] 
  ;; fdat[*,2] = 1-sigma uncertainty; inverse-variance weighted error
  ;;              propogation (very small) if weighted stack or
  ;;              standard error on the median if median stack
  ;; fdat[*,3] = weight (or median weight) [not mask array]
  ;; fdat[*,4] = continuum [not in SDSS spSpec]
  ;; fdat[*,5] = lower percentile if median stack [not in SDSS spSpec]
  ;; fdat[*,6] = upper percentile if median stack [not in SDSS spSpec]

  ngpix = (size(gstrct.gwave,/dim))[0] > 1

  if keyword_set(median) then begin
     if not keyword_set(percentile) then begin
        ;; Consider these like "errors" but just alternate percentiles
        ;; to save (where median is 50%)
;        percentile = replicate(gauss_pdf(1.),2) ; ~84.1%; +1 Gaussian sigma
;        percentile[0] = 1. - percentile[0]      ; ~15.9%; -1 Gaussian sigma
        percentile = [0.25,0.75]
     endif
     fdat = fltarr(ngpix,7,/nozero)
     fdat[*,1] = total(gstrct.gnspec,2)

     ;; Estimate error straight from the statistics
     for pp=0L,ngpix-1 do begin
        gd = where(finite(gstrct.gflux[pp,*]),ngd) ; and finite(gstrct.gweight[pp,*]),ngd)

        if ngd eq 0 then begin
           fdat[pp,[0,2,3,5,6]] = 0. ; must instantiate flux
           ;; 4 is continuum
        endif else begin

           if ngd eq 1 then begin
              fdat[pp,0] = gstrct.gflux[pp,gd]

              if keyword_set(fonly) then continue; save time

              fdat[pp,2] = sqrt(gstrct.gvariance[pp,gd]) ; lack of better to do
              fdat[pp,3] = gstrct.gweight[pp,gd]
              ;; 4 is continuum
              fdat[pp,[5,6]] = 0 ; percentiles not calculable
           endif else begin

              fdat[pp,0] = sdss_medianw(reform(gstrct.gflux[pp,gd]),$
                                        reform(gstrct.gweight[pp,gd]),/even,$
                                        /silent)

              if keyword_set(fonly) then continue ; save time


              ;; median weight per pix... hope to be something related to
              ;; completeness or 1 
              fdat[pp,3] = median(reform(gstrct.gweight[pp,gd]),/even) 

;              ;; http://davidmlane.com/hyperstat/A106993.html
;              ;; sigma_med = 1.253 * stddev() / sqrt(N)
;              ;; And standard deviation on a weighted quantity is:
;              ;; http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
;              ;; sd_w = sqrt( ( SUM( w_i*(x_i - X_w)^2 ) / ( (N-1)/N *
;              ;;                                SUM( w_i ) ) )
;              gstrct.gflux_w = gstrct.gflux[pp,gd]*gstrct.gweight[pp,gd]
;              twgt = total(gstrct.gweight[pp,gd],/double)
;              mean_w = total(gstrct.gflux_w,/double)/twgt
;              stddev_w_numer = total(gstrct.gweight[pp,gd]*(gstrct.gflux[pp,gd] - mean_w)^2,$
;                                     /double)
;              stddev_w_denom = (ngd-1)/float(ngd) * twgt
;              stddev_w = sqrt( stddev_w_numer / stddev_w_denom )
;              fdat[pp,2] = 1.253 * stddev_w / sqrt(ngd)


              ;; The standard error on a median isn't robust to outliers
              ;; and non-normal distributions (as the website declares),
              ;; let's try something like the median-absolute
              ;; deviation which *is* robust but have to handle the
              ;; weighted business. I figure that I want to pull the
              ;; MAD statistics towards the absolute deviations based
              ;; on weight... have to think about this
              ;; MAD = < | x_i - <x_i> | >, which actually is pretty
              ;; close to the 25th and 75th percentiles
;              fdat[pp,2] = sdss_medianw(reform(abs(gstrct.gflux[pp,gd]-fdat[pp,0])), $
;                                        reform(gstrct.gweight[pp,gd]),/even,/silent)
;              ;; The MAD on the weighted median flux is very large
;              ;; (nearly the quartiles)
;              ;; Let's assess on the variance
;              var_w = sdss_medianw(reform(gstrct.gvariance[pp,gd]),$
;                                   reform(gstrct.gweight[pp,gd]),/even,/silent)
;              fdat[pp,2] = sdss_medianw(reform(abs(gstrct.gvariance[pp,gd]-$
;                                                   var_w)),$
;                                        reform(gstrct.gweight[pp,gd]),/even,/silent)
;              fdat[pp,2] = sqrt(fdat[pp,2])
              ;; That's really larger, worse than the MAD of
              ;; the weighted median flux.
              ;; Since the fundamental problem is the large dispersion
              ;; of absolute flux of un-normalized sightlines,
              ;; let's pseudo-normalize them by their median
              ;; flux (included in the stack) then take the weighted
              ;; MAD (and yes, this doesn't have units of flux
              ;; but it has properties of interest and MC/bootstrap
              ;; will supercede).
              ;; KLC changed 11 Jul 2017 to just be MAD; now work with
              ;; /conti,/extra
;              if keyword_set(wvnrm) then $ ; b/c won't be ~1
                 fdat[pp,2] = sdss_medianw($
                              reform(abs(gstrct.gflux[pp,gd]-fdat[pp,0])), $
                              reform(gstrct.gweight[pp,gd]),/even,/silent) ;$
;              else $
;                 fdat[pp,2] = sdss_medianw($
;                              reform(abs(gstrct.gflux[pp,gd]/gstrct.medsnr_spec[gd,0]-1.)), $
;                              reform(gstrct.gweight[pp,gd]),/even,/silent)

              ;; For testing various error estimates method
;              print,gstrct.gwave[pp],fdat[pp,0],fdat[pp,2],$
;                    sdss_medianw(reform(abs(gstrct.gflux[pp,gd]-fdat[pp,0])), $
;                           reform(gstrct.gweight[pp,gd]),/even,/silent),$
;                    format='(f9.5,2x,3(f9.5,1x))'

              ;; Percentiles for error
              fdat[pp,5] = sdss_medianw(reform(gstrct.gflux[pp,gd]),$
                                        reform(gstrct.gweight[pp,gd]),$
                                        /even,percent=percentile[0],$
                                        /silent)
              fdat[pp,6] = sdss_medianw(reform(gstrct.gflux[pp,gd]),$
                                        reform(gstrct.gweight[pp,gd]),$
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
     ;; If /ivarwgt wasn't set in sdss_stackciv, this reduces to
     ;; just mean.
     fdat = dblarr(ngpix,5,/nozero)
     twgt = 1. / total(gstrct.gweight,2)
     fdat[*,0] = total(gstrct.gflux*gstrct.gweight,2,/nan) * twgt
     
     if not keyword_set(fonly) then begin ; save time
        fdat[*,1] = total(gstrct.gnspec,2) ; number of spectra per pixel
        fdat[*,2] = total(gstrct.gvariance*gstrct.gweight^2,2,/nan) * twgt^2
        fdat[*,2] = sqrt(fdat[*,2])
        fdat[*,3] = total(gstrct.gweight,2) ; weight
     endif

  endelse 
                             
  ;; fdat[*,4] = conti in sdss_stackciv (main)

  return, fdat
end                             ; sdss_stackciv_stack()



function sdss_stackciv_errmc, fdat, gstrct0, fexcl=fexcl, $
                              niter=niter, seed=seed, oseed=oseed
  if n_params() ne 2 then begin
     print,'Syntax -- sdss_stackciv_errmc(fdat, gstrct0, [fexcl=, '
     print,'                              niter=, seed=, oseed=])'
     return, -1
  endif

  ;; fexcl= fraction to exclude (if < 1) else total number

  if not keyword_set(niter) then niter = 1000
  if not keyword_set(fexcl) then fexcl = 0.25 ; 25%
  if keyword_set(seed) then oseed = seed
  
  ngpix = n_elements(fdat[*,0])
  if ngpix ne n_elements(gstrct0.gwave) then $
     stop,'sdss_stackciv_esterr() stop: input have different number of pixels'

  nspec = max(fdat[*,1])         ; largest number of spectra per pixel
  if fexcl lt 1 then nexcl = round(fexcl*nspec) $
  else nexcl = fexcl

  ;; Need to save the values of the flux from the stack
  fx_resmpl = fltarr(ngpix,niter,/nozero)

  for ii=0L,niter-1 do begin

     gstrct = gstrct0           ; reset

     ;; Exclude and re-sample with replacement
     iexcl = sort(randomu(oseed,nspec))
     iresmpl = sort(randomu(oseed,nspec))

     ;; Replace
     for rr=0L,nexcl-1 do begin
        gstrct.gflux[*,iexcl[rr]] = gstrct0.gflux[*,iresmpl[rr]]
        gstrct.gvariance[*,iexcl[rr]] = gstrct0.gvariance[*,iresmpl[rr]]
        gstrct.gweight[*,iexcl[rr]] = gstrct0.gweight[*,iresmpl[rr]]
        gstrct.gnspec[*,iexcl[rr]] = gstrct0.gnspec[*,iresmpl[rr]]
        gstrct.medsnr_spec[iexcl[rr],*] = gstrct0.medsnr_spec[iresmpl[rr],*]
     endfor                     ; loop rr=nloop

     new_fdat = sdss_stackciv_stack(gstrct, /fonly, $
                                    median=gstrct.median)

     fx_resmpl[*,ii] = new_fdat[*,0] ; save just the flux
     
  endfor                        ; loop ii=niter

  ;; Analyze just the distribution of the flux
  ;; now use the MAD
  tmp = reform(fdat[*,0])
  flux = rebin(tmp,ngpix,niter) ; dupliclate to get right dimensionality
  error = median( abs(fx_resmpl - flux), dim=2, /even)
  
  oseed = oseed[0]

  gstrct0 = create_struct(gstrct0,'FLUX_RESAMPLE',fx_resmpl)
  
;  x_splot,gstrct.gwave,fdat[*,0],ytwo=fdat[*,2],psym1=10,psym2=10,ythr=error,psym3=10,/block
;
;  idx = lindgen(10)
;  x_splot,gstrct.gwave,fdat[*,0],psym1=10,ytwo=fx_resmpl[*,idx[0]],psym2=10,ythr=fx_resmpl[*,idx[1]],psym3=10,yfou=fx_resmpl[*,idx[2]],psym4=10,yfiv=fx_resmpl[*,idx[3]],psym5=10,ysix=fx_resmpl[*,idx[4]],psym6=10,ysev=fx_resmpl[*,idx[5]],psym7=10,yeig=fx_resmpl[*,idx[6]],psym8=10

;   stop
  return, error

end                             ; sdss_stackciv_errmc()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_stackciv, civstrct_fil, outfil, debug=debug, clobber=clobber, $
                   gwave=gwave, wvmnx=wvmnx, wvnrm=wvnrm, final=final,$
                   wvmsk=wvmsk, cmplt_fil=cmplt_fil, litwgt=litwgt, $
                   ivarwgt=ivarwgt,$
                   civobs_corr=civobs_corr, conti=conti, cflg=cflg, $
                   median=median, percentile=percentile, refit=refit, $
                   ndblt=ndblt, reerr=reerr, qerr=qerr, _extra=extra

  if N_params() LT 2 then begin 
     print,'Syntax - sdss_stackciv, civstrct_fil, outfil, [/debug, /clobber, '
     print,'                   gwave=, wvmnx=, wvnrm=, /final, wvmsk=, '
     print,'                   cmplt_fil=, /civobs_corr, /conti, cflg=, /litwgt, '
     print,'                   /ivarwgt, /median, percentile=, /refit, /reerr, ndblt=, /qerr, _extra=]' 
     return
  endif 

  ;; Check file and clobber option
  test = file_search(outfil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     if keyword_set(refit) then goto, begin_fit $ ; SKIP
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
                    civcorr[test].zabs_orig[0],civcorr[test].ew_orig[0],$
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
  
  ;; This is memory intensive; try floats for now
  gstrct = {$
           median:keyword_set(median),$ ; whatever sdss_stackciv_errmc() needs
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
                 gstrct.gwave*(1.+civstr[ff].zabs_orig[0])/(1.+civstr[ff].z_qso) $
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
                             percentile=percentile, $ 
                             wvnrm=wvnrm) ; fudging

  ;; Trim leading and trailing (fdat[*,1] is number of spectra per pixel)
  gd = where(fdat[*,1] ne 0. and finite(fdat[*,2])) ; spectra added to pixels
  if gd[0] ne -1 then begin
     gapstrt = where(gd ne shift(gd,1)+1, ngap)
     gapstop = where(gd ne shift(gd,-1)-1)
     istrt = gd[gapstrt[0]]
     istop = gd[gapstop[ngap-1]]
     ngpix = istop - istrt + 1
     fdat = fdat[istrt:istop,*]
     ;; save space by triming
     tmp = {$
           median:gstrct.median,$ ; whatever sdss_stackciv_errmc() needs
           gwave:gstrct.gwave[istrt:istop,*],$
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
  if keyword_set(median) then begin
     sxaddpar,header,'PERLOW',percentile[0],'Lower percentile in fifth dimen'
     sxaddpar,header,'PERHIGH',percentile[1],'Upper percentile in sixth dimen'
     sxaddhist,'Fifth dimen is lower percentile',header,/comment
     sxaddhist,'Sixth dimen is upper percentile',header,/comment
  endif 


  begin_fit: 
  if keyword_set(refit) then begin
     ;; Read in
     fdat = xmrdfits(outfil,0,header,/silent)
;     civstr = xmrdfits(outfil,2,/silent)
     gstrct = xmrdfits(outfil,2,/silent) 
     print,'sdss_stackciv: re-fitting and finding in ',outfil
  endif

  if not keyword_set(qerr) or keyword_set(reerr) then begin
     if keyword_set(reerr) then $
        print,'sdss_stackciv: Monte Carlo error re-estimate'
     if not keyword_set(niter) then niter = 1000
     if not keyword_set(fexcl) then fexcl = 0.25 ; 25%
     ;; _extra= includes seed=, oseed=
     error = sdss_stackciv_errmc(fdat, gstrct, niter=niter, $
                                 fexcl=fexcl, _extra=extra)
     fdat[*,2] = error

     sxaddpar,header,'NITER',niter,'Number of MC iterations'
     sxaddpar,header,'FEXCL',fexcl,'Frac or num exclude each iter'
  endif                         ; /mcerr


  ;; Conti
  ;; _extra= includes lin_fil=, dvlin=, and lots of other stuff
  cstrct = sdss_stackciv_fitconti( fdat, wave=gstrct.gwave, debug=debug, $
                                   _extra=extra)
  fdat[*,4] = cstrct.conti[0:cstrct.npix-1,fix(alog(cstrct.cflg)/alog(2))]

  ;; Write
  mwrfits,fdat,outfil,header,/create,/silent ; ext = 0
  mwrfits,cstrct,outfil,/silent              ; ext = 1
;  mwrfits,civstr,outfil,/silent              ; ext = 2
  mwrfits,gstrct,outfil,/silent              ; ext = 3
  spawn,'gzip -f '+outfil
  print,'sdss_stackciv: created ',outfil

  ;; plot
  if keyword_set(debug) then begin
     x_specplot, outfil, ytwo=fdat[*,4], inflg=5, /lls, zin=1.e-6,/block
     stop
  endif 


end
