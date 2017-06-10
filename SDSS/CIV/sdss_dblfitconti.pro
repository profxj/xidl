;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_dblfitconti.pro               
; Author: Kathy Cooksey                      Date: 25 Jul 2011
; Project: 
; Description: 
; Input: 
; Optional Input:
; Output: 
; Optional Output:
; Example:
; History:
;   25 Jul 2011  Created by KLC
;      Aug 2011  Jointly updated by RAS, KLC
;   27 Sep 2011  Major revamp to have more functions, KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@sdss_fndlin                    ; resolve sdss_fndlin_fitspline()

function sdss_dblfitconti_premask, flux, error, eigconti,nbad_orig=nbad_orig
  if n_params() ne 3 then begin
     print,'Syntax - sdss_dblfitconti_premask(flux, error, eigconti, [nbad_orig=])'
     return,-1
  endif

  npix = (size(flux,/dim))[0]

  residual = flux-eigconti[*,0]
  premask = replicate(1,npix)
  bd = where(eigconti[*,1] eq 1 and residual lt 0.,nbad_orig) ; 1 b/c eigconti mask flipped

  if nbad_orig ne 0 then begin
     ;; Find span of masked out regions
     gapstrt0 = bd[where(bd ne shift(bd,1)+1)]
     gapstop0 = bd[where(bd ne shift(bd,-1)-1)]
     gapsz0 = gapstop0 - gapstrt0 + 1
     ;; If >= 2 pixels, go for more
     bd = where(gapsz0 ge 2,ngap)
     for gg=0,ngap-1 do begin
        gapstrt = gapstrt0[bd[gg]]
        gapstop = gapstop0[bd[gg]]
        gapsz = gapsz0[bd[gg]]
        
        ;; Grow with buffer of 2 pixels (see shift statement)
        ;; Grow to left until flux no longer growing (towards continuum)
        if gapstrt ne 0 then begin
           gd = where(flux[0:gapstrt] gt shift(flux[0:gapstrt-1],2) $
                      and error[0:gapstrt-1] ne 0.,ngd)
           if ngd gt 0 then gapstrt = gd[ngd-1]
        endif 
        ;; Grow to the right
        if gapstop ne npix-1 then begin
           gd = where(flux[gapstop:npix-1] gt $
                      shift(flux[gapstop:npix-1],-2) $
                      and error[gapstop:npix-1] ne 0.,ngd)
           if ngd gt 0 then gapstop = gapstop + gd[0]
        endif 
        gd = where(error[gapstrt:gapstop] ne 0.,ngd)
        if ngd ge 2 then $
           medresid = [median(residual[gd],/even),$
                       median(error[gd],/even)] $
        else medresid = [0.,-9999.]
        
        ;; Check size again and if these changed, mask them out
        ;; Also keep the gaps that are 2 <= N <= 10 pixels with median
        ;; flux decrement
        ;; Just don't let it be too crazy
        gsz = gapstop-gapstrt+1
        if ((gsz gt gapsz and gsz le 20) or $
            (gapsz le 10)) and medresid[0] lt -2.*medresid[1] then begin
           ;; Consider preventing a totally crazy region from being all
           ;; masked out; bias high
           premask[gapstrt + gd] = 0 ; bad
        endif 
        
     endfor                     ; loop gg=ngap

     ;; Final check to make sure no region too large (likely the
     ;; side of a bad emission line)
     bd = where(premask eq 0)
     if bd[0] ne 0 then begin
        ;; Find span of masked out regions
        gapstrt0 = bd[where(bd ne shift(bd,1)+1)]
        gapstop0 = bd[where(bd ne shift(bd,-1)-1)]
        bd = where(gapstop0 - gapstrt0 + 1 gt 30,ngap)
        for gg=0,ngap-1 do $
           premask[gapstrt0[bd[gg]]:gapstop0[bd[gg]]] = 1 ; good
     endif 
  endif                         ; ngap                     

  return, premask
end                             ; sdss_dblfitconti_premask()



function sdss_dblfitconti_fiteigspl, wave, flux, sigma, snr_spec, dblt_name, $
                                     eig_fil, cstrct_fil, sigrej, $
                                     splconti=splconti, contihdr=contihdr, $
                                     premask=premask, bsplmask=bsplmask, $
                                     qualflg=qualflg, chi_sqr=chi_sqr, $
                                     eigbasis=eigbasis, debug=debug, $
                                     _extra=extra
  ;; 
  if n_params() ne 8 then begin
     print,'Syntax - sdss_dblfitconti_fiteigspl(wave, flux, sigma, snr_spec, '
     print,'                                    dblt_name, eig_fil, cstrct_fil, sigrej, '
     print,'                                    [splconti=, contihdr=, '
     print,'                                    premask=, bsplmask=, qualflg=, /debug, '
     print,'                                    chi_sqr=, eigbasis=, _extra=]'
     return, -1
  endif 
  if keyword_set(debug) then $
     print,'sdss_dblfitconti_fiteigspl() debug: S/N = '+$
           string(snr_spec,format='(f6.2)')+'; sigrej = '+$
           string(sigrej,format='(f5.2)')

  ;; Hopefully don't have to read this; checks may slow us down
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)   
  if not keyword_set(eigbasis) then $
     eigbasis = xmrdfits(getenv("SDSSPATH")+$
                         "/eigenspectra/eigSpec_qso_all.fit",0, /silent)
  if size(eig_fil,/type) eq 7 then $
     eigconti = xmrdfits(sdssdir+eig_fil,0,contihdr,/silent)  $
  else begin
     eigconti = eig_fil
     if not keyword_set(contihdr) then $
        stop,'sdss_dblfitconti_fiteigspl(): contihdr must be set'
  endelse 
  if size(cstrct_fil,/type) eq 8 then cstrct = cstrct_fil $
  else cstrct = xmrdfits(sdssdir+cstrct_fil,1, /silent)

  ;; Check
  npix = (size(wave,/dim))[0]   ; use this value latter
  if cstrct.npix ne npix then $
     stop,'sdss_dblfitconti_fiteigspl(): npix != cstrct.npix'


  ;; Automatically-detected centroid masking
  absmask = replicate(1,cstrct.npix) ; all good
  nskip = 0
  for ilin=0, cstrct.ncent[sdss_getcflg(/eig,/index)]-1 do begin
     lam_cent = cstrct.centroid[ilin,sdss_getcflg(/eig,/index)]
     mn = min(lam_cent-wave,imn,/abs)
     if abs(flux[imn]-eigconti[imn,0]) ge sigrej*sigma[imn] then begin
        dv = abs((wave - lam_cent)/lam_cent * 299792.)
        absmask[where(dv LT 600.)] = 0
     endif else nskip++
  endfor                        ; loop ilin=cstrct.ncent[]
  if keyword_set(debug) then $
     print,'sdss_dblfitconti_fiteigspl() debug: Number of '+$
           strtrim(cstrct.ncent[sdss_getcflg(/eig,/index)],2)+$
           ' centroids skipped: ',nskip

  premask = absmask * sdss_dblfitconti_premask(flux, sigma, eigconti)
  if cstrct.ipix0 gt 0 then begin
     premask[0:cstrct.ipix0-1] = 0    ; bad; Lya forest
     eigconti[0:cstrct.ipix0-1,*] = 0 ; not NaN
     eigconti[0:cstrct.ipix0-1,1] = 1 ; bad; Lya forest
  endif 
  
  ;; Refit 
  neweig = eigqsoconti(wave[cstrct.ipix0:*]/(1d + cstrct.z_qso),$
                       flux[cstrct.ipix0:*],$
                       sigma[cstrct.ipix0:*],eigbasis,$
                       fitmask=(1-premask[cstrct.ipix0:*]),$
                       /silent,finalmask=finalmask,$
                       header=contihdr) ; will be modified
  eigconti[cstrct.ipix0:*,0] = neweig[*,0]
  eigconti[cstrct.ipix0:*,1] = finalmask ; good place to store it
  eigconti[cstrct.ipix0:*,2] = neweig[*,1]
  
  ;; Do a better job with masking using eigconti mask and finding
  ;; likely absorption regions
  premask = absmask * sdss_dblfitconti_premask(flux, sigma, eigconti)
  if cstrct.ipix0 gt 0 then $
     premask[0:cstrct.ipix0-1] = 0 ; bad; Lya forest

  ;; Normalize and remove masked out pixels (eigcont[*,1] = 1)
  ;; Calculate new flux sigma based on continuum fit
  fxnrm = flux/eigconti[*,0]
  ;; ernrm = sqrt(sigma^2*eigconti[*,0]^2 + eigconti[*,2]^2*flux^2)*ivconti^2
  ernrm = sdss_calcnormerr(flux,sigma,eigconti)
  bdpix = where(sigma eq 0.)
  if bdpix[0] ne -1 then begin
     fxnrm[bdpix] = 0.
     ernrm[bdpix] = 0.
  endif 

  ;; Fit spline to normalizd spectrum
  ;; _extra includes: pixlim=, flg_gap=, fx=, sig=, title=, /cchk,
  ;; qualflg=, everyn=[25], maxrej=[20], lower=[2.0], upper=[3.0],
  ;; nord=[3]  
  cstrct.npix = cstrct.npix - cstrct.ipix0
  chi_sqr = 1                   ; force return
  ;; Don't enter the redshift as the 4th parameter because don't
  ;; want to do sdss_fndlin_fitspline() funny-business with
  ;; the it's own eigenconti fit if it's a high-z QSO
  splstrct = sdss_fndlin_fitspline(wave[cstrct.ipix0:*],$
                                   fxnrm[cstrct.ipix0:*],$
                                   ernrm[cstrct.ipix0:*], 0.,$
                                   /groupbadpix, /sticky, $ ; typical
                                   ;; Pass back out
                                   sset=sset, qualflg=qualflg, $
                                   chi_sqr=chi_sqr, bsplmask=bsplmask, $
                                   ;; Pass in (some saves time)
                                   premask=premask[cstrct.ipix0:*],$
                                   pca_fil=pca, pca_head=pca_head, $ 
                                   cstrct_fil=cstrct,debug=debug, $
                                   silent=silent, _extra=extra) 
  cstrct.npix = npix            ; reset; npix saved from size(wave) above

  ;; Set up correct structure
  splconti = dblarr(cstrct.npix,3,/nozero) ; MUST set values
  splconti[cstrct.ipix0:*,0] = splstrct.conti[0:splstrct.npix-1,$
                                              sdss_getcflg(/spl,/index)]
  splconti[cstrct.ipix0:*,2] = splstrct.sigconti[0:splstrct.npix-1,$
                                                 sdss_getcflg(/spl,/index)]
  if cstrct.ipix0 gt 0 then begin
     splconti[0:cstrct.ipix0-1,*] = 0
     bsplmask = [replicate(0,cstrct.ipix0),bsplmask]
  endif 

  
  ;; Write header information (instead of passing it all around)
  sxaddpar,contihdr,'DBLTNAME',dblt.ion,'doublet setting S/N region'
  sxaddpar,contihdr,'MEDSNR',snr_spec,'median S/N in region of interest'
  sxaddpar,contihdr,'USESIGREJ',sigrej,'sigma rejection for abslin masking'  
  sxaddpar,contihdr,'SPLQUAL',qualflg,'Spline quality flag'
  sxaddpar,contihdr,'NSPLORD',sset.nord,'Number of orders in bspline fit'
  sxaddpar,contihdr,'NBKPTS',(size(sset.fullbkpt,/dim))[0],'Number of breakpoints in bspline'
  sxaddpar,contihdr,'CHI2SPL', chi_sqr[0], 'chisqr of splconti fit to norm flux'
  sxaddpar,contihdr,'CHI2SDOF', chi_sqr[1], 'spl-chi2 degrees of freedom'
  sxaddpar,contihdr,'PCHI2SPL', chi_sqr[2], 'spl-chi2 probability'
  
  ;; Also passing back out splstrct, absmask=, qualflg=, bsplmask=,
  return, eigconti
end                             ; sdss_dblfitconti_fiteigen()



function sdss_dblfitconti_fithybrid, wave, flux, sigma, $
                                     snr_spec, dblt_name, eig_fil, $
                                     cstrct_fil, eigconti=eigconti, $ 
                                     sigrej=sigrej, snr_cut=snr_cut, $
                                     silent=silent, debug=debug, plot=plot, $
                                     contihdr=contihdr, $
                                     extrap=extrap, _extra=extra
  ;; Does the actually fitting (easier to call by other routines,
  ;; e.g. sdss_completeness), modeled after sdss_fndlin_fitspline()
  ;; _extra includes what's passed to sdss_fndlin_fitspline()
  if n_params() ne 7 then begin
     print,'Syntax - sdss_dblfitconti_fithybrid, wave, flux, sigma, '
     print,'                                     snr_spec, dblt_name, eig_fil,'
     print,'                                     cstrct_fil, [eigconti=, '
     print,'                                     sigrej=, snr_cut=, /plot, '
     print,'                                     contihdr=, /silent, /debug, '
     print,'                                     _extra]'
     return,-1
  endif
  sdssdir = sdss_getsdssdir()
  wvlya_cut = 1230.             ; match sdss_fndlin_fitspline(), runeigqsoconti
  npix = (size(wave,/dim))[0]

  if not keyword_set(snr_cut) then snr_cut = 8
  if keyword_set(sigrej) then usesigrej = sigrej $
  else begin
     ;; S/N cut
     if snr_spec gt 2*snr_cut then usesigrej = 5. $
     else usesigrej = 1.
  endelse 


  ;; Read these items ONCE
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)   
  if size(eig_fil,/type) eq 7 then $
     eigconti = xmrdfits(sdssdir+eig_fil,0,contihdr,/silent)  $
  else eigconti = eig_fil
  if size(cstrct_fil,/type) eq 8 then cstrct = cstrct_fil $
  else cstrct = xmrdfits(sdssdir+cstrct_fil,1, /silent)

  ;; exclude Lya forest
  mn = min( abs(wave - wvlya_cut*(1d + cstrct.z_qso)), ipix0 ) 
  if cstrct.ipix0 ne ipix0 then $
     stop,'sdss_dblfitconti_fithybrid(): ipix0 != cstrct.ipix0'

  ;; Set new eigencontinuum which depends on fititng a new eigen
  ;; _extra includes:
  ;;   for sdss_dblfitconti_fiteigspl(): eigbasis= 
  ;;   for eigqsoconti(): gapmax=, weight=, topmarg=, lowmarg=,
  ;;   maxiter=, niter=, growmarg1=, chunksz=, ctol=, fail=, status=,
  ;;   chi_sqr=, stat_ctol=, _extra [...]
  ;;   for sdss_fndlin_fitspline(): 
  ;;   pixlim=, flg_gap=, fx=, sig=, title=, /cchk,
  ;;   qualflg=, everyn=[25], maxrej=[20], lower=[2.0], upper=[3.0],
  ;;   nord=[3]  

  eigconti = sdss_dblfitconti_fiteigspl(wave, flux, sigma, snr_spec, $
                                        dblt, eigconti, cstrct, usesigrej, $
                                        ;; Things passed out
                                        splconti=splconti, premask=premask, $
                                        qualflg=qualflg, chi_sqr=chi_sqr, $
                                        bsplmask=bsplmask, $
                                        ;; Things passed in (way in)
                                        contihdr=contihdr, $ ; will change
                                        debug=debug, _extra=extra)
  

  if ((qualflg and 2) eq 2 and snr_spec ge snr_cut) then begin
     ;; Try something else b/c this is high S/N
     ;; Mask out previously found absorption features (via
     ;; sdss_fndlin) for the old-eigen-normalized spectrum
     if not keyword_set(silent) then $
        print,'sdss_dblfitconti_fithybrid(): spline fail; '+$
              'force high-sigma absmask...'
     if usesigrej ge 5 then usesigrej = usesigrej + 1 $
     else usesigrej = 5
     
     eigconti = sdss_dblfitconti_fiteigspl(wave, flux, sigma, snr_spec, $
                                           dblt, eigconti, cstrct, usesigrej, $
                                           ;; Things passed out
                                           splconti=splconti, premask=premask,$
                                           qualflg=qualflg, chi_sqr=chi_sqr, $
                                           bsplmask=bsplmask, $
                                           ;; Things passed in (way
                                           ;; in)
                                           contihdr=contihdr, $ ; will change
                                           debug=debug, _extra=extra)

     if (qualflg and 2) eq 2 then begin
        if not keyword_set(silent) then $
           print,'... still failed'
        fail_list = [fail_list,hyb_fil[ss]] 
     endif 
  endif                         ; try better absmasking


  hybconti = dblarr(npix,3,/nozero) ; MUST set values
  if cstrct.ipix0 gt 0 then $
     hybconti[0:cstrct.ipix0-1,*] = 0


  ;; Time to actually set the hybrid continuum
  if ((qualflg and 2) eq 2) or (snr_spec lt snr_cut) then begin
     if (qualflg and 2) eq 2 then begin
        if not keyword_set(silent) then $
           print,'sdss_dblfitconti_fithybrid(): spline fail; '+$
              'hybrid = eigconti: ',qso_name[ss] 
        sxaddpar,contihdr,'CONTITYP','EIGENX','spline fail; new eigen'
     endif else begin
        if not keyword_set(silent) then $
           print,'sdss_dblfitconti_fithybrid(): low-S/N spectrum; '+$
                 'hybrid = eigconti: ',cstrct.qso_name
        sxaddpar,contihdr,'CONTITYP','EIGEN','low S/N; new eigen'
     endelse 
     hybconti = eigconti 
  endif else begin
     ;; Make combined spectrum and estimate errors, save the
     ;; new-eigen+spline error in hybrid[*,2]
     hybconti[*,0] = eigconti[*,0] * splconti[*,0]
     hybconti[*,1] = bsplmask eq 0 ; invert from 1 = good; 0 = bad
     hybconti[*,2] = sqrt( eigconti[*,2]^2 + splconti[*,2]^2 )
     sxaddpar,contihdr,'CONTITYP','HYBRID','high S/N; new eigen+spline'
  endelse

  if keyword_set(extrap) then begin
     ;; Use whatever is the extrapolated spectrum
     if size(extrap,/type) ne 8 then $
        stop,'sdss_dblfitconti_fithybrid() stop: extrap= should be structure of {wave, flux, [error, source]}'

     ;; Rename pertinent template parts
     wvr_qso = extrap.wave
     wvobs_qso = wvr_qso*(1+cstrct.z_qso)
     gd = where(wvobs_qso ge wave[cstrct.ipix0] and $
                wvobs_qso le wave[npix-1],ngd)
     if ngd eq 0 then $
        stop,'sdss_dblfitconit_fithybrid() stop: template does not span QSO'
     npix_qso = (size(wvr_qso,/dim))[0]
     fx_qso = extrap.flux
     tags = tag_names(extrap)
     test = where(stregex(tags,'error',/boolean,/fold_case),ntest)
     if ntest eq 1 then $
        er_qso = extrap.error $
     else er_qso = 0            ; keyword_set(er_qso) == 0
     med_qso = median(fx_qso[gd],/even)

     ;; Interpolate template onto SDSS scale (assumes template
     ;; higher-res than SDSS); default is linear; _extra= includes
     ;; /lsquadratic, /nan, /quadratic, /spline
     fx_qso_interp = interpol(fx_qso,wvobs_qso,wave,_extra=extra)
     ;; Check:
     ;; x_splot,wvobs_qso,extrap.flux,xtwo=wave[cstrct.ipix0:*],ytwo=fx_qso_interp

     ;; Establish normalization range to test
     ncorr = 1000               ; 1000 tests
     rslt_corr = fltarr(ncorr,2,/nozero) ; [scale, RMS]
     mn_conti = min(hybconti[cstrct.ipix0:*,0],max=mx_conti) 
     rslt_corr[*,0] = (mn_conti + $
                       (mx_conti-mn_conti)/(ncorr-1.)*findgen(ncorr))/$
                      med_qso
;     med_conti = median(hybconti[cstrct.ipix0:*,0],/even)
     

     ;; Determine "optimal" normalization (highest correlation coeff)
     for cc=0,ncorr-1 do begin
        ;; Chi^2 
        rslt_corr[cc,1] = sqrt(mean((fx_qso_interp[cstrct.ipix0:*]*$
                                     rslt_corr[cc,0] - $
                                     flux[cstrct.ipix0:*])^2))
        ;; Should do a check that I've actually passed through chi^2 min
     endfor                     ; loop cc=ncorr
     mn = min(rslt_corr[*,1],imn) ; add to header
     if imn eq 0 or imn eq ncorr-1 then $
        stop,'sdss_dblfitconti_fithybrid() stop: RMS min outside of range'
        ;; Check:
        ;; x_splot,rslt_corr[*,0],rslt_corr[*,1],psym1=4
     fx_qso_interp = fx_qso_interp * rslt_corr[imn,0]

     ;; Now for stitching the two; first find where closest in 20 Ang
     ;; window
     dpix = 20. ; Ang
     mn = min(wave[cstrct.ipix0]+dpix-wave[cstrct.ipix0:*],ilim,/abs)
     ilim = ilim + cstrct.ipix0
     mn = min((fx_qso_interp-hybconti)[cstrct.ipix0:ilim],iclosest,/abs)
     hybconti[0:cstrct.ipix0-1] = fx_qso_interp[0:cstrct.ipix0-1]
     if iclosest ge 1 then begin ; 0 or 1 pixel is a line
        ;; Connect linearly; fix the observed flux point as redward
        ;; (iclosest) and the QSO template as appropriate bluward of
        ;; cstrct.ipix0
        iclosest = iclosest + cstrct.ipix0
        slope = (hybconti[iclosest]-fx_qso_interp[cstrct.ipix0])/$
                (wave[iclosest]-wave[cstrct.ipix0])
        intercept = fx_qso_interp[cstrct.ipix0] - slope*wave[cstrct.ipix0]
        hybconti[cstrct.ipix0:iclosest] = intercept + $
                                          slope*wave[cstrct.ipix0:iclosest]
     endif

     x_splot,wave,flux,psym1=10,ytwo=hybconti,psym2=-3,$
             ythr=fx_qso_interp,psym3=-3,$
             lgnd=['Spectrum','Hybconti','Tmplt. Scaled'],/block

     test = where(stregex(tags,'source',/boolean,/fold_case),ntest)
     if ntest eq 1 then $
        source = extrap.source $
     else source = 0            ; keyword_set(source) == 0

     stop
  endif                         ; /extrap

  ;; Shift values around to store new ones
  sxaddpar,contihdr,'SNRCUT',snr_cut,'S/N threshold for hybrid conti'
  sxaddpar,contihdr,'CHI2EIG', sxpar(contihdr,'CHISQR'), 'chisqr of eigconti to flux'
  eigdof = sxpar(contihdr,'CHISQDOF')
  sxaddpar,contihdr,'CHI2EDOF', eigdof, 'eig-chi2 degrees of freedom'
  sxaddpar,contihdr,'PCHI2EIG', sxpar(contihdr,'PCHISQR'), 'eig-chi2 probability'
  ;; Calc new and replace old eig-values
  gd = where(hybconti[*,1] eq 0,ngd)
  chisqr = total((flux[gd]-hybconti[gd,0])^2/sigma[gd]^2) 
  dof = eigdof + chi_sqr[1]
  sxaddpar,contihdr,'CHISQR', chisqr, 'chisqr of continuum fit to flux'
  sxaddpar,contihdr,'CHISQDOF', dof, 'chisqr degrees of freedom'
  sxaddpar,contihdr,'PCHISQR', chisqr_pdf(chisqr,dof), 'chisqr probability'
  

  ;; Want S/N exactly as used otherwise
  if keyword_set(plot) or keyword_set(debug) then begin
     clr = getcolor(/load)
     ttl = strtrim(cstrct.qso_name,2)+': zqso = '+ $
           string(cstrct.z_qso,format='(f7.5)') + $
           ': S/N = '+strtrim(snr_spec,2)

     ;; Plot
     mask = replicate(!VALUES.F_NAN,npix)
     mask2 = mask
     bd = where(premask eq 0)
     if bd[0] ne -1 then mask2[bd] = flux[bd]
     bd = where(hybconti[*,1] eq 1)
     if bd[0] ne -1 then mask[bd] = flux[bd]
     
     ;; All plots
     sdss_chkconti,wave,flux,psym1=10, title=ttl, $
                   ytwo=sigma,psym2=10, $
                   color2=clr.red,$
                   ythr=hybconti[*,0],psym3=-3, color3=clr.limegreen, $
                   /block, ymnx=[min([sigma,flux],/nan,max=mx),mx], $
                   yfou=mask,psym4=4,color4=clr.limegreen,$
                   yfiv=mask2,psym5=1,color5=clr.purple,$
                   lgnd=['flux','sigma',$
                         'hybconti','newmask','premask'], $
                   _extra=extra


     ;; new_err = sqrt(sigma^2*hybconti[*,0]^2 + hybconti[*,2]^2*flux^2)/hybconti[*,0]^2
     new_err = sdss_calcnormerr(flux,sigma,hybconti)

     sdss_chkconti,wave,flux/hybconti[*,0],psym1=10, title=ttl, $ 
                   ytwo=new_err,psym2=10, $
                   color2=clr.red,$
                   /block, ymnx=[min([new_err,flux/hybconti[*,0]],/nan,max=mx),mx], $
                   ythr=mask/hybconti[*,0],psym3=4,color3=clr.limegreen,$
                   yfou=mask2/hybconti[*,0],psym4=1,color4=clr.purple,$
                   lgnd=['norm flux','total error','newmask','premask'], $
                   _extra=extra

     if keyword_set(debug) then $
        stop,'sdss_dblfitconti_fithybrid() debug: stopping before next iteration'

  endif                         ; /plot or /debug

  return, hybconti
end                             ; sdss_dblfitconti_fithybrid()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; _extra includes e.g. xsize=, ysize= for sdss_chkconti but also
;; e.g. dblt_name= for sdss_measuresnr()
pro sdss_dblfitconti, sdss_list, sdsssum, pca_fil=pca_fil, $
                      clobber=clobber, istrt=istrt, silent=silent, $
                      processor=processor, snr_cut=snr_cut, $
                      snrstrct_fil=snrstrct_fil,dblt_name=dblt_name,$
                      extrap=extrap,_extra=extra
  if n_params() ne 2 then begin
     print,'Syntax - sdss_dblfitconti, sdss_list, sdsssum, [pca_fil=, /plot, /clobber,'
     print,'                           /debug, /silent, _extra=]'
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Read in and do some checks
  readcol,sdss_list,spec_fil,format='a',/silent
  outdir = spec_fil[0]
  spec_fil = spec_fil[1:*]
  nfil = (size(spec_fil,/dim))[0]
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  sdsstab = xmrdfits(sdsssum, 1, /silent)
  if nfil NE (size(sdsstab,/dim))[0] then $
     stop,'sdss_dblfitconti: list and QSO structure must be same size'
  if keyword_set(processor) then begin
     ;; Takes precedence over istrt=
     sub = sdss_calcparalleljob(sdsstab, processor)
     istrt = sub[0]
     nfil = sub[1] + 1

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif else if not keyword_set(istrt) then istrt = 0L ; default

  ;; Create appropriate names
  eig_fil = sdss_getname(spec_fil,/spec,/eig,dir=cdir,plate=plate,$
                         root=qso_name)
  eig_fil = cdir + eig_fil
  hyb_fil = sdss_getname(spec_fil,/spec,/hyb,extrap=extrap) ; output name
  if outdir ne 'abslin/' and outdir ne 'conti/' then $
     outdir = outdir + '1d_26/' + plate + '/1d/' $
  else outdir = cdir
  hyb_fil = outdir + hyb_fil
  abslin_fil = sdss_getname(spec_fil,/spec,/abslin, dir=adir)      
  abslin_fil = adir + abslin_fil


  ;; Measure S/N and return the bounds used; keep consistent with
  ;; other codes
  ;; _extra includes: dvgal=, dvqso=
  if keyword_set(snrstrct_fil) then begin
     ;; Use structure to reconstruct
     if size(snrstrct_fil,/type) eq 7 then $
        snr_strct = xmrdfits(snrstrct_fil,1,/silent) $
     else snr_strct = snrstrct_fil

     subsnr_strct = 1
     sdss_getqsoinlist,sdss_list,subsnr_strct,snrstrct=snr_strct,/snr,/silent
     tags = tag_names(subsnr_strct[0])
     snrtag = where(tags eq 'SNR_'+strupcase(strtrim(dblt.ion,2)))
     wvtag = where(tags eq 'WVOBS_'+strupcase(strtrim(dblt.ion,2)))

     snr_arr = transpose(subsnr_strct.(snrtag))
     wvobs_lim = transpose(subsnr_strct.(wvtag))
  endif else $ 
     ;; Compute and this takes time
     snr_arr = sdss_measuresnr(spec_fil,wvlim_obs=wvobs_lim,sdsssum=sdsstab,$
                               dblt_name=dblt_name,_extra=extra)

  ;; Defaults
  eigbasis = xmrdfits(getenv("SDSSPATH")+"/eigenspectra/eigSpec_qso_all.fit",$
                      0, /silent)
  if not keyword_set(pca_fil) then $
     pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
  if size(pca_fil,/type) eq 7 then $
     pca = xmrdfits(pca_fil, 0, pca_head, /silent) $
  else begin  
     pca = pca_fil 
     if not keyword_set(pca_head) then $
        stop,'sdss_dblfitconti: must set pca_head' 
  endelse 

  fail_list = ''                ; print out at end

  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  for ss=istrt,nfil-1 do begin
     ;; Fast check
     test = file_search(sdssdir+hyb_fil[ss]+'*',count=ntest)
     if ntest ne 0 and not keyword_set(clobber) then begin
        print,'sdss_dblfitconti: skipping ',hyb_fil[ss]
        continue
     endif

     ;; Read data
     parse_sdss, sdssdir+spec_fil[ss], flux, wave, sig=sigma ;, head=head


     ;; _extra includes sigrej=
     hybconti = $
        sdss_dblfitconti_fithybrid(wave, flux, sigma, snr_arr[ss,2], $
                                   dblt, eig_fil[ss], abslin_fil[ss],  $
                                   ;; Options
                                   silent=silent, $
                                   ;; Things passed out
                                   eigconti=eigconti, $ ; new one
                                   ;; Things passed in (and deep)
                                   eigbasis=eigbasis, pca_fil=pca, $
                                   pca_head=pca_head, $
                                   contihdr=contihdr, $ ; to be modified
                                   xmnx=wvobs_lim[ss,*], $
                                   extrap=extrap, _extra=extra)


     ;; Write out and already ran checks earlier for file existing
     test = file_search(sdssdir+outdir[ss],count=ntest)
     if ntest eq 0 then $
        ;; Directory doesn't exist so file doesn't exist
        spawn,'mkdir -p '+sdssdir+outdir[ss]
     mwrfits,hybconti,sdssdir+hyb_fil[ss],contihdr,/create,/silent
     mwrfits,eigconti,sdssdir+hyb_fil[ss],/silent
     spawn,'gzip -f '+sdssdir+hyb_fil[ss]
     if not keyword_set(silent) then $
        print,'sdss_dblfitconti: created ',sdssdir+hyb_fil[ss]

  endfor                        ; loop ss=nfil

  nfail = n_elements(fail_list) ; have to use n_elements()
  if nfail gt 1 then begin
     print,''
     print,'sdss_dblfitconti: spline failed on the following'
     print,fail_list[1:nfail-1],format='(a)'
     print,''
  endif 

  ;; Summaries
  print,'sdss_dblfitconti: all done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_dblfitconti: Elapsed time for '+strtrim(nfil,2)+$
        ' QSOs (m) = ',dt/60.
  print,'sdss_dblfitconti: Average time per QSO (s) = ',dt/nfil


  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end
