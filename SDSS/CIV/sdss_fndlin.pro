;+ 
; NAME:
; sdss_fndlin   
;   Version 4.0
;
; PURPOSE:
;    Fits a continuum to SDSS QSO spectroscopic data interactively
;    using a PCA as a guess for the Lya and CIV emission lines.  The
;    code then automoatically searches for absorption line features
;    redward of Lya.
;
; CALLING SEQUENCE:
;   sdss_fndlin, sdss_list
;
; INPUTS:
;  sdss_list -- List of QSOs to analyze
;  sdsssum -- Summary file of QSO properties [required]; in
;              same order as filename lists!
;
; RETURNS:
;
; OUTPUTS:
;  One FITS file per QSO which contains the continuum and a list of
;  absorption lines detected. The 0th extension is structure of
;  form {sdsscontistrct} and 1st extension is detected lines.
;
; OPTIONAL KEYWORDS:
;  lsnr= -- Statistical significance for an absorption line feature
;           [default: 3.5]
;  /clobber  -- Redo the whole search even if the file exists
;  /redo  -- Redo the actual search even if the file exists (skip
;            conti)
;  processor= -- two element array of [ith, total] processors for a
;                parallel run with sdss_runparallelsrch.sh script,
;                where 1 <= ith <= (total <= 4). Overrides istrt=.
;  istrt= -- index at which to start
;  /nocompress -- Do not gzip suprress the output
;  /eig -- use eigencontinuum fit
;  /silent -- turn off print statements (should be faster)
;  /debug -- print some stuff
;
; OPTIONAL OUTPUTS:
;  /tcpu -- print information about CPU time used
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_fndlin,getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list',
;               getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit'
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   09-Sep-2003 Written by JXP
;   25-May-2011 Adopted from M Kao's sdss_qsolin, KLC
;   08-Jun-2011 Added processor=, KLC
;   11-Jun-2011 Centroids can overflow; set error in *srch(), KLC
;   25-Jul-2011 Check that normalizing leaves finite numbers, KLC
;   30-Jul-2011 Do three searches in one on 3 conti's
;   16-Aug-2011 Conti errors stored as just due to conti, KLC
;   30-Jul-2014 Expand abslin structure if necessary, KLC
;   15-Jul-2015 Explicitly round gpk; don't let x_findgauss "tune", KLC
;   10-Aug-2016 Make nfind optional input, KLC
;   21-Jul-2017 Sanity check to have no negative/zero centroids, KLC
;-
;------------------------------------------------------------------------------

function sdss_fndlin_srch, contistrct_fil, cflg, wave=wave, flux=flux, $
                           sigma=sigma, mask=mask, lsnr=lsnr, debug=debug, $
                           nfind=nfind
  if n_params() ne 2 then begin
     print,'Syntax - sdss_fndlin_srch(contistrct_fil, cflg, [wave=, flux=,'
     print,'                          sigma=, mask=, lsnr=, /debug, nfind=]'
     return, -1
  endif 
  cindx = fix(alog(cflg)/alog(2))
  
  ;; Actually find the lines in the region dictacted by information in 
  ;; contistrct (to be written) or wave, flux, sigma
  if not keyword_set( LSNR ) then lsnr = 3.5 ; stored in cstrct.snr_conv[npix:*,cindx]
  if not keyword_set( nfind ) then nfind = 1000L ; orig. default

  if size(contistrct_fil,/type) eq 8 then cstrct = contistrct_fil $ ; input is structure
  else cstrct = xmrdfits(contistrct_fil,1,/silent)                  ; read in

  ;; Other information
  if keyword_set(wave) or keyword_set(flux) or keyword_set(sigma) then begin
     nwv = (size(wave,/dim))[0]
     nfx =  (size(flux,/dim))[0]
     nsg = (size(sigma,/dim))[0]
     if nwv ne nfx and nwv ne nsg then $
           stop,'sdss_fndlin_srch(): must have equal wave, flux, sigma arrays'
  endif else begin
     ;; Must instantiate
     stop,'sdss_fndlin_srch(): currently must set wave, flux, sigma'
  endelse 
  npix = (size(wave,/dim))[0]
  if npix eq cstrct.npix then begin
     ipix0_orig = cstrct.ipix0
     cstrct.ipix0 = 0                                    ; for now 
     if keyword_set(debug) then $
        print,'sdss_fndlin_srch() debug: setting ipix0 == 0 because given full spectrum.'
  endif 

  if not keyword_set(mask) then mask = replicate(1,npix) ; all good

  
  ;; Cut around sky emission 5579 and 6302
  skylinwv = sdss_getskylinwave(dpix=dpix)  ; dpix = [5,2] or about 1.5 Ang
  mn = min(abs(wave-skylinwv[0]), i55)
  mn = min(abs(wave-skylinwv[1]), i63) ; s to match sdss_fndciv...
  bdbracket = [i55-lindgen(dpix[0]+1),i55+lindgen(dpix[0])+1, $
               i63-lindgen(dpix[1]+1),i63+lindgen(dpix[1])]
  keep = WHERE(bdbracket GE 0 OR bdbracket LE npix-1, nkeep)
  IF nkeep NE 0 THEN BEGIN
     bdbracket = bdbracket[keep]
     mask[bdbracket] = 0
  ENDIF
  if keyword_set(debug) then begin
     gd = where(mask gt 0,ngd)
     print,'sdss_fndlin_srch() debug: keeping '+strtrim(ngd,2)+$
           ' pixels outside sky emission'
  endif 
  
  !except = 0                   ; suppress messages from x_findnpeaks() and zfitmax()
  x_findgauss, 1.-flux, 1./sigma^2, xpeak=xpeak, gsn=gsn, /raw, $
               sn=sn, nfind=nfind ; don't want to run out of space
  !except = 1                     ; default
   
  if keyword_set(debug) then $
     print,'sdss_fndlin_srch() debug: x_findgauss found '+$
           strtrim((size(xpeak,/dim))[0],2)+' lines'
  cstrct.snr_conv[cstrct.ipix0:cstrct.npix-1,cindx] = gsn
  cstrct.snr_conv[cstrct.npix:*,cindx] = lsnr ; GOOD PLACE TO STORE

  gd = where(gsn gt lsnr)
  if gd[0] ne -1 then begin
     dwv = wave - shift(wave,1)
     dwv[0] = dwv[1]
     cstrct.deltaz_orig[cindx] = total(dwv[gd]) ; can't measure z b/c no wrest
  endif 

  ;; Map mask on original array to sn array
  tmpxpeak = xpeak
  tmpsn = sn
  maskpix = WHERE(mask EQ 0, nmaskpix)
  breakstart = WHERE(maskpix NE SHIFT(maskpix, 1)+1,countbreaks)
  breakstop  = WHERE(maskpix NE SHIFT(maskpix,-1)-1,ntest)
  if countbreaks ne ntest then $
     stop,'sdss_fndlin_srch(): did not find matching breakstart/stop'
  xarr = LINDGEN(npix)
  
  IF countbreaks NE 0 THEN BEGIN
     breakstart = maskpix[breakstart]
     breakstop = maskpix[breakstop]
     FOR ibreaks = 0L, countbreaks-1 DO BEGIN
        reject = WHERE( xpeak GT xarr[breakstart[ibreaks]] AND $
                        xpeak LT xarr[breakstop[ibreaks]], nrej )
        IF nrej NE 0 THEN BEGIN
           tmpxpeak[reject] = !VALUES.D_NAN
           tmpsn[reject] = !VALUES.D_NAN
        END
     ENDFOR                     ; loop ibreaks=countbreaks
  ENDIF

  ;; Trim ends (with sigma cut)
  good = WHERE(tmpsn GT lsnr AND FINITE(tmpsn), ngood)

  if keyword_set(debug) then $
     print,'sdss_fndlin_srch() debug: keeping '+strtrim(ngood,2)+$
           ' lines with S/N > '+strtrim(lsnr,2)

  if ngood NE 0 then begin
       
     ;; since gpk aren't integers, round 
     gpk = round(xpeak[good])
     
     ;; Output!
     ;; Just the unique ones! 
     indx = lindgen(npix)
     wv = wave[gpk]   
     unq = uniq( wv, sort(wv) ) 
     nunq = n_elements(unq)     ; cannot be size() if unq only 1
     cstrct.ncent[cindx] = nunq ;< (size(cstrct.centroid[*,cindx],/dim))[0] ; available in sdsscontistrct
     ncent_space = (size(cstrct.centroid[*,cindx],/dim))[0]
     while nunq gt ncent_space do begin
        ;; Make space instead of losing them
        print,'sdss_fndlin_srch(): WARNING! Cannot store all centroids. Expanding.'
        npix_space = (size(cstrct.conti[*,cindx],/dim))[0]
        ncent_space = 2*ncent_space ; double it (overkill?)
        tmp_cstrct = sdss_expandcontistrct(npix_space,nlin=ncent_space)
        cstrct = sdss_cpstrct(cstrct, tmp_cstrct)
     endwhile                   ; now cstrct has enough space
     if keyword_set(debug) then $
        print,'sdss_fndlin_srch() debug: keeping '+strtrim(cstrct.ncent[cindx],2)+$
              ' unique lines'

     ;; Trim
     gpk = gpk[unq]
     
     ;; Interpolate fractional pixels to wavelength array and
     ;; figure out the average uncertainty (knowing that SDSS
     ;; pixels aren't equal in Angstroms)
     wvcent = interpol(wave, indx, gpk)

     ;; Last sanity check (induced by error with a sdss_stackciv spec)
     test = where(wvcent gt 0,ntest)
     cstrct.ncent[cindx] = ntest
     cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx] = wvcent[test]
  endif                         ; ngood ne 0 for S/N cut
  
  if keyword_set(ipix0_orig) then cstrct.ipix0 = ipix0_orig ; restore
  
  return,cstrct

end                             ; sdss_fndlin_srch()


function sdss_fndlin_fitspline, wave, flux, sigma, zqso, premask=premask, chi_sqr=chi_sqr, $
                                cstrct_fil=cstrct_fil, everyn=everyn, sset=sset, $
                                maxrej=maxrej, lower=lower, upper=upper, $
                                nord=nord, groupbadpix=groupbadpix, sticky=sticky,$
                                flg_gap=flg_gap, $
                                debug=debug, title=title, bsplmask=bsplmask, $
                                pca_fil=pca_fil, pca_head=pca_head, $
                                qualflg=qualflg, silent=silent, $
                                nopca=nopca, _extra=extra
  ;; Fit the spline continuum
  if n_params() ne 4 then begin
     print,'Syntax - sdss_fndlin_fitspline(wave, flux, sigma, zqso, [premask=, chi_sqr='
     print,'                  cstrct_fil=, flg_gap=, /debug=, title=, '
     print,'                  pca_fil=, pca_head=, qualflg=, /nopca, '
     printk,'                 sset=, bsplmask=, /groupbadpix, /sticky, _extra=]) '
     return,-1
  endif 
  qualflg = 0

  ;; Defaults
  if not keyword_set(everyn) then everyn = 25
  if not keyword_set(maxrej) then maxrej = 20
  if not keyword_set(lower) then lower = 2.0
  if not keyword_set(upper) then upper = 3.0
  if not keyword_set(nord) then nord = 3

  ;; Take input, partially or fully insantiated sdsscontistrct
  if keyword_set(cstrct_fil) then begin
     if size(cstrct_fil,/type) eq 8 then cstrct = cstrct_fil $
     else cstrct = xmrdfits(cstrct_fil,1,/silent)
  endif else begin
     cstrct = { sdsscontistrct }
  endelse 
  cstrct.npix = (size(wave,/dim))[0]

  ;; Mask
  bdpix = where(sigma eq 0.,complement=gdiv) ; bad pixels; typically edges
  ivar = sigma * 0.
  ivar[gdiv] = 1./sigma[gdiv]^2

  ;; Logarithmic wavelength in emitted frame 
  rlam = alog10( wave / (1.+zqso) )
  flg_gap = 0
  
  if keyword_set(nopca) then begin
     ;; Skip all this PCA tilt business, likely because a stacked
     ;; spectrum continuum fit
     cstrct.ipix0 = 0 ; fit to all
     norm = flux[*]
     nrm_ivar = ivar[*]
     fconti = replicate(1.,  cstrct.npix)     
     title = 'No PCA for QSO, zqso = '+string(zqso,format='(f8.5)')  
     if keyword_set(debug) then $
        print,'sdss_fndlin_fitspline() debug: '+title
  endif else begin
     ;; Defaults
     wv_pca = alog10(1230.)
     if not keyword_set(pca_fil) then $
        pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
     if size(pca_fil,/type) eq 7 then $
        pca = xmrdfits(pca_fil, 0, pca_head, /silent) $
     else begin  
        pca = pca_fil 
        if not keyword_set(pca_head) then $
           stop,'sdss_fndlin_splconti(): must set pca_head' 
     endelse 
     
     ;; Some checks
     sz = size(pca,/dimensions)
     if sz[0] NE 2000L then stop,'sdss_fndlin_fitspline(): PCA not right dimensions'
     c0 = sxpar(pca_head, 'COEFF0')
     if c0 NE 3.0 then stop,'sdss_fndlin_fitspline(): PCA COEFF0 NE 3'
     
     xmed = dindgen(sz[0])*0.0001 + c0 ; 1000 Ang to 1584.5283 Ang
     ymed = pca[*,0]
     
     ;; Create PCA arrays
     delta_def = fltarr(sz[0])  ; 2000L
     d_ivar_def = fltarr(sz[0]) 
  
     ;; Call qso_tilt
     pca_qsotilt, rlam, flux, ivar, xmed, ymed, tilt, tiltivar, $
                  ANSWER=tilt_ans, FLG_TILT=flg_tilt, _EXTRA=extra ; pixlim=
  
  
     ;; Find wavelength edges for PCA analysis
     mn = min( abs(rlam - WV_PCA), ipstrt) ; wv/(1-zqso) = 1230 redward
     cstrct.ipix0 = ipstrt
     mn = min( abs(rlam - xmed[sz[0]-1]), ipend) ; wv/(1+zqso) = 1584.5 Ang
     
     ;; Check for bad flux at beginning
     if gdiv[0] GT ipend then begin
        ipend = 0L              ; gdiv set previously
     endif else begin
        ;; Fit when Lya is at lambda > 3800A but only if there's
        ;; not a large gap in these first 200 pixels; otherwise get
        ;; failure at choldc step
        test = where(ivar[0:ipend] eq 0., ntest)
        if ipend gt 200 and float(ntest)/ipend gt 0.9 then begin
           if not keyword_set(silent) then $
              print,'sdss_fndlin_fitspline(): large gap in beginning of high-z QSO, zqso = '+$
                    string(zqso,format='(f8.5)')
           qualflg = qualflg + 1
           ipend = 0L
           flg_gap = 1
           if keyword_set(debug) then $
              print,'sdss_fndlin_fitsplconit() debug: large gap in beginning of high-z QSO'
        endif
     endelse 
     
     if ipend GT 200 then begin
        if keyword_set(debug) then $
           print,'sdss_fndlin_fitspline() debug: this is a high-z QSO, zqso = '+$
                 string(zqso,format='(f8.5)')
        title = 'Fit when 1230(1+zqso) is at lambda > 3800A, zqso = '+$
                string(zqso,format='(f8.5)')
        
        ;; Find corresponding spots in xmed
        mn0 = min(abs(rlam[cstrct.ipix0]-xmed), ipc0)
        mn1 = min(abs(rlam[ipend]-xmed), ipc1)
        
        if (ipc1-ipc0) NE (ipend-cstrct.ipix0) then begin
           ;; Sub-arrays are not the same size but may be due to the
           ;; small rounding error in the min() call; try bumping one (in
           ;; the right direction, hence the sign sign)
           if (ipc1-ipc0) gt (ipend-cstrct.ipix0) then sign = -1 $
           else sign = 1.
           
           if mn1 gt mn0 then begin
              ;; Upper bound is more off; check that not maxing out
              if ipc1 + sign lt sz[0] then ipc1 = ipc1 + sign $
              else ipc0 = ipc0 - sign
           endif else begin
           ;;; Lower bound is more off; check that not running out of
           ;;; bounds 
              if ipc0-sign ge 0 then ipc0 = ipc0 - sign $
              else ipc1 = ipc1 + sign
           endelse 
           if (ipc1-ipc0) NE (ipend-cstrct.ipix0) then $
              stop,'sdss_fndlin_fitspline(): mismatch!'
           
           if keyword_set(debug) then $
              print,'sdss_fndlin_fitspline(): adjusted bounds in PCA'
        endif
        
        ;; Iterate 3 times to fix tilt
        for mm=0,2 do begin
           
           delta = delta_def
           delta[ipc0:ipc1] = tilt[cstrct.ipix0:ipend]/ymed[ipc0:ipc1] - 1.
           d_ivar = d_ivar_def
           d_ivar[ipc0:ipc1] = tiltivar[cstrct.ipix0:ipend] * ymed[ipc0:ipc1]^2
           
           eigenspec = double(pca[*,1:*])
           
           ;; pca_weight
           pca_weight = transpose( eigenspec ) * $
                        (d_ivar ## replicate(1,20) )
           alpha = pca_weight # eigenspec
           beta = (pca_weight # delta)[*]
           
           ;; Linear algebra
           choldc, alpha, p
           result = cholsol(alpha, p, beta)
           pca_fit = eigenspec # result
           
           
           ;; PCA Continuum
           conti = (1.+pca_fit[ipc0:ipc1]) * $
                   ymed[ipc0:ipc1] * tilt_ans[cstrct.ipix0:ipend]
           
           ;; Smooth
           conti = smooth(conti,7)
           
           ;; Reject 
           bad = where((conti-flux[cstrct.ipix0:ipend]) GT $
                       2.*sigma[cstrct.ipix0:ipend],nb)
           if nb NE 0 then begin
              nip = ipend-cstrct.ipix0+1
              grow = [bad, [(bad+1) < (nip-1)], [(bad-1) > 0]]
              ugrow = grow[uniq(grow,sort(grow))]
              tiltivar[cstrct.ipix0+ugrow] = 0.
           endif else begin
              if keyword_set(debug) then $
                 print,'sdss_fndlin_fitspline() debug: no change in tilt; break'
              break             ; quit loop because no change
           endelse
        endfor                  ; loop mm=0, 1, 2
        
     
        ;; Save
        fconti = tilt_ans
        fconti[cstrct.ipix0:ipend] = conti
     
     ;;;;;;;;;
        ;; Normalize
        norm = tilt[*,0]        
        norm[cstrct.ipix0:ipend] = flux[cstrct.ipix0:ipend]/conti
        
        nrm_ivar = tiltivar[*,0]
        nrm_ivar[cstrct.ipix0:ipend] = ivar[cstrct.ipix0:ipend]*conti^2
        
     endif else begin           ; Low z QSO
        norm = flux[*]
        nrm_ivar = ivar[*]
        fconti = replicate(1., cstrct.npix)
        if flg_gap eq 1 then $
           title = 'This is a high-z QSO treated like low-z, zqso = '+$
                   string(zqso,format='(f8.5)') $
        else title = 'This is a low-z QSO, zqso = '+string(zqso,format='(f8.5)')
        if keyword_set(debug) then $
           print,'sdss_fndlin_fitspline() debug: '+title
     endelse
  endelse ; PCA business

     
  ;; BSPLINE the red stuff; don't use premask previously
  ;; because just let PCA handle itself
  cut = norm[cstrct.ipix0:*]
  cutwv = rlam[cstrct.ipix0:*]
  cutiv = nrm_ivar[cstrct.ipix0:*] 
  if keyword_set(premask) then cutiv = cutiv * premask[cstrct.ipix0:*] ; 0 = bad, 1 = good
  ncut = cstrct.npix - cstrct.ipix0 

  sset = bspline_iterfit(cutwv, cut, invvar=cutiv, everyn=everyn, yfit=yfit, $
                         groupbadpix=groupbadpix, maxrej=maxrej, /silent, $
                         sticky=sticky, $
                         lower=lower, upper=upper, outmask=bsplmask, nord=nord)
  if cstrct.ipix0 ne 0 then $
     bsplmask = [replicate(0,cstrct.ipix0),bsplmask] ; 0=bad, 1=good

  ;; This may fail horribly for low S/N spectra; do a check
  bd = where(yfit eq 0.,nbd)
  if float(nbd)/ncut gt 0.9 then begin
     ;; 90% zero!   Remove inverse-variance weighting but mask out
     ;; bad portions
     gd = where(cutiv ne 0.)
;     tt = bspline_iterfit(cutwv[cstrct.ipix0+gd], cut[cstrct.ipix0+gd], $
     tt = bspline_iterfit(cutwv[gd], cut[gd], $
                          everyn=everyn, yfit=yfit2, outmask=bsplmask, $
                          groupbadpix=groupbadpix, maxrej=maxrej, /silent, $
                          sticky=sticky, nord=nord)
     tmp = replicate(0,cstrct.npix) ; 0 = bad
     tmp[cstrct.ipix0+gd] = bsplmask ; 1 = good
     bsplmask = tmp

     ;; Fix size
;     yfit[*] = 0. 
;     yfit[gd] = yfit2
     yfit = interpol(yfit2,cutwv[gd],cutwv)
     if not keyword_set(silent) then $
        print,'sdss_fndlin_fitspline(): WARNING! Spline failure; return unweighted '+$
              'spline'
     qualflg = qualflg + 2
     ;; Check again
     bd = where(yfit eq 0.,nbd)
     if float(nbd)/ncut gt 0.9 or $
        keyword_set(debug) then begin
        print,'sdss_fndlin_fitspline(): WARNING! Unweighted failed; blank spline.'
        qualflg = qualflg + 4 
     endif 
  endif 
  
  ;; CHK
  fconti[cstrct.ipix0:*] = fconti[cstrct.ipix0:*]*yfit ; bspline on top of PCA

  
  
  ;; Load continuum structure
  cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)] = fconti

  ;; Estimate error in the bspline region as:
  ;; abs(<F> - s) where <F> is the median flux in a window of everyn
  ;;           centered around given pixel, compared to the spline
  ;;           value s at that pixel, excluding masked out regions
  ;;           (which are interpolated later)
  sig_est = dblarr(cstrct.npix)
  if (qualflg and 4) eq 0 then begin
     ;; Can actually figure out error
     cut = flux
     bd = where(bsplmask eq 0)                ; includes any premask
     if bd[0] ne -1 then cut[bd] = !values.f_nan ; out of way, not medianed

     everyhalfn = fix(0.5*everyn)
;  if cstrct.ipix0 gt 1 then begin
;     ;; One big bin
;     sig_est[0:cstrct.ipix0-1] = abs(median(cut[0:cstrct.ipix0-1],/even,/double) - $
;                                     median(fconti[0:cstrct.ipix0-1],/even,/double))
;  endif 
     for ee=0,cstrct.npix-1 do begin
        ilo = 0 > (ee - everyhalfn)
        ihi = (cstrct.npix - 1) < (ee + everyhalfn + 1)
        sig_est[ee] = abs(median(cut[ilo:ihi],/even,/double) - fconti[ee])
     endfor 
     ;; any bin centered on masked-out portion, just interpolate over
     ;; because likely stuck in an absorption line which will go goofy
     ;; and just say it's unphysical to be greater than input error
     bd = where(finite(sig_est) eq 0 or bsplmask eq 0 or sig_est gt sigma,$
                complement=gd)
     if bd[0] ne -1 then begin
        tmp = interpol(sig_est[gd],wave[gd],wave)
        sig_est[bd] = tmp[bd]
     endif
  endif 
  cstrct.sigconti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)] = sig_est

  ;; Chi^2
  ;; probably should have the propogated error instead of just sigma
  if keyword_set(chi_sqr) then begin
     chi_sqr = fltarr(3)
     gd = where(sigma ne 0. and bsplmask eq 1,ngd)
     chi_sqr[0] = total((flux[gd]-cstrct.conti[gd,sdss_getcflg(/spl,/index)])^2/sigma[gd]^2)
     chi_sqr[1] = ngd - n_elements(sset.fullbkpt) - nord ; not size() b/c may be 1 element
     chi_sqr[2] = chisqr_pdf(chi_sqr[0],chi_sqr[1])
  endif 

  return, cstrct

end                             ; sdss_fndlin_fitspline


pro sdss_fndlin_fitallspl, list_fil, sdsssum, clobber=clobber,silent=silent, $
                           processor=processor, _extra=extra
  ;; Fit bspline continuum. Calls sdss_fndlin_fitspline()
  ;; _extra includes pixlim= for sdss_fndlin_fitspline(), recommended
  ;; to set equal to 7
  if n_params() ne 2 then begin
     print,'Syntax - sdss_fndlin_fitallspl, list_fil, sdsssum, [/clobber,processor=,_extra]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Read in
  readcol,list_fil,fil_spec,format='(a)',skip=1
  nfil = (size(fil_spec,/dim))[0]

  ;; Summary table
  sdsstab = xmrdfits(sdsssum, 1, /silent)
  if nfil NE (size(sdsstab,/dim))[0] then $
     stop,'sdss_fndlin_fitallspl: list and QSO structure must be same size'
  if keyword_set(processor) then begin
     ;; Takes precedence over istrt=
     sub = sdss_calcparalleljob(sdsstab, processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     if keyword_set(debug) then $
        print,'sdss_fndlin_fitallspl debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif else istrt = 0L ; default

  ;; Set continuum name by convention
  cfil = sdss_getname(fil_spec,/spec,/spl,dir=contidir) 
  test_fil = sdss_getname(sdsstab,/spl)
  bd = where(test_fil ne cfil)
  if bd[0] ne -1 then $
     stop,'sdss_fndlin_fitallspl: list and QSO structure not in same order'

  ;; Fil_Specs sdss_fndlin_fitspline() require; it could read them in but
  ;; better to pass them around this way (for now)
  wv_pca = alog10(1230.)
  pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
  pca = xmrdfits(pca_fil, 0, pca_head, /silent) 

  ;; Quality flag
  qualflg = intarr(nfil)

  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  ;; Loop
  for qq=istrt,nfil-1 do begin 
     
     ;; Overwrite?
     test = file_search(sdssdir+contidir[qq]+cfil[qq]+'*',count=ntest) 
     if ntest ne 0 and not keyword_set(clobber) then continue 

     if not keyword_set(silent) then $
        print,'sdss_fndlin_fitallspl: Fitting qq ',qq,' ',fil_spec[qq]
  
     ;; Read spectrum
     parse_sdss, sdssdir+fil_spec[qq], fx, wv, head=hdr, $
                 npix=npix, sig=er 

     ;; Create structure
     dat = dblarr(npix,3) 

     ;; Run fit
     chi_sqr = 1                ; force return
     cstrct = sdss_fndlin_fitspline(wv,fx,er,sdsstab[qq].z,bsplmask=bsplmask, $
                                    sset=sset, pca_fil=pca,pca_head=pca_head,$
                                    qualflg=qflg,chi_sqr=chi_sqr,$
                                    /groupbadpix,/sticky,_extra=extra) 
     qualflg[qq] = qflg
     ;; Create header and add keywords
     sxaddpar,hdr,'SPLQUAL',qflg,'Spline quality flag'
     sxaddpar,hdr,'NSPLORD',sset.nord,'Number of orders in bspline fit'
     sxaddpar,hdr,'NBKPTS',n_elements(sset.fullbkpt),'Number of breakpoints in bspline'
     sxaddpar,hdr,'CHISQR', chi_sqr[0], 'chisqr of continuum fit to flux'
     sxaddpar,hdr,'CHISQDOF', chi_sqr[1], 'chisqr degrees of freedom'
     sxaddpar,hdr,'PCHISQR', chi_sqr[2], 'chisqr probability'

     
     ;; Follow save convention
     dat[*,0] = cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)]
     dat[*,1] = bsplmask eq 0                          ; invert from 1 = good; 0 = bad
     dat[*,2] = cstrct.sigconti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)] 

     ;; Save and compress 
     test = file_search(sdssdir+contidir[qq],count=ntest)
     if ntest eq 0 then spawn,'mkdir -p '+sdssdir+contidir[qq]
     mwrfits,dat,sdssdir+contidir[qq]+cfil[qq],hdr,/create,/silent 
     mwrfits,sset,sdssdir+contidir[qq]+cfil[qq],/silent
     spawn,'gzip -f '+sdssdir+contidir[qq]+cfil[qq] 

  endfor                        ; for qq=nfil

  ;; Summaries
  if not keyword_set(silent) then $
     print, 'sdss_fndlin_fitallspl: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_fndlin_fitallspl: Elapsed time for '+strtrim(nfil,2)+$
        ' QSOs (m) = ',dt/60.
  print,'sdss_fndlin_fitallspl: Average time per QSO (s) = ',dt/nfil

  ;; Failures
  igap = where(qualflg eq 1,ngap) ; Gap towards blue
  iwgt = where(qualflg eq 2,nwgt) ; inverse variance weighting fail
  ibth = where(qualflg eq 3,nbth) ; Both
  ifail = where((qualflg and 4) eq 4,nfail) ; unweighted fit fail
  if ngap ne 0 then begin
     print,''
     print,'sdss_fndlin_fitallspl: Spectra with gap problem'
     print,fil_spec[igap],format='(a)'
  endif 
  if nwgt ne 0 then begin
     print,''
     print,'sdss_fndlin_fitallspl: Spectra with weighting problem'
     print,fil_spec[iwgt],format='(a)'
  endif 
  if nfail ne 0 then begin
     print,''
     print,'sdss_fndlin_fitallspl: Spectra with spline-fit fail'
     print,fil_spec[ifail],format='(a)'
  endif 
  if nbth ne 0 then begin
     print,''
     print,'sdss_fndlin_fitallspl: Spectra with gap and weighting problems'
     print,fil_spec[ibth],format='(a)'
  endif 

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end                             ;  sdss_fndlin_fitallspl


function sdss_fndlin_calcew, contistrct_fil, use_cflg=use_cflg, wave=wave, $
                             flux=flux, sigma=sigma, _extra=extra
  ;; Measure observed EW and bounds for all centroids
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fndlin_calcew(contistrct_fil [use_cflg=, wave=, flux=,'
     print,'                             sigma=, _extra=])'
     return, -1
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Read in file and set up params
  if size(contistrct_fil,/type) eq 8 then cstrct = contistrct_fil $ ; input is structure
  else cstrct = xmrdfits(sdssdir+contistrct_fil,1,/silent)                  ; read in
  if keyword_set(use_cflg) then cstrct.cflg = use_cflg
  cindx = fix(alog(cstrct.cflg)/alog(2))

  ;; Major check to see if anything needs to be done
  if cstrct.ncent[cindx] eq 0 then begin
     print,'sdss_fndlin_calcew(): no centroids to measure'
     return,cstrct
  endif 
 
  ;; Read in spectrum or check that everything passed in matches
  if not keyword_set(flux) then $
     parse_sdss, sdssdir+cstrct.sdss_obs[0], flux, wave, npix=npix, sig=sigma $
  else begin
     npix = (size(flux,/dim))[0]
     if not keyword_set(wave) or not keyword_set(sigma) then $
        stop,'sdss_fndlin_calcew(): must set wave and sigma arrays with flux'

     if npix ne cstrct.npix or npix ne (size(wave,/dim))[0] or $
        npix ne (size(sigma,/dim))[0] then $
           stop,'sdss_fndlin_calcew(): spectrum arrays not same size'
  endelse 

  ;; Set up continuum and fold in its errors
  conti = cstrct.conti[0:cstrct.npix-1,cindx]
  sig = sdss_calcnormerr(flux, sigma, cstrct, /unnorm, cflg=cstrct.cflg)

  ;; _extra includes /keepwvlim, /debug, /plot
  cstrct = sdss_ewciv(wave, flux, sig, conti, 'CIV', $
                      cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx], $
                      snr_conv=cstrct.snr_conv[0:cstrct.npix-1,cindx],$
                      istrct=cstrct, /generic, _extra=extra)
  
  return, cstrct
end                             ; sdss_fndlin_calcew()


pro sdss_fndlin_setallew, list_fil, clobber=clobber, silent=silent, $
                          processor=processor, _extra=extra
  ;; Loop through and do all the EW instantiation
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fndlin_setallew, list_fil, [/clobber, /silent, processor=, _extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Read in
  readcol,list_fil,fil_spec,format='(a)'
  outdir = fil_spec[0]
  fil_spec = fil_spec[1:*]
  nfil = (size(fil_spec,/dim))[0]

  if keyword_set(processor) then begin
     ;; Takes precedence over istrt=
     sub = sdss_calcparalleljob(fil_spec, processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     if keyword_set(debug) then $
        print,'sdss_fndlin_setallew debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif else istrt = 0L ; default

  contistrct_fil = sdss_getname(fil_spec,/spec,/abslin,dir=cdir)
  ;; Test output
  len = strpos(cdir[0],'/')
  if outdir ne strmid(cdir[0],0,len) then $
     outdir = outdir + strmid(cdir,len+1)
  outstrct_fil = outdir+contistrct_fil
  contistrct_fil = cdir+contistrct_fil

  for qq=istrt,nfil-1 do begin
     
     test = file_search(outdir[qq],count=ntest)
     if ntest eq 0 then $
        spawn,'mkdir -p ' + outdir[qq] $
     else $
        test = file_search(sdssdir+outstrct_fil[qq]+'*',count=ntest)

     if ntest eq 0 or keyword_set(clobber) then begin
        ;; _extra includes use_cflg, and things passed to sdss_ewciv()
        cstrct = sdss_fndlin_calcew(contistrct_fil[qq],_extra=extra)

        ;; Quick checks
        cindx = fix(alog(cstrct.cflg)/alog(2))
        if cstrct.ncent[cindx] eq 0 then continue ; 

        bd = where(cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,0] eq 0. or $
                   cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,1] eq 0. or $
                   cstrct.ew_orig[0:cstrct.ncent[cindx]-1] eq 0. or $
                   cstrct.sigew_orig[0:cstrct.ncent[cindx]-1] eq 0.)
        if bd[0] ne -1 then $
           stop,'sdss_fndlin_setallew: all centroid EWs not intantiated'

        ;; Write out
        mwrfits,cstrct,sdssdir+outstrct_fil[qq],/create,/silent
        spawn,'gzip -f '+sdssdir+outstrct_fil[qq]
        if not keyword_set(silent) then $
           print,'sdss_fndlin_setallew: wrote to ',outstrct_fil[qq]
     endif else $
        print,'sdss_fndlin_setallew: skipping ',outstrct_fil[qq]

  endfor                        ; loop qq=nfil

  print,'sdss_fndlin_setallew: All done!'

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu
end                             ; sdss_fndlin_setallew


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_fndlin, sdss_list, sdsssum, tcpu=tcpu, clobber=clobber, $
                 redo=redo, istrt=istrt, usesplfil=usesplfil, only_cflg=only_cflg, $
                 nocompress=nocompress, debug=debug, silent=silent, $
                 processor=processor, noBAL=noBAL, set_ew=set_ew, _extra=extra

  if  N_params() LT 2  then begin 
     print,'Syntax - ' + $
           'sdss_fndlin, sdss_list, sdsssum, [/CLOBBER, /REDO, ISTRT=, /debug,'
     print,'             /NOCOMPR, /TCPU, /usesplfil, /silent, only_cflg=,'
     print,'             /nocompress, processor=, /noBAL, /set_ew, _EXTRA=]'
     return
  endif 

  ;; can't use sdss_getsdssdir() because sdss_functions depend
  ;; on sdss_fndlin_fitspline() 
  sdssdir = getenv('SDSSPATH')+'/'+getenv('SDSSDR')+'/' 

  ;; Defaults
  pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'

  ;; Open pca just once and pass to sdss_fndlin_fitspline()
  pca = xmrdfits(pca_fil, 0, pca_head, /silent)
  
  ;; Read SDSS files
  ;; Parse out first line
  readcol, sdss_list, outdir, format='a', /silent, numline=1
  if not keyword_set(silent) then $
     print,'sdss_fndlin: Output directory set to $SDSSPATH/$SDSSDR/'+outdir[0]
  readcol, sdss_list, fil_sdss, format='A', /silent, skip=1
  nfil = (size(fil_sdss,/dim))[0]

  ;; Summary table
  sdsstab = xmrdfits(sdsssum, 1, /silent)
  if nfil NE (size(sdsstab,/dim))[0] then $
     stop,'sdss_fndlin: list and QSO structure must be same size'
  ;; Set BAL flags (_extra= includes /full); some IDL's can't return
  ;; the array balflg from sdss_fndbalqso() so do where() here
  ;; BALs
  if keyword_set(noBAL) then balflg = replicate(0,nfil) $
  else balflg = sdss_fndbalqso(sdsstab,/flgonly,_extra=extra)
  if not keyword_set(silent) then begin
     ibal = where(balflg gt 0,nbal)
     print,'sdss_fndlin: Number of BAL QSOs',nbal
  endif

  ;; Make names and set to whatever's written to sdss_list
  fil_abslin = sdss_getname(sdsstab,/abslin,root=qso_name,dir=absdir)
  test_fil = sdss_getname(fil_sdss,/spec,/abslin)
  bd = where(fil_abslin ne test_fil)
  if bd[0] ne -1 then $
     stop,'sdss_fndlin: list and QSO structure not in same order.'
  absdir = outdir[0] + strmid(absdir,strpos(absdir[0],'/')+1)
  cstrct_def = { sdsscontistrct }
  
  ;; Create other file names (faster to do bulk?)
  fil_eig = sdss_getname(sdsstab,/eig,dir=contidir)
  fil_hyb = sdss_getname(sdsstab,/hyb)
  if keyword_set(usesplfil) then $
     fil_spl = sdss_getname(sdsstab,/spl)
  
  if keyword_set(processor) then begin
     ;; Takes precedence over istrt=
     sub = sdss_calcparalleljob(sdsstab, processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     if keyword_set(debug) then $
        print,'sdss_fndlin debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif 
  if not keyword_set(istrt) then istrt = 0L ; default

  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  ;; LOOP
  for qq=istrt,nfil-1 do begin
     ;; Load output continuum structure
     cstrct = cstrct_def

     ;; Clobber or just REDO?
     fil_abslin[qq] = absdir[qq] + fil_abslin[qq]

     test = file_search(sdssdir+fil_abslin[qq]+'*',count=ntest) 
     if ntest NE 0 then begin
        if not keyword_set( clobber ) and keyword_set(redo) then begin
           ;; Do not overwrite but do read in and skip to finding
           ;; lines
           ;; Message
           if not keyword_set(silent) then $
              print, 'sdss_fndlin: Updating qq ', qq, ' ', fil_abslin[qq]

           ;; Read in (use test b/c has directory and .gz so
           ;; don't need xmrdfits)
           cstrct = mrdfits(test,1,/silent)
           
           ;; Parse SDSS 
           parse_sdss, sdssdir+fil_sdss[qq], flux, wave, NPIX=npix, SIG=sigma
        endif 
        if not keyword_set(clobber) and not keyword_set(redo) then begin
           ;; Do not overwrite and do not search again
           if not keyword_set(silent) then $
              print, 'Skipping ', qq,' '+fil_abslin[qq]
           continue
        endif 
        ;; Implicit else: if keyword_set(clobber) then just go ahead!
     endif
     if keyword_set(clobber) or ntest eq 0 then begin
        ;; Check if need to make directories
        test = file_search(sdssdir+absdir[qq],count=ntest)
        if ntest eq 0 then spawn,'mkdir -p '+sdssdir+absdir[qq]

        ;; Parse SDSS
        parse_sdss, sdssdir+fil_sdss[qq], flux, wave, NPIX=npix, SIG=sigma
        cstrct.sdss_obs[0] = fil_sdss[qq]
        cstrct.qso_name = qso_name[qq] ; JJJJJ-PPPP-FFF, sdss_getname()

        cstrct.z_qso = sdsstab[qq].z
        cstrct.ra = sdsstab[qq].ra 
        cstrct.dec = sdsstab[qq].dec 
        cstrct.mjd = sdsstab[qq].smjd 
        cstrct.plate = sdsstab[qq].plate 
        cstrct.fiber = sdsstab[qq].fiber 
        cstrct.rmag = sdsstab[qq].rtmag 

        cstrct.balflg = balflg[qq]
        cstrct.npix = npix

        ;; Message
        if not keyword_set(silent) then $
           print, 'sdss_fndlin: Begin qq ', qq, ' ', fil_sdss[qq]

     endif

     ;;;;;;;;;
     ;; SPLINE CONTINUUM
     ;;;;;;;;;
     if keyword_set(usesplfil) then begin
        ;; Read in and assume sdss_getname() returns the right
        ;; directory/ data-release
        if keyword_set(debug) then begin
           test = file_search(sdssdir+contidir[qq]+fil_spl[qq]+$
                              '*',count=ntest) ; x_chkfil() slow
           if ntest eq 0 then $
              stop,'sdss_fndlin: spline continuum file DNE'
        endif                   ; else assum it exits
        splconti = xmrdfits(sdssdir+contidir[qq]+fil_spl[qq],0,/silent)
        
        ;; Load continuum structure
        cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)] = splconti[*,0]
        cstrct.sigconti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)] = splconti[*,2]
        flg_gap = 0
     endif else $
        ;; _extra includes: pixlim=[10], everyn=[25], maxrej=[20],
                                   ;; lower=[2.0], upper=[3.0], nord=[3] 
        cstrct = sdss_fndlin_fitspline(wave, flux, sigma, cstrct.z_qso, $
                                       cstrct_fil=cstrct, $ ;inputs
                                       pca_fil=pca, pca_head=pca_head, $
                                       flg_gap=flg_gap, title=title, $ ; outputs
                                       silent=silent, debug=debug, /groupbadpix, /sticky, $
                                       _extra=extra) 
     
     ;;;;;;;;;
     ;; EIGENSPECTRA CONTINUUM
     ;;;;;;;;;
     ;; Get file and assume sdss_getname() returns the right
     ;; directory/data-release 
     test = file_search(sdssdir+contidir[qq]+fil_eig[qq]+'*',$
                        count=ntest) ; x_chkfil() too slow
     if ntest ne 0 then begin
        ;; Read in
        eigconti = xmrdfits(sdssdir+contidir[qq]+fil_eig[qq], 0,/silent) 
        cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/eig,/index)] = eigconti[*,0]
        cstrct.sigconti[0:cstrct.npix-1,sdss_getcflg(/eig,/index)] = eigconti[*,2]
     endif else begin
        if keyword_set(debug) then $
           print,'sdss_fndlin debug: eigen continuum file DNE'
        skip_eigen = 1
        if keyword_set(only_cflg) then $
           if only_cflg eq sdss_getcflg(/eig) then $
              stop,'sdss_fndlin: eigen continuum file DNE though targeting.'
     endelse
     
     ;;;;;;;;;
     ;; HYBRID CONTINUUM
     ;;;;;;;;;
     ;; Get file and assume sdss_getname() returns the right
     ;; directory/data-release 
     ;; Does not have to exist because may be trying to create it
     test = file_search(sdssdir+contidir[qq]+fil_hyb[qq]+'*',$
                        count=ntest) ; x_chkfil() too slow
     if ntest ne 0 then begin
        ;; Read in
        hybconti = xmrdfits(sdssdir+contidir[qq]+fil_hyb[qq], 0,/silent) 
        cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/hyb,/index)] = hybconti[*,0]
        cstrct.sigconti[0:cstrct.npix-1,sdss_getcflg(/hyb,/index)] = hybconti[*,2]
     endif else begin
        if keyword_set(debug) then $
           print,'sdss_fndlin debug: hybrid continuum file DNE'
        skip_hybrid = 1
        if keyword_set(only_cflg) then $
           if only_cflg eq sdss_getcflg(/hyb) then $
              stop,'sdss_fndlin: hybrid continuum file DNE though targeting.'
     endelse 
     
     ;;;;;;;;;
     ;; Find lines
     ;;;;;;;;;
     bdpix = where(sigma eq 0.) ; bad pixels; typically edges
     mask = replicate(1,cstrct.npix)
     if bdpix[0] ne -1 then mask[bdpix] = 0 ; mask out
     mn = min(wave-1230.*(1+cstrct.z_qso),ipstrt,/abs) ; replicate *fitspline() code
     cstrct.ipix0 = ipstrt                             ; may not get set otherwise
     if keyword_set(eig) and keyword_set(debug) then $
        print,'sdss_fndlin debug: using eigen conti'
     if keyword_set(spl) and keyword_set(debug) then $
        print,'sdss_fndlin debug: using spline conti'
     if keyword_set(hyb) and keyword_set(debug) then $
        print,'sdss_fndlin debug: using hybrid conti'

     for cc=0,2 do begin
        if keyword_set(only_cflg) then $
           if only_cflg ne 2^cc then continue
        if keyword_set(skip_eigen) and cc^2 eq sdss_getcflg(/eig) then $
           continue             ; don't want to search a blank conti
        if keyword_set(skip_hybrid) and cc^2 eq sdss_getcflg(/hyb) then $
           continue             

        ;; Normalize just redward of Lya
        ;; propogate errors
        fx = flux/cstrct.conti[0:cstrct.npix-1,cc]
        sig = sdss_calcnormerr(flux,sigma,cstrct,cflg=2^cc,baderrval=9.e9)

        bd = where(finite(fx) eq 0 or finite(sig) eq 0.,nbd) ; b/c conti bad
        if bd[0] ne -1 then begin
           fx[bd] = 0.
           sig[bd] = 9.e9
           if keyword_set(debug) then $
              print,'sdss_fndlin debug: normalizing made non-finite fx and sig: ',$
                    strtrim(nbd,2)+' pixels'
        endif
        cstrct.cflg = 2^cc ; unclear how this will affect /redo
        cstrct = sdss_fndlin_srch(cstrct,cstrct.cflg,$
                                  wave=wave[cstrct.ipix0:cstrct.npix-1],$
                                  flux=fx[cstrct.ipix0:*],sigma=sig[cstrct.ipix0:*], $
                                  mask=mask[cstrct.ipix0:cstrct.npix-1],$
                                  debug=debug,_extra=extra) ; lsnr=
     endfor                                                 ;loop cc=2


     ;; Measure EW and bounds for whatever cflg is set
     ;; Use the original error because function will do the extra
     ;; error itself
     if keyword_set(set_ew) then $
        cstrct = sdss_fndlin_calcew(cstrct,wave=wave,flux=flux,$
                                    sigma=sigma,debug=debug,_extra=extra)


     
     ;;;;;;;;;
     ;; Write all to file
     ;;;;;;;;;
     mwrfits, cstrct, sdssdir+fil_abslin[qq], /create, /silent

     if not keyword_set(NOCOMPRESS) then spawn, 'gzip -f '+sdssdir+fil_abslin[qq]
     
     if keyword_set(tcpu) and ((qq+1) mod 1000L) eq 0 and not keyword_set(silent) $
     then begin
        tt = systime(/seconds)
        dt = (tt-tlast)         ; seconds
        print,'sdss_fndlin: Elapsed time (m) = ',dt/60.
        print,'sdss_fndlin: Average time per qq (s) = ',dt/1000L
        tlast = tt
     endif 

     if keyword_set(debug) then print,' ' ; easier to read
  endfor                                  ; loop qq=nfil


  ;; Final messages
  if not keyword_set(silent) then $
     print, 'sdss_fndlin: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_fndlin: Elapsed time for '+strtrim(nfil-istrt,2)+$
        ' QSOs (m) = ',dt/60.
  print,'sdss_fndlin: Average time per QSO (s) = ',dt/(nfil-istrt)


  if keyword_set(processor) then $
     cpu, restore=save_cpu

end
