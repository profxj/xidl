;+ 
; NAME:
; sdss_cleanspec
;    Version 1.0
;
; PURPOSE:
;   "Clean" SDSS spectrum by taking out automatically detected
;   features and random sampling back in error (onto continuum). 
;
; CALLING SEQUENCE:
;   
;   sdss_cleanspec( spec_fil, frac_rm, seed=, ewobs_lim=, cflg=, /debug, 
;                   save_fil, /clobber, /silent)
;
; INPUTS: 
;   spec_fil -- spectrum name
;   frac_rm -- Fraction (0 to 1) of sdss_fndlin centroids to clean out
;
; RETURNS: 
;   cleaned spectrum in SDSS format [npix,3] where flux is 0th, mask
;   is 1st, and sigma is 2nd 2D elements
;
; OUTPUTS: 
;
; OPTIONAL KEYWORDS:
;   ewobs_lim= -- observed EW limit for sdss_fndlin centroids  
;                 (default: 0.0 Ang)
;   cflg= -- conti-type to use (default: sdsscontistrct.cflg); if
;            input cflg and structure cflg do not m
;   seed= -- (default: IDL choose)
;   /save_fil -- (default file name from sdss_getname())
;   /clobber -- obey /save_fil even if cleanspec file exists
;   /silent -- minimize spurious messages
;   /debug -- print extra messages and plot
;
; OPTIONAL OUTPUTS:
;   cent_mask= -- sdsscontistrct.centroid masked out (of size [300])
;   oseed= -- last random seed used, can be used for next process
;
; COMMENTS: 
;
; EXAMPLES:
;    cleanspec = sdss_cleanspec(spec_fil,0.3,/debug,/save_fil,/clobber)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15 Jul 2011  Created by KLC
;   17 Jul 2015  Modify how regions combined/excluded, KLC
;-
;------------------------------------------------------------------------------

function sdss_cleanspec, wave, flux, sigma, abslin_fil, frac_rm, $
                         ewobs_lim=ewobs_lim, header=header, $
                         cflg=cflg, seed=seed, oseed=oseed, cent_mask=cent_mask, $
                         cleanspec_fil=cleanspec_fil,clobber=clobber, $
                         silent=silent, debug=debug

  if n_params() ne 5 then begin
     print,'Syntax - sdss_cleanspec( wave, flux, sigma, abslin_fil, frac_rm, [ewobs_lim=, cflg=, '
     print,'                         cleanspec_fil=, seed=, /clobber, /silent, /debug])'
     return, -1
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Defaults
  if not keyword_set(ewobs_lim) then ewobs_lim = 0. ; keep all positive
  if keyword_set(seed) then oseed = seed            ; else IDL handles

  if keyword_set(cleanspec_fil) then begin
     dir = strmid(cleanspec_fil,0,strpos(cleanspec_fil,'/',/reverse_search))
     test = file_search(sdssdir+dir,count=ntest)
     if ntest eq 0 then spawn,'mkdir -p '+sdssdir+dir
  endif 
     
  ;; Read in files
  if size(abslin_fil,/type) eq 8 then cstrct = abslin_fil $
  else cstrct = xmrdfits(sdssdir+abslin_fil,1,/silent)
  npix = (size(wave,/dim))[0]
  cleanspec = dblarr(npix,3,/nozero)
  if size(abslin_fil,/type) eq 8 then cstrct = abslin_fil $
  else cstrct = xmrdfits(sdssdir+abslin_fil,1,/silent)
  if npix ne cstrct.npix then $
     stop,'sdss_cleanspec() stop: npix != abslin.npix'
  if keyword_set(cflg) then begin
     if cflg ne cstrct.cflg and keyword_set(debug) then $
        print,'sdss_cleanspec() debug: input cflg != abslin cflg so boundaries may be off.'
     cstrct.cflg = cflg
  endif 
  cindx = fix(alog(cstrct.cflg)/alog(2))

  ;; Make complicated mask with 1 = good, -1 = (sigma=0), 0 = masked line
  mask = replicate(1,npix)
  bdpix = where(sigma eq 0.) ; prevents use of gaps
  if bdpix[0] ne -1 then mask[bdpix] = -1
  ;; Set up return mask 
  cent_mask = round(0*cstrct.centroid[*,cindx]) ; make integer


  ;; Abort if no reason to be here so do three tests:
  nrmcent = cstrct.ncent[cindx] ; 1) anything to clean out?
  if nrmcent ne 0 then begin
     ;; 2) anything strong enough to clean out?
     rmcent = where(cstrct.ew_orig[0:nrmcent-1] ge ewobs_lim,nrmcent)
     if nrmcent ne 0 then begin
        ;; 3) Remove a random fraction by drawing a uniformly random
        ;; number and keep those that have num_rand <= frac_rm
        ncentgewvobs = nrmcent
        num = randomu(oseed,nrmcent) ; uniform
        gd = where(num le frac_rm,nrmcent)
        if nrmcent ne 0 then rmcent = rmcent[gd] $ ; ones to keep
        else rmmsg = 'none <= frac_rm'
     endif else rmmsg = 'none >= ewobs_lim'
  endif else rmmsg = 'no centroids'

  if nrmcent eq 0 then begin
     if not keyword_set(silent) then $
        print,'sdss_cleanspec(): '+rmmsg+'; return original spectrum'
     cleanspec[*,0] = flux
     cleanspec[*,1] = mask eq 0 ; invert from 1 = good; 0 = bad
     cleanspec[*,2] = sigma

     if keyword_set(cleanspec_fil) then begin
        test = file_search(sdssdir+cleanspec_fil+'*',count=ntest)
        if ntest eq 0 or keyword_set(clobber) then begin
           ;; What if header is not set?
           mwrfits,cleanspec,sdssdir+cleanspec_fil,header,/create,/silent
           spawn,'gzip -f '+sdssdir+cleanspec_fil
        endif         
     endif

     return, cleanspec          ; EXIT
  endif 
  
  ;; Set up return mask 
  cent_mask[rmcent] = 1         ; removed = 1; kept = 0
  if keyword_set(debug) then begin
     print,'sdss_cleanspec() debug: '+strtrim(nrmcent,2)+' to clean of '+$
           strtrim(ncentgewvobs,2)+' with EW >= '+$
           string(ewobs_lim,format='(f5.2)')+' of '+$
           strtrim(cstrct.ncent[cindx],2)+' total'
     print,cstrct.centroid[rmcent,cindx]
  endif 

  ;; Set continuum and normalize
  conti = cstrct.conti[0:cstrct.npix-1,cindx]
  normfx = fltarr(cstrct.npix,/nozero)
  normsig = normfx
  normfx[cstrct.ipix0:*] = flux[cstrct.ipix0:*]/conti[cstrct.ipix0:*]
;  normsig = sdss_calcnormerr(flux,sigma,cstrct)
  normsig[cstrct.ipix0:*] = sigma[cstrct.ipix0:*]/conti[cstrct.ipix0:*] ; don't want conti error included
  if cstrct.ipix0 ne 0 then begin
     normfx[0:cstrct.ipix0-1] = 0.
     normsig[0:cstrct.ipix0-1] = 0.
  endif

  ;; Set up output
  cleanspec[*,0] = normfx
  cleanspec[*,2] = normsig

  ;; Use bounds of line defined in abslin structure (cstrct) via
  ;; sdss_ewciv; have to mask out all lines
  pixlim = intarr(cstrct.ncent[cindx],2,/nozero)
  for ll=0L,cstrct.ncent[cindx]-1 do begin
     dum = min(cstrct.wvlim_orig[ll,0]-wave,pixmin,/abs)
     dum = min(cstrct.wvlim_orig[ll,1]-wave,pixmax,/abs)
     mask[pixmin:pixmax] = 0    ; mask out
     pixlim[ll,*] = [pixmin,pixmax]
  endfor                        ; loop ll=cstrct.ncent

  ;; Find ranges of regions to be filled in
  bd = where(mask eq 0,nbd)
  istrt = bd[where(bd ne shift(bd,1)+1, nstrt)]
  istop = bd[where(bd ne shift(bd,-1)-1, nstop)]
  if nstrt ne nstop then $
     stop,'sdss_cleanspec() stop: did not find matching gap start/stop'

  ;; Randomly sample in regions surrounding known absorption features
  ;; so have to mask out all regions even if only scrubbing a fraction
  nrng0 = 5L                    ; just a magic number
  rng0 = lindgen(nrng0)+1
  ;; Check that there's enough gap between features; last one won't count
  bd = where(istop+nrng0+1 ge shift(istrt,-1),nbd,$
             complement=gd,ncomplement=ngd)
  if nbd gt 1 then begin
     nbd0 = nbd
     if nbd eq nstrt then begin
        ;; One big bin
        istrt = istrt[0]
        istop = istop[nstrt-1]
     endif else begin
        ;; Combine by rearranging bounds but need to find which gd is
        ;; going to get modified b/c not as simple as bd+1 since bd can
        ;; have consecutive numbers
        nwstrt = istrt
        nwstop = istop
        bd = bd[0:nbd-2]
        nwstrt[bd+1] = istrt[bd]
        nwstop[bd] = istop[bd+1]
        istrt = [nwstrt[gd],nwstrt[nstrt-1]] ; should never throw exception
        istop = [nwstop[gd],nwstop[nstrt-1]]
        bd = where(istop+nrng0+1 ge shift(istrt,-1),nbd,$
                   complement=gd,ncomplement=ngd)
        while nbd gt 1 do begin ; since last will always trigger
           stop,'sdss_cleanspec() stop: this part of revised code not tested as of 15 Jul 2015; proceed with caution'
           nwstrt = istrt[gd]   ; should never throw exception
           nwstop = istop[gd]
           bd = bd[0:nbd-2]
           nwstrt[bd+1] = istrt[bd]
           nwstop[bd+1] = istop[bd]
           istrt = nwstrt
           istop = nwstop
           bd = where(istop+nrng0+1 ge shift(istrt,-1),nbd,$
                      complement=gd,ncomplement=ngd)           
        endwhile                ; nbd > 1
        
        if keyword_set(debug) then $
           print,'sdss_cleanspec() debug: combined '+strtrim(nbd0-1,2)+' regions'
     endelse                    ; nbd < nstrt
  endif 
  
  
  for ii=0L,nrmcent-1 do begin
     ;; Find appropriate bounds given istrt and istop (and may recover
     ;; the same excised bounds because of the combining above)
     ibnd = (where(cstrct.centroid[rmcent[ii],cindx] ge wave[istrt] and $
                   cstrct.centroid[rmcent[ii],cindx] le wave[istop]))[0]
     if keyword_set(debug) then begin
        print,''
        print,'sdss_cleanspec() debug: starting region '+strtrim(ii,2)+': '+$
              string(cstrct.wvlim_orig[rmcent[ii],*],format='(2(f7.2,1x))')
     endif

     if ibnd eq -1 then begin
        ;; The centroid is actually at the edge of a border
        mn = min(cstrct.centroid[rmcent[ii],cindx] - wave[istrt],ibnd,/abs)
        print,'sdss_cleanspec(): assume can find bound with min()',$
              cstrct.centroid[rmcent[ii],cindx],wave[istrt[ibnd]]
     endif 

     ;; Reset
     cnt = 0L
     done = 0
     rng = rng0
     nrng = nrng0
     ilo = pixlim[rmcent[ii],0]

     ;; Loop until satisfied
     while not done do begin
        lcl = [istrt[ibnd]-rng,istop[ibnd]+rng] ; sample only neighboring pixels
        num = randomu(oseed,nrng*2)         ; uniform
        nsig = randomn(oseed,nrng*2)        ; Gaussian
        lcl = lcl[sort(num)]    ; random order

        ;; Check within spectrum bounds
        reg = where(lcl ge 0 and lcl lt npix,nreg)
        if nreg eq 0 then stop,'sdss_cleanspec() stop: all out of bounds'
        lcl = lcl[reg]

        ;; Check not sampling masked out pixels
        reg = where(mask[lcl] eq 1,nreg)
        if nreg eq 0 then begin
           ;; Increase boundaries but don't increment count
           if keyword_set(debug) then $
              stop,'sdss_cleanspec() debug stop: not enough suitable pixels in nrng = '$
                    +strtrim(nrng,2)+'; doubling size'
           nrng = 2*nrng
           rng = lindgen(nrng) + 1
           if nrng gt npix then begin
              stop,'sdss_cleanspec() stop: not enough suitable pixels in spectrum'
           endif else $
              continue          ; loop again
        endif 
        lcl = lcl[reg]

        ;; Add noise to normalize, masked out region
        ihi = (ilo+nreg-1) < pixlim[rmcent[ii],1]
        lcl = lcl[0:(ihi-ilo)]

        ;; Using the error like this has obvious (and probably
        ;; detrimental) effects around emission lines; consider
        ;; interpolating. 
        ;; However, interpolating will cause the error in bad sky
        ;; subtraction to be wrong.
        ;; Maybe will have to do some local characterization of the
        ;; error and sample from that. 
        cleanspec[ilo:ihi, 0] = 1. + nsig*normsig[lcl]
        cleanspec[ilo:ihi, 2] = normsig[lcl]

        ;; Set up for next loop
        ilo = ihi + 1
        if ilo gt pixlim[rmcent[ii],1] then done = 1 ; finished!

        cnt = cnt + 1
        if (cnt mod 3) eq 0 then nrng = nrng * 2 ; may need to go farther
        if keyword_set(debug) then $
           print,'         ...while pass #',cnt

     endwhile                                    ; loop while not done

  endfor                                         ; loop ii=nstrt
  

  ;; Output
  cleanspec[*,0] = cleanspec[*,0]*conti
  cleanspec[*,1] = mask eq 0    ; invert from 1 = good; 0 = bad
  cleanspec[*,2] = cleanspec[*,2]*conti
  if cstrct.ipix0 ne 0 then begin
     ;; Make regular flux
     cleanspec[0:cstrct.ipix0-1,0] = flux[0:cstrct.ipix0-1]
     cleanspec[0:cstrct.ipix0-1,1] = 0
     cleanspec[0:cstrct.ipix0-1,2] = sigma[0:cstrct.ipix0-1]
  endif 

  test = where(finite(cleanspec) eq 0)
  if test[0] ne -1 then stop,'sdss_cleanspec() stop: cleanspec has non-finite values'

  if keyword_set(cleanspec_fil) then begin
     ;; Check
     test = file_search(sdssdir+cleanspec_fil+'*',count=ntest)
     if ntest eq 0 or keyword_set(clobber) then begin
        mwrfits,cleanspec,sdssdir+cleanspec_fil,header,/create,/silent
        spawn,'gzip -f '+sdssdir+cleanspec_fil
        if not keyword_set(silent) then $
        print,'sdss_cleanspec(): wrote ',cleanspec_fil
     endif else print,'sdss_cleanspec(): will not clobber ',cleanspec_fil
  endif 

  if keyword_set(debug) then begin
     ;; Show output
     clr = getcolor(/load)
     gap = where(mask eq 0)
     x_splot,wave,flux,psym1=10,color1=clr.black,$
             ytwo=sigma,psym2=10,$
             ythr=cleanspec[*,0],psym3=10,color3=clr.blue,$
             yfou=cleanspec[*,2],psym4=10,color4=clr.magenta,$
             xfiv=wave[gap],yfiv=cleanspec[gap,0],psym5=4,color5=clr.purple,$
             ysix=conti,psym6=10,color6=clr.limegreen,$
             title='Original Spectrum and Clean',/block,$
             lgnd=['Flux','Error','Clean flux','Clean error','Gaps','Conti']
  endif 

  oseed = oseed[0]   ; not the array
  return, cleanspec
end
