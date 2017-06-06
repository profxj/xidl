;+ 
; NAME:
;   EIGQSOCONTI 
;     Version 2.0
;
; AUTHOR:
;   Melodie M. Kao
;   Massachusetts Institute of Technology
;   Kavli Institute for Astrophysics and Space Research
;   77 Massachusetts Avenue, Building 37-287, Cambridge, MA 02139 
;   melodie.kao@alum.mit.edu
;   mkao@caltech.edu
;
; PURPOSE:
;    Fits a continuum to spectroscopic data using an iterative PCA
;    method.  You must supply your own eigenspectra.  A good resource
;    is C.W. Yip et al., 2004, AJ, 128, 2603.
;
; CALLING SEQUENCE:
;   eigflux = EIGQSOCONTI( wave, flux, error, eigspecfil, [gapmax=,     
;                FITMASK=, /DEBUG, WEIGHT=, CHI_SQR=,
;                TOPMARG=, LOWMARG=, MAXITER=, NITER=, STAT_TOL=
;                GROWMARG1=, CHUNKSZ=, CTOL=, FAIL=, STATUS=, _EXTRA=])
;
; INPUTS:
;  wave     -- Array of wavelength values
;  flux     -- Array of flux values   
;  error    -- Array of flux errors   
;  eigspecfil -- Eigenspectra as either array [npix,neig+1], fits file
;                of array, list to make into array.
;
; RETURNS: 
;  n x 2 array, where n = # pixels in original wave/flux arrays.
;  finaleig[*,0] = eigencontinuum
;  finaleig[*,1] = errors from eigencontinuum fit ONLY (*not* total
;                  errors including flux errors)
;
; OPTIONAL KEYWORDS:
;  /debug    -- Toggles on debugging plots
;  /plot -- shows final 
;  gapmax= -- number of contiguous pixels allowed to be above conti;
;             default 3
;  fitmask=  -- fixed mask to always apply
;  
;  growmarg1= --Scalar double: any pixels with flux less than 
;               -1*(mean + growmarg1*error) will be rejected PLUS the
;               pixels on either side of that pixel will also be
;               rejected; defualt = 2.5
;  topmarg=  -- Scalar double: In all iterations AFTER the first
;               fitting iteration, any pixels with flux greater than 
;               +1*(mean + topmarg1*error) will be rejected and not
;               fitted for consquent iterations. Default = 3.0
;  lowmarg=  -- Scalar double: In all iterations AFTER the first
;               fitting iteration, any pixels with flux less than 
;               -1*(mean + lowmarg1*error) will be rejected and not
;               fitted for consquent iterations.  Default = 2.
;  maxiter=  -- Scalar integer: maximum # of allowed fitting
;               iterations.  In general, 3 or 4 is enough.  Iterations
;               rarely exceed 10.  Default = 10
;  chunksz=  -- Scalar double: # of pixels per bin over which lowest
;               FCHUNK of flux values will be rejected for initial fit.
;               Gets converted into an integer with a FLOOR function
;               later.  This is used to prevent the first iteration
;               from fitting absorption lines and underfitting the
;               actual continuum, a problem which never gets resolved
;               without a good first fit.  Value should be larger than
;               broadest absorption features of interest.  Default =
;               50.0
;  fchunk=   -- Scalar float: fraction of lowest-flux pixels in each
;               "chunk" to mask out in the first interation. Default = 0.3
;  chforce=  -- Boolean (i.e., just set it): force instantiating the
;               prefitmask with the "chunk" masking even if FITMASK is
;               also set.
;  wv_full=  -- Wavelength range for final eigenconti (extrapolated
;               from wave)
;  idx_sub=  -- Indices in wv_full to which wave applies
;
; OPTIONAL OUTPUTS:
;  CHI_SQR= -- [chi^2, DOF, P(chi^2,/DOF)]
;  WEIGHT= -- final weights used
;  NITER= -- number of iterations until ending criterion met
;  STAT_TOL= -- history of median(flux - conti)
;  FAIL= -- error matrix status
;  STATUS= -- exit status
;  finalmask= -- final mask used
;
; COMMENTS:
;  PLOT INFORMATION:
;      
;
; EXAMPLES:
;   eigarr = eigqsoconti(waveFit, fluxFit, fluxerrFit, eigfile, 
;                        [FITMASK=fitmask, TOPMARG=5.0,
;                        LOWMARG=3.5,maxiter=3, GROWMARG1=3.5,
;                        CHUNKSZ=50.0, /DEBUG])
;
; PROCEDURES/FUNCTIONS CALLED:
;   eigenrecstrct.pro             
;   eigenrecstrct_mij.pro         
;   eigenrecstrct_mij_interpspec.pro     
;   getcolor.pro
;   x_splot.pro
;   legend.pro (IDLUTILS)
;
; REVISION HISTORY:
;   10 Jul 2011 Written by M. Kao 
;   28 Jul 2011 Major reorg and overhaul, KLC
;   31 Dec 2014  Enable FITMASK and "chunk" mask, KLC
;   15 Jul 2015  Correct growing prefitmask for absorption, KLC
;    6 Jun 2017  Enable extrapolation to larger wavelength range, KLC
;
;-
; Copyright (C) 2011, Melodie Kao
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;;;;
;------------------------------------------------------------------------------
@eigenrecstrct                  ; compile b/c function and want eigentrecstrct_mij

function eigqsoconti, wave, flux, error, eigenArrFull, GAPMAX=gapmax, $
                      FITMASK=fitmask, DEBUG=debug, HEADER=header, $
                      FINALMASK=finalmask, WEIGHT=weight, TOPMARG=topmarg,    $
                      LOWMARG=lowmarg, MAXITER=maxiter, NITER=niter, GROWMARG1=growmarg1,  $
                      CHUNKSZ=chunksz, FCHUNK=fchunk, CHFORCE=chforce, $
                      CTOL=ctol,  FAIL=fail, WV_FULL=wv_full, IDX_SUB=idx_sub, $
                      STATUS=status, CHI_SQR=chi_sqr, STAT_CTOL=stat_ctol, _EXTRA=extra
  
  IF  N_PARAMS() LT 4  THEN BEGIN 
     ;; KLC: Show that this is a function with the params in the Syntax
     print,'Syntax - ' + $
           'eigqsoconti( wave, flux, error, eigspecfil, [GAPMAX=, FITMASK=, WEIGHT=,'
     print,'               HEADER=, FINALMASK=, TOPMARG=, LOWMARG=, MAXITER=, '
     print,'               GROWMARG1=, CHUNKSZ=, FCHUNK=, CHFORCE=, CTOL=, FAIL=, STAT_CTOL=,'
     PRINT,'               NITER=, CHI_SQR=, STATUS=, STAT_CTOL=, /DEBUG, _EXTRA=   ] ) '
     return, -1
  ENDIF
  sdssdir = sdss_getsdssdir()
  npix = (size(wave,/dim))[0]

  ;;; Default values
  if not keyword_set(gapmax)    then gapmax    = 3 ; # contiguous pixels allowed to be high
  IF NOT KEYWORD_SET(topmarg)   THEN topmarg   =  3.0 ; high sigma clip
  IF NOT KEYWORD_SET(lowmarg)   THEN lowmarg   =  2.0 ; low sigma clip
  IF NOT KEYWORD_SET(maxiter)   THEN maxiter   =  10
  IF NOT KEYWORD_SET(growmarg1) THEN growmarg1 =  2.5 ; increase bounds on this low sigma clip
  IF NOT KEYWORD_SET(chunksz)   THEN chunksz   =  50.0
  IF NOT KEYWORD_SET(fchunk)    THEN fchunk    =  0.3 ; 30%
  if NOT KEYWORD_SET(ctol)      THEN ctol      =  1.e-4


  ;; _extra= includes: modes=, zQSO=, zin=
  allEigflux = eigenrecstrct_mij_interpspec(eigenArrFull, wave, _extra=extra)
  sz = size(allEigflux,/dim)
  nmodes = sz[0]
  
  ;;; Check to make sure that wave, flux, and error arrays are same
  ;;; length
  IF npix NE (size(flux,/dim))[0] or npix ne (size(error,/dim))[0] THEN $
     stop, "eigqsoconti: wave, flux, and error must be same length"
  IF npix ne sz[1] then $
     stop,'eigqsoconti: interpolating eigenspectra failed',npix,sz[1]

  ;;; Set values
  chi_sqr = fltarr(3)           ; [chi^2, DOF, P(chi^2,DOF)]
  stat_ctol = fltarr(maxiter)   ; store all values of convergence
  nanarr = replicate(!VALUES.D_NAN,npix)
  maskcount = intarr(npix)      ; number of times same pixel masked out

  sav_flux = flux               ; because want to modify below
  sav_error = error
  bd = where(error eq 0.)
  if bd[0] ne -1 then error[bd] = !VALUES.F_NAN

  ;;; Mask out pixels if given a fitting mask
  IF KEYWORD_SET(fitmask) THEN BEGIN
     if npix ne (size(fitmask,/dim))[0] then $
        stop,'eigqsoconti: spectrum and fitmask must be same length' ; Just do check here
     reject = WHERE( fitmask EQ 1 , nreject )
     IF nreject NE 0 THEN BEGIN
        flux[reject] = !VALUES.F_NAN
        error[reject] = !VALUES.F_NAN
     ENDIF
  ENDIF 
  
  ;;;;;;;;;; FIRST ITERATION: Fit only the top 70% of flux for every
  ;;;;;;;;;; chunksz pixels.  This prevents fitting absorption lines for
  ;;;;;;;;;; the first iteration and prevents rejecting MgII, Lya,
  ;;;;;;;;;; CaII bumps due to underfitting 
  tempnewflux    = flux
  tempnewfluxerr = error
 
  ;; KLC: don't want to be doing logical steps in halves, which
  ;; means the logic flow is:
  ;; 1) do the first fit after masking out bottom 30%
  ;;       a) calc chi^2 and save pre-fit mask and post-fit mask
  ;; 2) Enter while loop [save old conti and old mask, which is post-fit]
  ;; 3) Make new mask
  ;; 4) Fit
  ;;       a) calc chi^2 and save post-fit mask
  ;;       b) check if we're done


  ;;; Define top seven pixels with highest flux in every
  ;;; successive chunksz pixels.
  if not keyword_set(fitmask) or keyword_set(chforce) then begin
     nbinchnk      = LONG(FLOOR(npix/chunksz)) ; # of bins of chunksz pix; KLC: missing last bin
     bracket       = LINDGEN(chunksz)
     firstmask = INTARR(nbinchnk*(fchunk*chunksz)) ; initial rejected pix
     nbottompct = fix(fchunk*chunksz)            ; integer
     bottompct   = LINDGEN(nbottompct)
     FOR ii=0L, nbinchnk[0]-1 DO BEGIN
        indices         = bracket+ii*chunksz
        sortedFluxChunk =  SORT( flux[indices] ) 
        firstmask[bottompct+ii*nbottompct] = $
           indices[sortedFluxChunk[bottompct]]
     ENDFOR
     tempnewflux[firstmask] = !VALUES.D_NAN
     tempnewfluxerr[firstmask] = !VALUES.D_NAN
  endif
  prefitmask = where(finite(tempnewflux,/nan),nprefitmask)
  
  if keyword_set(debug) then begin
     print,''
     ymnx = [MIN(flux-lowmarg*error,/NAN), MAX(flux,/NAN)]
  endif 

  ;; ;;;;;;;;;
  ;; Initial fit; doesn't count as an iteration
  ;; _extra includes: /silent
  eigflux = EIGENRECSTRCT(wave, tempnewflux, tempnewfluxerr, allEigflux, $
                          rejectedflux=rejectedflux, WEIGHT=weight, $
                          _extra=extra) 
  
  postfitmask = WHERE( FINITE(rejectedflux) OR weight EQ 0.0 , npostfitmask, $
                       complement=valid, ncomplement=nvalid )     

  chi_sqr[0] = total((flux[valid]-eigflux[valid])^2/$
                    tempnewfluxerr[valid]^2, /NAN)
  chi_sqr[1] = nvalid - nmodes - 1.d
  chi_sqr[2] = chisqr_pdf(chi_sqr[0], chi_sqr[1])
  chisq = chi_sqr[0] / chi_sqr[1]


  
  ;;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;;; %%%%%%%%%%%%%%%%%%%% BEGIN debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;;; %%%%%%                                                            %%%%%%
  IF KEYWORD_SET(debug) THEN BEGIN 
     clr = GETCOLOR(/load)

     residuals        = tempnewflux - eigflux ; new residuals
     maskedflux       = nanarr
     rejectedresid    = nanarr
     IF npostfitmask NE 0 THEN $
        rejectedresid[postfitmask] = flux[postfitmask] ; masked in eigenrecstrct
     IF prefitmask[0] NE -1 THEN maskedflux[prefitmask] = flux[prefitmask] $
     ELSE maskedflux = 0        ; KLC: then just won't be plotted
     print, "# valid pixels = ", nvalid
     print, "mean residuals = ", mean(residuals[valid])
     print, "nprefitmask    = ", nprefitmask,$
            "        npostfitmask   = ", npostfitmask
     x_splot, wave,                                                       $
              flux,                        PSYM1=10, COLOR1=clr.black,    $
              YTWO=maskedflux,             PSYM2= 2, COLOR2=clr.seagreen, $
              YTHR=residuals,              PSYM3=10, COLOR3=clr.cyan,     $
              YFOU=topmarg*tempnewfluxerr, PSYM4=10, COLOR4=clr.blue,  $ ; show next cut
              YFIV=rejectedresid,          PSYM5=-6, COLOR5=clr.magenta,  $
              YSIX=-1.*lowmarg*tempnewfluxerr,PSYM6=10, COLOR6=clr.charcoal,$ ; next cut
              YSEV=eigflux,                PSYM7=10, COLOR7=clr.red,      $
              YMNX=ymnx,                      $
              TITLE = "Initial: " + $
              string(chisq,chi_sqr[2],format='("CHISQR=",f6.3," (P=",f6.3,")")'),  $
              /BLOCK,                                                     $
              LGND = ['flux','prev. mask','residuals','lowmarg',           $
                      'new rejections','topmarg','curr. conti']
     
  ENDIF             ; /debug
  ;;; %%%%%%                                                            %%%%%%
  ;;; %%%%%%%%%%%%%%%%%%%% END debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;   BEGIN Iterative Approach    ;;;;;;;;;;;;;;;;;;;;;;;;; 
  ;;;;;;;                                                              ;;;;;;;
  status = 'NOTDONE'
  niter = 0
  WHILE status EQ 'NOTDONE' DO BEGIN 
     if keyword_set(debug) then print,''

     ;; ;;;;;;;;;
     ;; Buffer old values
     oldeigflux = eigflux
     oldpostfitmask = postfitmask        
     noldpostfitmask = npostfitmask

     ;; ;;;;;;;;;
     ;; New masking, allowing previously rejected pixels back in, and
     ;; based on the results of previous conti fit
     tempnewflux    = flux
     tempnewfluxerr = error

     residuals = tempnewflux - oldeigflux

     prefitmask = WHERE( finite(tempnewflux,/nan)  or finite(tempnewfluxerr,/nan) or $
                         residuals LT -1.*lowmarg*tempnewfluxerr OR $
                         residuals GT topmarg*tempnewfluxerr OR $
                         maskcount GT 5, nprefitmask, $        ; yes, 5 is magic #
                         complement=valid, ncomplement=nvalid) ; for new masking
     
     IF nprefitmask NE 0 THEN BEGIN
        ;; Take the deepest absorption (presumably) and buffer those pixels
        absorp = WHERE( residuals[prefitmask] LT -1.*growmarg1*tempnewfluxerr[prefitmask], $
                        nabs )
        IF nabs NE 0 THEN BEGIN
           absorp = prefitmask[absorp]
           ;; KLC, 15 Jul 2015: added absorp to the extra mask and,
           ;; since maskxtr maps directly to e.g., flux array,
           ;; prefitmask is equal to a uniq subset of maskstr
           ;; (previously had prefitmask[uniq(...)])
           maskxtr = [absorp - 1, absorp, prefitmask, absorp + 1]
           prefitmask = maskxtr[uniq(maskxtr,sort(maskxtr))]
           prefitmask = prefitmask[where(prefitmask ge 0 and prefitmask le npix-1,$
                                         nprefitmask)] ; should always have pixels
        ENDIF

        tempnewflux[prefitmask]    = !VALUES.D_NAN
        tempnewfluxerr[prefitmask] = !VALUES.D_NAN

        ;; ;;;;;;;;;
        ;; Now find any large region where the conti is too low
        emis = where(residuals[prefitmask] gt 0. and $
                     finite(error[prefitmask]),nemis)
        if nemis ne 0 then begin
           emis = prefitmask[emis]
           gapstrt = where(emis ne shift(emis,1)+1, nemis) ; NaN !> 0
           gapstop = where(emis ne shift(emis,-1)-1,ntest)
           if nemis ne ntest then $
              stop,'eigqsoconti: did not find matching emission start/stop'
           gapsz = gapstop - gapstrt + 1
           gap = where(gapsz ge gapmax,nemis)
           if nemis ne 0 then begin
              gapstrt = emis[gapstrt[gap]]
              gapstop = emis[gapstop[gap]]
              for kk=0,nemis-1 do begin
                 tempnewflux[gapstrt[kk]:gapstop[kk]] = flux[gapstrt[kk]:gapstop[kk]]
                 tempnewfluxerr[gapstrt[kk]:gapstop[kk]] = error[gapstrt[kk]:gapstop[kk]]
              endfor 
              ;; Reset mask
              prefitmask = where(finite(tempnewflux,/nan),nprefitmask)
           endif                ; ngap != 0
        endif else nemis = 0    ; nemis != 0

     ENDIF ELSE BEGIN
        nabs = 0                ;  for printing purposes
        nemis = 0
        if keyword_set(debug) then $
           print,'eigqsoconti debug: troubled by nothing to mask out'
     endelse 


     ;; :::::::::
     ;; New fit and error
     ;; _extra includes: /silent
     eigflux = EIGENRECSTRCT(wave, tempnewflux, tempnewfluxerr, allEigflux, $
                             rejectedflux=rejectedflux, WEIGHT=weight, _extra=extra) 
     
     postfitmask = WHERE( FINITE(rejectedflux) OR weight EQ 0.0, npostfitmask, $
                          complement=valid, ncomplement=nvalid )

     ;; ;;;;;;;;;
     ;; Check exit criterion
     if nvalid eq 0 then begin
        ;; This will be bad chi^2
        chi_sqr[0] = total((flux-eigflux)^2/$
                           tempnewfluxerr^2, /NAN)
        chi_sqr[1] = npix - nmodes - 1.d
        chi_sqr[2] = chisqr_pdf(chi_sqr[0], chi_sqr[1])
        chisq = chi_sqr[0] / chi_sqr[1]

        if status ne 'NOTDONE' then status = 'NOPIX'
        if keyword_set(debug) then $
           print,'eigqsoconti debug: exiting fitting because no valid pixels left to fit.'
     endif else begin
        chi_sqr[0] = total((flux[valid]-eigflux[valid])^2/$
                           tempnewfluxerr[valid]^2, /NAN)
        chi_sqr[1] = nvalid - nmodes - 1.d
        chi_sqr[2] = chisqr_pdf(chi_sqr[0], chi_sqr[1])
        chisq = chi_sqr[0] / chi_sqr[1]

        ;; All comparisons to the region of the spectrum used in the
        ;; actual fit (excludes absorption lines).
        ;; Fractional change; just a metric
        test = where(abs(eigflux[valid] - oldeigflux[valid])/$
                     (0.5*(eigflux[valid] + oldeigflux[valid])) gt 1.e-5, ntest)

        ;; Average change across the board 
        stat_ctol[niter] = median(flux[valid]-eigflux[valid],/even)

        if keyword_set(debug) then $
           print,'eigqsoconti debug: niter = '+strtrim(niter,2) + $
                 '; fraction changed = '+strtrim(float(ntest)/nvalid,2) + $
                 '; median residual = '+strtrim(stat_ctol[niter],2)

        ;; Done Test 1
        if ntest eq 0 then begin
           if status eq 'NOTDONE' then status = 'NOCHANGE'
           if keyword_set(debug) then $
              print,'eigqsoconti debug: exiting fitting because no change.'
        endif 
        ;; Done Test 2
        if abs(stat_ctol[niter]) le ctol then begin
           if status eq 'NOTDONE' then status = 'HITCTOL'
           if keyword_set(debug) then $
              print,'eigqsoconti debug: exiting fitting because residual is small: '+$
                    strtrim(stat_ctol[niter],2)
        endif 
     endelse                    ;  done with goodness-of-fit checks

     ;; Done Test 3
     if npostfitmask ne 0 then begin
        maskcount[postfitmask] = maskcount[postfitmask] + 1
        test1 = nanarr
        test2 = nanarr
        if npostfitmask ne 0 then test1[postfitmask] = 1
        if noldpostfitmask ne 0 then test2[oldpostfitmask] = 1
        test = where(test1 eq test2, ntest) ; NaN != NaN, conveniently
        if ntest eq npostfitmask and ntest eq noldpostfitmask then begin
           if keyword_set(debug) then $
              print,'eigqsoconti debug: previous post-fit mask and current post-fit mask same.'
           if status eq 'NOTDONE' then status = 'EQFITMASK' ; these are obsure
        endif 
     endif

     niter = niter + 1        
     if niter eq maxiter then $
        if status eq 'NOTDONE' then status = 'MAXITER' ; others take precedence


     ;;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;;; %%%%%%%%%%%%%%%%%%%% BEGIN debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;;; %%%%%%                                                         %%%%%%
     IF KEYWORD_SET(debug) THEN BEGIN
        residuals        = tempnewflux - eigflux ; new residuals
        maskedflux       = nanarr
        rejectedresid    = nanarr
        IF npostfitmask       NE 0 THEN rejectedresid[postfitmask] = flux[postfitmask]
        IF prefitmask[0] NE -1 THEN $
           maskedflux[prefitmask] = flux[prefitmask] $
        else maskedflux = 0
        print, "# valid pixels = ", nvalid
        print, "mean residuals = ", mean(residuals[valid])
        print, "nprefitmask    = ", nprefitmask,$
               "; npostfitmask = ", npostfitmask, "; nabs = ", nabs, $
               "; nemis = ", nemis
        print, 'median residuals = ',stat_ctol[niter-1], $
               '        exit status = ',status
        x_splot, wave,                                                     $
                 flux,                       PSYM1=10, COLOR1=clr.black,   $
                 YTWO=oldeigflux,            PSYM2=10, COLOR2=clr.orange,  $
                 YTHR=maskedflux,            PSYM3= 2, COLOR3=clr.seagreen,$
                 YFOU=residuals,             PSYM4=10, COLOR4=clr.cyan,    $
                 YFIV=topmarg*tempnewfluxerr, PSYM5=10, COLOR5=clr.blue,     $
                 YSIX=rejectedresid,         PSYM6=-6, COLOR6=clr.magenta, $
                 YSEV=-1*lowmarg*tempnewfluxerr,PSYM7=10, COLOR7=clr.charcoal,     $
                 YEIG=eigflux,               PSYM8=10, COLOR8=clr.red,     $
                 YMNX=ymnx,                    $
                 TITLE = "Iter #"+strtrim(niter,2)+': '+$
                 string(chisq,chi_sqr[2],format='("CHISQR=",f6.3," (P=",f6.3,")")')+$
                 '---'+status,$
                 /BLOCK,                                                   $
                 LGND = ['flux','prev. conti','prev. mask','residuals',   $
                         'topmarg','new rejections','lowmarg','curr. conti']
        if status ne 'NOTDONE' then $
           stop,'eigqsoconti debug stop: end interation'
     ENDIF                      ; /debug
     ;;; %%%%%%                                                         %%%%%%
     ;;; %%%%%%%%%%%%%%%%%%%% END debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDWHILE                      ; status == 'NOTDONE'


  ;; Calculate error: and technically should compute chi^2 again...
  eigerror = 1                  ; trigger return
  eigflux2 = EIGENRECSTRCT(wave, tempnewflux, tempnewfluxerr, allEigflux, FAIL=fail, $
                          eigerror=eigerror, _extra=extra) 
  bd = where(eigflux2 ne eigflux,nbd)
  if nbd ne 0 then $
     stop,'eigqsoconti stop: WARNING! re-calculated fit with error != original',$
          nbd,median(eigflux2[bd]/eigflux[bd]-1,/even) ; frac diff

  

  IF keyword_set(wv_full) THEN BEGIN
     ;; Return eigenconti extrapolated to wavelength not included in fit
     ;; (e.g., Lya forest)
     npix_full = (size(wv_full,/dim))[0]
     IF npix_full LT npix THEN $
        stop,'eigqsoconti stop: wv_full smaller than wave'


     IF NOT keyword_set(idx_sub) THEN BEGIN
        print,'eigqsoconti: WARNING!!! Assuming extrapolating eigenconti blueward'
        idx_sub = lindgen(npix_full-npix) ; [0:npix_new-npix-1]
     ENDIF

     IF n_elements(idx_sub) NE npix THEN $
        stop,'eigqsoconti stop: idx_sub not same size as wave'

     ;; re-interpolate PCA spectra; _extra= includes: modes=, zQSO=, zin=
     newEigflux = eigenrecstrct_mij_interpspec(eigenArrFull, wv_full, _extra=extra)

     ;; expand arrays
     tempfullflux = replicate(!VALUES.D_NAN,npix_full)
     tempfullfluxerr = dblarr(npix_full) ; err = 0 is sign to make weight = 0

     tempfullflux[idx_sub] = tempnewflux
     tempfullfluxerr[idx_sub] = tempnewfluxerr

     ;; New "fit" (see caveates above)
     eigerror_full =  1         ; trigger return
     eigflux_full = EIGENRECSTRCT(wv_full, tempfullflux, tempfullfluxerr, $
                                  newEigflux, FAIL=fail_full, $
                                  eigerror=eigerror_full, _extra=extra)

     ;; Sanity check
     bd = where(eigflux_full[idx_sub] ne eigflux,nbd)
     if nbd ne 0 then $
        print,'eigqsoconti stop: WARNING! extrapolated fit with error != original',$
             nbd,median(eigflux_full[idx_sub[bd]]/eigflux[bd]-1.,/even) ; frac diff

     ;; Return values; chi_sqr, status, and fail set previously
     finaleig      = DBLARR(npix_full, 2, /nozero)
     finaleig[*,0] = eigflux_full
     finaleig[*,1] = eigerror_full
     tmpmask       = intarr(npix)
     if npostfitmask ne 0 then tmpmask[postfitmask] = 1 ; masked!
     finalmask     = replicate(1,npix_full)
     finalmask[idx_sub] = tmpmask
  ENDIF ELSE BEGIN
     
     ;; Return values; chi_sqr, status, and fail set previously
     finaleig      = DBLARR(npix, 2, /nozero)
     finaleig[*,0] = eigflux
     finaleig[*,1] = eigerror
     finalmask     = intarr(npix)
     if npostfitmask ne 0 then finalmask[postfitmask] = 1 ; masked!
  ENDELSE                                                 ; wv_full = 0

  
  if keyword_set(sav_flux) then begin
     ;; Restore so the external call doesn't get messed up
     flux = sav_flux
     error = sav_error
  endif



  ;; Header
  if keyword_set(header) then begin
     ;; Save information to header
     ;;  fail      -- Integer for each spectrum fitted:
     ;;               -1 = Small pivot element used to invert matrix,
     ;;                    significant accuracy probably lost
     ;;                0 = Matrix inversion successful, accuracy pretty good
     ;;                1 = Singular array, inversion is invalid, no
     ;;                continuum
     sxaddpar, header, 'NITER', niter, 'Number of eigqsoconti interations'
     sxaddpar, header, 'FAILFLG', fail, 'eigenrecstrct_mij FAIL status'
     sxaddpar, header, 'EXTSTAT', status, 'eigqsoconti EXIT status'
     sxaddpar, header, 'EXTCTOL', stat_ctol[niter-1], 'final median flux-conti residual'
     sxaddpar, header, 'CHISQR', chi_sqr[0], 'chisqr of continuum fit to flux'
     sxaddpar, header, 'CHISQDOF', chi_sqr[1], 'chisqr degrees of freedom'
     sxaddpar, header, 'PCHISQR', chi_sqr[2], 'chisqr probability'
     if keyword_set(wv_full) then $
        sxaddpar, header, 'WV_FULL', npix, 'number of fit pixels; extrapolated'
  endif 

  RETURN, finaleig
  ;;;;;;;                                                              ;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;   END Iterative Approach   ;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end
