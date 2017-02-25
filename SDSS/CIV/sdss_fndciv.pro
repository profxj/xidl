 ;+ 
 ; NAME:
 ; sdss_fndciv   
 ;    Version 1.0
 ;
 ; PURPOSE:
 ;    Brute force algorithm which examines the absorption lines
 ;    detectedin a SDSS spectrum and searches for matches with CIV
 ;
 ; CALLING SEQUENCE:
 ;  qalstrct = sdss_fndciv(abswav,'CIV',zem)
 ;
 ; INPUTS:
 ;   abswave - Array of observed wavelengths of detected absorption
 ;            lines from sdss_qsolin or sdss_fndlin
 ;   dblt_name - Doublet name or structure (compatible with
 ;               dblt_retrieve())
 ;   zqso - redshift of QSO 
 ;
 ; RETURNS:
 ;    qstrct -- {qalcharstrct} formatted structure of candidate
 ;               doublets; qstrct.ndla2 = 0 if no candidates.
 ;               Tag zlim2 added.
 ;   
 ; OUTPUTS:
 ;
 ; OPTIONAL KEYWORDS:
 ;    orig -- use original search mechanism (loop over z subarry)
 ;            instead of matching up by centroids
 ;    istrct -- input {qalcharstrct} formatted structure for candidate
 ;              doublets
 ;    zmin= -- sets min wavelength searched (default: careful
 ;             consideration of Lya or OI/SiII forest, rest doublet
 ;             wavelength and SDSS spectral coverage)
 ;    zmax= -- sets max wavelength searched (default: similar to zmin)
 ;    dvtol= -- allowed velocity offset to define candidate doublet
 ;              (default; 150 km/s)
 ;    dvgal= -- minimum redshift to search, used to exclude Galactic
 ;              (default: +5000 km/s)
 ;    dvqso= -- set max redshift search relative to QSO z (default:
 ;              -3,000 km/s)
 ;    dvbal= -- set max redshift search relative to BAL QSO z (default:
 ;              -10,000 km/s)
 ;    flg_bal -- if set, then use dvbal instead of dvqso
 ;    dz= -- redshift increments to search (default: 0.0001)
 ;    lsnr2= -- for un-paired candidate lines, use this statistical
 ;              significance to look for partner; requires that
 ;              wave, flux, sigma, cstrct_fil be set but won't
 ;              happen unless lsnr also set
 ;    wave= -- wavelenth array
 ;    flux= -- normalized flux
 ;    sigma= -- normalized error of spectrum
 ;    cstrct_fil= -- continuum structure
 ;    /silent -- turn off print statements (should be faster)
 ;    /debug -- print some stuff
 ;
 ; OPTIONAL OUTPUTS:
 ;
 ; COMMENTS:
 ;
 ; EXAMPLES:
 ;   see sdss_civsearch.pro
 ;
 ;PROCEDURES/FUNCTIONS CALLED:
 ;
 ; REVISION HISTORY:
 ;   07-May-2002 Written by JXP
 ;   28-Apr-2011 Modified min search wavelength, KLC
 ;   30-May-2011 Made generic doublet search, KLC
 ;   06-Jun-2011 Turned into function and allowed input, KLC
 ;   09-Jun-2011 Enable searching for lower S/N weaker doublet line, KLC
 ;-
 ;------------------------------------------------------------------------------

@sdss_fndlin                    ; resolve sdss_fndlin_srch()
function sdss_fndciv_prslin, zarr, score, sc, istrct, grade=grade, $
                             gapsz=gapsz, silent=silent, $
                             debug=debug
  ;; Take the main (/orig) sdss_fndciv result and dish it out with the
  ;; returned number of lines 
  if n_params() ne 4 then begin
     print,'Syntax - sdss_fndciv_prslin(zarr, score, sc, istrct, [grade='
     print,'                            gapsz=, /silent, /debug]'
     return,-1
  endif 
  if not keyword_set(grade) then grade = 1.0   ; perfect
  if not keyword_set(gapsz) then gapsz = 0.01  ; in  redshift space
  qstrct = istrct

  ;; Cap number of candidates
  ;; Save time by not setting things to zero; just be careful to not
  ;; access things beyond what's instantiated
  nmx = 100
  z = dblarr(nmx,/nozero)
  zlim = dblarr(nmx,2,/nozero)  ; added to {qalcharstrct}
  dlasc = fltarr(nmx,/nozero)
  hits = fltarr(nmx,/nozero)

  iscore = where(score EQ grade, ngscore)

  if iscore[0] NE -1 then begin

     ghits = sc[iscore]
     gscore = score[iscore]
     gz = zarr[iscore]
     dgz = gz-shift(gz,1)
     dgz[0] = 0.
     gapzindex = where(dgz GT gapsz, ngap) 

     ;; Use mean instead of median because may have even number of
     ;; pixels; don't use max anymore because that skews the result
     case ngap of
        0: begin
           z[0] = mean(gz)
           zlim[0,*] = [min(gz,max=mx),mx]
           dlasc[0] = max(gscore)
           hits[0] = max(ghits)
        end
        1: begin
           gapz =  gz[gapzindex]

           window1 = where(gz LT gapz[0])
           z[0] = mean(gz[window1])
           zlim[0,*] = [min(gz[window1],max=mx),mx]
           dlasc[0] = max(gscore[window1])
           hits[0] = max(ghits[window1])

           window2 = where(gz GE gapz[0])
           z[1] = mean(gz[window2])
           zlim[1,*] = [min(gz[window2],max=mx),mx]
           dlasc[1] = max(gscore[window2])
           hits[1] = max(ghits[window2])
        end
        else: begin
           gapz =  gz[gapzindex]

           for j=0, ngap-1 do begin
              case j of 
                 0: begin
                    window = where(gz LT gapz[0])
                    z[j] = mean(gz[window])
                    zlim[j,*] = [min(gz[window],max=mx),mx]
                    dlasc[j] = max(gscore[window])
                    hits[j] = max(ghits[window])
                 end
                 ngap-1: begin
                    window1 = where(gz LT gapz[j] AND gz GE gapz[j-1])
                    z[ngap-1] = mean(gz[window1])
                    zlim[ngap-1,*] = [min(gz[window1],max=mx),mx]
                    dlasc[ngap-1] = max(gscore[window1])
                    hits[ngap-1] = max(ghits[window1])

                    window2 = where(gz GE gapz[j])
                    z[ngap] = mean(gz[window2])
                    zlim[ngap,*] = [min(gz[window2],max=mx),mx]
                    dlasc[ngap] = max(gscore[window2])
                    hits[ngap] = max(ghits[window2])
                 end
                 else: begin
                    window = where(gz LT gapz[j] and gz GE gapz[j-1])
                    z[j] = mean(gz[window])
                    zlim[j,*] = [min(gz[window],max=mx),mx]
                    dlasc[j] = max(gscore[window])
                    hits[j] = max(ghits[window])
                 end
              endcase           ; case j
           endfor               ; loop j=ngap
        end
     endcase                    ; ncase gap >= 0

     qstrct.nDLA2 = ngap + 1
     if qstrct.nDLA2 gt nmx then $
        stop,'sdss_fndciv_prslin: too many candidates found'

     ;; Must trim out uninstantiated parts
     srt = sort(z[0:qstrct.nDLA2-1])
     qstrct.DLA_zabs2[0:qstrct.nDLA2-1] = z[srt]
     qstrct.DLA_score2[0:qstrct.nDLA2-1] = dlasc[srt]
     qstrct.DLA_hits[0:qstrct.nDLA2-1] = hits[srt]
     qstrct.zlim2[0:qstrct.nDLA2-1,*] = zlim[srt,*]

     if not keyword_set(silent) or keyword_set(debug) then begin
        print, 'sdss_fndciv_prslin: z with score = '+string(grade,format='(f3.1)')+':'
        print, qstrct.DLA_zabs2[0:qstrct.nDLA2-1]
     endif  
  endif else qstrct.nDLA2 = 0

  return, qstrct

end                             ;  sdss_fndciv_prslin


function sdss_fndciv_setdblt, abswave, dblt_name, zlim, istrct, dvtol=dvtol, $
                              mask=mask, silent=silent, debug=debug
  ;; Take main (default) sdss_fndciv result and dish it out with the
  ;; number of lines 
  if n_params() ne 4 then begin
     print,'Syntax - sdss_fndciv_setdblt(abswave, dblt_name, zlim, istrct, [dvtol=,'
     print,'                             mask=, /silent, /debug]'
     return,-1
  endif 
  nmx = (size(istrct.dla_zabs1,/dim))[0] ; space limit in qstrct structure
  zcand = dblarr(nmx*2)         ; will be accepting duplicates
  nzcand = 0
  nabswave = (size(abswave,/dim))[0]
  mask = replicate(-1,nabswave)               ; has to be available to return
  cinv = 1./3.e5                               ; km^-1 s
  if not keyword_set(dvtol) then dvtol = 150.  ; km/s; matches absorbers

  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  qstrct = istrct


  ;; Find lines by dvtol and input abswave
  ;; note that sdss_fndlin used 6300 Ang which is a bit annoying
  ;; Assume everything stronger doublet line; look for weaker and visa
  ;; versa
  ;; Make the redshift boundaries a little fuzzy because /orig version
  ;; would find these edge things
  skylinwv = sdss_getskylinwave(dwv=dwv)
  sub = where(abs(abswave-skylinwv[0]) gt dwv[0] and $
              abs(abswave-skylinwv[1]) gt dwv[1] and $
              abswave ge (1+zlim[0]-dvtol*cinv*(1+zlim[0]))*dblt.wvI and $
              abswave le (1+zlim[1]+dvtol*cinv*(1+zlim[1]))*dblt.wvII, nsub)
  if nsub eq 0 then begin
     if keyword_set(debug) then $
        print,'sdss_fndciv_setdblt() debug: no lines detected away from sky lines and in zlim.'
     qstrct.nDLA2 = 0
     ;; Just abort
     return, qstrct
  endif

  subabswave = abswave[sub]
  submask = replicate(0,nsub)
  wvobsII = dblt.wvII * subabswave/dblt.wvI 
  ;; Try opposite search; picks up blends
  wvobsI = dblt.wvI * subabswave/dblt.wvII

  for ii=0,nsub-1 do begin
     ;; Going to use 2*dvtol to better match /orig, which essentially
     ;; allowed +/- dvtol from both centroids 
     mnII = min(abs(wvobsII[ii] - subabswave),imnII) ; closest 1550
     if mnII lt wvobsII[ii]*2.*dvtol*cinv then begin
        ;; flag 1548 b/c assumed 1550 has 1548
        submask[ii] = submask[ii] or 1
        ;; Take average for redshift
        zcand[nzcand] = 0.5*(wvobsII[ii]+subabswave[imnII])/dblt.wvII - 1.
        nzcand = nzcand + 1
        if nzcand eq (size(zcand,/dim))[0] then begin
           print,'sdss_fndciv_setdblt(): Cannot store all zcand. Expanding.'
           tmp = dblarr((size(zcand,/dim))[0]*2)
           tmp[0:nzcand-1] = zcand
           zcand = tmp
        endif
     endif 
     mnI = min(abs(wvobsI[ii] - subabswave),imnI) ; closest 1548
     if mnI lt wvobsI[ii]*2.*dvtol*cinv then begin
        ;; this might be superfluous except blend
        ;; flag 1550 b/c assumed 1548 has 1550
        submask[ii] = submask[ii] or 2 ; so 3 indicates blend
        ;; Take average for redshift
        zcand[nzcand] = 0.5*(wvobsI[ii]+subabswave[imnI])/dblt.wvI - 1.
        nzcand = nzcand + 1
        if nzcand eq (size(zcand,/dim))[0] then begin
           print,'sdss_fndciv_setdblt(): Cannot store all zcand. Expanding.'
           tmp = dblarr((size(zcand,/dim))[0]*2)
           tmp[0:nzcand-1] = zcand
           zcand = tmp
        endif
     endif

  endfor                        ; loop ii=nsub

  ;; Flesh out mask; -1 will be outside of bounds
  mask[sub] = submask

  if nzcand ne 0 then begin
     gd = where(zcand gt 0.,nz)
     zstr = string(zcand[gd],format='(f9.7)') ; makes uniq easier
     unq = uniq(zstr,sort(zstr))
     z = zcand[gd[unq]]

     ;; Output
     qstrct.nDLA2 = n_elements(z) ; cannot use size(/dim) b/c doesn't allow for 1 element
     while qstrct.nDLA2 gt nmx do begin
        ;; Make space instead of losing time
        print,'sdss_fndciv_setdblt(): WARNING! Cannot store all candidates. Expanding.'
        nmx = 2*(size(istrct.dla_zabs1,/dim))[0] ; double it (overkill?)
        tmp_qstrct = sdss_expandqalstrct(nmx,qstrct0=qstrct) ; preserve ZLIM2
        qstrct = sdss_cpstrct(qstrct, tmp_qstrct)
     endwhile                   ; now qstrct has enough space
     srt = sort(z)
     qstrct.DLA_zabs2[0:qstrct.nDLA2-1] = z[srt]
     qstrct.zlim2[0:qstrct.nDLA2-1,0] = z[srt] - dvtol*cinv*(1+z[srt]) ; default width
     qstrct.zlim2[0:qstrct.nDLA2-1,1] = z[srt] + dvtol*cinv*(1+z[srt]) ; default width
     
     if not keyword_set(silent) or keyword_set(debug) then begin
        print, 'sdss_fndciv_setdblt(): z of candidates found by centroids:'
        print, qstrct.DLA_zabs2[0:qstrct.nDLA2-1]
     endif  
  endif else qstrct.nDLA2 = 0


  return, qstrct
end                             ; sdss_fndciv_setdblt()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sdss_fndciv, abswave, dblt_name, zqso, cflg, istrct=istrct, dvtol=dvtol, $
                      dvgal=dvgal, dvqso=dvqso, dvbal=dvbal, flg_bal=flg_bal, $
                      zmin=zmin, zmax=zmax, dz=dz, lsnr2=lsnr2, wave=wave, $
                      flux=flux, sigma=sigma, cstrct_fil=cstrct_fil, $
                      orig=orig, debug=debug, silent=silent, _extra=extra

  if n_params() ne 4 then begin
     print,'Syntax - sdss_fndciv(abswave, dblt_name, zqso, cflg, [istrct=,dvtol=,dvgal=,'
     print,'                     dvqso=, dvbal=, /flg_bal, zmin=, zmax=, dz=, '
     print,'                     lsnr2=, wave=, flux=, sigma=, cstrct_fil=, '
     print,'                     /silent, /debug])'
     return,-1
  endif 

  ;; Defaults
  if not keyword_set(dz) then dz = 0.0001        ; redshift bins
  if not keyword_set(dvtol) then dvtol = 150.    ; km/s; matches absorbers
  if not keyword_set(dvgal) then dvgal = 5000.   ; km/s; outside Galaxy
  if not keyword_set(dvqso) then dvqso = -3000.  ; km/s; relative to QSO
  if not keyword_set(dvbal) then dvbal = -10000. ; km/s; relative to QSO if BAL
  ;; Faster to use already instantiated file
  if keyword_set(istrct) then qstrct = istrct $
  else qstrct = create_struct({qalcharstrct},'zlim2',dblarr(100,2))

  ;; Other parameters
  cinv = 1./3.e5                ; km^-1 s
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else  dblt = dblt_retrieve(dblt_name)

  ;; Figure out bounds and don't use defaults in sdss_measuresnr()
  tmp = sdss_measuresnr('blank',wvlim_obs=wvlim_obs,/no_snr,dvgal=dvgal,$
                        dvqso=dvqso, $
                        dblt_name=dblt, zqso=zqso)

  ;; Avoid sky lines if at all possible; if no candidates, don't go there
  if max(abswave) LE 8000. then wvlim_obs[1] = 8000. < wvlim_obs[1]

  ;; Make search (loop over little dz's) efficient by careful
  ;; consideration of limits
  if not keyword_set( ZMIN ) then begin
     ;; Set by SDSS spectral coverage, Lya or OI/SiII "forest", and
     ;; HVCs of Galaxy
     qstrct.start_wave = wvlim_obs[0]
     zmin = qstrct.start_wave/dblt.wvI - 1.
  endif 
  if not keyword_set( ZMAX ) then begin
     ;; Set by SDSS spectral coverage, emission line of QSO (modulo
     ;; dvqso)
     if keyword_set(flg_bal) then begin
        if keyword_set(debug) then $
           print,'sdss_fndciv debug: BAL flagged and dvbal adopted'
        qstrct.flg_bal = flg_bal ; see sdss_getbalflg()
     endif 
     qstrct.DLA_quality[0] = wvlim_obs[1]
     zmax = qstrct.DLA_quality[0]/dblt.wvII - 1.
  endif 

  if zmin gt zmax then begin
     ;; With the dv cuts, this can happen and just want to exit
     qstrct.nDLA2 = 0
     if not keyword_set(silent) then $
        print,'sdss_fndciv: QSO has no searchable spectrum after dv cuts'
     return, qstrct
  endif 
  qstrct.DLA_z[0] = zmax - zmin ; just to save to get to civstr

  if keyword_set(debug) then $
     print,'sdss_fndciv: z limits ',string(zmin,zmax,format='(f8.5,"--",f8.5)')
  skylinwv = sdss_getskylinwave(dwv=dwv)

  ;; ;;;;;;;;;
  ;; Find doublets; two flavors
  if keyword_set(orig) then begin
     ;; ;;;;;;;;;
     ;; Find lines by small z divisions
     ;; loop over z values from zmin to zmax in .0001 increments, at each z
     ;; calculate array of shifted wavelengths (abswavearr) corresponding
     ;; to each rest wavelength 
     nsearch = long(ceil((zmax-zmin)/dz))
     zarr = zmin + dz*dindgen(nsearch)
     score = fltarr(nsearch,/nozero)
     sc = replicate(0.,nsearch)
     
     if keyword_set(debug) then $
        print,'sdss_fndciv debug: number of dz bins = ',strtrim(nsearch)

     for j=0L, nsearch-1 do begin 
        abswavearr = [dblt.wvI,dblt.wvII]*(1.+zarr[j])
        
        ;; Check not in sky emission 5579 or 6302 +/- 5 Ang
        goodlines = where((abs(abswavearr-skylinwv[0]) GT dwv[0]) AND $ 
                          (abs(abswavearr-skylinwv[1]) GT dwv[1]), nposs) 
        if nposs NE 2 then continue ; skip
        
        ;;at each z, search for a match of each shifted rest wavelength with a
        ;;wavelength in the absorption wavelength array, within a certain
        ;;tolerance, and count how many matches are found 

        for i=0, 1 do begin     ; both lines of doublet
           x = where(abs(abswavearr[i] - abswave) $
                     LT abswavearr[i]*dvtol*cinv, count) 
           if (count NE 0) then sc[j] = sc[j]+1
        endfor
        
        ;;define score array as scores for each dz increment in z
        ;;define z arr as z values for each increment
        score[j] = sc[j]/float(nposs)
     endfor                     ; loop j=nsearch

     ;; Find discrete lines that represent stronger line of doublet
     ;; instantiate tags: nDLA2, DLA_zabs2, DLA_score, DLA_hits
     ;; _extra includes gapsz= 
     qstrct = sdss_fndciv_prslin(zarr, score, sc, qstrct, grade=1.0, $
                                 silent=silent, debug=debug, _extra=extra)

  endif else begin
     ;; ;;;;;;;;;
     ;; Find lines by dvtol and input abswave
     ;; note that sdss_fndlin used 6300 Ang even which is a bit annoying
     qstrct = sdss_fndciv_setdblt(abswave, dblt, [zmin,zmax], qstrct, dvtol=dvtol, $
                                  mask=cand_mask, $ ; return
                                  silent=silent,debug=debug)
  endelse 

  if keyword_set(lsnr2) then begin
     ;; Try alternate S/N cut on convolved spectrum, for lines without
     ;; partners. BEWARE: there are a lot of exit points in this
     ;; section when it doesn't find requisite stuff.
     ;; Find redshifts that pass muster above but don't have partner (yet)
     if keyword_set(orig) then begin
        ;; Find redshifts that pass muster above but don't have partner (yet)
        halfqstrct = sdss_fndciv_prslin(zarr, score, sc, qstrct, grade=0.5, $
                                        /silent, debug=debug, _extra=extra)

        ;; Just abort; no "half" doublets to check
        if halfqstrct.nDLA2 eq 0 then begin
           if keyword_set(debug) then $
              print,'sdss_fndciv debug: no "half" doublets to check'
           return, qstrct       ; EXIT
        endif 

        nwzmin = min(halfqstrct.DLA_zabs2[0:halfqstrct.nDLA2-1],max=mx)
        nwzmax = mx+dvtol*cinv < zmax
     endif else begin
        ;; Check if any sdss_fndlin line not 1548
        lftover = where((cand_mask and 1) ne 1,nlft)
        if nlft eq 0 then begin
           if keyword_set(debug) then $
              print,'sdss_fndciv debug: no un-paired lines to check'
           return, qstrct       ; EXIT
        endif 

        ;; Maximum range
        nwzmin = min(abswave[lftover]/dblt.wvII - 1.) > zmin
        nwzmax = max(abswave[lftover]/dblt.wvI - 1.)+dvtol*cinv < zmax
     endelse 
     
     ;; Do the search with bad regions masked out as done in
     ;; sdss_fndlin exactly.  Will return same lines given in call
     ;; with abswave plus the new LSNR2 lines.
     mask = replicate(1,(size(wave,/dim))[0])
     sig = sigma
     bdpix = where(sigma eq 0.)
     if bdpix[0] ne -1 then begin
        ;; You do need both in sdss_fndlin_srch(), as used previously 
        mask[bdpix] = 0
        sig[bdpix] = 9.e9       ; keep ivar real
     endif 
     ;; So occasionally going to pass in full spectrum when actually
     ;; want just redward of Lya but just going to let
     ;; sdss_fndlin_srch() handle it.
     cstrct = sdss_fndlin_srch(cstrct_fil, cflg, lsnr=lsnr2, debug=debug,$
                               mask=mask, wave=wave, flux=flux, sigma=sig)
     ;; don't bother setting or using cstrct.cflg; keep
     ;; consistent in the here and now
     cindx = fix(alog(cflg)/alog(2))
     if keyword_set(debug) then $
        print,'sdss_fndciv debug: number of old vs new centroids:',$
              (size(abswave,/dim))[0], cstrct.ncent[cindx]
     
     ;; Remove the centroids for anything already used
     for ii=0,qstrct.nDLA2-1 do begin
        mn = min(cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx]-$
                 dblt.wvI*(1+qstrct.DLA_zabs2[ii]),imn,/abs)
        cstrct.centroid[imn,cindx] = 0.
     endfor                     ; loop ii=qstrct.nDLA2
     gd = where(cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx] gt 0.)

     ;; Nothing new to check!
     if gd[0] eq -1 then begin
        if keyword_set(debug) then $
           print,'sdss_fndciv debug: no new centroids to check with LSNR2 > ',$
                 strtrim(lsnr2,2)
        return, qstrct          ; EXIT
     endif 

     ;; Search again for doublets, parroting the parent call but
     ;; avoiding calling this section (e.g. no lsnr2 set) and limiting
     ;; to be just the redshifts range possible to find dblt.ionI S/N >
     ;; LNSR from sdss_fndlin and dblt.ionII S/N > LSNR2 here.
     ;; Should return what's in qstrct
     nwqstrct = $
        sdss_fndciv(cstrct.centroid[gd,cindx], $
                    dblt_name, zqso, cflg, istrct=istrct, dvtol=dvtol, $
                    dvgal=dvgal, dvqso=dvqso, dvbal=dvbal, flg_bal=flg_bal,$
                    zmin=nwzmin,zmax=nwzmax, dz=dz, orig=orig, $
                    debug=debug, /silent)

     if nwqstrct.nDLA2 eq 0 then begin
        if keyword_set(debug) then $
           print,'sdss_fndciv debug: no new doublets with LSNR2 > ',$
                 strtrim(lsnr2,2)
        return, qstrct          ; EXIT
     endif 

     ;; Make sure any new doublets pulled in have dblt.ionI stronger
     ;; (LSNR) than dblt.ionII (found now with LSNR2) by comparing
     ;; redshifts in halfqstrct.DLA_zabs2 or leftover abswave elements
     for ii=0,nwqstrct.nDLA2-1 do begin
        if keyword_set(orig) then $
           mn = min(abs(nwqstrct.DLA_zabs2[ii]-$
                        halfqstrct.DLA_zabs2[0:halfqstrct.nDLA2-1]),imn) $
        else mn = min(abs(nwqstrct.DLA_zabs2[ii]-abswave[lftover]/dblt.wvI-1),$ ;?
                      imn)
        if mn lt dvtol*cinv*(1+nwqstrct.DLA_zabs2[ii]) then begin
           ;; Load it up
           qstrct.nDLA2 = qstrct.nDLA2 + 1
           qstrct.DLA_zabs2[qstrct.nDLA2-1] = nwqstrct.DLA_zabs2[ii]
           qstrct.DLA_score2[qstrct.nDLA2-1] = nwqstrct.DLA_score2[ii]
           qstrct.DLA_hits[qstrct.nDLA2-1] = nwqstrct.DLA_hits[ii]
           qstrct.zlim2[qstrct.nDLA2-1,*] = nwqstrct.zlim2[ii,*]
           qstrct.snr[qstrct.nDLA2-1] = lsnr2 ; flag this as a secondary thing
        endif else $
           if keyword_set(debug) then $
              print,'sdss_fndciv debug: new candidate does not have right LSNR ratio: ',$
                    nwqstrct.DLA_zabs2[ii]
     endfor                     ; loop ii=nwqstrct.nDLA2-1 

  endif                         ; lsnr2=

  return,qstrct                 ; EXIT

end 


