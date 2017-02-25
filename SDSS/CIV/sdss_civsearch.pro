;+ 
; NAME:
; sdss_civsearch   
;    Version 2.0
;
; PURPOSE:
;    Given a list of SDSS quasars, search these spectra for CIV
;    absorbers and output a structure of detections and EWs
;
; CALLING SEQUENCE:
;  sdss_civsearch, filename,sdsssum,  cand_fil, zmin=, $
;                      rmax=
;
; INPUTS:
;   filename - File containing a (long) list of SDSS quasar data files
;              (1d spectra, .fit)
;   sdsssum - Summary file of QSO properties [required]; in
;              same order as filename lists!
;   cand_fil - Name for outputted structure of CIV absorbers
;
; RETURNS:
;
; OUTPUTS:
;  cand_fil - Name for outputted structure of CIV absorbers, the
;             first extension has the full sorted structure array
;             and the second, just the broad lines
;
; OPTIONAL KEYWORDS:
;  zmin= -- Minimum redshift of the doublet [default: 3820/wvI-1]
;  zmax= -- Maximum redshift of the doublet [default: 9200/wvII-1]
;  dblt_name= -- Doublet name or structure (compatible with
;               dblt_retrieve())
;  rmax= -- limiting R magnitude for QSOs [default: 23.]
;  processor= -- two element array of [ith, total] processors for a
;                parallel run with sdss_runparallelsrch.sh script,
;                where 1 <= ith <= (total <= 4). 
;  /plot -- show some stuff
;  /debug -- print some stuff
;  /silent -- turn off print statements (should be faster)
;  /excl_broad -- do not search for broad lines as likely candidates
;  _extra= -- passed to sdss_fndciv; see code
;
; OPTIONAL OUTPUTS:
;  /tcpu -- print information about CPU time used
;
; COMMENTS:
;
; EXAMPLES:
;
;PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;   30-May-2011 Major revamp to be generic and run with new
;               structures, KLC
;   08-Jun-2011 Added processor=, KLC
;   22-Nov-2011 Fix when use conti-including-errors, KLC
;   02-Dec-2011 Force gwidth to be snr_conv at zabs, KLC
;   25-Jan-2013 Fix bug where-by sightlines with zqso > zmax exluded, KLC
;-
;------------------------------------------------------------------------------
@sdss_fndlin                    ; resolve sdss_fndlin_calcew()

pro sdss_civsearch_setgwidth, civ_fil, out_fil, clobber=clobber, final=final
  ;; New feature to let gwidth tag in sdsscivstrct (originally for the
  ;; Gaussian sigma from an EW fit) be the convolved S/N at zabs_orig,
  ;; where the convolved S/N is measured in sdss_fndlin (so stored in
  ;; the abslin file)
  if n_params() ne 2 then begin
  endif 

  test = file_search(out_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_civsearch_setgwidth: exiting; file exists and will not clobber ',$
           out_fil
     return
  endif 
  sdssdir = sdss_getsdssdir()

  if size(civ_fil,/type) eq 7 then civstr = xmrdfits(civ_fil,1,/silent) $
  else civstr = civ_fil
  nciv = (size(civstr,/dim))[0] > 1
  tags = tag_names(civstr[0])
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
  endelse 

  ;; Avoid reading in too much
  prev_absfil = ''
  for iciv=0L,nciv-1 do begin
     if civstr[iciv].abslin_fil ne prev_absfil then begin
        ;; Read spectrum (need wavelength)
        parse_sdss,sdssdir+civstr[iciv].sdss_obs[0],flux,wave
        ;; Read conti
        cstrct = xmrdfits(sdssdir+civstr[iciv].abslin_fil,1,/silent)
        cindx = fix(alog(civstr[iciv].cflg)/alog(2)) ; assume all same cflg in LOS
     endif ; read in

     for gg=0,1 do begin
        ;; Loop over doublet
        mn = min(wave-civstr[iciv].wrest[gg]*(1+civstr[iciv].(ztag)[gg]),$
                 imn,/abs)
        civstr[iciv].gwidth[gg] = cstrct.snr_conv[imn,cindx]
     endfor                     ; loop gg=0,1
     
     ;; print progress
     if iciv ne 0 and (iciv mod 500) eq 0 then $
        print,'sdss_civsearch_setgwidth: iciv = ',iciv

  endfor                        ; loop iciv=nciv

  ;; Write out (clobber test made at beginning)
  mwrfits,civstr,out_fil,/create,/silent
  spawn,'gzip -f '+out_fil
  print,'sdss_civsearch_setgwidth: created ',out_fil

end                             ; sdss_civsearch_setgwidth


function sdss_civsearch_broad, contistrct_fil, dblt_name, count=count, $
                               set_ew=set_ew, wave=wave, flux=flux, $
                               sigma=sigma, silent=silent, zmin=zmin, $
                               zmax=zmax, dvgal=dvgal, dvqso=dvqso, $
                               dvbal=dvbal, flg_bal=flg_bal, $
                               civstr_tmplt=civstr_tmplt, _extra=extra
  ;; Find broad lines that should be candidates
  if n_params() ne 2 then begin
     print,'Syntax - sdss_civsearch_broad(contistrct_fil, dblt_name, [/set_ew,'
     print,'                              count=, wave=, flux=, sigma=, /silent,'
     print,'                              zmin=, zmax=, dvgal=, dvqso=,'
     print,'                              dvbal=, /flg_bal, _extra=])'
     return,-1
  endif
  sdssdir = sdss_getsdssdir()
  if not keyword_set(dvgal) then dvgal = 5000.   ; km/s; outside Galaxy
  if not keyword_set(dvqso) then dvqso = -3000.  ; km/s; relative to QSO
  if not keyword_set(dvbal) then dvbal = -10000. ; km/s; relative to QSO if BAL

  ;; Can save time if using pre-instantiated structure
  if not keyword_set(civstr_tmplt) then civstr_tmplt = {sdsscivstrct}

  ;; Read in params
  if size(contistrct_fil,/type) eq 8 then cstrct = contistrct_fil $ ; input is structure
  else cstrct = xmrdfits(sdssdir+contistrct_fil,1,/silent)                  
  
  ;; cflg better be right; would be affected by use_cflg= in _extra=
  ;; for sdss_fndlin_calcew()
  cindx = fix(alog(cstrct.cflg)/alog(2))

  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  
  ;; Figure out bounds and don't use defaults in sdss_measuresnr()
  tmp = sdss_measuresnr('blank',wvlim_obs=wvlim_obs,/no_snr,dvgal=dvgal,$
                        dvqso=dvqso, $
                        dblt_name=dblt, zqso=cstrct.z_qso)
  ;; Avoid sky lines if at all possible; if no candidates, don't go there
  if max(cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx]) LE 8000. then $
     wvlim_obs[1] = 8000. < wvlim_obs[1]

  if not keyword_set( zmin ) then $
     zmin = wvlim_obs[0]/dblt.wvI - 1.
  if not keyword_set( zmax ) then begin
     ;; Set by SDSS spectral coverage, emission line of QSO (modulo
     ;; dvqso)
     if keyword_set(flg_bal) and not keyword_set(silent) then $
           print,'sdss_civsearch_broad() debug: BAL flagged and dvbal adopted'
     zmax = wvlim_obs[1]/dblt.wvII - 1.
  endif 
  if not keyword_set(silent) then $
     print,'sdss_civsearch_broad(): z limits ',$
           string(zmin,zmax,format='(f8.5,"--",f8.5)')

  wvc_rest = 0.5*(dblt.wvI + dblt.wvII)
  dvdblt_rest = 2.998e5 * (dblt.wvII-dblt.wvI) / wvc_rest
  
  ;; Perhaps need to set the bounds (but lots of things need to be set
  ;; up correctly do to this, so better if didn't have to; or
  ;; should do check this is the case)
  ;; _extra includes things for sdss_fndlin_calcew():
  ;; use_cflg=, _extra [/keepwvlim, /debug, /plot]
  if keyword_set(set_ew) then $
     cstrct = sdss_fndlin_calcew(cstrct,wave=wave,flux=flux,sigma=sigma,$
                                 _extra=extra) ; will call *calcnormerr()
  
  wvc_lin =  0.5 * (cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,0] + $
                    cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,1])
  dvlin = 2.998e5 * (cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,1] - $
                     cstrct.wvlim_orig[0:cstrct.ncent[cindx]-1,0]) / $
          wvc_lin
  
  ;; Check
  cand = where(dvlin ge 1.5*dvdblt_rest,count) ; magic number = fudge factor
  if count eq 0 then return,civstr_tmplt  ; EXIT point

  ;; Instantiate appropriate structure
  civstr = replicate(civstr_tmplt,count)
;  civstr.gwidth[0:1] = cstrct.snr_conv[cstrct.npix,cindx] ; LSNR from sdss_fndlin
  civstr.rating[9] = sdss_getblendflg(/self)              ; self-blend

  mask = replicate(1,count)     ; 1 = good/keep, 0 = bad/reject

  ;; Read in spectrum
  if not keyword_set(flux) then $
     parse_sdss, sdssdir+cstrct.sdss_obs[0], flux, wave, npix=npix, sig=sigma $
  else begin
     npix = (size(flux,/dim))[0]
     if not keyword_set(wave) or not keyword_set(sigma) then $
        stop,'sdss_civsearch_broad(): must set wave and sigma arrays with flux'

     if npix ne cstrct.npix or npix ne (size(wave,/dim))[0] or $
        npix ne (size(sigma,/dim))[0] then $
           stop,'sdss_civsearch_broad(): spectrum arrays not same size'
  endelse 

  ;; Set up continuum; it's errors should already be included!!!!
  conti = cstrct.conti[0:cstrct.npix-1,cindx]
  sig = sdss_calcnormerr(flux,sigma,cstrct,/unnorm)


  ;; Basic info
  civstr.qso_name = cstrct.qso_name
  civstr.mjd = cstrct.mjd   ;  can't pull from fil_sdss header
  civstr.plate = cstrct.plate
  civstr.fiber = cstrct.fiber
  civstr.z_qso = cstrct.z_qso
  civstr.cflg = cstrct.cflg
  civstr.balflg = cstrct.balflg     
  civstr.sdss_obs = cstrct.sdss_obs 
  if size(contistrct_fil,/type) eq 7 then $
     civstr.abslin_fil = contistrct_fil $
  else begin
     civstr.abslin_fil = sdss_getname(cstrct.sdss_obs[0],/spec,/abslin,$
                                      dir=cdir)
     civstr.abslin_fil = cdir[0]+civstr.abslin_fil
  endelse 
  civstr.zabs_orig[0] = wvc_lin[cand]/wvc_rest - 1. ; estimate

  for ii=0L,count-1 do begin
     ;; Set doublet info, as best as possible
     civstr[ii] = sdss_ewciv(wave[cstrct.ipix0:cstrct.npix-1], $
                             flux[cstrct.ipix0:cstrct.npix-1], $
                             sig[cstrct.ipix0:cstrct.npix-1], $
                             conti[cstrct.ipix0:cstrct.npix-1], $
                             dblt, civstr[ii].zabs_orig[0], istrct=civstr[ii], $
                             snr_conv=cstrct.snr_conv[cstrct.ipix0:cstrct.npix-1,cindx],$
                             _extra=extra)

     if civstr[ii].zabs_orig[0] lt zmin or $
        civstr[ii].zabs_orig[0] gt zmax then begin
        mask[ii] = 0
        continue
     endif 
     
     for gg=0,1 do begin
        ;; Set meta data
        mn = min(civstr[ii].wrest[gg]*(1+civstr[ii].zabs_orig[gg])-wave,imn,/abs)
        civstr[ii].gwidth[gg] = cstrct.snr_conv[imn,cindx]
     endfor                     ; loop gg=0,1

     ;; Since many lines may make the initial velocity cut, have to do
     ;; a sanity check because lots of narrow lines are getting
     ;; chopped up. 
     ;; First, inspect EW because don't accept junk here; true
     ;; broad guys should result in positive EW in both lines
     if civstr[ii].ew_orig[0] lt 0. or civstr[ii].ew_orig[1] lt 0. then $
        mask[ii] = 0 $
     else begin
        ;; Do shape checks, namely how the flux looks around key portions
        mn = min(wave-civstr[ii].wvlim_orig[0,1],ihiI,/abs)
        mn = min(wave-civstr[ii].wvlim_orig[1,0],iloII,/abs)
        mn = min(wave-dblt.wvI*(1.+civstr[ii].zabs_orig[0]),iwvI,/abs)
        mn = min(wave-dblt.wvII*(1.+civstr[ii].zabs_orig[1]),iwvII,/abs)
        idx = [-1,0,1]          ; going to take mean of range
        mean_fxI = mean(flux[iwvI+idx])
        mean_fxII = mean(flux[iwvII+idx])
        ;; If the mid-point is lower than either other flux-weighted centroid,
        ;; then it's one line that's been divided
        if mean(flux[ihiI+idx]) lt mean_fxI and $
           mean(flux[iloII+idx]) lt mean_fxII then $
              mask[ii] = 0 $
        else begin
           ;; Check strength at center of each bounds, secondary check on above
           mn = min(wave-mean(civstr[ii].wvlim_orig[0,*]),icI,/abs)
           mn = min(wave-mean(civstr[ii].wvlim_orig[1,*]),icII,/abs)
           if mean(flux[icI+idx]) gt mean_fxI and $
              mean(flux[icII+idx]) gt mean_fxII then $
                 mask[ii] = 0 $
           else begin
              if (abs(ihiI-iwvI) le 1 and abs(iloII-iwvII) le 1) or $
                 (ihiI eq iwvI or iloII eq iwvII) then $
                    mask[ii] = 0 $
              else begin
                 ;; Don't want min flux in range to be at boundary
                 mn = min(wave-civstr[ii].wvlim_orig[0,0],iloI,/abs)
                 mn = min(wave-civstr[ii].wvlim_orig[1,1],ihiII,/abs)
                 mn = min(flux[iloI:ihiII],imn)
                 if abs(imn-ihiI) le 1 or abs(imn-iloII) le 1 then $
                    mask[ii] = 0
              endelse           ; absolute min
           endelse              ; center flux > mean fxI/II

           ;; If well-aligned in boxes, keep (actually pretty good for
           ;; picking off MgII); undoes ill effects of previous
           if abs(icI-iwvI) le 1 and abs(icII-iwvII) le 1 then $
              mask[ii] = 1

        endelse                 ; boundary flux min
        
     endelse                    ; EW test
  endfor                        ; loop ii=count
                    
  
  gd = where(mask eq 1,count)
  if count eq 0 then begin
     if not keyword_set(silent) then $
        print,'sdss_civsearch_broad(): no broad candidates found'
     return,civstr_tmplt        ; EXIT point
  endif else begin
     civstr = civstr[gd]
     if not keyword_set(silent) then begin
        print,'sdss_civsearch_broad(): z of candidates found in broad lines:'
        print,civstr.zabs_orig[0]
     endif 
     return, civstr
  endelse 

end                             ; sdss_civsearch_broad()



function sdss_civsearch_srch, wave, flux, sigma, dblt_name, abslin_fil, $
                              cflg=cflg, qstrct_tmplt=qstrct_tmplt, $
                              count=count, broad_count=broad_count, $
                              broadcivstr=broadcivstr, civstr_tmplt=civstr_tmplt,$
                              excl_broad=excl_broad, debug=debug, silent=silent,$
                              _extra=extra
  if n_params() ne 5 then begin
     print,'Syntax - sdss_civsearch_srch(wave, flux, sigma, dblt_name, abslin_fil,'
     print,'                             [cflg=, qstrct_tmplt=, count=, '
     print,'                             count_broad=, broadcivstr=, '
     print,'                             civstr_tmplt=, /excl_broad, /debug, '
     print,'                             /silent, _extra=])'
     return,-1
  endif 
  sdssdir = sdss_getsdssdir()
  count = 0
  broad_count = 0

  if not keyword_set(qstrct_tmplt) then $
     qstrct_tmplt = create_struct({qalcharstrct},'zlim2',dblarr(100,2))
  if not keyword_set(civstr_tmplt) then $
     civstr_tmplt = {sdsscivstrct}
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)


  ;; Search 
  ;; Redshift limits will make a more efficient search
  ;; sdss_fncdiv: _EXTRA include dvtol=, dvgal=, dvqso=, dz=
  ;; read test b/c has directory and *.gz
  if size(abslin_fil,/type) eq 8 then cstrct = abslin_fil $
  else begin
     if keyword_set(debug) then begin
        test = file_search(sdssdir+abslin_fil+'*',count=ntest) ; x_chkfil() slow
        if ntest EQ 0 then $
           stop,'sdss_civsearch_srch(): abslin file DNE ',abslin_fil
     endif 
     cstrct = xmrdfits(sdssdir+abslin_fil,1,/silent) ; candidate lines
  endelse 

  if keyword_set(cflg) then cstrct.cflg = cflg    ; overwrite
  cindx = fix(alog(cstrct.cflg)/alog(2))
  if cstrct.ncent[cindx] eq 0 then begin
     if not keyword_set(silent) then begin
        if size(abslin_fil,/type) eq 7 then $
           print,'sdss_civsearch_srch(): skipping; no lines ',abslin_fil $
        else print,'sdss_civsearch_srch(): skipping; no lines ',$
                   cstrct.qso_name,cstrct.z_qso
     endif
     if keyword_set(debug) then print,''

     return, civstr_tmplt       ;  EXIT
  endif 
  
  ;; Use conti+spectrum errors but sdss_ewciv wants unnormalized
  ;; error (though it's OK if it just multiplies/divides by
  ;; conti as it deems right)
  conti = cstrct.conti[0:cstrct.npix-1,cindx] 
  sig = sdss_calcnormerr(flux,sigma,cstrct,/unnorm) 
  if keyword_set(debug) then begin
     if cstrct.cflg eq sdss_getcflg(/hyb)  then $
        print,'sdss_civsearch_srch() debug: using hybrid conti' $
     else begin
        if cstrct.cflg eq sdss_getcflg(/eig)  then $
           print,'sdss_civsearch_srch() debug: using eigen conti' $
        else print,'sdss_civsearch_srch() debug: using spline conti'
     endelse 
  endif
  

  if not keyword_set(excl_broad) then begin
     ;; Check if exceptionally broad lines should be candidates 
     ;; Use z range returned by sdss_Fndciv
     if keyword_set(debug) then begin
        print,''
        print,'sdss_civsearch_srch() debug: '+$
              'Searching for broad lines that can be candidates'
     endif 
     ;; _extra= includes zmin=, zmax=, /set_ew
     civstr = sdss_civsearch_broad(cstrct, dblt, count=nbroad, $
                                   wave=wave, flux=flux, sigma=sigma, $ ; calls *calcnormerr()
                                   flg_bal=cstrct.balflg, $
                                   debug=debug, silent=silent, $
                                   civstr_tmplt=civstr_tmplt,_extra=extra)
     if nbroad gt 0 then begin
        ;; Save and prepare for normally-selected guys
        if count eq 0 then allciv = civstr $
        else allciv = [allciv, civstr]
        count = count + nbroad
        ;; Save for secondary structure
        if broad_count eq 0 then broadcivstr = civstr $
        else broadcivstr = [broadcivstr, civstr] 
        broad_count = broad_count + nbroad
        if n_elements(broadcivstr) ne broad_count then stop
     endif 
  endif                         ; search broad lines

  
  ;; For looping purposes must make sure zmn/x_fnd must not change
  ;; and sdss_fndciv should be allowed to change them, just not
  ;; globally
  ;; _extra= includes dvtol=, dvgal=, dvqso=, dvbal=, dz=
  ;;         lsnr2=, gapsz=, /orig
  ;; wave=, flux=, sigma= and cstrct_fil= will only be used if
  ;; lsnr2= passed in through _extra.
  ;; Passing in cstrct should be safe because it's not modified.
  qstrct_tmplt.snr = cstrct.snr_conv[cstrct.npix,cindx] ; default for WVII search
  qstrct = sdss_fndciv(cstrct.centroid[0:cstrct.ncent[cindx]-1,cindx], $
                       dblt, cstrct.z_qso, cstrct.cflg, $
                       istrct=qstrct_tmplt, $
                       flg_bal=cstrct.balflg, cstrct_fil=cstrct, $
                       wave=wave, flux=flux/conti, sigma=sig/conti, $
                       silent=silent, debug=debug, _extra=extra)

  
  if qstrct.nDLA2 gt 0 then begin
     ;; Further measurements
     ;; Fill up EW in new {sdsscivstrct}
     civstr = replicate(civstr_tmplt,qstrct.ndla2)
     ;; Copy over some information; civstr is instantiated
     ;; with: zabs_orig[0:1], wrest, ew_orig, sigew_orig, ewflg,
     ;; ncolm_orig, signcolm_orig, ncolmflg, wvlim_orig
     civstr.qso_name = cstrct.qso_name
     civstr.mjd = cstrct.mjd
     civstr.plate = cstrct.plate
     civstr.fiber = cstrct.fiber
     civstr.ra = cstrct.ra
     civstr.dec = cstrct.dec
     civstr.Rmag = cstrct.Rmag 
     civstr.z_qso = cstrct.z_qso
     civstr.cflg = cstrct.cflg
     civstr.balflg = cstrct.balflg     
     civstr.sdss_obs = cstrct.sdss_obs 
     if size(abslin_fil,/type) eq 7 then civstr.abslin_fil = abslin_fil
     civstr.deltaz_cand = qstrct.DLA_z[0] ; actually path
;     civstr.gwidth[0] = cstrct.snr_conv[cstrct.npix,cindx] ; default
;     civstr.gwidth[1] = qstrct.snr                         ; LSNR2 if used

     if civstr[0].deltaz_cand eq 0 then $
        stop,'sdss_civsearch_srch(): no path length?!'
     
     for kk=0L,qstrct.ndla2-1 do begin
        ;; sdss_ewciv automatically measures for both elements
        ;; of doublet assuming same redshift and some width
        ;; based on the doublet seperation (which will be revised)
        civstr[kk] = sdss_ewciv(wave[cstrct.ipix0:cstrct.npix-1], $
                                flux[cstrct.ipix0:cstrct.npix-1], $
                                sig[cstrct.ipix0:cstrct.npix-1], $
                                conti[cstrct.ipix0:cstrct.npix-1], $
                                dblt, qstrct.dla_zabs2[kk], $
                                zlim=qstrct.zlim2[kk,*],$
                                snr_conv=cstrct.snr_conv[cstrct.ipix0:cstrct.npix-1,cindx],$
                                istrct=civstr[kk], debug=debug)
        for gg=0,1 do begin
           ;; Set meta data
           mn = min(civstr[kk].wrest[gg]*(1+civstr[kk].zabs_orig[gg])-wave,imn,/abs)
           civstr[kk].gwidth[gg] = cstrct.snr_conv[imn,cindx]
        endfor                  ; loop gg=0,1
     endfor                     ; loop kk=qstrct.ndla2

     if count eq 0 then allciv = civstr $
     else allciv = [allciv, civstr]
     count = count + qstrct.ndla2
        
  endif else $                  ; strct.ndla2 ne 0
     if keyword_set(debug) then $
        print,'sdss_civsearch_srch() debug: no candidate lines found'

;  if keyword_set(debug) then $
;     print,'sdss_civsearch_srch() debug: CIV count = '+strtrim(count,2)+$
;           '; broad count = '+strtrim(broad_count,2)

  if count gt 0 then return, allciv $
  else return, civstr_tmplt
  
end                             ; sdss_civsearch_srch()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_civsearch, filename, sdsssum, cand_fil, zmin=zmin, $
                    zmax=zmax, dblt_name=dblt_name, rmax=rmax, $
                    debug=debug, tcpu=tcpu, silent=silent, $
                    processor=processor, _extra=extra

  if  N_params() LT 3 then begin 
     print,'Syntax - ' + $
           'sdss_civsearch, filename, sdsssum, cand_fil, [zmin=, zmax=, dblt_name='
     print,'            rmax=, /debug, /tcpu, /silent,'
     print,'            processor=, _extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  
  ;; Defaults
  if not keyword_set( dblt_name ) then dblt_name = 'CIV'
  dblt = dblt_retrieve(dblt_name) 
  wvspecmnx = sdss_getspecwave()
  if not keyword_set( ZMAX ) then begin
     zmax = wvspecmnx[1]/dblt.wvII - 1.
     zmx_fnd = 0                ; let sdss_fndciv do it's own cut
  endif else zmx_fnd = zmax     ; let user set sdss_fndciv
  if not keyword_set( ZMIN ) then begin
     zmin = wvspecmnx[0]/dblt.wvI - 1.
     zmn_fnd = 0
  endif else zmn_fnd = zmin
  if keyword_set(debug) then $
     print,'sdss_civsearch: z limits: ',string(zmin,zmax,format='(f8.5,"--",f8.5)')
  if not keyword_set( RMAX ) then rmax = 23.0


  ;; List
  readcol, filename, fil_spec, FORMAT='A', /silent
  nfil = n_elements(fil_spec)
  absdir = fil_spec[0]
  fil_spec = fil_spec[1:nfil-1]
  nfil = nfil-1

  ;; Summary table
  sdsstab = xmrdfits(sdsssum, 1, /silent)
  if nfil NE n_elements(sdsstab) then $
     stop,'sdss_civsearch: list and QSO structure must be same size'

  ;; Truncate
  ;; Increase efficiency by not even looping over QSOs that
  ;; don't stand a chance.
  ;; However, sdss_fndciv actually allows zabs to be > zem.
  ;; 25 Jan 2013 --- oh, this looks like a terrible bug... cuts off
  ;;                 high redshift systems, KLC
  gd = where(sdsstab.z GE zmin AND sdsstab.rtmag LT Rmax, nfil); AND $
;             sdsstab.z LE zmax, nfil) 

  if nfil EQ 0 then begin
     print,'sdss_civsearch: no QSOs with magnitude and redshift necessary; exiting' 
     return
  endif else begin
     fil_spec = fil_spec[gd]
     sdsstab = sdsstab[gd]
     if not keyword_set(silent) then $
        print,'sdss_civsearch: Number of QSOs in z cut: ',nfil
  endelse

  ;; generate abslin names
  fil_abslin = sdss_getname(sdsstab,/abslin,root=qso_name,dir=contidir)
  test_fil = sdss_getname(fil_spec,/spec,/abslin)
  bd = where(fil_abslin ne test_fil)
  if bd[0] ne -1 then $
     stop,'sdss_civsearch: list and QSO structure not in same order.'
  contidir = absdir + strmid(contidir,strpos(contidir[0],'/')+1)
  fil_abslin = contidir + fil_abslin

  ;;  Make the structure (and want to pass back boundaries
  qstrct_tmplt = create_struct({qalcharstrct},'zlim2',dblarr(100,2))
  qstrct_tmplt.nDLA2 = -1           ; flags hasn't been processed

  ;; more effcient to load this once and copy it around
  civstr_def = {sdsscivstrct}
  
  flg_civ = 0
  flg_broadciv = 0
  istrt = 0L
  if keyword_set(processor) then begin
     sub = sdss_calcparalleljob(sdsstab, processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     outfil = cand_fil + '.' + strtrim(processor[0],2)
     if keyword_set(debug) then $
        print,'sdss_civsearch debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif else $
     outfil = cand_fil
  
  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  ;; LOOP
  for qq=istrt, nfil-1 do begin
     
     ;; instantiate some information but only when actually looped on
     if not keyword_set(silent) then $
        print, 'sdss_civsearch: Analysing file qq ',qq,' ', fil_spec[qq]
     
     ;; Grab names
     parse_sdss, sdssdir+fil_spec[qq], flux, wave, sig=sigma, head=head

     ;; Structure
     ;; _extra= includes cflg=, excl_broad=, /orig
     zmn0 = zmn_fnd
     zmx0 = zmx_fnd

     civstr = sdss_civsearch_srch(wave, flux, sigma, dblt, fil_abslin[qq],$ 
                                  debug=debug, silent=silent, $
                                  qstrct_tmplt=qstrct_tmplt, $
                                  civstr_tmplt=civstr_def, $
                                  zmin=zmn_fnd, zmax=zmx_fnd, $
                                  count=count, broad_count=broad_count, $
                                  broadcivstr=broadcivstr, $
                                  _extra=extra)

     ;; Reset
     zmn_fnd = zmn0
     zmx_fnd = zmx0
     
     if count ne 0 then begin
        if flg_civ eq 0 then begin
           allciv = civstr 
           flg_civ = 1          ; now check for duplicates
        endif else allciv = [allciv, civstr]
;        if keyword_set(debug) then $
;           print,'sdss_civsearch debug: CIV count = '+strtrim(n_elements(allciv),2)
     endif 
     if broad_count ne 0 then begin
        if flg_broadciv eq 0 then begin
           allbroadciv = broadcivstr 
           flg_broadciv = 1     ; now check for duplicates
        endif else allbroadciv = [allbroadciv, broadcivstr]
        if keyword_set(debug) then $
           print,'sdss_civsearch debug: broad CIV count = '+$
                 n_elements(allroadciv)
     endif 

     ;; Output (backup)
     if ((qq + 1) MOD 1000L) EQ 0 and flg_civ eq 1 then begin
        mwrfits, allciv, outfil, /create, /silent
        if flg_broadciv eq 1 then $
           mwrfits, allbroadciv, outfil, /silent ; ext=2
        spawn,'gzip -f '+outfil
        if keyword_set(tcpu) and not keyword_set(silent) then begin
           tt = systime(/seconds)
           dt = (tt-tlast)      ; seconds
           print,'ssss_civsearch: Elapsed time (m) = ',dt/60.
           print,'ssss_civsearch: Average time per qq (s) = ',dt/1000L
           tlast = tt
        endif                   ; /tcpu output
     endif                      ; backup file

     if keyword_set(debug) then print,'' ; make space
  endfor                                 ; loop qq=istrt,nfil
  
  ;; Final Output
  ;; Sort by QSO because makes sdss_chkciv run faster
  if flg_civ eq 1 then begin
     allciv = sdss_srtcivstrct(allciv)
     mwrfits, allciv, outfil, /create, /silent 
  endif else begin
     ;; print dummy "structure"
     mwrfits, -1, outfil, /create, /silent
     print,'sdss_civsearch: no candidates found'
  endelse 
  if flg_broadciv eq 1 then begin
     allbroadciv = sdss_srtcivstrct(allbroadciv)
     mwrfits, allbroadciv, outfil, /silent ; ext=2
  endif else $
     ;; print dummy "structure"
     mwrfits, -1, outfil, /silent ; ext=2
  spawn,'gzip -f '+outfil
  
  print, "sdss_civsearch: created ", outfil

  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_civsearch: Elapsed time for '+strtrim(nfil-istrt,2)+$
        ' QSOs (m) = ',dt/60.
  print,'ssss_civsearch: Average time per QSO (s) = ',dt/(nfil-istrt)

  print,''
  print,'sdss_civsearch: Summary'
  print,'                z limits: ',string(zmin,zmax,format='(f8.5,"--",f8.5)')
  print,'                Number of QSOs in z cut: ',nfil
  if flg_civ eq 1 then $
     print,'                Number of '+dblt.ion+' candidates found: ',$
           strtrim(n_elements(allciv),2)
  if flg_broadciv eq 1 then $
     print,'                Number of '+dblt.ion+' broad candidates found: ',$
           strtrim(n_elements(allbroadciv),2)

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end
