;+ 
; NAME:
; sdss_completeness
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
; INPUTS: 
;
; RETURNS: 
;
; OUTPUTS: 
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: 
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21 Sep 2011  Created by KLC
;-
;------------------------------------------------------------------------------
@sdss_fndlin                    ; resolve sdss_fndlin_srch()
@sdss_dblfitconti               ; resolve sdss_dblfitconti_fithybrid()
@sdss_civsearch                 ; resolve sdss_civsearch_srch()
@sdss_genprof                   ; resolve loads of things

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_cpciv2mcstrct, mcstrct, mcindx, civstrct, civindx
  ;; Copy over pertinent values from civstrct to mcstrct in the
  ;; locations specified
  if n_params() ne 4 then begin
     print,'sdss_completeness_cpciv2mcstrct, mcstrct, mcindx, civstrct, civindx'
     return
  endif 

  ;; Copy over values
  if civstrct[civindx].gwidth[1] eq 0. then stop             ; check
  mcstrct.lsnr[mcindx,0:1] = civstrct[civindx].gwidth[0:1]    ; wvI, wvII
  
  mcstrct.wvlim_rec[mcindx,0,*] = civstrct[civindx].wvlim_orig[0,*]    ; wvI
  mcstrct.wvlim_rec[mcindx,1,*] = civstrct[civindx].wvlim_orig[1,*]    ; wvII
  
  mcstrct.zabs_rec[mcindx,*] = civstrct[civindx].zabs_orig[0:1] ; wvI, wvII
  mcstrct.sigzabs_rec[mcindx,*] = civstrct[civindx].sigzabs_orig[0:1] 
  
  mcstrct.ncolm_rec[mcindx,*] = civstrct[civindx].ncolm_orig[0:1] ; wvI, wvII
  mcstrct.signcolm_rec[mcindx,*] = civstrct[civindx].signcolm_orig[0:1] 
  mcstrct.ncolmflg_rec[mcindx,*] = civstrct[civindx].ncolmflg[0:1]
  
  mcstrct.ew_rec[mcindx,*] = civstrct[civindx].ew_orig[0:1]          ; wvI, wvII
  mcstrct.sigew_rec[mcindx,*] = civstrct[civindx].sigew_orig[0:1]    ; wvI, wvII
     
end                             ; sdss_completeness_cpciv2mc



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_completeness_compare, mcstrct, civstrct, civmask=civmask, $
                                    ignore_flg=ignore_flg, allmtch=allmtch, $
                                    srch_both=srch_both, debug=debug
  if n_params() ne 2 then begin
     print,'Syntax -  sdss_completeness_compare(mcstrct, civstrct, [/ignore_flg,'
     print,'                                  civmask=, allmtch=, /srch_both, /debug]'
     return,-1
  endif 
  newmcstrct = mcstrct
  civmask = civstrct.rating[0] * 0  ; 0 = not recovered, >=1 = recovered

  ;; Match on maximal bounds; assume mcstrct is smaller
  ;; this *MUST* be only one mcstrct structure
  wvobs = mcstrct.wrest[*,0]*(1+mcstrct.zabs_input[*,0])
  allmtch = transpose(civstrct.rating[0:1]) * 0 - 1 ; [ndim,2]

  if keyword_set(srch_both) then $
     wvobs_civ = civstrct.wrest[0]*(1+civstrct.zabs_orig[0])

  for ss=0,mcstrct.nsys-1 do begin
     if keyword_set(srch_both) then $
        mtch = where(mcstrct.qso_name eq civstrct.qso_name and $
                     ((wvobs[ss] ge civstrct.wvlim_orig[0,0] and $
                       wvobs[ss] le civstrct.wvlim_orig[0,1]) or $
                      (mcstrct.wvlim_input[ss,0,0] le wvobs_civ and $
                       mcstrct.wvlim_input[ss,0,1] ge wvobs_civ)),nmtch) $
     else $
        mtch = where(mcstrct.qso_name eq civstrct.qso_name and $
                     wvobs[ss] ge civstrct.wvlim_orig[0,0] and $
                     wvobs[ss] le civstrct.wvlim_orig[0,1],nmtch)
     if nmtch eq 0 then begin
        if keyword_set(debug) then $
           print,'sdss_completeness_compare() debug: ',$
                 'input doublet not recovered ',mcstrct.qso_name,$
                 string(mcstrct.zabs_input[ss,0],mcstrct.ew_input[ss,0,0],$
                       format='(f7.5,2x,f6.2)')
        continue                 ; everything left blank
     endif 

     allmtch[mtch,0] = mtch
     allmtch[mtch,1] = ss

     if nmtch gt 1 then begin
        ;; Take closest
        mn = min(mcstrct.zabs_input[ss,0]-civstrct[mtch].zabs_orig[0],imn)
        mtch = mtch[imn]
     endif else mtch = mtch[0]
     civmask[mtch]++            ; matched

     if newmcstrct.flg_rec[ss] eq 1 and not keyword_set(ignore_flg) then $
        stop,'sdss_completeness_compare() stop: double matching mcstrct system'
     newmcstrct.flg_rec[ss] = 1 ; recovered

     ;; Copy values over
     sdss_completeness_cpciv2mcstrct, newmcstrct, ss, civstrct, mtch

  endfor                        ; loop ss=mcstrct.nsys

  if keyword_set(debug) and not keyword_set(notrec) then $
     print,'sdss_completeness_compare() debug: Input and recovered doublets',$
           newmcstrct.nsys,total(newmcstrct.flg_rec)
  
  return,newmcstrct
end                             ; sdss_completeness_compare()



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_los, dblt_name, spec_fil, abslin_fil, $
                           config_fil, nabstot, mcstrct_fil, mchdr=mchdr,$
                           vpstrct_fil=vpstrct_fil, mcstrct_tmplt=mcstrct_tmplt,$
                           eigbasis=eigbasis, pca_fil=pca_fil, $
                           pca_head=pca_head, debug=debug, clobber=clobber,$
                           seed=seed, oseed=oseed, quick=quick, user=user,$
                           _extra=extra
  if n_params() ne 6 then begin
     print,'Syntax - sdss_completeness_los, dblt_name, spec_fil, abslin_fil, config_fil, nabstot, '
     print,'                           mcstrct_fil, [mchdr=, vpstrct_fil=, mcstrct_tmplt='
     print,'                           eigbasis=, pca_fil=, pca_head=, /debug, /clobbber,'
     print,'                           seed=, oseed=, /quick, /user, _extra=]'
     return
  endif 

  sdssdir = sdss_getsdssdir()
  ;; Real workhorse on a given sightline and parameters

  ;; Params
  if keyword_set(seed) then oseed = seed ; will return to same variable
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  if size(config_fil,/type) eq 8 then config_strct = config_fil $
  else config_strct = sdss_genprof_config(config_fil,header=mchdr)
  if not keyword_set(mcstrct_tmplt) then $
     mcstrct_tmplt = {sdssmcstrct} 
  if not keyword_set(mchdr) then begin
     ;; Create properly formatted header by brute force
     tmp_fil = 'sdss_completeness_los_mchdr.fit'
     mwrfits,mcstrct_tmplt,tmp_fil,/create,/silent
     mwrfits,config_strct,tmp_fil,/silent
     dum = xmrdfits(tmp_fil,1,mchdr,/silent)
;     spawn,'rm '+tmp_fil
  endif 

  if keyword_set(quick) then begin
     ncleanloop = 10            ; make more diversity in spectra
  endif else begin
     if keyword_set(user) then begin
        ncleanloop = 1 
     endif else begin
        ncleanloop = 50         ; a lot of fits
        
        if not keyword_set(eigbasis) then $ ; better not to read in every time
           eigbasis = xmrdfits(getenv("SDSSPATH")+"/eigenspectra/eigSpec_qso_all.fit",$
                               0, /silent)
        if not keyword_set(pca_fil) then $
           pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
        if size(pca_fil,/type) eq 7 then $
           pca = xmrdfits(pca_fil, 0, pca_head, /silent) $
        else begin  
           pca = pca_fil 
           if not keyword_set(pca_head) then $
              stop,'sdss_completeness_los stop: must set pca_head' 
        endelse 
     endelse                    ; conti
  endelse                       ; user or conti

  ;; Useful numbers
  cflg_eig = sdss_getcflg(/eig)
  cindx_eig = sdss_getcflg(/eig,/index)
  cflg_hyb = sdss_getcflg(/hyb)
  cindx_hyb = sdss_getcflg(/hyb,/index)


  cstrct0 = xmrdfits(sdssdir+abslin_fil,1,/silent)
  if cstrct0.cflg ne cflg_hyb then $
     stop,'sdss_completeness_los stop: completeness test geared toawards hybrid conti only'
  
  if keyword_set(debug) then $
     print,'sdss_completeness_los debug: processing ',cstrct0.qso_name,$
           string(cstrct0.z_qso,format='(f8.5)')

  tmp = sdss_measuresnr('blank', wvlim_obs=wvlim_obs, /no_snr, dblt_name=dblt, $
                        zqso=cstrct0.z_qso) 
  
  ;; Instantiate instantiate LOS info tags:
  ;; QSO_NAME, MJD, PLATE, FIBER, RA, DEC, RMAG, Z_QSO, BALFLG,
  ;; CFLG 
  mcstrct_sub0 = sdss_cpstrct(cstrct0,mcstrct_tmplt)

  ;; Read in spectrum, original abslin structure; prepare to pass around 
  parse_sdss, sdssdir+spec_fil, flux, wave, head=hdr, npix=npix, sig=sigma
  wave_rest = wave/(1.+cstrct0.z_qso)
  mask0 = replicate(1,npix)     ; should never change
  bdpix = where(sigma eq 0.)
  if bdpix[0] ne -1 then mask0[bdpix] = 0


  ;; Figure out median resolution in region sensitive to doublet
  resstrct = xmrdfits(sdssdir+spec_fil,6,/silent) 
  gd = where(wave ge wvlim_obs[0] and wave le wvlim_obs[1], ngd) 
  if size(resstrct,/type) eq 8 then begin
     if ngd eq 1 then $
        config_strct.dvelo = resstrct[gd].dispersion * config_strct.pixscale $
     else $
        config_strct.dvelo = median(resstrct[gd].dispersion) * config_strct.pixscale 
  endif else begin
     config_strct.dvelo = config_strct.pixscale
     print,'sdss_completeness_los: resolution ext=6 DNE ',spec_fil
  endelse 

  ;; ;;;;;;;
  ;; Set the number of absorbers based on wavelength coverage of spectrum
  config_strct.zlim = wvlim_obs/dblt.wvI - 1. 
  config_strct.nabs = round(config_strct.dndz*$
                            (config_strct.zlim[1]-config_strct.zlim[0])) > 1
  pixmin0 = gd[0]
  pixmin = pixmin0 - 10 > cstrct0.ipix0 ; buffer search zone
  dum = min(wave-dblt.wvII*(1+config_strct.zlim[1]), pixmax0, /abs)
  pixmax = pixmax0 + 10 < npix - 1

  ;; Read or generate a whole lot of profiles
  if keyword_set(vpstrct_fil) then begin
     ;; Scramble systemic redshift and order of systems; this is SLOW!!!
     ;; but after allowing nabstot= to truncate things, it's less slow
     if keyword_set(debug) then $
        print,'sdss_completeness_los debug: scrambling redshifts for LOS zlim=',$
              string(config_strct.zlim[0],'--',config_strct.zlim[1],$
                     format='(f8.5,a3,f8.5)')
     vpstrct = sdss_genprof_setnewz(vpstrct_fil,config_strct.zlim,seed=oseed,$
                                    oseed=oseed,/reorder,$
                                    nabstot=(nabstot > config_strct.nabs),$
                                    newindex=index_sys)

     nvpsys = (size(index_sys,/dim))[0]
  endif else begin

     ;; This is very slow!!!
     if keyword_set(debug) then $
        print,'sdss_completeness_los debug: generating '+strtrim(nabstot,2)+$
              ' profiles; this may take awhile...'
     nabs = config_strct.nabs
     config_strct.nabs = nabstot
     vpstrct = sdss_genprof_mc(dblt, config_strct, seed=oseed, oseed=oseed, $
                               debug=debug)
     config_strct.nabs = nabs
     ;; Grep out the systems (id_lin --> wvI)
     index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nvpsys)

     if keyword_set(debug) then $
        print,'sdss_completeness_los debug:                         ... done'
  endelse 

  ;; Set up for big loop 
  uber_count = 0L
  istart = 0L
  loop_count = 0L
  clean_count = 0L
  draw_count = 0L

  ;; Instantiate mcstrct header
  if sxpar(mchdr,'NABSTOT') eq 0 then begin
     sxaddpar,mchdr,'NABSTOT',nabstot,'Number requested for input profiles'
     sxaddpar,mchdr,'NCLEAN',ncleanloop,'Number loops using same cleaned spectrum'
     sxaddpar,mchdr,'ILOOP',loop_count,'Number of total loops passed'
     sxaddpar,mchdr,'IINPUT',uber_count,'Number of total input profiles'
     sxaddpar,mchdr,'ICLEAN',clean_count,'Number of cleaned spectrum versions'
     sxaddpar,mchdr,'IDRAWN',draw_count,'Number of re-drawn profiles'
     sxaddpar,mchdr,'OSEED',0,'Random seed last used'
  endif

  ;; ;;;;;;;
  ;; Loop over redshift
  while uber_count lt nabstot do begin
     ;; Select centroids to scrub based on observed EW; only do once 
     ;; in awhile (every 100 loops below)
     ;; _extra includes ewobs_lim=
     cleanspec = sdss_cleanspec(wave, flux, sigma, cstrct0, config_strct.frac_rm, $
                                header=hdr, seed=oseed, oseed=oseed, debug=0, $
                                cent_mask=cent_mask,_extra=extra) 
     mcstrct_sub0.cent_mask = cent_mask ; removed = 1; kept = 0
     clean_count++

     ;; ;;;;;;;
     ;; Randomly generate many profiles to be input into the same
     ;; cleaned spectrum
     for vv=0, ncleanloop-1 do begin
        mcstrct_sub = mcstrct_sub0 ; preserve info and reset

        ;; Take substrcture based on the fact that the systems in
        ;; vpstrct should be ordered consecutively
        if istart+config_strct.nabs ge nvpsys then begin
           if keyword_set(vpstrct_fil) then begin
              ;; Scramble systemic redshift; this is SLOW!!!
              ;; but less slow now that passing in nabstot= to limit
              ;; looping there but MUST pass in original file to have
              ;; chance of getting new profiles
              if keyword_set(debug) then $
                 print,'sdss_completeness_los debug: scrambling redshifts again'
              vpstrct = sdss_genprof_setnewz(vpstrct_fil,config_strct.zlim,$
                                             /reorder,$
                                             nabstot=(nabstot > config_strct.nabs),$
                                             newindex=index_sys,$
                                             seed=oseed,oseed=oseed)
              nvpsys = (size(index_sys,/dim))[0] 
           endif else begin
              ;; New pool of many profiles; this is very slow!!!
              ;; but at least the redshifts are already in the right range
              ;; really never want to get here
              if keyword_set(debug) then $
                 print,'sdss_completeness_los debug: regenerating '+$
                       strtrim(nabstot,2)+' profiles...'
              nabs = config_strct.nabs
              config_strct.nabs = nabstot
              vpstrct = sdss_genprof_mc(dblt, config_strct, seed=oseed, $
                                        oseed=oseed,debug=debug)
              config_strct.nabs = nabs

              ;; Grep out the systems (id_lin --> wvI)
              index_sys = where(vpstrct.id_comp eq 0 and $
                                vpstrct.id_lin eq 1,nvpsys)
              if keyword_set(debug) then $
                 print,'                                   ... done'
           endelse 

           istart = 0L          ; reset
           draw_count++         ; increment
        endif                   ; out of profiles

        ;; Grab out the right number of systems
        for ss=0,config_strct.nabs-1 do begin
           gd = where(vpstrct.id_sys eq vpstrct[index_sys[istart]].id_sys)
           if ss eq 0 then $
              sub = gd $
           else sub = [sub,gd]
           istart++
        endfor                  ; loop ss=config_strct.nabs
        vpstrct_sub = vpstrct[sub]
        uber_count = uber_count + config_strct.nabs
        
        if keyword_set(debug) then begin
           gd = where(vpstrct_sub.id_comp eq 0 and vpstrct_sub.id_lin eq 1,ngd)
           print,'sdss_completeness_los debug: input profiles:'
           printcol,vpstrct_sub[gd].id_sys,vpstrct_sub[gd].z_sys,$
                    vpstrct_sub[gd].ew_obs,vpstrct_sub[gd].n_sys,$
                    vpstrct_sub[gd].ncomp,$
                    vpstrct_sub[gd].wrest*(1+vpstrct_sub[gd].z_sys),$
                    format='(i10,2x,f7.5,2x,f5.2,2x,f5.2,2x,i2,2x,f9.4)'
           test = uniq(vpstrct_sub[gd].id_sys,sort(vpstrct_sub[gd].id_sys))
           if n_elements(test) ne ngd then $
              stop,'sdss_completeness_los debug: mismatch systems in sub'
        endif 

        
        ;; ;;;;;;;
        ;; Once in a great while throw in some huge absorber near
        ;; systemic of QSO 
        ;; TO BE IMPLEMENTED!!!
        

        ;; ;;;;;;;
        ;; Insert into spectrum 
        ;; Save spectrum in SDSS format [npix, 5] so
        ;; that 0th = flux, 1st = pixel mask (scrubbed lines); 2nd =
        ;; error, 3rd = voigt profile
        newspec = sdss_genprof(wave, dblt, config_strct, cleanspec, $
                               seed=oseed, oseed=oseed, $
                               conti=cstrct0.conti[0:cstrct0.npix-1, cindx_hyb], $
                               vpstrct=vpstrct_sub, debug=0)
        mcstrct_sub.snr = median(newspec[pixmin0:pixmax0, 0]/$
                                 newspec[pixmin0:pixmax0, 2])


        mask = mask0            ; sdss_fndlin_srch() modifies
        if keyword_set(quick) then begin
                                ;: just jump to the end; won't have cstrcte
           cstrcth = cstrct0

           ;; Trim the centroids out of bounds
           gd = where(cstrct0.centroid[*,cindx_hyb] ge pixmin0 and $
                      cstrct0.centroid[*,cindx_hyb] le pixmax0,ngd)

           if keyword_set(debug) then $
              print,'sdss_completeness_los debug: /quick mode; ',$
                    'trim centroids from/to',cstrct0.ncent[cindx_hyb],ngd

           cstrcth.ncent[cindx_hyb] = ngd
           if ngd gt 0 then $
              cstrcth.centroid[0:ngd-1,cindx_hyb] = $
              cstrct0.centroid[gd,cindx_hyb]

        endif else begin

           ;; ;;;;;;;
           ;; Fit eigen-conti and find centroids in just the region
           ;; sensitive to doublet
           if keyword_set(debug) then $
              print,'sdss_completeness_los debug: fiting new eigen conti'
           contihdr = hdr
           neweig = eigqsoconti(wave_rest[cstrct0.ipix0:*], $
                                flux[cstrct0.ipix0:*], $
                                sigma[cstrct0.ipix0:*], eigbasis, $
                                finalmask=finalmask, header=contihdr, $
                                /silent)
           eigconti0 = dblarr(cstrct0.npix, 3, /nozero) ; MUST set values
           eigconti0[cstrct0.ipix0:*, 0] = neweig[*, 0]
           eigconti0[cstrct0.ipix0:*, 1] = finalmask
           eigconti0[cstrct0.ipix0:*, 2] = neweig[*, 1]
           if cstrct0.ipix0 gt 0 then begin
              eigconti0[0:cstrct0.ipix0-1,*] = 0
           endif 
           ;; Update structure and find lines
           cstrcte = cstrct0
           cstrcte.cflg = cflg_eig
           cstrcte.conti[0:cstrct0.npix-1, cindx_eig] = eigconti0[*, 0]
           cstrcte.sigconti[0:cstrct0.npix-1, cindx_eig] = eigconti0[*, 2]


           ;; ;;;;;;;
           ;; FNDLIN #1: has to be for full spectrum for proper mask
           if keyword_set(debug) then $
              print,'sdss_completeness_los debug: finding new centroids for hybrid mask'
           fx = newspec[*, 0] / cstrcte.conti[0:cstrct0.npix-1,cindx_eig]
           sig = sdss_calcnormerr(newspec[*,0],newspec[*,2],cstrcte,$
                                  baderrval=9.e9)
           cstrcte = sdss_fndlin_srch(cstrcte, cstrcte.cflg, wave=wave, $
                                      flux=fx, sigma=sig, mask=mask, $
                                      debug=debug) 


           ;; ;;;;;;;
           ;; Fit hybrid-conti and find centroids in just the region
           ;; sensitive to doublet; mimic the call in sdss_dblfitconti
           if keyword_set(debug) then $
              print,'sdss_completeness_los debug: fiting new hybrid conti'
           hybconti = $
              sdss_dblfitconti_fithybrid(wave, newspec[*, 0], newspec[*, 2], $
                                         mcstrct_sub.snr, dblt, eigconti0, cstrcte, $
                                         ;; Options
                                         /silent, $
                                         ;; Things passed out
                                         eigconti=eigconti, $ ; new one
                                         ;; Things passed in (and deep)
                                         eigbasis=eigbasis, pca_fil=pca, $
                                         pca_head=pca_head, $
                                         contihdr=contihdr) ; to be modified
           ;; Update structure and find lines
           cstrcth = cstrcte
           cstrcth.cflg = cflg_hyb
           cstrcth.conti[0:cstrct0.npix-1, cindx_hyb] = hybconti[*, 0]
           cstrcth.sigconti[0:cstrct0.npix-1, cindx_hyb] = hybconti[*, 2]
        endelse                 ; full blown test with conti fit


        ;; ;;;;;;;
        ;; FNDLIN #2
        if keyword_set(debug) then $
           print,'sdss_completeness_los debug: finding new centroids'
        cstrctf = cstrcth
        cstrctf.ipix0 = 0       
        cstrctf.npix = pixmax - pixmin + 1
        cstrctf.conti[0:cstrctf.npix-1,*] = cstrcth.conti[pixmin:pixmax,*]
        cstrctf.sigconti[0:cstrctf.npix-1,*] = cstrcth.sigconti[pixmin:pixmax,*]
        msk = mask[pixmin:pixmax]
        fx = newspec[pixmin:pixmax, 0] / $
             cstrctf.conti[0:cstrctf.npix-1,cindx_hyb]
        sig = sdss_calcnormerr(newspec[pixmin:pixmax,0],newspec[pixmin:pixmax,2],$
                               cstrctf,baderrval=9.e9)
        cstrctf = sdss_fndlin_srch(cstrctf, cstrctf.cflg, $
                                   wave=wave[pixmin:pixmax], flux=fx, $
                                   sigma=sig, mask=msk, debug=debug) 
        mask[pixmin:pixmax] = msk ; put back
        ;; Put back for storage (if /user)
        cstrcth.ncent[cindx_hyb] = cstrctf.ncent[cindx_hyb]
        cstrcth.centroid[*,cindx_hyb] = cstrctf.centroid[*,cindx_hyb]
        cstrcth.snr_conv[pixmin:pixmax,cindx_hyb] = $ ; just the snippet
           cstrctf.snr_conv[0:cstrctf.npix-1,cindx_hyb]

        cstrcth.wvlim_orig[*,0] = cstrctf.wvlim_orig[*,0]
        cstrcth.wvlim_orig[*,1] = cstrctf.wvlim_orig[*,1]
        cstrcth.ew_orig = cstrctf.ew_orig
        cstrcth.sigew_orig = cstrctf.sigew_orig


        ;; ;;;;;;;
        ;; Find candidate doublets; mimic call in sdss_civsearch
        ;; Make sure cflg set appropriately before
        ;; sdss_civsearch_srch() calls sdss_calcnormerr()
        if keyword_set(debug) then $
           print,'sdss_completeness_los debug: finding new candidate doublets'
        civstr = sdss_civsearch_srch(wave[pixmin:pixmax], $
                                     newspec[pixmin:pixmax, 0], $
                                     newspec[pixmin:pixmax, 2], dblt, cstrctf, $
                                     debug=0, /silent, $
                                     qstrct_tmplt=qstrct_tmplt, $
                                     civstr_tmplt=civstr_tmplt, $
                                     /set_ew, count=count, _extra=extra) 

        if keyword_set(user) then begin
           ;; ;;;;;;;
           ;; Need to save the new spectrum, write its name to the
           ;; candidate structure, save the continuum structure and
           ;; the candidate structure
           ;; Put them in the same dir (e.g., mcfkiv/)
           tmp = sdss_getname(spec_fil,/spec,user=loop_count+1,mc=dblt.ion,dir=nwdir,$
                              quick=quick)
           newspec_fil = sdss_getname(spec_fil,/spec,user=loop_count+1) ; avoid 0
           newspec_fil = nwdir+newspec_fil
           test = file_search(newspec_fil+'*',count=ntest)
           if ntest eq 0 or keyword_set(clobber) then begin
              mwrfits,newspec,sdssdir+newspec_fil,hdr,/create,/silent
              spawn,'gzip -f '+sdssdir+newspec_fil
              if keyword_set(debug) then $
                 print,'sdss_completeness_los debug: created ',newspec_fil
           endif else begin
              stop,'sdss_completeness_los stop: file exists ',newspec_fil
           endelse

           nwcstrct_fil = sdss_getname(spec_fil,/spec,user=loop_count+1,$ ; avoid 0
                                       /abslin)
           nwcstrct_fil = nwdir+nwcstrct_fil ; put in 
           test = file_search(nwcstrct_fil+'*',count=ntest)
           if ntest eq 0 or keyword_set(clobber) then begin
              mwrfits,cstrcth,sdssdir+nwcstrct_fil,/create,/silent ; copied back info
              spawn,'gzip -f '+sdssdir+nwcstrct_fil
              if keyword_set(debug) then $
                 print,'sdss_completeness_los debug: created ',nwcstrct_fil
           endif else begin
              stop,'sdss_completeness_los stop: file exists ',nwcstrct_fil
           endelse

           if count gt 0 then begin
              ;; File names
              civstr.sdss_obs[0] = newspec_fil
              civstr.abslin_fil = nwcstrct_fil
              civstr.cmplt_fil = mcstrct_fil ; for once useful
              
              ;; Save data (not recovered ones just for kicks saved below)
              if keyword_set(civstr_all) then $
                 civstr_all = [civstr_all,civstr] $
              else civstr_all = civstr
           endif                ; count > 0
        endif                   ; /user


        ;; ;;;;;;;
        ;; Compare to input doublets
        ;; sdss_cpvp2mcstrct() instantiates tags: ID_SYS, NSYS, WREST,
        ;; WVLIM_INPUT, ZABS_INPUT, NCOLM_INPUT, B_INPUT, EW_INPUT
        mcstrct_sub = sdss_cpvp2mcstrct(vpstrct_sub,mcstrct_sub) ; copy data
        if count gt 0 then begin

           if keyword_set(debug) then $
              print,'sdss_completeness_los debug: matching input to recovered doublets'
           ;; sdss_completeness_compare() matches and instanties: FLG_REC,
           ;; WVLIM_REC, ZABS_REC, SIGZABS_REC, NCOLM_REC,
           ;; SIGNCOLM_REC, EW_REC, SIGEW_REC
           mcstrct_sub = sdss_completeness_compare(mcstrct_sub,civstr,debug=debug,$
                                                  civmask=civmask)
           if keyword_set(user) then begin
              nciv = (size(civstr_all,/dim))[0]
              civstr_all[nciv-count:*].rating[8] = civmask
           endif 
;           printcol,mcstrct_sub.zabs_input[0,0],mcstrct_sub.ew_input[0,0],$
;                    mcstrct_sub.zabs_rec[0,0],mcstrct_sub.ew_rec[0,0],$
;                    mcstrct_sub.flg_rec[0],$
;                    mcstrct_sub.wrest[*,0]*(1+mcstrct_sub.zabs_input[*,0])
;           printcol,civstr.zabs_orig[0],civstr.wvlim_orig[0,0],$
;                    civstr.wvlim_orig[0,1]


        endif 


        ;; ;;;;;;;
        ;; Measure EW and values for not recovered doublets
        lost = where(mcstrct_sub.flg_rec[0:mcstrct_sub.nsys-1] eq 0,nlost)
        for ll=0,nlost-1 do begin
           if ll eq 0 and keyword_set(debug) then $
              print,'sdss_completeness_los debug: measuring values for ',$
                    'not-recovered doublets'
           ;; sig made from [pixmin:pixmax], cstrctf
           ;; cstrctf is made from cstrcth but having everything
           ;; shifted to [0:npix-1]
           tmpciv = sdss_ewciv(wave[pixmin:pixmax],newspec[pixmin:pixmax, 0], $
                               sig*cstrctf.conti[0:cstrctf.npix-1,cindx_hyb], $
                               cstrctf.conti[0:cstrctf.npix-1,cindx_hyb], $
                               dblt, mcstrct_sub.zabs_input[lost[ll]], $
                               istrct=civstr_tmplt, $
                               snr_conv=cstrctf.snr_conv[0:cstrctf.npix-1,cindx_hyb],$
                               debug=0)
           for gg=0,1 do begin  ; must instantiate properly
              mn = min(tmpciv.wrest[gg]*(1+tmpciv.zabs_orig[gg])-$
                       wave[pixmin:pixmax],imn,/abs)
              tmpciv.gwidth[gg] = cstrctf.snr_conv[imn,cindx_hyb]
           endfor               ; loop gg=0,1
           
           sdss_completeness_cpciv2mcstrct, mcstrct_sub, lost[ll], tmpciv, 0

           if keyword_set(user) then begin
              ;; Must instantiate values and will save
              tmpciv.qso_name = cstrct0.qso_name
              tmpciv.z_qso = cstrct0.z_qso
              tmpciv.sdss_obs[0] = newspec_fil
              tmpciv.abslin_fil = nwcstrct_fil
              tmpciv.cmplt_fil = mcstrct_fil ; for once useful

              if keyword_set(civstr_allrej) then $
                 civstr_allrej = [civstr_allrej,tmpciv] $
              else civstr_allrej = tmpciv
           endif                ; /user
        endfor                  ; loop ll=nlost

        ;; ;;;;;;;
        ;; Save 
        if keyword_set(mcstrct) then mcstrct = [mcstrct, mcstrct_sub] $
        else mcstrct = mcstrct_sub
        
        if keyword_set(debug) then begin
           printcol,mcstrct.zabs_input[0,0],mcstrct.ew_input[0,0],$
                    mcstrct.zabs_rec[0,0],mcstrct.ew_rec[0,0],mcstrct.flg_rec[0]
           stop,'sdss_completeness_los debug stop: end of loop'
        endif

     endfor                     ; loop vv=ncleanloop


     loop_count++ 
     if (loop_count mod 500) eq 0 then begin
        ;; Update header (comments included above) and save
        sxaddpar,mchdr,'ILOOP',loop_count
        sxaddpar,mchdr,'IINPUT',uber_count
        sxaddpar,mchdr,'ICLEAN',clean_count
        sxaddpar,mchdr,'IDRAWN',draw_count
        sxaddpar,mchdr,'OSEED',oseed[0]
        if keyword_set(user) then begin
           ;; Back up recovered candidate structure (ext = 1), not
           ;; recovered candidate structure (ext = 2), MC structure
           ;; (ext = 3), and configuration file (ext = 4)
           if keyword_set(civstr_all) then $
              mwrfits,civstr_all,sdssdir+mcstrct_fil,mchdr,/create,/silent $
           else mwrfits,civstr_tmplt,sdssdir+mcstrct_fil,mchdr,/create,/silent 
           if keyword_set(civstr_allrej) then $
              mwrfits,civstr_allrej,sdssdir+mcstrct_fil,/silent $
           else mwrfits,civstr_tmplt,sdssdir+mcstrct_fil,/silent
           mwrfits,mcstrct,sdssdir+mcstrct_fil,/silent
        endif else $
           mwrfits,mcstrct,sdssdir+mcstrct_fil,mchdr,/create,/silent
        mwrfits,config_strct,sdssdir+mcstrct_fil,/silent
        if keyword_set(debug) then $
           print,'sdss_completeness_los debug: saving intermediate file, loop =',$
                 loop_count
     endif 
  endwhile                      ; uber_count < nabstot

  
  ;; Update header and save
  sxaddpar,mchdr,'ILOOP',loop_count
  sxaddpar,mchdr,'IINPUT',uber_count
  sxaddpar,mchdr,'ICLEAN',clean_count
  sxaddpar,mchdr,'IDRAWN',draw_count
  sxaddpar,mchdr,'OSEED',oseed[0]
  if keyword_set(user) then begin
     ;; Save recovered candidate structure (ext = 1), not recovered
     ;; candidate structure (ext = 2), MC structure (ext = 3), and
     ;; configuration file (ext = 4)
     if keyword_set(civstr_all) then $
        mwrfits,civstr_all,sdssdir+mcstrct_fil,mchdr,/create,/silent $
     else mwrfits,civstr_tmplt,sdssdir+mcstrct_fil,mchdr,/create,/silent 
     if keyword_set(civstr_allrej) then $
        mwrfits,civstr_allrej,sdssdir+mcstrct_fil,/silent $
     else mwrfits,civstr_tmplt,sdssdir+mcstrct_fil,/silent
     mwrfits,mcstrct,sdssdir+mcstrct_fil,/silent
  endif else $
     mwrfits,mcstrct,sdssdir+mcstrct_fil,mchdr,/create,/silent
  mwrfits,config_strct,sdssdir+mcstrct_fil,/silent
  spawn,'gzip -f '+sdssdir+mcstrct_fil
  print,'sdss_completeness_los: created ',mcstrct_fil
  
end                             ; sdss_completeness_los


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_completeness_dz, list_fil, snrstrct_fil, list=list, $
                               ewbinsize=ewbinsize, zbinsize=zbinsize, $
                               ewlim=ewlim, zlim=zilm, dblt_name=dblt_name, $
                               dx_2darr=dx_2darr, rz_2darr=rz_2darr,$
                               debug=debug, silent=silent, _extra=extra
  ;; Compute 2D grid of total unblocked redshift and pathlength
  ;; WARNING!!! This code needs to match what's called in
  ;; sdss_completeness_czw 
  ;; only used when calling fulllist_fil= in sdss_completeness_czw
  if n_params() ne 2 then begin
     print,'Syntax - sdss_completeness_dz( list_fil, snrstrct_fil, [list=, '
     print,'                               ewbinsize=, zbinsize=, '
     print,'                               ewlim=, zlim=, dblt_name=, '
     print,'                               dx_2darr=, rz_2darr=,'
     print,'                               /debug, /silent, _extra=)'
     return,-1
  endif 

  ;; Default parameters
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  if not keyword_set(zbinsize) then zbinsize = 0.05 ;0.005
  if not keyword_set(zlim) then zlim = [1., 6.]
  if not keyword_set(ewbinsize) then ewbinsize = 0.1 ;0.05 ; Ang
  if not keyword_set(ewlim) then ewlim = [0.05, 5.]  ; Ang

  sdssdir = sdss_getsdssdir()
  cosmology = sdss_setcosmology(_extra=extra) ; _extra = cosmology[3]

  ;; Read in data
  if keyword_set(list) then spec_fil = list_fil $
  else readcol,list_fil,spec_fil,format='a',skip=1,/silent
  nfil = (size(spec_fil,/dim))[0]
  
  if size(snrstrct_fil,/type) eq 7 then $
     snrstrct = xmrdfits(snrstrct_fil,1,/silent) $ 
  else snrstrct = snrstrct_fil

  ;; Set up arrays
  ew_global = sdss_mkewarr(ewlim, ewbinsize, newbin=newbin)
  z_global = sdss_mkzarr(zlim, zbinsize, nzbin=nzbin)
  dz_global = replicate(zbinsize,nzbin)
  x_global = cosm_xz(z_global,/silent,/exact,/noinit)
  dx_global = cosm_xz(z_global+zbinsize,/silent,/exact,/noinit) - x_global

  rz_2darr = lonarr(nzbin*newbin) ; holds info over all LOSs; normalize by newbin
  dz_2darr = fltarr(nzbin*newbin)
  dx_2darr = dz_2darr
  rzlos_2darr = rz_2darr        ; holds info for LOS; reset ever LOS loop
  dzlos_2darr = dz_2darr
  dxlos_2darr = dz_2darr

  ;; Main loop
  iew_global = lindgen(newbin)
  for ff=0L,nfil-1 do begin
     ;; ;;;;;;;
     ;; Fractional contribution to dz(z), and dX(z) bins by LOS
     ;; Excludes Sky Lines (5579, 6302) as sdss_fndciv
     ;; Calls sdss_measuresnr() to get proper bounds
     ;; _extra= includes dvem=
     dz_los = sdss_calcdzlos( snrstrct[ff].z_qso, dblt, zlim, z_global, $
                              sdssdir+spec_fil[ff],$
                              zlim_los=zlim_los, zbinsize=zbinsize, $
                              nzpix=nztmp, iz_los=iz_los, dx_los=dx_los, $
                              x_global=x_global, dx_global=dx_global, $
                              xlim_los=xlim_los, rz_los=rz_los, debug=0, $
                              csomology=cosmology, _extra=extra)

     if nztmp eq 0 then begin
        if not keyword_set(silent) then $
           print,'sdss_completeness_dz: LOS outside of desired redshift range; skipping ',$
              snrstrct[ff].qso_name,string(snrstrct[ff].z_qso,format='(f8.5)')
        continue                ; skip
     endif 
     
     ;; Contribution of just this LOS but to all EW bins in grid
     for zz=0,nztmp-1 do begin  ; small array (likely) but bad stride
        indx = zz + iz_los[0] + iew_global * nzbin
        dzlos_2darr[indx] = dz_los[zz]
        dxlos_2darr[indx] = dx_los[zz]
        rzlos_2darr[indx] = rz_los[zz] ; could have gap
        if zz eq 0 then index_los = indx $
        else index_los = [index_los, indx]
     endfor                     ; loop zz=nztmp


     ;; Append this one LOS contribution; every other rzlos_2darr
     ;; element should be zero
     rz_2darr = rz_2darr + rzlos_2darr
     dz_2darr = dz_2darr + dzlos_2darr
     dx_2darr = dx_2darr + dxlos_2darr

     ;; ;;;;;;;
     ;; Reset for next loop
     rzlos_2darr[index_los] = 0. 
     dzlos_2darr[index_los] = 0.
     dxlos_2darr[index_los] = 0.

  endfor                        ; loop ff=nfil

  return, dz_2darr
end                             ; sdss_completeness_dz()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_czw, list_fil, snrstrct_fil, cmplt_fil, $
                           dblt_name=dblt_name, rec_param=rec_param, $
                           ewbinsize=ewbinsize, zbinsize=zbinsize, $ ; global
                           ewlim=ewlim, zlim=zlim, mccum_fil=mccum_fil, $
                           czn=czn, snrmin=snrmin, snrbinsize=snrbinsize, $
                           zmin=zmin, z2binsize=z2binsize, $ ; averaging
                           civobs_corr=civobs_corr, fulllist_fil=fulllist_fil, $
                           fullsnrstrct_fil=fullsnrstrct_fil, user=user,$
                           clobber=clobber, debug=debug, silent=silent, $
                           biasuser_fil=biasuser_fil, dvem=dvem, _extra=extra
  ;; Process sdssmcstrct from sdss_completeness_los
  if n_params() ne 3 then begin
     print,'Syntax - sdss_completeness_czw, list_fil, snrstrct_fil, cmplt_fil, '
     print,'                          [dblt_name=, /rec_param, ewbinsize=, zbinsize=,'
     print,'                          ewlim=, zlim=, cosmology=, /czn, mccum_fil=, '
     print,'                          snrmin=, snrbinsize=, zmin=, z2binsize=, '
     print,'                          civobs_corr=, fulllist_fil=, fullsnrstrct_fil=, '
     print,'                          user=, /clobber, /debug, /silent, biasuser_fil=, '
     print,'                          dvem=, _extra=]'
     return
  end
  sdssdir = sdss_getsdssdir()

  test = file_search(cmplt_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_completeness_czw: file exists and will not clobber; exiting'
     return                     ; EXIT
  endif 

  ;; Default parameters
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  if not keyword_set(zbinsize) then zbinsize = 0.05 ;0.005
  if not keyword_set(zlim) then zlim = [1., 6.]
  if not keyword_set(ewbinsize) then ewbinsize = 0.1 ;0.05 ; Ang
  if not keyword_set(ewlim) then ewlim = [0.05, 5.]  ; Ang
  ;; Following should match sdss_geetqsosubset() call
  if not keyword_set(snrbinsize) then snrbinsize = 0.25
  if not keyword_set(snrmin) then snrmin = 4
  if not keyword_set(z2binsize) then z2binsize = 0.25
  if not keyword_set(zmin) then zmin = zlim[0]

  cosmology = sdss_setcosmology(_extra=extra) ; _extra = cosmology[3]

  ;; Read in meta data
  readcol,list_fil,spec_fil,format='a',skip=1,/silent
  nfil = (size(spec_fil,/dim))[0]
  
  if size(snrstrct_fil,/type) eq 7 then $
     snrstrct = xmrdfits(snrstrct_fil,1,/silent) $ 
  else snrstrct = snrstrct_fil
  tags = tag_names(snrstrct)
  snrtag = (where(tags eq 'SNR_'+strupcase(dblt.ion)))[0]  

  ;; _extra includes /quick
  mcstrct_fil = sdss_getname(spec_fil,/spec,mc=dblt.ion,dir=mcdir,$
                             user=user,_extra=extra)
  mcstrct_fil = mcdir + mcstrct_fil

  ew_global = sdss_mkewarr(ewlim, ewbinsize, newbin=newbin)
  z_global = sdss_mkzarr(zlim, zbinsize, nzbin=nzbin)
  dz_global = replicate(zbinsize,nzbin)
  x_global = cosm_xz(z_global,/silent,/exact,/noinit)
  dx_global = cosm_xz(z_global+zbinsize,/silent,/exact,/noinit) - x_global

  ;; Template structure
  ;; R(z) varies from 0 to 1 and represents the fraction of bin probed
  ;; by LOS
  f2darr = fltarr(nzbin*newbin)
  l2darr = lonarr(nzbin*newbin)
  losstr_tmplt = {qso_name:'', $
                  z_qso:0., $
                  zlim:[0.,0.], $
                  snr:0., $
                  iz_qso:0L, $
                  isnr:0L, $
                  iz_snr:0L, $          ; = iz + nzbin * isnr
                  flg_mc:0, $           ; MC exists? 0 = no; 1 = yes
                  qso_name_used:'', $   ; which LOS actually used
                  index_sys:l2darr-1, $ ; all -1
                  numer_2darr:f2darr, $ ; = rzlos_2darr * nrec_2darr
                  denom_2darr:f2darr, $ 
                  nrec_2darr:f2darr, $   ; want floats even though whole #s
                  ninput_2darr:f2darr, $ ; for error analysis
                  rzlos_2darr:l2darr, $  ; holds info for LOS
                  dzlos_2darr:f2darr, $
                  dxlos_2darr:f2darr $
                 }
  los_arr = replicate(losstr_tmplt,nfil)
  ;; Load up generic information
  los_arr.qso_name = snrstrct.qso_name
  los_arr.z_qso = snrstrct.z_qso
  los_arr.snr = snrstrct.(snrtag)[2]
  los_arr.iz_qso = floor((los_arr.z_qso - zmin)/z2binsize)
  los_arr.isnr = floor((los_arr.snr - snrmin)/snrbinsize)
  nz2bin = max(los_arr.iz_qso)
  los_arr.iz_snr = los_arr.iz_qso + los_arr.isnr * nz2bin

  ;; ;;;;;;;
  ;; Main Loop
  ;; Read in files and process one by one
  iew_global = lindgen(newbin)
  for ff=0L,nfil-1 do begin 
     ;; ;;;;;;;
     ;; WARNING!!! This code needs to match sdss_completeness_dz()
     ;; Fractional contribution to dz(z), and dX(z) bins by LOS
     ;; Excludes Sky Lines (5579, 6302) as sdss_fndciv
     ;; Calls sdss_measuresnr() to get proper bounds
     dz_los = sdss_calcdzlos( snrstrct[ff].z_qso, dblt, zlim, z_global, $
                              sdssdir+spec_fil[ff], $ ; account for blanks in spectra
                              zlim_los=zlim_los, zbinsize=zbinsize, $
                              nzpix=nztmp, iz_los=iz_los, dx_los=dx_los, $
                              x_global=x_global, dx_global=dx_global, $
                              xlim_los=xlim_los, rz_los=rz_los, debug=0, $
                              csomology=cosmology, dvem=dvem, _extra=extra)
     if size(zlim_los,/n_dim) gt 1 then $
        los_arr[ff].zlim = [min(zlim_los,max=mx),mx] $ ; could have gap
     else begin
        los_arr[ff].zlim = zlim_los
        zlim_los = transpose(zlim_los) ; [1,2]
     endelse 

     if nztmp eq 0 then begin
        ;; Skip instantiating anything so nothing will pass through
        if not keyword_set(silent) then $
           print,'sdss_completeness_czw: LOS outside of desired redshift range; skipping ',$
                 snrstrct[ff].qso_name,string(snrstrct[ff].z_qso,format='(f8.5)')
        continue                ; skip
     endif 
     
     ;; Contribution of just this LOS but to all EW bins in grid
     for zz=0,nztmp-1 do begin  ; small array (likely) but bad stride
        indx = zz + iz_los[0] + iew_global * nzbin
        los_arr[ff].dzlos_2darr[indx] = dz_los[zz]
        los_arr[ff].dxlos_2darr[indx] = dx_los[zz]
        los_arr[ff].rzlos_2darr[indx] = rz_los[zz] ; could have gap
        if zz eq 0 then index_los = indx $
        else index_los = [index_los, indx]
     endfor                     ; loop zz=nztmp


     ;; ;;;;;;;
     ;; Now for the part where only some information exists
     ;; mcstr.nsys should be same for every one 
;     if keyword_set(debug) then $
;        test = file_search('/scratch/'+mcstrct_fil[ff]+'*',count=ntest) $
;     else $
     test = file_search(sdssdir+mcstrct_fil[ff]+'*',count=ntest)
     if ntest ne 0 then begin
        los_arr[ff].flg_mc = 1  ; exists
        if keyword_set(user) then mcstr = xmrdfits(test,3,/silent) $
        else mcstr = xmrdfits(test,1,/silent)

        ;; Select parameter to bin over
        if not keyword_set(ztag) then begin
           tags = tag_names(mcstr[0])
           if keyword_set(rec_param) then begin
              ;; Use *_REC tags
              ztag = (where(tags eq 'ZABS_REC'))[0]
              ncolmtag = (where(tags eq 'NCOLM_REC'))[0]
              signcolmtag = (where(tags eq 'SIGNCOLM_REC'))[0]
              flgtag = (where(tags eq 'NCOLMFLG_REC'))[0] ; might be -1
              ewtag = (where(tags eq 'EW_REC'))[0]
              sigewtag = (where(tags eq 'SIGEW_REC'))[0]
           endif else begin
              ;; Use *_INPUT tags
              ztag = (where(tags eq 'ZABS_INPUT'))[0]
              ncolmtag = (where(tags eq 'NCOLM_INPUT'))[0]
              flgtag = -1       ; by fiat
              ewtag = (where(tags eq 'EW_INPUT'))[0]
           endelse 

           ;; Select the tag over which the binning is occurring
           if keyword_set(czn) then begin
              vartag = ncolmtag
              if keyword_set(rec_param) then sigvartag = signcolmtag
           endif else begin
              vartag = ewtag
              if keyword_set(rec_param) then sigvartag = sigewtag
           endelse 
        
        endif                   ; ff = 0


        ;; ;;;;;;;
        ;; Contribution to completeness by input systems
        unq = uniq(mcstr.nsys)
        if n_elements(unq) ne 1 then stop
        for ss=0,mcstr[0].nsys-1 do begin
           
           ;; Mask out absobers in gaps, e.g., dvem, if som exist
           mask = mcstr.nsys*0 ; 1 = good, 0 = bad
           for ll=0,n_elements(zlim_los[*,0])-1 do begin
              gd = where(mcstr.(ztag)[ss,0] ge zlim_los[ll,0] and $
                         mcstr.(ztag)[ss,0] lt zlim_los[ll,1]) ; should it be [ss,1]?
              if gd[0] ne -1 then mask[gd]++
           endfor               ; loop ll = n_elements(zlim_los[*,0])
           gdmc = where(mask gt 0)

           if keyword_set(mccum_fil) then begin
              ;; Analyze meta information
              ;; save before changing ncolm
              if keyword_set(mcew_arr) then begin
                 mcew_arr = [mcew_arr,mcstr.(ewtag)[ss,0]]
                 mcncolm_arrI = [mcncolm_arrI,mcstr.(ncolmtag)[ss,0]]
                 mcz_arr = [mcz_arr,mcstr.(ztag)[ss,0]]
                 mclsnr_arr = [mclsnr_arr,mcstr.lsnr[ss,0]]
                 mcmask_arr = [mcmask_arr,mask]

                 mcew_arrII = [mcew_arrII,mcstr.(ewtag)[ss,1]]
                 mcncolm_arrII = [mcncolm_arrII,mcstr.(ncolmtag)[ss,1]]
                 mcz_arrII = [mcz_arrII,mcstr.(ztag)[ss,1]]
                 mclsnr_arrII = [mclsnr_arrII,mcstr.lsnr[ss,1]]
                 
                 if keyword_set(rec_param) then begin
                    mcsigew_arr = [mcsigew_arr,mcstr.(sigewtag)[ss,0]]
                    mcsigncolm_arrI = [mcsigncolm_arrI,mcstr.(signcolmtag)[ss,0]]
                    mcsigew_arrII = [mcsigew_arrII,mcstr.(sigewtag)[ss,1]]
                    mcsigncolm_arrII = [mcsigncolm_arrII,mcstr.(signcolmtag)[ss,1]]
                 endif
                 
                 mcflg_arr = [mcflg_arr,mcstr.flg_rec[ss]]
              endif else begin
                 mcew_arr = mcstr.(ewtag)[ss,0]
                 mcncolm_arrI = mcstr.(ncolmtag)[ss,0]
                 mcz_arr = mcstr.(ztag)[ss,0]
                 mclsnr_arr = mcstr.lsnr[ss,0]
                 mcmask_arr = mask
                 
                 mcew_arrII = mcstr.(ewtag)[ss,1]
                 mcncolm_arrII = mcstr.(ncolmtag)[ss,1]
                 mcz_arrII = mcstr.(ztag)[ss,1]
                 mclsnr_arrII = mcstr.lsnr[ss,1]
                 
                 if keyword_set(rec_param) then begin
                    mcsigew_arr = mcstr.(sigewtag)[ss,0]
                    mcsigncolm_arrI = mcstr.(signcolmtag)[ss,0]
                    mcsigew_arrII = mcstr.(sigewtag)[ss,1]
                    mcsigncolm_arrII = mcstr.(signcolmtag)[ss,1]
                 endif
                 
                 mcflg_arr = mcstr.flg_rec[ss]
              endelse 
              
           endif                ; mccum_fil=


           if keyword_set(czn) then begin
              ;; ;;;;;;;
              ;; Binning by column density, need to do extra stuff
              if keyword_set(rec_param) then begin
                 ;; Easier to copy over to right format
                 mccivstr = sdss_cpmc2civstrct(mcstr,ss)
                 ;; Want to use the same metric as used for the real
                 ;; absorbers, which means combo leverage on the AOD
                 ;; column density
                 ;; _extra= includes /estflg
                 ncolm = sdss_calcncolmdblt(mccivstr,dblt.ion,/log,$
                                            signcolm=signcolm,ncolmflg=ncolmflg,$
                                            silent=(not keyword_set(debug)),$
                                            _extra=extra)
                 mcstr.(ncolmtag)[ss,0] = ncolm ; overwrite, but saved
                 mcstr.(signcolmtag)[ss,0] = signcolm

                 mcstr.(ncolmtag)[ss,1] = mcstr.(ncolmtag)[ss,0]
                 mcstr.(signcolmtag)[ss,1] = mcstr.(signcolmtag)[ss,0]

                 if flgtag ne -1 then begin
                    mcstr.(flgtag)[ss,0] = ncolmflg
                    mcstr.(flgtag)[ss,1] = mcstr.(flgtag)[ss,0]
                 endif 
              endif

              if keyword_set(mccum_fil) then begin
                 ;; Save the new column densities (for system)
                 if keyword_set(mcncolm_arr) then begin
                    ;; concatenate
                    mcncolm_arr = [mcncolm_arr,mcstr.(ncolmtag)[ss,0]]
                    if keyword_set(rec_param) then $
                       mcsigncolm_arr = [mcsigncolm_arr,$
                                         mcstr.(signcolmtag)[ss,0]]
                    if flgtag ne -1 then $
                       mcsigncolmflg_arr = [mcsigncolm_arr,$
                                            mcstr.(flgtag)[ss,0]]
                 endif else begin
                    ;; instantiate
                    mcncolm_arr = mcstr.(ncolmtag)[ss,0]
                    if keyword_set(rec_param) then $
                       mcsigncolm_arr = mcstr.(signcolmtag)[ss,0]
                    if flgtag ne -1 then $
                       mcsigncolmflg_arr = mcstr.(flgtag)[ss,0]
                 endelse 
              endif             ; mccum_fil=
           endif                ; /czn


           ;; ;;;;;;;
           ;; Bins are defined on the LHS, so limit to everything from
           ;; lower limit to upper limit (previously had
           ;; +zbinsize/ewbinsize but that's wrong)
           if gdmc[0] eq -1 then begin
              print,'sdss_completeness_czw: no injected doublets in un-masked spectrum ',$
                    snrstrct[ff].qso_name,string(snrstrct[ff].z_qso,format='(f8.5)')
              continue          ; skip
           endif 

           gd = where(mcstr[gdmc].(ztag)[ss,0] ge zlim[0] and $
                      mcstr[gdmc].(ztag)[ss,0] lt zlim[1] and $
                      mcstr[gdmc].(vartag)[ss,0] ge ewlim[0] and $
                      mcstr[gdmc].(vartag)[ss,0] lt ewlim[1],ngd)

           if ngd eq 0 then begin
              print,'sdss_completeness_czw: no system ss='+strtrim(ss,2)+$
                    ' MC profiles in z and EW bounds ',$
                    snrstrct[ff].qso_name,string(snrstrct[ff].z_qso,format='(f8.5)'),$
                    zlim[0],zlim[1],zbinsize,$
                    ewlim[0],ewlim[1],ewbinsize
              continue          ; skip
           endif else gd = gdmc[gd]

           iz_sys = floor((mcstr[gd].(ztag)[ss,0] - zlim[0]) / zbinsize) 
           iew_sys = floor((mcstr[gd].(vartag)[ss,0] - ewlim[0]) / ewbinsize) ; wvI
           index_sys = iz_sys + iew_sys * nzbin                               ; linear indices
           los_arr[ff].index_sys[index_sys] = index_sys


           if keyword_set(debug) then begin
              ;; Curious why these can wander out of bounds...
              zmn = min(iz_sys,max=zmx)
              if zmn lt iz_los[0] or zmx gt iz_los[1] then $
                 print,'sdss_completeness_czw debug: input system z range '+$
                       'larger than LOS z limits, ss =',$
                       ss,' '+mcstr[0].qso_name
           endif                ; /debug

           ;; "Hit" in grid only if mcstr[gd].flg_rec[ss] = 1
           los_arr[ff].nrec_2darr[index_sys] = $
              los_arr[ff].nrec_2darr[index_sys] + mcstr[gd].flg_rec[ss] 
           los_arr[ff].ninput_2darr[index_sys]++

           los_arr[ff].numer_2darr[index_sys] = $
              los_arr[ff].numer_2darr[index_sys] + $
              mcstr[gd].flg_rec[ss] * los_arr[ff].rzlos_2darr[index_sys]
           los_arr[ff].denom_2darr[index_sys] = $
              los_arr[ff].denom_2darr[index_sys] + $
              los_arr[ff].rzlos_2darr[index_sys]

           if keyword_set(debug) then begin
              ;; Check that rzlos_2darr[index_sys] all instantiated;
              ;; should never have a case where there's an input
              ;; system with no R(z) value
              test = where(los_arr[ff].rzlos_2darr[index_sys] ne 1,ntest)
              if test[0] ne -1 then $
                 print,'sdss_completeness_czw debug: R(z) subset != 1',ntest
           endif                ; /debug

        endfor                  ; loop ss=mcstr[0].nsys

     endif                      ; mcstrct_fil exists for LOS

  endfor                        ; loop ff=nfil
  

  ;; ;;;;;;;;
  ;; Since may be working on subsample, have to duplicate sampled LOSs
  ;; in lieu of unsampled LOSs (according to zqso and S/N)
;  rng = where((los_arr.zlim[0] ge zlim[0] and los_arr.zlim[0] lt zlim[1]) or $
;              (los_arr.zlim[1] gt zlim[0] and los_arr.zlim[1] le zlim[1]))
  ;; Trying to fix MgII dN/dX problems
  rng = where(los_arr.zlim[0] lt zlim[1] or los_arr.zlim[1] gt zlim[0])
  if rng[0] ne -1 then begin
     sub = where(los_arr[rng].flg_mc eq 0,complement=smpl)
     if sub[0] eq -1 then begin
        unsmpl = [-1]           ; skip next part
     endif else begin
        ;; Have to duplicate
        unsmpl = rng[sub]
        smpl = rng[smpl]
     endelse 

     if unsmpl[0] ne -1 then begin
        ;; Loop and whittle down by figuring out what needs to be
        ;; duplicated 
        unq = uniq(los_arr[unsmpl].iz_snr,sort(los_arr[unsmpl].iz_snr))
        unsmpl_indx = los_arr[unsmpl[unq]].iz_snr
        nfil = (size(unsmpl_indx,/dim))[0] > 1 ; foil singularity
        unq = uniq(los_arr[smpl].iz_snr,sort(los_arr[smpl].iz_snr))
        smpl_indx = los_arr[smpl[unq]].iz_snr

        for ff=0L,nfil-1 do begin ; loop over *UN*sampled LOS (should be smaller)
           unsmplsub = unsmpl[where(los_arr[unsmpl].iz_snr eq $
                                    unsmpl_indx[ff],nunsmplsub)]        ; always match
           gd = where(los_arr[smpl].iz_snr eq unsmpl_indx[ff],nsmplsub) ; not so
           if gd[0] eq -1 then begin
              ;; The binning may be different and should just find
              ;; closest in zqso and S/N
              iz = unsmpl_indx[ff] mod nz2bin
              isnr = (unsmpl_indx[ff] - iz) / nz2bin
              if not keyword_set(iz_smpl) then begin
                 iz_smpl = smpl_indx mod nz2bin
                 isnr_smpl = (smpl_indx - iz_smpl) / nz2bin
              endif
              mn = min(sqrt((iz-iz_smpl)^2 + (isnr-isnr_smpl)^2),imn)
              indx = iz_smpl[imn] + isnr_smpl[imn] * nz2bin
              gd = where(los_arr[smpl].iz_snr eq indx,nsmplsub)

              if keyword_set(debug) then $
                 print,'sdss_completeness_czw: desired zqso and S/N and closest',$
                       zmin+iz*z2binsize, snrmin+isnr*snrbinsize, $
                       zmin+iz_smpl[imn]*z2binsize, snrmin+isnr_smpl[imn]*snrbinsize
           endif 
           smplsub = smpl[gd]
           
           norm = 1./float(nsmplsub) ; for averaging
           los_arr[unsmplsub].qso_name_used = $
              strjoin(los_arr[smplsub].qso_name,'; ') ; save info
           for gg=0,nsmplsub-1 do begin
              ;; For sampled LOS(s) (where mcstrct_fil exists), which EW
              ;; bins were actually tested
              index_sys = where(los_arr[smplsub[gg]].index_sys ge 0,nindex)

              if nindex eq 0 then begin
                 ;; This results when the representative sightline
                 ;; didn't manage to have coverage... 
                 print,'sdss_completeness_czw: representative LOS does not have coverage ', $
                       los_arr[smplsub[gg]].qso_name,los_arr[smplsub[gg]].z_qso
                 continue
              endif 

              ;; For unsampled LOS(s), need to take average of any
              ;; sampled and contributing LOS(s)
              if nunsmplsub gt 1 then begin
                 ;; Need to inflate to match unsmplsub dimensions
                 add_numer = rebin(los_arr[smplsub[gg]].numer_2darr[index_sys],$
                                   nindex, nunsmplsub)
                 add_denom = rebin(los_arr[smplsub[gg]].denom_2darr[index_sys],$
                                   nindex, nunsmplsub)
              endif else begin
                 add_numer = los_arr[smplsub[gg]].numer_2darr[index_sys]
                 add_denom = los_arr[smplsub[gg]].denom_2darr[index_sys]
              endelse 

              los_arr[unsmplsub].numer_2darr[index_sys] = $
                 los_arr[unsmplsub].numer_2darr[index_sys] + add_numer * norm
              los_arr[unsmplsub].denom_2darr[index_sys] = $
                 los_arr[unsmplsub].denom_2darr[index_sys] + add_denom * norm
           endfor               ; loop gg=nsmplsub

        endfor                  ; loop ff=nfil

     endif                      ; nunsmpl != 0
  endif else begin
     stop,'sdss_completeness_czw stop: no LOS covering given zlim; null results.'
  endelse 


  ;; ;;;;;;;
  ;; Combine information
  ;; Collpase
  nrec_2darr = total(los_arr.nrec_2darr,2)
  ninput_2darr = total(los_arr.ninput_2darr,2)
  numer_2darr = total(los_arr.numer_2darr,2)
  denom_2darr = total(los_arr.denom_2darr,2)

  rz_2darr = total(los_arr.rzlos_2darr,2)
  dz_2darr = total(los_arr.dzlos_2darr,2)
  dx_2darr = total(los_arr.dxlos_2darr,2)
  

  ;; M(z,W) = N_rec / N_in (i.e. number recovered automatically over
  ;; number input)
  ;; If perfect search, M(z,W) = 1
  mzw_2darr = f2darr            ; can reform to [nzbin,newbin]
  gd = where(denom_2darr ne 0.,complement=bd)
  if gd[0] ne -1 then $
     mzw_2darr[gd] = numer_2darr[gd] / denom_2darr[gd]
  
  ;; A(z,W) = N_accept / N_rec (i.e. number interactively accepted
  ;; over number recovered)
  ;; False positives handled later, for now this can be 1
  azw_2darr = mzw_2darr * 0. + 1. ; for now

  ;; F(z,W) = continuum correction
  fzw_2darr = azw_2darr

  ;; Sensitivity is the sum of the product of the completeness of each
  ;; QSO and 1 (if (z, W) bin coverd by QSO) or 0 (if not). This
  ;; latter part is denoted R(z).
  ;; Cq(z,W) = M(z,W)*A(z,W)
  ;; g(z,W) = SUM( Cq(z,W) * R(z) , QSOs)
  ;; let h(z) = SUM( R(z), QSOs )
  czw_2darr = mzw_2darr * azw_2darr ; mzw_2darr folds in R(z) actually

  ;; _extra includes cl=, sigma=
  sigczw_2darr = sdss_calcsigbinom(ninput_2darr, nrec_2darr, $
                                   verbose=(not keyword_set(silent)),_extra=extra)
  

  ;; ;;;;;;;
  ;; Output
  ;; Have to create structures in situ because they will have variable
  ;; lengths 
  cmpltstr = { $
             ;; General parameters, so can re-create this
             list_fil:list_fil, $
             dblt_name:dblt.ion, $
             rec_param:keyword_set(rec_param), $
             civobs_corr:keyword_set(civobs_corr), $
             userbias:keyword_set(biasuser_fil), $
             czn:keyword_set(czn), $
             nzbin:nzbin, $
             zlim:zlim, $
             zbinsize:zbinsize, $
             newbin:newbin, $
             ewlim:ewlim, $
             ewbinsize:ewbinsize, $
             cosmology:cosmology, $
             ;; Now for actual results in this bin
             ninput_2darr:ninput_2darr, $
             nrec_2darr:nrec_2darr, $
             rz_2darr:rz_2darr, $ ; can reform(2darr, nzbin, newbin)
             dz_2darr:dz_2darr, $
             dx_2darr:dx_2darr, $
             czw_2darr:czw_2darr, $
             sigczw_2darr:sigczw_2darr $ ; [*,0] = lower, [*,1] = upper
             }
  if keyword_set(dvem) then $   ; Record
     tmp = create_struct(cmpltstr,'dvem',dvem) $
  else tmp = create_struct(cmpltstr,'dvem',0) 
  cmpltstr = tmp
  

  if keyword_set(fulllist_fil) then begin
     ;; Since may be working on subsample, have to calculate new
     ;; dz_2darr, dX_2darr, rz_2darr
     ;; SHOUlD BE OBSOLETE SINCE NOW INCLUDE EVERYTHING
     ;; (and this may now be broken)
     fulldz_2darr = sdss_completeness_dz(fulllist_fil, fullsnrstrct_fil, $
                                         ewbinsize=ewbinsize, zbinsize=zbinsize, $
                                         ewlim=ewlim, zlim=zlim, dblt_name=dblt, $
                                         dx_2darr=fulldx_2darr, $
                                         rz_2darr=fullrz_2darr, debug=0, $
                                         cosmology=cosmology, silent=silent, $
                                         dvem=dvem)

     ;; Append the information and copy over
     newstrct = create_struct(cmpltstr, $
                              'fulllist_fil',fulllist_fil, $
                              'fullrz_2darr',fullrz_2darr, $ 
                              'fulldz_2darr',fulldz_2darr, $
                              'fulldx_2darr',fulldx_2darr )
     cmpltstr = sdss_cpstrct(cmpltstr,newstrct)
  endif 
  
  
  ;;;;;;;;;
  ;; Factor in the effect of pathlength blocked by actual CIV systems 
  ;; Since did *not* remove any/all during completeness tests, real
  ;; CIV systems will artificially allow any weaker simulated profile
  ;; laid on top
  if keyword_set(civobs_corr) then begin
     print,'SDSS_COMPLETENESS_CZW: APPLYING CIVOBS_CORR!!!'

     ;; _extra includes /default, rating=, dvqso=, /noBAL, /final
     cmpltstr = sdss_applycivobscorr(cmpltstr, civobs_corr, dvem=dvem, $
                                     rec_param=rec_param, _extra=extra)
  endif                         ; /civobs_corr


  ;; ;;;;;;;
  ;; Account for user bias
  ;; C(W) = Nrec / Ninput * Naccept / Nrec
  if keyword_set(biasuser_fil) then begin
     
     print,'SDSS_COMPLETENESS_CZW: APPLYING BIASUSER_FIL!!!'
     cmpltstr = sdss_applyuserbias(cmpltstr, biasuser_fil, $
                                   civobs_corr=0, _extra=extra)
     
  endif                         ; biasuser_fil

  
  ;; ;;;;;;;
  ;; Write to file
  mwrfits,cmpltstr,cmplt_fil,/create,/silent
  if keyword_set(debug) then $
     mwrfits,los_arr,cmplt_fil,/silent
  spawn,'gzip -f '+cmplt_fil
  print,'sdss_completeness_czw: created ',cmplt_fil

  if keyword_set(mccum_fil) then begin
     ;; Just the facts, Ma'am
     nmc = (size(mcew_arr,/dim))[0]
     dum_arr = fltarr(nmc,2,/nozero)
     mcinfostrct = {rec_param:keyword_set(rec_param),$
                    zabs:dum_arr, $
                    ew:dum_arr, $
                    ncolm:dum_arr, $
                    lsnr:dum_arr, $
                    mask:mcmask_arr, $
                    flg_rec:mcflg_arr}
     mcinfostrct.zabs[*,0] = mcz_arr
     mcinfostrct.ew[*,0] = mcew_arr
     mcinfostrct.ncolm[*,0] = mcncolm_arrI
     mcinfostrct.lsnr[*,0] = mclsnr_arr

     mcinfostrct.zabs[*,1] = mcz_arrII
     mcinfostrct.ew[*,1] = mcew_arrII
     mcinfostrct.ncolm[*,1] = mcncolm_arrII
     mcinfostrct.lsnr[*,1] = mclsnr_arrII

     if keyword_set(rec_param) then begin
        ;; Append the error information
        mcinfostrct = create_struct(mcinfostrct,$
                                    'sigew',dum_arr,'signcolm',dum_arr)
        mcinfostrct.sigew[*,0] = mcsigew_arr
        mcinfostrct.signcolm[*,0] = mcsigncolm_arrI
        mcinfostrct.sigew[*,1] = mcsigew_arrII
        mcinfostrct.signcolm[*,1] = mcsigncolm_arrII
     endif

     if keyword_set(czn) then begin
        ;; Append the total column density adopted
        mcinfostrct = create_struct(mcinfostrct,'ncolmdblt',mcncolm_arr)
        if keyword_set(rec_param) then $
           mcinfostrct = create_struct(mcinfostrct,'signcolmdblt',mcsigncolm_arr)
        if flgtag ne -1 then $
           mcinfostrct = create_struct(mcinfostrct,'ncolmflgdblt',mcncolmflg_arr)
     endif 

     test = file_search(mccum_fil+'*',count=ntest)
     if ntest eq 0 or keyword_set(clobber) then begin
        mwrfits,mcinfostrct,mccum_fil,/silent
        spawn,'gzip -f '+mccum_fil
        print,'sdss_completeness_czw: created ',mccum_fil
     endif else $
        stop,'sdss_completeness_czw stop: will not overwrite ',mccum_fil
  endif                         ; /debug append output


  ;; ;;;;;;;
  ;; Jumbles of plots for debugging and sanity checking
  if keyword_set(debug) then begin
     czw_grid = reform(czw_2darr,nzbin,newbin)
     rz_grid = reform(rz_2darr,nzbin,newbin) 
     rz_arr = total(rz_grid,2)/float(newbin)
     
     dz_grid = reform(dz_2darr,nzbin,newbin) 
     dz_arr = total(dz_grid,2)/float(newbin) 
     
     dx_grid = reform(dx_2darr,nzbin,newbin)
     dx_arr = total(dx_grid,2)/float(newbin)

     ;; x_splot,z_global,rz_arr,psym1=10,xtitle='global z',ytitle='R(z)',/block
     ;; x_splot,z_global,rz_arr*dz_global,psym1=10,ytwo=dz_arr,psym2=10,xtitle='global z',ytitle='R(z)*global dz',/block
     ;; x_splot,z_global,rz_arr*dx_global,psym1=10,ytwo=dx_arr,psym2=10,xtitle='global z',ytitle='R(z)*global dX',/block
     ;; x_splot,z_global,dx_arr/(rz_arr*dx_global),psym1=10,xtitle='global z',ytitle='dX / (R(z)*global dX)',ytwo=dz_arr/(rz_arr*dz_global),psym2=10,ymnx=[0.,1.],/block

     ;; x_splot,z_global,dz_arr,psym1=10,ytwo=dx_arr,xtitle='global z',psym2=10,ytitle='pathlength',lgnd=['dz','dX'],ymnx=[0.,max(dx_arr)],/block
     ;; x_splot,denom_2darr,numer_2darr,psym1=4,xtitle='Number input',ytitle='Number recovered',xtwo=[0,1e6],ytwo=[0,1e6],/block  
     ee=4 & de=4 & 
     lgnd='EW='+[string(ew_global[ee],format='(f5.2)'),$
                 string(ew_global[ee+de],format='(f5.2)'),$
                 string(ew_global[ee+2*de],format='(f5.2)')]
     x_splot,z_global,czw_grid[*,ee],psym1=10,ytwo=czw_grid[*,ee+de],$
             psym2=10,ythr=czw_grid[*,ee+2*de],psym3=10,xtitle='global z',$
             ytitle='C(z,W)',ymnx=[0,max(czw_grid[*,ee+2*de])],lgnd=lgnd,/block

     zz=10 & dz=6 
     lgnd='z='+[string(z_global[zz],format='(f5.2)'),$
                string(z_global[zz+dz],format='(f5.2)'),$
                string(z_global[zz+2*dz],format='(f5.2)')]
     ;; x_splot,ew_global,rz_grid[zz,*],psym1=10,ytwo=rz_grid[zz+dz,*],psym2=10,ythr=rz_grid[zz+2*dz,*],psym3=10,xtitle='global EW (Ang)',ytitle='R(z)',ymnx=[0,max(rz_grid[zz+2*dz,*])],lgnd=lgnd,/block

     ;; gd=where(mcflg_arr eq 1 and mcz_arr lt 2) & bd=where(mcflg_arr eq 0 and mcz_arr lt 2)
     gd=where(mcflg_arr eq 1,complement=bd) 
     ;; x_splot,mcz_arr[gd],mcew_arr[gd],xtwo=mcz_arr[bd],ytwo=mcew_arr[bd],psym1=3,psym2=3,xtitle='z',ytitle='EW',/block

     ;; sub = where(mcz_arr lt 2) 
     if keyword_set(czn) then begin
        histew=histogram(mcncolm_arr,binsize=ewbinsize,min=ewlim[0],max=ewlim[1],loc=locew) 
        histewgd=histogram(mcncolm_arr[gd],binsize=ewbinsize,min=ewlim[0],max=ewlim[1])
        histewbd=histogram(mcncolm_arr[bd],binsize=ewbinsize,min=ewlim[0],max=ewlim[1])
     endif else begin
        histew=histogram(mcew_arr,binsize=ewbinsize,min=ewlim[0],max=ewlim[1],loc=locew) 
        histewgd=histogram(mcew_arr[gd],binsize=ewbinsize,min=ewlim[0],max=ewlim[1])
        histewbd=histogram(mcew_arr[bd],binsize=ewbinsize,min=ewlim[0],max=ewlim[1])
     endelse 
     ;; x_splot,locew,histew,ytwo=histewgd,ythr=histewbd,psym1=10,psym2=10,psym3=10,xtitle='EW (Ang)',ytitle='Number',lgnd=['Inp','Rec','Not Rec'],/block 

     ;; x_splot,locew,histewgd/float(histew),ytwo=histewbd/float(histew),psym1=10,psym2=10,xtitle='EW (Ang)',ytitle='Fraction',lgnd=['Rec','Not Rec'],/block


     clr = getcolor(/load)
     contour, czw_grid, z_global, ew_global, levels=[0.1,0.5,0.95], $
              yrange=[0.,2.], xrange=[1.4,5], /xstyle, /ystyle, xtitle='z', $
              ytitle='EW (Ang)',title='C(z,W) Contours: 10%, 50%, 95%',$
              background=clr.white,color=clr.black,thick=1,charsize=2,$
              c_colors=[clr.black,clr.red,clr.blue]

     numer_grid=reform(numer_2darr,nzbin,newbin)
     denom_grid=reform(denom_2darr,nzbin,newbin)
     numer=total(numer_grid,1) & denom=total(denom_grid,1) ; collapse over EW
     
     ;; x_splot,locew,histew,ytwo=histewgd,ythr=histewbd,psym1=10,psym2=10,psym3=10,xtitle='EW (Ang)',ytitle='Number',lgnd=['Inp','Rec','Not Rec','Numer','Denom'],xfou=ew_global,yfou=numer,psym4=10,xfiv=ew_global,yfiv=denom,psym5=10,/block 
     ;; 
     ;; x_splot,locew,histewgd/float(histew),ytwo=histewbd/float(histew),psym1=10,psym2=10,xtitle='EW (Ang)',ytitle='Fraction',lgnd=['Rec','Not Rec','Numer/Denom'],xthr=ew_global,ythr=numer/float(denom),psym3=10,/block
     
     ;; x_splot,locew,histewgd/float(numer[0:49]),psym1=10,ytwo=histew/float(denom[0:49]),psym2=10,ymnx=[0.,5.],xtitle='EW (Ang)',ytitle='Ratio',lgnd=['Rec/Numer','Total Input/Denom'],xthr=locew,ythr=histew/float(max(histew)),psym3=10,/block 
  endif                         ; /debug 
  
end                             ; sdss_completeness_czw


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_plot, cmplt_fil, zindx=zindx, ewindx=ewindx, $
                            dzindx=dzindx, dewindx=dewindx, gz=gz, $
                            gx=gx, contour=contour, full=full


  if size(cmplt_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil
  dblt = dblt_retrieve(cmpltstr.dblt_name)

  tags = tag_names(cmpltstr)
  if keyword_set(full) then begin
     dztag = (where(tags eq 'FULLDZ_2DARR'))[0]
     dxtag = (where(tags eq 'FULLDX_2DARR'))[0]
  endif else begin
     dztag = (where(tags eq 'DZ_2DARR'))[0]
     dxtag = (where(tags eq 'DX_2DARR'))[0]
  endelse 
  
  ;; Transform information to grid
  czw_grid = reform(cmpltstr.czw_2darr,cmpltstr.nzbin,cmpltstr.newbin)

  if not keyword_set(zindx) and not keyword_set(ewindx) then begin
     ;; plot both
     zindx = 10
     ewindx = 4
  endif 
  if not keyword_set(dzindx) then dzindx = 6
  if not keyword_set(dewindx) then dewindx = 4

  if keyword_set(gz) then begin
     yttl = 'g(z,W)' 
     dz_grid = reform(cmpltstr.(dztag),cmpltstr.nzbin,cmpltstr.newbin)
     gzw_grid = czw_grid * dz_grid
  endif else begin
     if keyword_set(gx) then begin
        yttl = 'g(X,W)' 
        dX_grid = reform(cmpltstr.(dxtag),cmpltstr.nzbin,cmpltstr.newbin)
        gzw_grid = czw_grid * dX_grid
     endif else begin
        yttl = 'C(z,W)'
        gzw_grid = czw_grid
     endelse 
  endelse 
  xttl_z = 'z'+strtrim(floor(dblt.wvI),2)
  xttl_ew = 'EW (Ang)'

  ;; Set up necessary arrays
  ew_global = sdss_mkewarr(cmpltstr.ewlim, cmpltstr.ewbinsize)
  if (size(cmpltstr.zlim,/dim))[0] gt 1 then $
     z_global = cmpltstr.zlim[*,0] $
  else $
     z_global = sdss_mkzarr(cmpltstr.zlim, cmpltstr.zbinsize)

  if keyword_set(ewindx) then begin
     lgnd = 'EW='+[string(ew_global[ewindx],format='(f5.2)'),$
                   string(ew_global[ewindx+dewindx],format='(f5.2)'),$
                   string(ew_global[ewindx+2*dewindx],format='(f5.2)')]
     x_splot,z_global,gzw_grid[*,ewindx],psym1=10,$
             ytwo=gzw_grid[*,ewindx+dewindx],$
             psym2=10,ythr=gzw_grid[*,ewindx+2*dewindx],psym3=10,$
             xtitle=xttl_z,$
             ytitle=yttl,ymnx=[0,max(gzw_grid[*,ewindx+2*dewindx])],$
             lgnd=lgnd,/block
  endif 

  if keyword_set(zindx) then begin
     lgnd = 'z='+[string(z_global[zindx],format='(f5.2)'),$
                  string(z_global[zindx+dzindx],format='(f5.2)'),$
                  string(z_global[zindx+2*dzindx],format='(f5.2)')]
     x_splot,ew_global,gzw_grid[zindx,*],psym1=10,$
             ytwo=gzw_grid[zindx+dzindx,*],psym2=10,$
             ythr=gzw_grid[zindx+2*dzindx,*],psym3=10,$
             xtitle=xttl_ew,ytitle=yttl,$
             ymnx=[0,max(gzw_grid[zindx+2*dzindx,*])],lgnd=lgnd,/block
  endif 

  if keyword_set(contour) then begin
     clr = getcolor(/load)
     contour, czw_grid, z_global, ew_global, levels=[0.1,0.5,0.95], $
              yrange=[0.,2.], xrange=[1.4,5], /xstyle, /ystyle, xtitle='z', $
              ytitle='EW (Ang)',title='g(z,W) Contours: 10%, 50%, 95%',$
              background=clr.white,color=clr.black,thick=1,charsize=2,$
              c_colors=[clr.black,clr.red,clr.blue]
  endif 

end                             ; sdss_completeness_plot


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_cum, n_per_bin, froot, civ_fil, round_scale=round_scale, $
                           dblt_name=dblt_name,final=final, z=z, $
                           proot=proot, debug=debug, nociv=nociv, $
                           out_fil=out_fil,ewciv_lim=ewciv_lim,_extra=extra
  ;; Hardcode the various subsets desired
  if n_params() lt 2 then begin
     print,'Syntax - sdss_completeness_cum, n_per_bin, froot, [civ_fil, round_scale=, '
     print,'                           dblt_name=, /final, /z, proot=, /debug, '
     print,'                           /nociv, out_fil=,ewciv_lim=, _extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  zfudge = 1e-5                ; rounding errors in single bins

  ;; Parameters for setting binning
  if not keyword_set(round_scale) then round_scale = 100.
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  ;; Contruct other file names
  if not keyword_set(proot) then $
     proot = sdssdir+'inputs/SNRge4/dr7qso_'+dblt.ion+'_noBAL' ; param root
  list_fil = proot+'.list'
  snrstrct_fil = proot+'_SNR.fit'

  if not keyword_set(nociv) then begin
     ;; Read in data
     ;; _extra includes rating=, zlim=, ewlim=, nlim=, dvgal=, dvqso=, devm=,
     ;; /noBAL, /unblend, /default, /dropbox, civstrct_fil=
     civstr = sdss_getcivstrct(civ_fil,/default,final=final,count=nciv,$
                               ewlim=ewciv_lim,_extra=extra)
     tags = tag_names(civstr[0])
     if keyword_set(final) then $
        ztag = (where(tags eq 'ZABS_FINAL'))[0] $
     else ztag = (where(tags eq 'ZABS_ORIG'))[0] 
     
     ;; Work up in redshift
     zsrt = civstr[sort(civstr.(ztag)[0])].(ztag)[0]
  endif else nciv = 0

  if keyword_set(z) then begin
     ;; n_per_bin is actually redshift bins
     zlim = n_per_bin
     nzbin = (size(zlim[*,0],/dim))[0] > 1 ; foil singularity
  endif else begin
     ;; Find max number of bins to put everything in
     hist = sdss_histogram(civstr.(ztag)[0],n_per_bin,location=loc,$
                           szloc=dloc, nbin=nzbin)
     zlim = fltarr(nzbin,2,/nozero)
     zlim[*,0] = round(loc*round_scale)/round_scale    
     zlim[0,0] = floor(loc[0]*round_scale)/round_scale 
     zlim[*,1] = shift(zlim[*,0],-1) ; make match
     zlim[nzbin-1,1] = ceil((loc[nzbin-1]+dloc[nzbin-1])*round_scale)/$
                       round_scale
  endelse 

  ;; Save infomation
  suffix = strarr(nzbin)
  n_per_bin_use = lonarr(nzbin)
  out_fil = strarr(nzbin+1)
  
  for ii=0,nzbin-1 do begin
     suffix[ii] = string(zlim[ii,0],zlim[ii,1],format='(f4.2,"z",f4.2)')
     if not keyword_set(nociv) then begin
        ;; Figure out number in bin
        gd = where(zsrt ge zlim[ii,0] and zsrt lt zlim[ii,1],ngd)
        n_per_bin_use[ii] = ngd
        
        ;; Informational string to append to froot (filename)
        num = strtrim(n_per_bin_use[ii],2)
        bd = where(strlen(num) lt max(strlen(num)))
        while bd[0] ne -1 do begin
           num[bd] = '0'+num[bd]
           bd = where(strlen(num) lt max(strlen(num)))
        endwhile
        suffix[ii] = suffix[ii] + 'n' + num

        ;; Set for next loop
        if ngd ne 0 then $
           istart = gd[ngd-1] + 1
     endif 
     
     ;; Print information
     if keyword_set(debug) then begin
        if ii eq 0 then $
           print,'zlo','zhi','dz','#',format='(3(a5,2x),a5)'
        print,zlim[ii,0],zlim[ii,1],zlim[ii,1]-zlim[ii,0],n_per_bin_use[ii],$
              format='(3(f5.3,2x),i5)'
     endif                      ; /debug

     ;; Generate completeness file for one redshift bin
     ;; _extra includes: /rec_param, ewbinsize=, ewlim=, /czn, 
     ;;                  snrmin=, snrbinsize=, zmin=, z2binsize=, 
     ;;                  /clobber, /silent, civobs_corr=, dvem=, [_extra=]
     ;; That _extra would include: cosmology=,  cl=, /quick
     ;; Inflate the bin zbinsize just to avoid making an extra zbin
     ;; due to rounding issues
     ztmp = reform(zlim[ii,*])  ; remove extra dimension of 1
     out_fil[ii] = froot+suffix[ii]+'.fit'
     sdss_completeness_czw, list_fil, snrstrct_fil, out_fil[ii], $
                            zlim=ztmp,dblt_name=dblt, debug=debug, $
                            zbinsize=zlim[ii,1]-zlim[ii,0]+zfudge,$
                            _extra=extra
     
     ;; Append the files into new completeness structure
     cmpltstr = xmrdfits(out_fil[ii],1,/silent)
     if ii eq 0 then begin
        ;; Instantiate
        iew_global = lindgen(cmpltstr.newbin)
        ninput_2darr = fltarr(nzbin*cmpltstr.newbin,/nozero)
        nrec_2darr = ninput_2darr
        rz_2darr = lonarr(nzbin*cmpltstr.newbin,/nozero)
        dz_2darr = ninput_2darr
        dx_2darr = dz_2darr
        czw_2darr = dz_2darr
        sigczw_2darr = fltarr(nzbin*cmpltstr.newbin,2,/nozero)
     endif 
     index = ii + iew_global*nzbin
     ninput_2darr[index] = cmpltstr.ninput_2darr
     nrec_2darr[index] = cmpltstr.nrec_2darr
     rz_2darr[index] = cmpltstr.rz_2darr
     dz_2darr[index] = cmpltstr.dz_2darr
     dx_2darr[index] = cmpltstr.dx_2darr
     czw_2darr[index] = cmpltstr.czw_2darr
     sigczw_2darr[index,0] = cmpltstr.sigczw_2darr[*,0]
     sigczw_2darr[index,1] = cmpltstr.sigczw_2darr[*,1]
  endfor                        ; loop ii=nzbin

  ntot = round(total(n_per_bin_use))
  if nciv ne ntot and not keyword_set(z) then $
     stop,'sdss_completeness_cum stop: did not bin correctly; missing data'

  ;; ;;;;;;;
  ;; Create uber file of irregular grids
  nwcmpltstr = { $
               ;; General parameters, so can re-create this
               list_fil:list_fil, $
               dblt_name:dblt.ion, $
               rec_param:cmpltstr.rec_param, $
               civobs_corr:cmpltstr.civobs_corr, $
               userbias:cmpltstr.userbias, $
               czn:cmpltstr.czn, $
               nzbin:nzbin, $
               zlim:zlim, $
               zbinsize:(zlim[*,1]-zlim[*,0])+zfudge, $
               newbin:cmpltstr.newbin, $
               ewlim:cmpltstr.ewlim, $
               ewbinsize:cmpltstr.ewbinsize, $
               cosmology:cmpltstr.cosmology, $
               ;; Now for actual results in this bin
               ninput_2darr:ninput_2darr, $
               nrec_2darr:nrec_2darr, $
               rz_2darr:rz_2darr, $ ; can reform(2darr, nzbin, newbin)
               dz_2darr:dz_2darr, $
               dx_2darr:dx_2darr, $
               czw_2darr:czw_2darr, $
               sigczw_2darr:sigczw_2darr, $ ; [*,0] = lower, [*,1] = upper
               dvem:cmpltstr.dvem $
               }
  
  ;; Informational string to append to froot (filename)
  nwsuffix = string(zlim[0,0],zlim[nzbin-1,1],format='(f4.2,"z",f4.2)')
  if not keyword_set(nociv) then begin
     if keyword_set(z) then $
        num = strtrim(ntot,2) $
     else num = strtrim(n_per_bin,2)
     nwsuffix = nwsuffix + 'n' + num
  endif 
  out_fil[nzbin] = froot+nwsuffix+'cum.fit'
  mwrfits,nwcmpltstr,out_fil[nzbin],/create,/silent
  spawn,'gzip -f '+out_fil[nzbin]
  print,'sdss_completeness_cum: created ',out_fil[nzbin]

end                             ;  sdss_completeness_cum


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_estncolm, cmplt_fil, mcinp_fil, mcrec_fil, ostrct_fil, $
                                civstrct_fil=civstrct_fil,check=check,$
                                incl_col=incl_col,$
                                clobber=clobber,silent=silent,_extra=extra

  if n_params() ne 4 then begin
     print,'Syntax - sdss_completeness_estncolm, cmplt_fil, mcinp_fil, mcrec_fil, ostrct_fil, '
     print,'                                /incl_col,civstrct_fil=,/check,/clobber,_extra='
     return
  endif

  test = file_search(ostrct_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_completeness_estncolm: file exists; will not clobber ',$
           ostrct_fil
     return
  endif 

  civstr = sdss_getcivstrct(civstrct_fil,count=nciv,/default,_extra=extra)

  if size(cmplt_fil,/type) eq 7 then $
     cmpltstr = xmrdfits(cmplt_fil,1,/silent) $
  else cmpltstr = cmplt_fil

  if size(mcinp_fil,/type) eq 7 then $
     mcinpstr = xmrdfits(mcinp_fil,1,/silent) $
  else mcinpstr = mcinp_fil
  nmc = (size(mcinpstr.zabs,/dim))[0]
  imcall = lindgen(nmc)         ; may need

  if size(mcrec_fil,/type) eq 7 then $
     mcrecstr = xmrdfits(mcrec_fil,1,/silent) $
  else mcrecstr = mcrec_fil
  nmcrec = (size(mcrecstr.zabs,/dim))[0] 
  if nmcrec ne nmc then $
     stop,'sdss_completeness_estncolm: MC input and recovered structures mismatch',$
          nmc,nmcrec
  if not keyword_set(mcrecstr.rec_param) then $
     print,'sdss_completeness_estncolm: WARNING! not using recovered MC params.'

  if keyword_set(check) then begin
     ;; Going to want to check that estimated column and measured
     ;; column somewhat in agreement
     dblt = dblt_retrieve(civstr[0].wrest[0])
     ncolm = sdss_calcncolmdblt(civstr,dblt.ion,signcolm=signcolm,/log,$
                                ncolmflg=ncolmflg,/silent)
     ;; log(N+/-sigN) = logN + log(1 +/- ln(10)*siglogN)
     ncolmlo = ncolm + alog10((1 - alog(10)*signcolm) > 1.e-15) > 0.
     ncolmhi = ncolm + alog10(1 + alog(10)*signcolm)
     
     check_str = {qso_name:civstr.qso_name,$
                  zabs_orig:civstr.zabs_orig[0],$
                  ew_orig:civstr.ew_orig[0],$
                  sigew_orig:civstr.sigew_orig[0],$
                  ncolm_orig:ncolm,$
                  signcolm_orig:signcolm,$
                  ncolmflg:ncolmflg,$
                  num_mc:lonarr(nciv,/nozero),$
                  flg_fail:lonarr(nciv),$
                  ew_rec:fltarr(nciv,/nozero),$
                  sigew_rec:fltarr(nciv,/nozero),$
                  ncolm_rec:fltarr(nciv,/nozero),$
                  signcolm_rec:fltarr(nciv,/nozero) $
                 }
  endif 

  for zz=0,cmpltstr.nzbin-1 do begin
     sub = where(civstr.zabs_orig[0] ge cmpltstr.zlim[zz,0] and $
                 civstr.zabs_orig[0] lt cmpltstr.zlim[zz,1],nsub)
     if nsub eq 0 then begin
        print,'sdss_completeness_estncolm: no doublets in z bin ',$
              cmpltstr.zlim[zz,*]
        continue
     endif 
     subciv = civstr[sub]

     gdmc = where(mcrecstr.zabs ge cmpltstr.zlim[zz,0] and $
                  mcrecstr.zabs lt cmpltstr.zlim[zz,1],ngdmc)
     if ngdmc eq 0 then begin
        print,'sdss_completeness_estncolm: no MC profiles in z bin ',$
              cmpltstr.zlim[zz,*]
        ;; Take all
        gdmc = imcall
        ngdmc = nmc
     endif 
     
     for cc=0L,nsub-1 do begin
        nsigma = 1
        flg_rec = 1
        if keyword_set(incl_col) then $
           nlim = subciv[cc].ncolm_orig[0] + $
                  alog10(1 + $
                         [-1,1]*alog(10.)*subciv[cc].signcolm_orig[0]) $
        else nlim = [1,!values.f_infinity]
        mtch = where(mcrecstr.ew[gdmc] ge $
                     subciv[cc].ew_orig[0]-subciv[cc].sigew_orig[0] $
                     and mcrecstr.ew[gdmc] le $
                     subciv[cc].ew_orig[0]+subciv[cc].sigew_orig[0] $
                     and mcrecstr.flg_rec[gdmc] eq 1 and $
                     mcrecstr.ncolm[gdmc] ge nlim[0] and $
                     mcrecstr.ncolm[gdmc] le nlim[1],nmtch)

        if nmtch eq 0 then begin
           ;; Check 1
           print,'sdss_completeness_estncolm: no recovered profiles matching doublet ',$
                 subciv[cc].qso_name,subciv[cc].zabs_orig[0],$
                 subciv[cc].ew_orig[0]
           ;; Don't required recovered
           flg_rec = 0
           mtch = where(mcrecstr.ew[gdmc] ge $
                        subciv[cc].ew_orig[0]-subciv[cc].sigew_orig[0] $
                        and mcrecstr.ew[gdmc] le $
                        subciv[cc].ew_orig[0]+subciv[cc].sigew_orig[0] and $
                        mcrecstr.ncolm[gdmc] ge nlim[0] and $
                        mcrecstr.ncolm[gdmc] le nlim[1],nmtch)
        endif                   ; no match w/ recovered flag

        if nmtch eq 0 then begin
           ;; Check 2
           print,'sdss_completeness_estncolm: none of the MC profiles match'
           nsigma = 2
           flg_rec = 1
           for ss=0,1 do begin
              ;; Back to require recovered, but step outward in sigEW
              if keyword_set(incl_col) then $
                 nlim = subciv[cc].ncolm_orig[0] + $
                        alog10(1 + $
                               [-1,1]*nsigma*alog(10.)*subciv[cc].signcolm_orig[0]) 
              mtch = where(mcrecstr.ew[gdmc] ge $
                           subciv[cc].ew_orig[0]-nsigma*subciv[cc].sigew_orig[0] $
                           and mcrecstr.ew[gdmc] le $
                           subciv[cc].ew_orig[0]+nsigma*subciv[cc].sigew_orig[0] $
                           and mcrecstr.flg_rec[gdmc] eq 1 and $
                           mcrecstr.ncolm[gdmc] ge nlim[0] and $
                           mcrecstr.ncolm[gdmc] le nlim[1], nmtch)
              nsigma++
              if nmtch ne 0 then break
           endfor               ; loop ss=0,1 
           nsigma--
           print,'sdss_completeness_estncolm: include out to nsigEW = ',$
                 nsigma
        endif                   ; no match for any input


        if nmtch eq 0 then begin
           ;; Check 3
           print,'sdss_completeness_estncolm: none of the MC profiles match in EW and column'
           flg_rec = 1
           mtch = where(mcrecstr.ew[gdmc] ge $
                        subciv[cc].ew_orig[0]-subciv[cc].sigew_orig[0] $
                        and mcrecstr.ew[gdmc] le $
                        subciv[cc].ew_orig[0]+subciv[cc].sigew_orig[0] $
                        and mcrecstr.flg_rec[gdmc] eq 1,nmtch)
        endif

        if nmtch eq 0 then begin
           ;; Check 4
           stop
        endif
        

        ;; Take median and standard deviation of INPUT values as real
        ;; column density
        ;; Error estimate question:
        ;; stddev(logN) != stddev(10.^ncolm) / (alog(10)*median(10.^ncolm))
        mtch = gdmc[mtch]
        subciv[cc].ncolm_final[0:1] = median(mcinpstr.ncolm[mtch],/even)
        subciv[cc].signcolm_final[0:1] = stddev(mcinpstr.ncolm[mtch])
        subciv[cc].ncolmflg[0:1] = sdss_getlimflg() ; perfect measurement
        subciv[cc].notes = 'Nmc='+strtrim(nmtch,2)+$
                           '; NsigEW='+strtrim(nsigma,2)+$
                           '; flg_rec='+strtrim(flg_rec,2)

        if keyword_set(check) then begin
           ;; Check that ncolm_rec[mtch] and ncolm is reasonable, +/-3 sigma
           med = median(mcrecstr.ncolm[mtch],/even)
           stddev = stddev(mcrecstr.ncolm[mtch])
           if med lt ncolmlo[sub[cc]] or med gt ncolmhi[sub[cc]] then begin
;              print,'sdss_completeness_estnolm: check fail ',zz,cc,med,$
;                    ncolmlo[sub[cc]],ncolmhi[sub[cc]]
              check_str.flg_fail[sub[cc]] = 1
           endif 
           check_str.num_mc[sub[cc]] = nmtch
           check_str.ew_rec[sub[cc]] = median(mcrecstr.ew[mtch],/even)
           check_str.sigew_rec[sub[cc]] = stddev(mcrecstr.ew[mtch])
           check_str.ncolm_rec[sub[cc]] = med
           check_str.signcolm_rec[sub[cc]] = stddev
        endif 

     endfor                     ; loop cc=nsub

     ;; Put back
     civstr[sub] = subciv
  endfor                        ; loop zz=cmpltstr.nzin

  ;; Copy all over
  civstr.zabs_final = civstr.zabs_orig
  civstr.sigzabs_final = civstr.sigzabs_orig
  civstr.ew_final = civstr.ew_orig
  civstr.sigew_final = civstr.sigew_orig

  ;; Write out; already checked
  mwrfits,civstr,ostrct_fil,/create,/silent
  if keyword_set(check) then $
     mwrfits,check_str,ostrct_fil,/silent
  spawn,'gzip -f '+ostrct_fil
  print,'sdss_completeness_estncolm: created ',ostrct_fil

end                             ; sdss_completeness_estncolm


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_completeness_diff, cmplt1_fil, cmplt2_fil, redshift=redshift,$
                                    psfil=psfil,yttl=yttl,csize=csize,psize=psize,$
                                    lthick=lthick,xrng=xrng,_extra=extra
  if n_params() ne 2 then begin
     print,'Syntax - sdss_completeness_diff(cmplt1_fil,cmplt2_fil,[/redshift])'
     return,-1
  endif
  if not keyword_set(csize) then csize = 2.
  if not keyword_set(psize) then psize = 1.5
  if not keyword_set(lthick) then lthick = 4

  cmpltstr1 = xmrdfits(cmplt1_fil,1,/silent)
  cmpltstr2 = xmrdfits(cmplt2_fil,1,/silent)

  czw_diff = cmpltstr1.czw_2darr / cmpltstr2.czw_2darr
  if not keyword_set(yttl) then yttl = cmplt1_fil+' / '+cmplt2_fil 
  print,yttl

  yrange = [0.,2.]
  if keyword_set(redshift) then begin
     xdat = sdss_mkzarr(cmpltstr1.zlim,cmpltstr1.zbinsize)
     zarr = sdss_mkewarr(cmpltstr1.ewlim,cmpltstr1.ewbinsize)
     index = lindgen(cmpltstr1.nzbin)
     title = ['!8z!X!D1548!N','Completeness','!8W!X!Dr,1548!N']
     xrange = cmpltstr1.zlim+[-1.,1.]*cmpltstr1.zbinsize
     print,'EW','dCmin','@z','dCmax','@z','Mean','Median','Stddev',$
           format='(a5,2x,2(a7,2x,a5,2x),3(a7,2x))'
     nbin = cmpltstr1.newbin
  endif else begin
     xdat = sdss_mkewarr(cmpltstr1.ewlim,cmpltstr1.ewbinsize)
     zarr = sdss_mkzarr(cmpltstr1.zlim,cmpltstr1.zbinsize)
     index = lindgen(cmpltstr1.newbin)
     title = ['!8W!X!Dr,1548!N','Completeness','!8z!X!D1548!N']
     xrange = cmpltstr1.ewlim+[-1.,1.]*cmpltstr1.ewbinsize
     print,'z','dCmin','@EW','dCmax','@EW','Mean','Median','Stddev',$
           format='(a5,2x,2(a7,2x,a5,2x),3(a7,2x))'
     nbin = cmpltstr1.nzbin
  endelse 

  clr = getcolor(/load)
  window,0
  wset,0
  plot,[0],[0],/nodata,color=clr.black,thick=lthick,$
       background=clr.white,charsize=charsize,/xstyle,/ystyle,$
       xtitle=title[0],ytitle=title[1],psym=10,$
       xrange=xrange,yrange=yrange
  stats = {cmplt1_fil:cmplt1_fil,$
           cmplt2_fil:cmplt2_fil, $
           redshift:keyword_set(redshift), $
           zdat:zarr, $
           ymnx:fltarr(nbin,2,/nozero), $
           xmnx:fltarr(nbin,2,/nozero), $
           ymed:fltarr(nbin,2,/nozero), $   ; value and abs(value)
           ystddev:fltarr(nbin,2,/nozero), $ ; value and abs(value)
           ymean:fltarr(nbin,2,/nozero), $
           fitpar:fltarr(2,/nozero), $
           covpar:fltarr(2,2,/nozero) $
          }
  for zz=0,nbin-1 do begin
     if keyword_set(redshift) then $
        ydat = czw_diff[zz*cmpltstr1.nzbin+index] $
     else ydat = czw_diff[zz+cmpltstr1.nzbin*index]

     oplot,xdat,ydat,thick=lthick,psym=10,color=clr.gray,linestyle=(zz mod 7)

     gd = where(finite(ydat) eq 1)
     if gd[0] eq -1 then begin
        stats.ymnx[zz,*] = !values.f_nan
        stats.xmnx[zz,*] = !values.f_nan
        stats.ymed[zz,*] = !values.f_nan
        stats.ystddev[zz,*] = !values.f_nan
        stats.ymean[zz] = !values.f_nan
     endif else begin
        ymn = min(ydat[gd],imn,max=ymx,subscript_max=imx)
        
        stats.ymnx[zz,*] = [ymn,ymx] 
        stats.xmnx[zz,*] = [xdat[gd[imn]],xdat[gd[imx]]]
        stats.ymed[zz,*] = [median(ydat[gd],/even),median(abs(ydat[gd]),/even)]
        stats.ystddev[zz,*] = [stddev(ydat[gd]),stddev(abs(ydat[gd]))]
        stats.ymean[zz] = mean(ydat[gd])
     endelse 

     print,stats.zdat[zz],stats.ymnx[zz,0],stats.xmnx[zz,0],stats.ymnx[zz,1],$
           stats.xmnx[zz,1],stats.ymean[zz],stats.ymed[zz,0],stats.ystddev[zz,0],$
           format='(f5.2,2x,2(f7.4,2x,f5.2,2x),3(f7.4,2x))'
  endfor                        ; zz=cmpltstr1.nzbin

  gd = where(finite(stats.ymed[*,0]) eq 1)
  stats.fitpar = linfit(stats.zdat[gd],stats.ymed[gd,0],$
                        measure_errors=stats.ystddev[gd,0],$
                        sigma=sigrslt,covar=covrslt)
  stats.covpar = covrslt
  printcol,stats.fitpar,sigrslt,stats.covpar[0,*],stats.covpar[1,*]
  
  if keyword_set(psfil) then begin
     ;; _extra= includes /encaps
     x_psopen,psfil,/maxs,/portrait,_extra=extra 
     !p.multi = [1,1,1]
     !x.margin = [6.8,0.5]
     !y.margin = [3.1,2.4]
     clr = getcolor(/load)
  endif else begin
     window, 1
     wset, 1
  endelse 
  if not keyword_set(xrng) then $
     xrng = [stats.zdat[0]*0.5,stats.zdat[nbin-1]*1.01]
  mnx = [min(stats.ymed[gd,0]-stats.ystddev[gd,0]),$
         max(stats.ymed[gd,0]+stats.ystddev[gd,0])]
  plot,[0],[0],/nodata,color=clr.black,background=clr.white,$
       thick=lthick,charsize=csize,/ystyle,/xstyle,xtitle=title[2],$
       ytitle=yttl,yrange=mnx,xrange=xrng
  oploterror,stats.zdat,stats.ymed[*,0],stats.ystddev[*,0],$
             psym=4,color=clr.black,errcolor=clr.black,symsize=psize
  fit = stats.fitpar[0]+stats.fitpar[1]*stats.zdat
  oplot,stats.zdat,fit,color=clr.red,thick=lthick
  varfit = stats.covpar[0,0] + stats.zdat^2*stats.covpar[1,1] + $
           2 * stats.zdat * stats.covpar[0,1]
  oplot,stats.zdat,fit+sqrt(varfit),thick=0.5*lthick,$
        linestyle=1,color=clr.red
  oplot,stats.zdat,fit-sqrt(varfit),thick=0.5*lthick,$
        linestyle=1,color=clr.red

  if keyword_set(psfil) then begin
     x_psclose
     print,'sdss_completeness_diff(): created ',psfil
  endif

  return,stats
end                             ; sdss_completeness_diff()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_user, list_fil, out_fil, dblt_name=dblt_name, clobber=clobber,$
                            _extra=extra
  if n_params() ne 2 then begin
     print,'Syntax - sdss_completeness_user, list_fil, out_fil, [dblt_name=, /clobber,'
     print,'                            _extra=]'
  endif

  sdssdir = sdss_getsdssdir()

  test = file_search(out_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_completeness_user: file exists; will not clobber ',out_fil
     return
  endif 

  readcol,list_fil,spec_fil,format='a',skip=1
  nfil = (size(spec_fil,/dim))[0]

  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  ;; _extra= includes /quick
  mcstrct_fil = sdss_getname(spec_fil,/spec,mc=dblt.ion,dir=mcdir,$
                             /user,_extra=extra) ; all will have *-001.fit
  mcstrct_fil = mcdir + mcstrct_fil

  for ff=0,nfil-1 do begin

     ;; Read in extensions and store
     ;; recovered candidate structure (ext = 1), not
     ;; recovered candidate structure (ext = 2), MC structure
     ;; (ext = 3)
     dat = xmrdfits(sdssdir+mcstrct_fil[ff],1,/silent) ; recovered
     if dat[0].wrest[0] ne 0. then begin
        if keyword_set(civstr_all) then civstr_all = [civstr_all,dat] $
        else civstr_all = dat
     endif 

     dat = xmrdfits(sdssdir+mcstrct_fil[ff],2,/silent) ; not recovered
     if dat[0].wrest[0] ne 0. then begin
        if keyword_set(civstr_rej) then civstr_rej = [civstr_rej,dat] $
        else civstr_rej = dat
     endif 

     dat = xmrdfits(sdssdir+mcstrct_fil[ff],3,/silent) ; MC structure
     if keyword_set(mcstrct) then mcstrct = [mcstrct,dat] $
     else mcstrct = dat

  endfor    ; loop ff=nfil

  ;; Everything must exist
  if not keyword_set(civstr_all) then civstr_all = [0]
  if not keyword_set(civstr_rej) then civstr_rej = [0]

  ;; Write output (clobber check happened at beginning of code)
  mwrfits,civstr_all,out_fil,/create,/silent ; ext = 1
  mwrfits,civstr_rej,out_fil,/silent         ; ext = 2
  mwrfits,mcstrct,out_fil,/silent            ; ext = 3
  spawn,'gzip -f '+out_fil
  print,'sdss_completeness_user: created ',out_fil

end                             ; sdss_completeness_user


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_completeness_userbias, user_fil, userrate_fil, out_fil, $
                                dblt_name=dblt_name, zlim=zlim,$
                                n_per_bin=n_per_bin, cmplt_fil=cmplt_fil, nuser=nuser, $
                                ewlim=ewlim, nrand=nrand, seed=seed,$
                                oseed=oseed, noplot=noplot, $
                                debug=debug, silent=silent, $
                                clobber=clobber, dvqso=dvqso, dvem=dvem, $
                                stats=stats, _extra=extra
  ;; Analyze results of User ratings of doublets
  if n_params() ne 3 then begin
     print,'Syntax - sdss_completeness_userbias, user_fil, userrate_fil, out_fil, '
     print,'                                [dblt_name=, zlim=,n_per_bin=, '
     print,'                                ewlim=, nrand=, seed=,oseed=, /noplot, '
     print,'                                /debug, /silent,/clobber, dvem=,'
     print,'                                cmplt_fil=, dvqso=, /stats, _extra=]'
     return
  endif 

  ;; Several other default parameters are set below (after
  ;; restore_outfil:) 
  if not keyword_set(dvqso) then dvqso = -3000. ; km/s default
  c = 299792.458
  if keyword_set(dvem) then begin
     n_em = size(dvem,/n_dim)
     if n_em ne 1 then $ ; [1,3]
        stop,'sdss_completeness_userbias stop: cannot handle dvem that is not [wrest, -dvem, +dvem]'
  endif 
  if not keyword_set(dblt_name) then dblt_name = 'FkIV'
  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  if not keyword_set(n_per_bin) then n_per_bin = 100 
  if keyword_set(zlim) then begin
     if n_elements(zlim) eq 2 then begin
        ;; want consistent form
        zbin = transpose(zlim)  ; size() = [1,2]
     endif else zbin = zlim
  endif else begin
     ;; Default is at 7000 Ang, where sky lines kick in and full rng
     tmp = 7000./dblt.wvI-1
     sdssswv = sdss_getspecwave()
     mnx = (sdsswv/dblt.wvI-1)
     zbin = transpose([mnx,[mnx[0],tmp],[tmp,mnx[1]]])
  endelse
  nzbin = n_elements(zbin[*,0])

  ;; For plotting
  clr = getcolor(/load)
  afp_clr = clr.red
  afp_psym = 4                  ; diamond
  fpr_clr = clr.limegreen
  fpr_psym = 5                  ; triangele
  tp_clr = clr.black
  tp_psym = 7                   ; Cross
  gdrtg = sdss_getrating(/good)
  wreststring = strtrim(floor(dblt.wvI),2)
  ewstring = '!8W!X!Dr,'+wreststring+'!N'

  ;; Test output file
  test = file_search(out_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) and not keyword_set(stats) then begin
     print,'sdss_completeness_userbias: will not clobber existing file; skip to analysis ',out_fil
     goto,restore_outfil
  endif 

  ;; Read in information
  civrej = xmrdfits(user_fil,2,/silent) ; rejected elements (sanity check)
  civrej.rating[0] = sdss_getrating(/unrated)   
  ncivrej = (size(civrej,/dim))[0]
  
  mcstr = xmrdfits(user_fil,3,/silent) ; raw mcstrct; apply dvqso later
  nmcstr = (size(mcstr,/dim))[0]
  nmcciv = total(mcstr.nsys)
  nsys_max = max(mcstr.nsys)

  
  ;; Combine duplicates in the rated file and re-measure values (by
  ;; rating groups)
  ;; dvqso, dvem cut applied later 
  civrate = xmrdfits(userrate_fil,1,/silent) ; rated elements

  civaccept = sdss_getcivstrct(civrate,rating=[gdrtg,sdss_getrating(/def)],$
                               count=ncivaccept)
  ;; _extra includes dvsys=
  nwcivaccept = sdss_combineciv(civaccept,count=nnwcivacccept,debug=0,$
                                dvqso=dvqso,_extra=extra)
  civrateno = $
     sdss_getcivstrct(civrate,$
                      rating=[sdss_getrating(/bad),sdss_getrating(/maybe)],$
                      count=ncivrateno)
  nwcivrateno = sdss_combineciv(civrateno,count=nnwcivrateno,debug=0,$
                                _extra=extra)
  sdss_ewciv_strct, [nwcivaccept,nwcivrateno], 'nwcivrate.fit', $ ; SLOW!
                    dblt_name=dblt, debug=0, _extra=extra
  civrate = sdss_srtcivstrct('nwcivrate.fit')
  ncivrate = (size(civrate,/dim))[0]
  spawn,'rm nwcivrate.fit*'

  ;; Set up brute force matching
  mask_civrate = replicate(-1,ncivrate,2)
  mask_civrej = replicate(-1,ncivrej,2)
  
  for ii=0L,nmcstr-1 do begin
     ;; Save information for analysis
     dvqso_mc = (mcstr[ii].zabs_rec[0:mcstr[ii].nsys-1,0]-$
                 mcstr[ii].z_qso)/(1+mcstr[ii].z_qso)*c
     if keyword_set(dvem) then begin
        zem_mc = dvem[0]/dblt.wvI*(1+mcstr[ii].z_qso) - 1.
        dvem_mc = c*(mcstr[ii].zabs_rec[0:mcstr[ii].nsys-1,0]-zem_mc)/(1+zem_mc)
     endif
     if ii eq 0 then begin
        mcflg_rec = mcstr[ii].flg_rec[0:mcstr[ii].nsys-1] 
        mcdvqso = dvqso_mc
        if keyword_set(dvem) then mcdvem = dvem_mc
     endif else begin
        mcflg_rec = [mcflg_rec,mcstr[ii].flg_rec[0:mcstr[ii].nsys-1]]
        mcdvqso = [mcdvqso, dvqso_mc]
        if keyword_set(dvem) then mcdvem = [mcdvem, dvem_mc]
     endelse

     ;; Loop over systems per sightline
     for jj=0,mcstr[ii].nsys-1 do begin

        if mcstr[ii].flg_rec[jj] eq 1 then begin
           ;; Search in the rated systems (allow
           ;; duplicates)... /srch_both uses both mcstr and civstr
           ;; wvlim for combining
           nwmcstr = $
              sdss_completeness_compare(mcstr[ii], civrate, $
                                        /ignore_flg, allmtch=allmtch,$
                                        /srch_both, debug=0)
           mtch = where(allmtch[*,0] ge 0 and allmtch[*,1] ge 0)
           if mtch[0] ne -1 then begin
              ;; May be overwriting other detections
              mask_civrate[mtch,0] = ii              ; mcstr element
              mask_civrate[mtch,1] = jj              ; system within that
              civrate[mtch].wrest[9] = mcstr[ii].snr ; just save <S/N> 
           endif else begin
              print,'sdss_completeness_userbias: no match for recovered mcstr ',$
                    mcstr[ii].qso_name,mcstr[ii].zabs_input[jj,0],mcstr[ii].ew_input[jj,0],$
                    format='(a,1x,a15,1x,f7.5,1x,f6.2)'
              if keyword_set(debug) then stop
           endelse
        endif else begin
           ;; Else search in the rejected systems (everything should be
           ;; found here; allows duplicates)
           nwmcstr = $
              sdss_completeness_compare(mcstr[ii], civrej, $
                                        /ignore_flg, allmtch=allmtch,$
                                        /srch_both, debug=0)
           mtch = where(allmtch[*,0] ge 0 and allmtch[*,1] ge 0)
           if mtch[0] ne -1 then begin
              ;; May be overwriting other detections
              mask_civrej[mtch,0] = ii              ; mcstr element
              mask_civrej[mtch,1] = jj              ; system within that
              civrej[mtch].wrest[9] = mcstr[ii].snr ; just save <S/N>
              if civrej[mtch[0]].z_qso eq 0. then $
                 civrej[mtch].z_qso = mcstr[ii].z_qso ; previous typo, may become obsolete
           endif else $
              ;; Must have matches
              stop,'sdss_completeness_userbias stop: no match for rejected mcstr ',ii,jj
        endelse 
     endfor                     ; loop jj=mcstr[ii].nsys
  endfor                        ; loop ii=nmcstr


  ;; ;;;;;;;
  ;; Apply dvqso cut
  mcsub = where(mcdvqso lt dvqso,nmcciv)
  mcflg_rec = mcflg_rec[mcsub]

  gdrate = where((civrate.zabs_orig[0]-civrate.z_qso)/(1+civrate.z_qso)*c lt dvqso, $
                 ncivrate)
  civrate = civrate[gdrate]
  mask_civrate = mask_civrate[gdrate,*]

  gdrej = where((civrej.zabs_orig[0]-civrej.z_qso)/(1+civrej.z_qso)*c lt dvqso,$
             ncivrej)
  civrej = civrej[gdrej]
  mask_civrej = mask_civrej[gdrej,*]

 
  ;; ;;;;;;;
  ;; Apply dvem cut
  if keyword_set(dvem) then begin
     mcsub = where(mcdvem lt dvem[1] or mcdvem ge dvem[2],nmcciv)
     mcflg_rec = mcflg_rec[mcsub]
     
     zem_ion = dvem[0]/dblt.wvI*(1+civrate.z_qso) - 1.
     dvem_ion = c*(civrate.zabs_orig[0]-zem_ion)/(1+zem_ion)
     gdrate = where(dvem_ion lt dvem[1] or dvem_ion ge dvem[2], ncivrate)
     civrate = civrate[gdrate]
     mask_civrate = mask_civrate[gdrate,*]
     
     zem_ion = dvem[0]/dblt.wvI*(1+civrej.z_qso) - 1.
     dvem_ion = c*(civrej.zabs_orig[0]-zem_ion)/(1+zem_ion)
     gdrej = where(dvem_ion lt dvem[1] or dvem_ion ge dvem[2], ncivrej)
     civrej = civrej[gdrej]
     mask_civrej = mask_civrej[gdrej,*]     
  endif                         ; dvem=


  ;; ;;;;;;;
  ;; Report
  rec = where(mcflg_rec eq 1,nrec,complement=rej,ncomplement=nrej)

  print,''
  print,''
  ;; The not-recovered ones are lost forever and taken care of by
  ;; the basic completeness tests
  print,'Total input profiles < dvqso: ',string(nmcciv,format='(i5)')
  print,'Total number automatically recovered (not): ',$
        string(nrec,nrej,format='(i5,1x,"(",i5,")")')

  ;; The candidates *not* associated with an input system are the
  ;; accepted false positives
  trurec_indx = where(mask_civrate[*,0] ge 0,ntrurec_indx,$
                      complement=fpr_indx,ncomplement=nfpr_indx)
  print,'Total '+dblt.ion+' candidates: ',string(ncivrate,format='(i5)')
  print,'Number of recovered '+dblt.ion+' matched up (not): ',$
        string(ntrurec_indx,nfpr_indx,format='(i5,1x,"(",i5,")")')
;  sdss_chkciv,civrate[fpr_indx],strlowcase(dblt.ion)+'_fpr_indx.fit',dblt_name=dblt.ion

  ;; Check ratings
  if ntrurec_indx ne 0 then begin
     ;; Which true ones were rated as not real? (false negatives)
     tp_indx = where(civrate[trurec_indx].rating[0] ge gdrtg,ntp_indx,$
                     complement=fn_indx,ncomplement=nfn_indx)
     print,'    True positives (false negatives): ',$
           string(ntp_indx,nfn_indx,format='(i5,1x,"(",i5,")")')
  endif
  if nfpr_indx ne 0 then begin
     ;; Which false ones were rated as real? (accepted false positives)
     afp_indx = where(civrate[fpr_indx].rating[0] ge gdrtg,nafp_indx,$
                      complement=tn_indx,ncomplement=ntn_indx)
     print,'    True negatives (accepted false positives): ',$
           string(ntn_indx,nafp_indx,format='(i5,1x,"(",i5,")")')
  endif

  ;; All the rejected ones should be matched up because the only
  ;; reason they exist is because we looked in there...
  more_input_indx = where(mask_civrej[*,0] ge 0,nmore_input_indx,$
                          complement=unknown_indx,$
                          ncomplement=nunknown_indx)
  print,'Total rejected '+dblt.ion+' profiles: ',$
        string(ncivrej,format='(i5)')
  print,'Number of rejected '+dblt.ion+' matched up (not): ',$
        string(nmore_input_indx,nunknown_indx,format='(i5,1x,"(",i5,")")')

  print,''
  print,''

  if keyword_set(stats) then return ; EXIT


  ;; ;;;;;;; 
  ;; Output: concatenate results and save only what have to so
  ;; input parameters can change
  bd = where(civrate.wrest[9] eq 0.,nbd) ; set S/N for everything
  if nbd ne 0 then begin
     unq = uniq(civrate[bd].qso_name,sort(civrate[bd].qso_name))
     unq = bd[unq]
     nunq = (size(unq,/dim))[0] > 1
     for bb=0,nunq-1 do begin
        sub = bd[where(civrate[unq[bb]].qso_name eq $
                       civrate[bd].qso_name)]
        mtch = where(civrate[unq[bb]].qso_name eq $
                     mcstr.qso_name,nmtch)
        if nmtch ne 0 then $
           civrate[sub].wrest[9] = mcstr[mtch[0]].snr $
        else stop
     endfor                     ; loop bb=nbd
  endif                         ; instantiate spec S/N in wrest[9]

  trurecstr = civrate[trurec_indx] 
  mask_trurecstr = mask_civrate[trurec_indx,*]

  falsestr = civrate[fpr_indx]
  mask_falsestr = mask_civrate[fpr_indx,*]


  ;; ;;;;;;;
  ;; Begin serious analysis
  restore_outfil:
  if ntest ne 0 and not keyword_set(clobber) then begin
     print,'sdss_completeness_userbias: restoring save file ',out_fil
     ostrct = xmrdfits(out_fil,1,/silent)
     trurecstr = xmrdfits(out_fil,2,/silent)
     falsestr = xmrdfits(out_fil,3,/silent)
     ;; Just restore the basics and let new input paramters overwrite
     ;; things 
     mask_trurecstr = ostrct.mask_trurecstr
     mask_falsestr = ostrct.mask_falsestr
     tp_indx = ostrct.tp_indx
     afp_indx = ostrct.afp_indx
  endif 

  if keyword_set(debug) then begin
     ;; ;;;;;;;
     ;; Make and plot a lot of diagnostic histograms
     ;; change the xttl2d, yttl2d commented lines to plot different
     ;; param space with contour
     xttl = ''
     yttl = ''
     tags = tag_names(trurecstr)
     for ii=0,5 do begin
        yttl = 'Fraction'
        if ii le 2 then begin
           case ii of 
              0: begin          ; EW
                 itag = (where(tags eq 'EW_ORIG'))[0]
                 rng = [0.,3.5,0.1] ; [min,max,delta]
                 xttl =  ewstring
              end
              1: begin          ; z1548
                 itag = (where(tags eq 'ZABS_ORIG'))[0]
                 rng = [1.4,4.6,0.1] ; [min,max,delta]
                 xttl = '!8z!X!D'+wreststring+'!N'
              end
              2: begin          ; local S/N
                 itag = (where(tags eq 'GWIDTH'))[0]
                 rng = [0.,20.,0.5] ; [min,max,delta]
                 xttl = 'Centroid S/N'
              end
           endcase 
           hist = sdss_histogram(trurecstr.(itag)[0],n_per_bin,loc=loc,$
                                 szloc=dloctp,min=rng[0],max=rng[1])
           histtp = sdss_bintoloc(trurecstr[tp_indx].(itag)[0],loc,dloctp)
           histfpr = sdss_histogram(falsestr.(itag)[0],n_per_bin,$
                                    loc=locfpr,szloc=dlocfpr,$
                                    min=rng[0],max=rng[1])
           histafp = sdss_bintoloc(falsestr[afp_indx].(itag)[0],locfpr,dlocfpr)
        endif else begin
           case ii of
              3: begin
                 ;; Wavelength bounds
                 dvlim_allin = 2.998e5*(trurecstr.wvlim_orig[0,1] - $
                                        trurecstr.wvlim_orig[0,0])/trurecstr.wrest[0]
                 rng = [600.,10000,250.]
                 xttl = '!9d!X!8v!X!Dlim!N (km s!E-1!N)'
                 hist = sdss_histogram(dvlim_allin,n_per_bin,loc=loc,$
                                       szloc=dloctp,min=rng[0],max=rng[1])
                 histtp = sdss_bintoloc(dvlim_allin[tp_indx],loc,dloctp)

                 dvlim_fp = 2.998e5*(falsestr.wvlim_orig[0,1] - $
                                     falsestr.wvlim_orig[0,0])/falsestr.wrest[0]
                 histfpr = sdss_histogram(dvlim_fp,n_per_bin,$
                                          loc=locfpr,szloc=dlocfpr,$
                                          min=rng[0],max=rng[1])
                 histafp = sdss_bintoloc(dvlim_fp[afp_indx],locfpr,dlocfpr)
              end
              4: begin
                 ;; Just stored this in a random place above b/c mcstr
                 ;; had the information
                 rng = [3.,40.,0.5]
                 xttl = '<S/N>'
                 hist = sdss_histogram(trurecstr.wrest[9],n_per_bin,loc=loc,$
                                       szloc=dloctp,min=rng[0],max=rng[1])
                 histtp = sdss_bintoloc(trurecstr[tp_indx].wrest[9],loc,dloctp)

                 histfpr = sdss_histogram(falsestr.wrest[9],n_per_bin,$
                                          loc=locfpr,szloc=dlocfpr,$
                                          min=rng[0],max=rng[1])
                 histafp = sdss_bintoloc(falsestr[afp_indx].wrest[9],locfpr,dlocfpr)
              end
              5: begin
                 ;; EW/sigEW
                 rng = [0.,30.,1.0]
                 xttl = ewstring+'/!9s!X!D!8W!X!N'
                 newsig_allin = trurecstr.ew_orig[0]/trurecstr.sigew_orig[0]
                 bd = where(trurecstr.sigew_orig[0] eq 0.)
                 if bd[0] ne -1 then newsig_allin[bd] = 0.
                 hist = sdss_histogram(newsig_allin,n_per_bin,loc=loc,$
                                       szloc=dloctp,min=rng[0],max=rng[1])
                 histtp = sdss_bintoloc(newsig_allin[tp_indx],loc,dloctp)

                 newsig_fp = falsestr.ew_orig[0]/falsestr.sigew_orig[0]
                 bd = where(falsestr.sigew_orig[0] eq 0.)
                 if bd[0] ne -1 then newsig_fp[bd] = 0.
                 histfpr = sdss_histogram(newsig_fp,n_per_bin,$
                                          loc=locfpr,szloc=dlocfpr,$
                                          min=rng[0],max=rng[1])
                 histafp = sdss_bintoloc(newsig_fp[afp_indx],locfpr,dlocfpr)
              end
           endcase
        endelse 

        ;; _extra includes sigma=, cl=, /verbose
        loctp = loc + 0.5*dloctp
        hist = float(hist)

        locfpr = locfpr + 0.5*dlocfpr
        histfpr = float(histfpr)

        fpr_frac = histfpr/dlocfpr ; normalize
        mx = max(fpr_frac)
        fpr_frac = fpr_frac/mx
        sigfpr_frac = sdss_calcsigpoiss(histfpr, _extra=extra)
        sigfpr_frac[*,0] = sigfpr_frac[*,0]/(mx*dlocfpr)
        sigfpr_frac[*,1] = sigfpr_frac[*,1]/(mx*dlocfpr)
        afp_frac = histafp/histfpr
        sigafp_frac = sdss_calcsigbinom(histfpr, histafp, _extra=extra)
        tp_frac = histtp/hist
        sigtp_frac = sdss_calcsigbinom(hist, histtp, _extra=extra)

        
        ;; Plot it up
        window, ii              ; new window
        wset, ii                ; plot to it
        
        plot,[0],[0],/nodata,/xsty,/ysty,color=clr.black,$
             background=clr.white,xrange=rng[0:1],yrange=[0.,1.],$
             xtitle=xttl,ytitle=yttl,charsize=2.,thick=2.

        ;; False positive rate (wrong things) = green triangles
        ;; shifted for plotting purposes
        oploterror,locfpr*1.15,fpr_frac,0.5*dlocfpr,sigfpr_frac[*,0],$
                   /lobar,color=fpr_clr,$
                   errcolor=fpr_clr,/nohat,psym=fpr_psym,symsize=2,thick=2
        oploterror,locfpr*1.15,fpr_frac,0.5*dlocfpr,sigfpr_frac[*,1],$
                   /hibar,color=fpr_clr,$
                   errcolor=fpr_clr,/nohat,psym=fpr_psym,symsize=2,thick=2

        ;; Accepted False positive (wrong things) = red diamonds
        ;; shifted for plotting purposes
        oploterror,locfpr*1.05,afp_frac,0.5*dlocfpr,sigafp_frac[*,0],$
                   /lobar,color=afp_clr,$
                   errcolor=afp_clr,/nohat,psym=afp_psym,symsize=2,thick=2
        oploterror,locfpr*1.05,afp_frac,0.5*dlocfpr,sigafp_frac[*,1],$
                   /hibar,color=afp_clr,$
                   errcolor=afp_clr,/nohat,psym=afp_psym,symsize=2,thick=2

        ;; False negatives (missed things) = black crosses
        oploterror,loctp,tp_frac,0.5*dloctp,sigtp_frac[*,0],$
                   /lobar,color=tp_clr,$
                   errcolor=tp_clr,/nohat,psym=tp_psym,symsize=2,thick=2
        oploterror,loctp,tp_frac,0.5*dloctp,sigtp_frac[*,1],$
                   /hibar,color=tp_clr,$
                   errcolor=tp_clr,/nohat,psym=tp_psym,symsize=2,thick=2

        legend,['True Positive','False Positive Distribution','Accepted False Positives'],$
               psym=[tp_psym,fpr_psym,afp_psym],$
               colors=[tp_clr,fpr_clr,afp_clr],textcolors=[tp_clr,fpr_clr,afp_clr],$
               /top

     endfor                     ; loop ii=0,3
  endif else ii = 0             ; end debug



  ;; ;;;;;;;
  ;; Try fitting exponential growth
  if not keyword_set(ewlim) then ewlim = [0.05,1.0] ; [lo, hi]
  if not keyword_set(nrand) then nrand = 1000L
  if keyword_set(seed) then oseed = seed
  nfitpar = 3

  ;; Temporary arrays
  coeff_arr = fltarr(nrand,nfitpar,2,/nozero) 
  coeff_false = fltarr(nfitpar,3,/nozero)
  
  coeff_tp = fltarr(nzbin,nfitpar,3,/nozero) ; [z,[a0,low,high]]
  coeff_afp = fltarr(nzbin,nfitpar,3,/nozero)
  for zz=0,nzbin-1 do begin

     for ss=0,1 do begin
        if ss eq 0 then civstr = trurecstr $ ; to test false negatives
        else civstr = falsestr                  ; to test false positives

        ;; Select out subsamples
        sub = where(civstr.zabs_orig[0] ge zbin[zz,0] and $
                    civstr.zabs_orig[0] lt zbin[zz,1] and $
                    civstr.ew_orig[0] ge ewlim[0] and $
                    civstr.ew_orig[0] lt ewlim[1],nsubcivstr)
        subcivstr = civstr[sub]
        
        ;; Find the false or true samples
        ;; ss=0: true positive: Matched up so real but given low
        ;; rating
        ;; mask should be >= 0
        false = where(subcivstr.rating[0] ge gdrtg,nfalse) ; not civrej
        ;; ss=1: Accepted False positive: Not matched up so not real but
        ;; given high rating 
        ;; Mask should be < 0


        ;; Monte Carlo errors by drawing from Gaussian distribution and
        ;; scrambling the fractions by their errors
        !quiet = 1              ; suppress gaussfit() messages
        for rr=0,nrand-1 do begin
           ;; Draw random numbers from normal distribution
           if rr eq 0 then num = fltarr(nsubcivstr) $ ; basic
           else num = randomn(oseed,nsubcivstr)
           ewdat = subcivstr.ew_orig[0] + num*subcivstr.sigew_orig[0]
           
           ;; bin things up
           ;; Reference is of all in EW and zabs bin
           hist = sdss_histogram(ewdat,n_per_bin,loc=loc,szloc=dloc,$
                                 min=ewlim[0],max=ewlim[1])

           hist_false = sdss_bintoloc(ewdat[false],loc,dloc)
           
           ;; Set up the fractions, including asymmetric, binomial
           ;; errors
           nloc = (size(loc,/dim))[0]
           hist = float(hist)      ; must be float
           false_frac = hist_false/hist
           sigfalse_frac = sdss_calcsigbinom(hist, hist_false, _extra=extra)
           
           loc = loc + 0.5*dloc

           ;; Draw a random sign to know which asymmetric error to use
           num = randomn(oseed,nloc)
           pos = where(num gt 0.)
           sig = sigfalse_frac[*,0]                             ; lower
           if pos[0] ne -1 then sig[pos] = sigfalse_frac[pos,1] ; upper

           ;; Initial guess
           coeff = [max(false_frac,subscript_min=imn), -5.]
           if nfitpar eq 3 then coeff = [coeff, loc[imn]]
           rslt = curvefit(loc, false_frac, 1./sig^2, coeff, sigcoeff, $
                           function_name='sdss_expgrowthfunc')

           coeff_arr[rr,*,0] = coeff
           coeff_arr[rr,*,1] = sigcoeff

           if rr eq 0 then begin
              ;; Save the basics
              if ss eq 0 then begin
                 tp_loc = loc
                 tp_dloc = dloc
                 tp_frac = false_frac
                 sigtp_frac = sigfalse_frac
              endif else begin
                 afp_loc = loc
                 afp_dloc = dloc
                 afp_frac = false_frac
                 sigafp_frac = sigfalse_frac                 
              endelse 
           endif 
           
        endfor                  ; loop rr=nrand
        !quiet = 0              ; restore standard IDL messages
        
        ;; Combine the random coefficients and estimate errors by the
        ;; cumulative distribution and the percentiles spanning 1
        ;; sigma
        ;; And added the standard deviation business because the
        ;; false-positive rate best-fit parameter space looks like the
        ;; exponential growth curve itself and just goes crazy
        for jj=0,nfitpar-1 do begin
           gd = where(coeff_arr[*,jj,1] ne 0. and $
                      finite(coeff_arr[*,jj,0]) and $
                      finite(coeff_arr[*,jj,1]),ngd)
           coeff_false[jj,0] = median(coeff_arr[gd,jj,0])

           srt = gd[sort(coeff_arr[gd,jj,0])]
           cumrand = findgen(ngd)/float(ngd)             ; for Gauss fit error est.
           lo = min(cumrand-(1.-gauss_pdf(1.)),ilo,/abs) ; ~15.9%
           hi = min(cumrand-gauss_pdf(1.),ihi,/abs)      ; ~84.1%
           md = min(cumrand-0.5,imd,/abs)                ; 50%
           coeff_false[jj,1] = stddev(coeff_arr[srt[0:imd],jj,0]) < $
                               (coeff_false[jj,0] - coeff_arr[srt[ilo],jj,0])
           coeff_false[jj,2] = stddev(coeff_arr[srt[imd:*],jj,0]) < $
                              (coeff_arr[srt[ihi],jj,0] - coeff_false[jj,0])    
        endfor                  ; loop jj=nfitpar

        if keyword_set(debug) then begin
           ;; Look at error
           x_splot,coeff_arr[*,0,0],coeff_arr[*,1,0],psym1=4,$
                   xtwo=coeff_false[0,0]+[-1.,1]*coeff_false[0,1:2],$
                   ytwo=replicate(coeff_false[1,0],2),psym2=10,$
                   xthr=replicate(coeff_false[0,0],2),$
                   ythr=coeff_false[1,0]+[-1.,1]*coeff_false[1,1:2],$
                   psym3=10,xmnx=[0.,4],ymnx=[-10.,5.],$
                   xtitle='C0',ytitle='C1',$
                   title=string(zbin[zz,0],zbin[zz,1],$
                                format='(f4.2,"<=z<",f4.2)'),/block
           if nfitpar eq 3 then begin
              x_splot,coeff_arr[*,1,0],coeff_arr[*,2,0],psym1=4,$
                      xtwo=coeff_false[1,0]+[-1.,1]*coeff_false[1,1:2],$
                      ytwo=replicate(coeff_false[2,0],2),psym2=10,$
                      xthr=replicate(coeff_false[1,0],2),$
                      ythr=coeff_false[2,0]+[-1.,1]*coeff_false[2,1:2],$
                      psym3=10,xmnx=[-10.,5.],ymnx=[-1.,1.],$
                      xtitle='C1',ytitle='C2',$
                      title=string(zbin[zz,0],zbin[zz,1],$
                                   format='(f4.2,"<=z<",f4.2)'), /block
              x_splot,coeff_arr[*,0,0],coeff_arr[*,2,0],psym1=4,$
                      xtwo=coeff_false[0,0]+[-1.,1]*coeff_false[0,1:2],$
                      ytwo=replicate(coeff_false[2,0],2),psym2=10,$
                      xthr=replicate(coeff_false[0,0],2),$
                      ythr=coeff_false[2,0]+[-1.,1]*coeff_false[2,1:2],$
                      psym3=10,xmnx=[0.,4.],ymnx=[-1.,1.],$
                      xtitle='C0',ytitle='C2',$
                      title=string(zbin[zz,0],zbin[zz,1],$
                                   format='(f4.2,"<=z<",f4.2)'), /block
           endif
        endif 

        ;; Save
        if ss eq 0 then coeff_tp[zz,*,*] = coeff_false $
        else coeff_afp[zz,*,*] = coeff_false

     endfor                     ; loop ss=0,1
     
     ;; Set up fit plot
     fit_afp = sdss_calcexpgrowthfunc(afp_loc,reform(coeff_afp[zz,*,*]), $
                                      sigma=sigfit_afp)
     fit_tp = sdss_calcexpgrowthfunc(tp_loc,reform(coeff_tp[zz,*,*]), $
                                     sigma=sigfit_tp)
     
     if not keyword_set(noplot) then begin
        ;; Plot
        window, ii
        wset, ii

        plot,[0],[0],/nodata,/xsty,/ysty,xrange=ewlim[0:1],yrange=[0.,1.0],$
             color=clr.black,background=clr.white,charsize=2,thick=2,$
             xtitle=ewstring,ytitle='Fraction',$
             title=string(zbin[zz,0],zbin[zz,1],$
                          format='(f4.2," <= !8z!X < ",f4.2)')

        oploterror,afp_loc,afp_frac,0.5*afp_dloc,sigafp_frac[*,0],$
                   /lobar,/nohat,$
                   color=afp_clr,errcolor=afp_clr,psym=afp_psym,thick=2,symsize=2
        oploterror,afp_loc,afp_frac,0.5*afp_dloc,sigafp_frac[*,1],$
                   /hibar,/nohat,$
                   color=afp_clr,errcolor=afp_clr,psym=afp_psym,thick=2,symsize=2
        oplot,afp_loc,fit_afp,linestyle=0,thick=1,color=afp_clr
        oplot,afp_loc,fit_afp-sigfit_afp[*,0],linestyle=2,thick=1,color=afp_clr
        oplot,afp_loc,fit_afp+sigfit_afp[*,1],linestyle=2,thick=1,color=afp_clr

        oploterror,tp_loc,tp_frac,0.5*tp_dloc,sigtp_frac[*,0],$
                   /lobar,/nohat,$
                   color=tp_clr,errcolor=tp_clr,psym=tp_psym,thick=2,symsize=2
        oploterror,tp_loc,tp_frac,0.5*tp_dloc,sigtp_frac[*,1],$
                   /hibar,/nohat,$
                   color=tp_clr,errcolor=tp_clr,psym=tp_psym,thick=2,symsize=2
        oplot,tp_loc,fit_tp,linestyle=0,thick=1,color=tp_clr
        oplot,tp_loc,fit_tp-sigfit_tp[*,0],linestyle=2,thick=1,color=tp_clr
        oplot,tp_loc,fit_tp+sigfit_tp[*,1],linestyle=2,thick=1,color=tp_clr

        legend,['True Positives?','Accepted False Positives'],$
               psym=[tp_psym,afp_psym],$
               colors=[tp_clr,afp_clr],textcolors=[tp_clr,afp_clr],$
               /top
        
        ii++

     endif                      ; noplot=0
  endfor                        ; zz=nzbin

  if keyword_set(cmplt_fil) then begin
     ;; ;;;;;;;
     ;; Make completeness files
     ;; I'm pretty sure the directory structure is confused so
     ;; best if proot is supplied for
     ;; inputs/SNRge4/dr7qso_FkIV_noBAL[.list/_SNR.fit]
     if not keyword_set(nuser) then nuser = 1 ; sdss_getname() append
     froot = strmid(cmplt_fil,0,strpos(cmplt_fil,'_',/reverse_search)+1) 
     ;; _extra includes /rec_param,proot=
     sdss_completeness_cum,zbin,froot,dblt_name=dblt,silent=silent,$
                           user=nuser,clobber=clobber,/z,$
                           quick=0,civobs_corr=0,biasuser_fil=0,/nociv,$
                           out_fil=cmpltcum_fil,dvqso=dvqso,dvem=dvem,_extra=extra
  endif else cmpltcum_fil = ''


  ;; ;;;;;;;
  ;; Save all the necessary information
  ;; ext 1: summary structure
  ;; ext 2: all true recovered CIVSTRCT
  ;; ext 3: false positive CIVSTRCT
  ;; tp_indx maps to trurec_indx so applies to trurecstr fine
  ;; afp_indx maps to fpr_indx so applies to falsestr fine  
  ostrct = { $
           dblt_name:dblt.ion, $
           dvqso:dvqso, $
           cmplt_fil:cmpltcum_fil[0:nzbin-1], $ ; same zbins
           n_per_bin:n_per_bin, $
           nfitpar:nfitpar, $
           nzbin:nzbin, $
           zlim:zbin, $
           ewlim:ewlim, $
           nrand:nrand, $
           mask_trurecstr:mask_trurecstr, $
           tp_indx:tp_indx, $
           mask_falsestr:mask_falsestr, $
           afp_indx:afp_indx, $
           coeff_tp:coeff_tp, $ ; [z,[a0,low,high]]
           coeff_afp:coeff_afp $
           }
  if keyword_set(dvem) then $   ; Record
     tmp = create_struct(ostrct,'dvem',dvem) $
  else tmp = create_struct(ostrct,'dvem',0)
  ostrct = tmp
  mwrfits,ostrct,out_fil,/create,/silent
  mwrfits,trurecstr,out_fil,/silent
  mwrfits,falsestr,out_fil,/silent
  spawn,'gzip -f '+out_fil
  print,'sdss_completeness_userbias: created ',out_fil

  oseed = oseed[0]              ; for output

  if keyword_set(debug) then $
     stop,'sdss_completeness_userbias debug stop: last chance'
end                             ; sdss_completeness_userbias



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_completeness, list_fil, sdsssum, nabstot, dblt_name=dblt_name, $
                       config_fil=config_fil, vpstrct_fil=vpstrct_fil, $
                       processor=processor, eigbasis=eigbasis, pca_fil=pca_fil, $
                       pca_head=pca_head, tcpu=tcpu, istrt=istrt, $
                       clobber=clobber, debug=debug, _extra=extra
  if n_params() ne 3 then begin
     print,'Syntax - sdss_completeness, list_fil, sdsssum, nabstot, [dblt_name=, '
     print,'                       config_fil=, vpstrct_fil=, /tcpu, istrt=,'
     print,'                       processor=, eigbasis=, pca_fil, pca_head=,'
     print,'                       /clobber, /debug, _extra=]'
     return
  endif
  if not keyword_set(dblt_name) then $
     stop,'sdss_completeness stop: must set dblt_name='

  sdssdir = sdss_getsdssdir()


  ;; Params
  if keyword_set(seed) then oseed = seed ; will return to same variable
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  if size(config_fil,/type) eq 8 then config_strct = config_fil $
  else config_strct = sdss_genprof_config(config_fil,header=mchdr)
  mcstrct_tmplt = {sdssmcstrct} 
  tmp_fil = 'sdss_completeness_mchdr.fit'
  mwrfits,mcstrct_tmplt,tmp_fil,/create,/silent
  mwrfits,config_strct,tmp_fil,/silent
  dum = xmrdfits(tmp_fil,1,mchdr,/silent)
  spawn,'rm '+tmp_fil
  
  ;; Read in
  readcol,list_fil,spec_fil,skip=1,format='a',/silent
  nfil = (size(spec_fil,/dim))[0] > 1 ; foil singularity

  if keyword_set(processor) then begin
     ;; Takes precedence over istrt=
     sub = sdss_calcparalleljob(sdsssum, processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     if keyword_set(debug) then $
        print,'sdss_completeness debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif 
  if not keyword_set(istrt) then istrt = 0L ; default


  if not keyword_set(eigbasis) then $ ; better not to read in every time
     eigbasis = xmrdfits(getenv("SDSSPATH")+"/eigenspectra/eigSpec_qso_all.fit",$
                         0, /silent)
  if not keyword_set(pca_fil) then $
     pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
  if size(pca_fil,/type) eq 7 then $
     pca = xmrdfits(pca_fil, 0, pca_head, /silent) $
  else begin  
     pca = pca_fil 
     if not keyword_set(pca_head) then $
        stop,'sdss_completeness stop: must set pca_head' 
  endelse 

  if keyword_set(vpstrct_fil) then begin
     if size(vpstrct_fil,/type) eq 7 then $
        vpstrct0 = xmrdfits(vpstrct_fil,1,/silent) $
     else vpstrct0 = vpstrct_fil
  endif else vpstrct0 = 0       ; let sdss_completeness_los handle it

  ;; Set output file names
  ;; _extra includes /quick, /user
  abslin_fil = sdss_getname(spec_fil,/spec,/abslin,dir=absdir)
  abslin_fil = absdir + abslin_fil
  mcstrct_fil = sdss_getname(spec_fil,/spec,mc=dblt.ion,dir=mcstrctdir,$
                             _extra=extra)
  mcstrct_fil = mcstrctdir + mcstrct_fil

  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  for ff=istrt,nfil-1 do begin
     ;; Check file status
     test = file_search(sdssdir+mcstrct_fil[ff]+'*',count=ntest)
     if ntest eq 0 or keyword_set(clobber) then begin
        if ntest eq 0 then begin
           test2 = file_search(sdssdir+mcstrctdir[ff],count=ntest2)
           if ntest2 eq 0 then $
              spawn,'mkdir -p '+sdssdir+mcstrctdir[ff]
        endif 
     endif else begin
        print,'sdss_completeness: will not clobber ',mcstrct_fil[ff]
        continue
     endelse 

     ;; Main function, handles everything, including adjusting and
     ;; trimming any given vpstrct0 to the redshift and nabstot
     ;; necessary 
     ;; _extra includes /quick

     if keyword_set(debug) then print,''
     sdss_completeness_los,dblt,spec_fil[ff],abslin_fil[ff],config_strct,nabstot,$
                           mcstrct_fil[ff],vpstrct_fil=vpstrct0,mchdr=mchdr,$
                           mcstrct_tmplt=mcstrct_tmplt,eigbasis=eigbasis,$
                           pca_head=pca_head,pca_fil=pca_fil,$
                           seed=oseed,oseed=oseed,$ ; iterate
                           debug=debug,clobber=clobber,_extra=extra

     if keyword_set(tcpu) and ((ff+1) mod 500) eq 0 then begin
        tt = systime(/seconds)
        dt = (tt-tlast)         ; seconds
        print,'sdss_completeness: Elapsed time (m) = ',dt/60.
        print,'sdss_completeness: Average time per ff (s) = ',dt/500
        tlast = tt
     endif 

  endfor                        ; loop ff=nfil

  print,'sdss_completeness: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_completeness: Elapsed time for '+strtrim(nfil-istrt,2)+$
        ' QSOs (m) = ',dt/60.
  print,'sdss_completeness: Average time per QSO (s) = ',dt/(nfil-istrt)

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end                             ; sdss_completeness
