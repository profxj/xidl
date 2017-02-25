;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_synthspec.pro               
; Author: Kathy Cooksey                      Date: 26 Aug 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Driver program for completeness tests
; Input: 
;   qso -- name of QSO as used to make directories
;   zqso -- redshift of QSO
; Optional:
;   outdir -- directory name to store output (default: synth)
;             $MLSS_DIR/analysis/qso/outdir
;   nabs -- number of profiles to generate (default: 5000)
;   zlim -- redshift range (default: -1000 km/s to zqso+5000 km/s)
;   config_fil -- driver for profile generate (e.g. nlim, blim...)
;                 (default: $MLSS_DIR/pro/civ_genprof.config'
;   /savfits -- whether or not to save the massive amounts of files
; Output: 
;   (in $MLSS_DIR/analysis/qso/outdir)
;   qso_genprof.fits -- full range of profiles created for QSO
;   clean_[spectrum name] -- spectrum with absorption lines removed
;         -- also produced are clean*instr.lst, autoclean*, recclean*
;   true[#]_[spectrum root]civ.fits -- input (true) synthetic lines
;   synth[#]_[spectrum name] -- synthetic spectrum
;   synth[#]_[spectrum_root]instr.lst -- temporary file; can be deleted
;   auto[#]_[spectrum root]abslin.fits -- automatically detected features
;   auto[#]_[spectrum_root]_[instrument].fits/.tab
;   auto[#]_[spectrum_root]_dopb.tab
;   rec[#]_[spectrum root]civ.fits -- output (recovered) synthetic lines
; Example:
; History:
;   26 Aug 2008  created by KLC
;   28 Aug 2008  added search of clean spectrum
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_genprof
@civ_complete
pro civ_synthspec_dndz,qso,zqso,siiv=siiv,outdir=outdir,nabs=nabs,zlim=zlim,$
                       config_fil=config_fil,prof_fil=prof_fil,savfits=savfits,$
                       _extra=extra
savfits = 1
if (N_params() lt 2) then begin
    print,'Syntax - '+$
      'civ_synthspec,qso,zqso,outdir=outdir,nabs=nabs,[zlim=zlim,'+$
      '              config_fil=config_fil,_extra=extra]'
    return
endif

if keyword_set(siiv) then dblt_name = 'SiIV' $
else dblt_name = 'CIV'

pwd,curdir                      ; to return

rootdir = getenv('MLSS_DIR')+'/analysis/'+qso+'/'
cd,rootdir

;; Default values
if not keyword_set(nabs) then nabs = 15000
if not keyword_set(zlim) then zlim = [-1000./3e5,5000./3e5*(1+zqso)+zqso]
if not keyword_set(outdir) then outdir = 'synth'+dblt_name
if not keyword_set(config_fil) then config_fil = $
  getenv('MLSS_DIR')+'/pro/civ_genprof.config'

;; Useful items
civ = dblt_retrieve(dblt_name)
config_strct = civ_genprof_config(config_fil)
dndzobs = config_strct.nabsobs/config_strct.dzunblck

;; File names
instr_list = 'lists/'+qso+'instr.lst'
spawn,'ls -d '+outdir,test,err
if strtrim(test,2) eq '' then spawn,'mkdir '+outdir

if not keyword_set(prof_fil) then begin
    prof_fil = outdir + '/'+qso+'_genprof.fits'
    ;; Extras include: /calcew, seed, debug
    civ_genprof,prof_fil,nabs,zlim,config_fil=config_strct,$
                siiv=siiv,_extra=extra 
endif 


;; Read instrument file
readcol,instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f' 
instr_fil = getenv('MLSS_DIR')+'/'+instr_fil
nlist = n_elements(instr_fil)
if nlist LT 7 then stop,'civ_synthspec: too few spectra'

instr_fil_nodir = instr_fil
for ii=0,nlist-1 do begin
    prs = strsplit(instr_fil[ii],'/',count=nprs,/extract)
    instr_fil_nodir[ii] = prs[nprs-1]
endfor                          ; loop instruments

for ispec=0,nlist-1 do begin
    ;; Verify spectrum exists
    test = file_search(instr_fil[ispec],count=ntest)
    if ntest eq 0 then continue

    ;; Read in spectrum
    if ispec le 6 then begin
        fx = x_readspec(instr_fil[ispec], SIG=sig, wav=wave, NPIX=npix, $
                        inflg=3,head=head) 
        spec_root = strmid(instr_fil_nodir[ispec],0,$
                           strpos(instr_fil_nodir[ispec],'.fit'))
        if not keyword_set(fwhm) then fwhm = 5
    endif else begin            ; HST
        spos = strpos(instr_fil[ispec], 'f.fits')
        sig_fil = strmid(instr_fil[ispec], 0, spos)+'e.fits'
        fx = x_readspec(instr_fil[ispec], SIG=sig, wav=wave, $
                        NPIX=npix,fil_sig=sig_fil, inflg=0, head=head)
        spec_root = strmid(instr_fil_nodir[ispec],0,$
                           strpos(instr_fil_nodir[ispec],'_f.fit'))
        if not keyword_set(fwhm) then fwhm = 2
    endelse
    ;; Shift wave
    wave = wave + inst_dw[ispec]*inst_w0[ispec]/wave

    ;; Check if CIV detectable
    if min(wave) lt civ.wvI*(1+zlim[0]) then continue

    dzspec = [min(wave)/civ.wvI-1 > zlim[0], max(wave)/civ.wvII-1 < zlim[1]]
    ;; Remove absorption lines 
    cleanspec_fil = outdir+'/clean_'+instr_fil_nodir[ispec]
    inspec = {wave:wave,flux:fx,error:sig}
    civ_cleanspec,inspec,cleanspec,varflr=varflr,_extra=extra ; seed
    if ispec le 6 then $
      mwrfits,cleanspec,cleanspec_fil,head,/create,/silent $
    else begin
        mwrfits,cleanspec.flux,cleanspec_fil,head,/create,/silent 
        spos = strpos(cleanspec_fil, 'f.fits')
        sig_fil = strmid(cleanspec_fil, 0, spos)+'e.fits'
        mwrfits,cleanspec.error,sig_fil,head,/create,/silent
    endelse

    ;; Create clean spectrum instrument file
    synthinstr_fil = outdir+'/clean_'+spec_root+'instr.lst'
    openw,1,synthinstr_fil
    len = strlen('analysis/'+qso+'/'+cleanspec_fil)
    fmt = '(a'+strtrim(len,2)+',2x,f6.4,2x,f6.4)'
    for ii=0,ispec-1 do printf,1,'dne.fits',0.0,1.0,format=fmt
    printf,1,'analysis/'+qso+'/'+cleanspec_fil,0.0,1.0,format=fmt
    if ispec ne nlist-1 then $
      for ii=ispec+1,nlist-1 do printf,1,'dne.fits',0.0,1.0,format=fmt
    close,1
        
    ;; Recover un-cleaned features
    auto_fil = outdir+'/autoclean_'+spec_root+'abslin.fits'
    civrec_root = outdir+'/recclean_'+spec_root+strlowcase(dblt_name)
    civrec_fil = civrec_root+'.fits'
    search_all,synthinstr_fil,auto_fil,3.,/calcb,name=qso,$
      root=outdir+'/autoclean_'+spec_root
    civ_find,auto_fil,siiv=siiv,root=civrec_root,zlim=dzspec,$
             instrfil=synthinstr_fil
    test = file_search(civrec_fil,count=ntest) 
    if ntest eq 0 then continue 

    ;; Other post-processing
    civstrct = xmrdfits(civrec_fil,1,/silent) 
    civstrct.qso = qso
    civstrct.zqso = zqso
    if keyword_set(savfits) then begin
        mwrfits,civstrct,civrec_fil,/create,/silent 
        civ_calcewn,civrec_fil,'analysis/'+qso+'/'+synthinstr_fil
        civ_aodm,civrec_fil,'analysis/'+qso+'/'+synthinstr_fil,/log,nbin=3
    endif else begin
        civ_calcewn,civstrct,'analysis/'+qso+'/'+synthinstr_fil
        civ_aodm,civstrct,'analysis/'+qso+'/'+synthinstr_fil,/log,nbin=3
    endelse 


    ;; Loop to well-sample the parameter space (including the
    ;; wavelength coverage of the spectra)
    nsub = round(dndzobs * (dzspec[1] - dzspec[0]))
    nreal = round(nabs/dndzobs) ; "nabs" = dN_MC/dz so nabs = dN_MC/dz*dzspec

    ;; Save names to pass to civ_complete
    trufits = strarr(nreal)
    recfits = strarr(nreal)

    for isynth=0,nreal-1 do begin
        num = strtrim(isynth+1,2)
        ;; _extra may be indx,ncomp,noblend,nlim,blim
        subprof = civ_genprof_sub(prof_fil,nsub,dzspec,_extra=extra) 
        trufits[isynth] = outdir+'/true'+num+'_'+spec_root+$
                          strlowcase(dblt_name)+'.fits'
        ;; Must write this file
        mwrfits,subprof,trufits[isynth],/create,/silent
        
        ;; Create vplin and synthesize noise
        synthspec_fil = outdir+'/synth'+num+'_'+instr_fil_nodir[ispec]
        civ_setvplin,subprof,vplin,_extra=extra ; /dvabs
        vpfx = civ_voigt(cleanspec.wave,vplin,fwhm=fwhm)
        civ_synthnoise,cleanspec,varflr,vpfx,synthspec,$
          _extra=extra          ;conti,seed,debug
        if ispec le 6 then $
          mwrfits,synth,synthspec_fil,head,/create,/silent $
        else begin
            mwrfits,synthspec.flux,synthspec_fil,head,/create,/silent
            spos = strpos(synthspec_fil, 'f.fits')
            sig_fil = strmid(synthspec_fil, 0, spos)+'e.fits'
            mwrfits,synthspec.error,sig_fil,head,/create,/silent
        endelse
        if keyword_set(savfits) then $
          mwrfits,vpfx,outdir+'/vp'+num+'_'+instr_fil_nodir[ispec],head,$
          /create,/silent       ; could append
        civ_synthew,vpfx,subprof,/totew,wave=cleanspec.wave

        ;; Create one-spectrum instrument file
        synthinstr_fil = outdir+'/synth'+num+'_'+spec_root+'instr.lst'
        openw,1,synthinstr_fil
        len = strlen('analysis/'+qso+'/'+synthspec_fil)
        fmt = '(a'+strtrim(len,2)+',2x,f6.4,2x,f6.4)'
        for ii=0,ispec-1 do printf,1,'dne.fits',0.0,1.0,format=fmt
        printf,1,'analysis/'+qso+'/'+synthspec_fil,0.0,1.0,format=fmt
        if ispec ne nlist-1 then $
          for ii=ispec+1,nlist-1 do printf,1,'dne.fits',0.0,1.0,format=fmt
        close,1
        
        ;; Recover synthetic lines
        auto_fil = outdir+'/auto'+num+'_'+spec_root+'abslin.fits'
        civrec_root = outdir+'/rec'+num+'_'+spec_root+strlowcase(dblt_name)
        recfits[isynth] = civrec_root+'.fits'
        ;; search_all creates root_INSTR.fits/.tab and root_dopb.info
        search_all,synthinstr_fil,auto_fil,3.,/calcb,name=qso,$
          root=outdir+'/auto'+num+'_'+spec_root
        civ_find,auto_fil,root=civrec_root,zlim=dzspec,$
                 instrfil=synthinstr_fil,siiv=siiv
        test = file_search(recfits[isynth],count=ntest) 
        if ntest eq 0 then continue 

        ;; Other post-processing
        civstrct = xmrdfits(recfits[isynth],1,/silent) 
        civstrct.qso = qso
        civstrct.zqso = zqso
        civ_calcewn,civstrct,'analysis/'+qso+'/'+synthinstr_fil
        civ_aodm,civstrct,'analysis/'+qso+'/'+synthinstr_fil,/log,nbin=3
        civ_flag,civstrct,dvlim=10.
        ;; Must write out this file
        mwrfits,civstrct,recfits[isynth],/create,/silent 

        ;; Clean up un-necessary byproducts
        spawn,'\rm '+outdir+'/*allID.*' ; search_all
        spawn,'\rm '+outdir+'/*dopb.info' ; search_all
        spawn,'\rm '+outdir+'/*'+strlowcase(dblt_name)+'.fits.bkp' ; civ_find
        spawn,'\rm '+outdir+'/*orig.fits'    ; civ_find
        if ispec le 6 then spawn,'\rm '+outdir+'/*_[LS]i[FC][12][ab].*' $
        else spawn,'\rm '+outdir+'/*_HST.*'
    endfor                      ; loop realizations

    ;; Completeness check
    ostrct_fil = 'cmplt_'+spec_root+'.fits'
    civ_complete,trufits,recfits,ostrct_fil=ostrct_fil, siiv=siiv,$
                 _extra=extra   ; bincolm,binz

    if not keyword_set(savfits) then begin
       ;; Delete everything
;       spawn,'\rm '+outdir+'/?*clean*.fits' ; for study of features in clean
;       spawn,'\rm '+outdir+'/vp[0-9]*.fits' ; w/clean*.fits, recreate
       spawn,'\rm '+outdir+'/synth[0-9]*.fits'
;       spawn,'\rm '+outdir+'/rec[0-9]*.fits' ; during testing phase, need
;       spawn,'\rm '+outdir+'/true[0-9]*.fits' ; during testing phase, need
    endif
endfor                          ; loop instrument file

;; Return
cd,curdir
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_synthspec,qso,zqso,siiv=siiv,outdir=outdir,nabs=nabs,$
                  zlim=zlim,ncomp=ncomp,$
                  config_fil=config_fil,savfits=savfits,$
                  resume=resume,_extra=extra
  if (N_params() lt 2) then begin
     print,'Syntax - '+$
           'civ_synthspec,qso,zqso,outdir=outdir,nabs=nabs,[zlim=zlim,'+$
           '              config_fil=config_fil,_extra=extra]'
     return
  endif
  pwd,curdir                    ; to return

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  rootdir = getenv('MLSS_DIR')+'/analysis/'+qso+'/'
  cd,rootdir

  ;; Default values
  if not keyword_set(nabs) then nabs = 100
  if not keyword_set(zlim) then zlim = [-1000./3e5,5000./3e5*(1+zqso)+zqso]
  if not keyword_set(outdir) then outdir = 'synth'+dblt_name
  if not keyword_set(config_fil) then config_fil = $
     getenv('MLSS_DIR')+'/pro/civ_genprof.config'

  ;; Useful items
  civ = dblt_retrieve(dblt_name)
  dvciv = round(2.998e5 * (civ.wvII-civ.wvI)/civ.wvI)
  config_strct = civ_genprof_config(config_fil)
  nbinncolm = ceil((config_strct.nlim[1]-config_strct.nlim[0])/$
                   config_strct.binncolm) 
  if (nbinncolm-1)*config_strct.binncolm+config_strct.nlim[0] ge $
     config_strct.nlim[1] then nbinncolm = nbinncolm-1

  ;; File names
  instr_list = 'lists/'+qso+'instr.lst'
  spawn,'ls -d '+outdir,test,err
  if strtrim(test,2) eq '' then spawn,'mkdir '+outdir


  ;; Read instrument file
  readcol,instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f' 
  instr_fil = getenv('MLSS_DIR')+'/'+instr_fil
  nlist = n_elements(instr_fil)
  if nlist LT 7 then stop,'civ_synthspec: too few spectra'

  instr_fil_nodir = instr_fil
  for ii=0,nlist-1 do begin
     prs = strsplit(instr_fil[ii],'/',count=nprs,/extract)
     instr_fil_nodir[ii] = prs[nprs-1]
  endfor                        ; loop instruments

  for ispec=0,nlist-1 do begin
     ;; Verify spectrum exists
     test = file_search(instr_fil[ispec],count=ntest)
     if ntest eq 0 then continue

     ;; Read in spectrum
     if ispec le 6 then begin
        fx = x_readspec(instr_fil[ispec], SIG=sig, wav=wave, NPIX=npix, $
                        inflg=3,head=head) 
        spec_root = strmid(instr_fil_nodir[ispec],0,$
                           strpos(instr_fil_nodir[ispec],'.fit'))
        if not keyword_set(fwhm) then fwhm = 5
     endif else begin           ; HST
        spos = strpos(instr_fil[ispec], 'f.fits')
        sig_fil = strmid(instr_fil[ispec], 0, spos)+'e.fits'
        fx = x_readspec(instr_fil[ispec], SIG=sig, wav=wave, $
                        NPIX=npix,fil_sig=sig_fil, inflg=0, head=head)
        spec_root = strmid(instr_fil_nodir[ispec],0,$
                           strpos(instr_fil_nodir[ispec],'_f.fit'))
        if not keyword_set(fwhm) then fwhm = 2
     endelse
     ;; Shift wave
     wave = wave + inst_dw[ispec]*inst_w0[ispec]/wave

     ;; Check if CIV detectable
     if (max(wave) lt civ.wvI*(1+zlim[0])) or $
        (min(wave) gt civ.wvII*(1+zlim[1])) then continue

     ;; Measure dopb to use with searchspec (based on resolution)
     dwv = wave - shift(wave,1)
     dwv[0] = dwv[1]
     if ispec le 6 then r = median(0.2*wave/dwv) $ ; 5 pix per R
     else r = median(0.5*wave/dwv)                 ; 2 pix per R
     dopb = round(2.998e5*2.354/r)

     dzspec = [min(wave)/civ.wvI-1 > zlim[0], max(wave)/civ.wvII-1 < zlim[1]]
     if dzspec[0] gt dzspec[1] then continue
     ;; Remove absorption lines 
     cleanspec_fil = outdir+'/clean_'+instr_fil_nodir[ispec]
     inspec = {wave:wave,flux:fx,error:sig}


     ;; New run 
     civ_cleanspec,inspec,cleanspec,varflr=varflr,instr_flg=2^ispec,$
                   _extra=extra ; search_fil,seed
     if keyword_set(savfits) then begin
        if ispec le 6 then $
           mwrfits,cleanspec,cleanspec_fil,head,/create,/silent $
        else begin
           mwrfits,cleanspec.flux,cleanspec_fil,head,/create,/silent 
           spos = strpos(cleanspec_fil, 'f.fits')
           sig_fil = strmid(cleanspec_fil, 0, spos)+'e.fits'
           mwrfits,cleanspec.error,sig_fil,head,/create,/silent
        endelse
     endif 
     
     ;; Recover un-cleaned features
     auto_fil = outdir+'/autoclean_'+spec_root+'abslin.fits'
     civrec_root = outdir+'/recclean_'+spec_root+strlowcase(dblt_name)
     civrec_fil = civrec_root+'.fits'

     ;; May want to save auto_fil to do false-positive tests
     if keyword_set(savfits) then begin
        searchspec,cleanspec.wave,cleanspec.flux,cleanspec.error,3.,dopb,$
                   instr=2^ispec,dw0=inst_dw[ispec], w0=inst_w0[ispec], /shift, $
                   ostrct=auto,header=auto_hd,/silent
        
        mwrfits,auto,auto_fil,auto_hd,/create,/silent
        
        if size(auto,/type) eq 8 then begin
           civ_find,auto,ostrct=civstrct,zlim=dzspec,siiv=siiv,$
                    inspec=cleanspec,instr_flg=2^ispec,_extra=extra ;/silent
           if size(civstrct,/type) eq 8 then begin
              ;; Other post-processing
              civstrct.qso = qso
              civstrct.zqso = zqso
              civ_calcewn,civstrct,inspec=cleanspec,instr_flg=2^ispec
              civ_aodm,civstrct,inspec=cleanspec,/log,nbin=3
              civ_flag,civstrct,dvlim=10.
              if keyword_set(savfits) then $
                 mwrfits,civstrct,civrec_fil,/create,/silent 
           endif                ; CIV found
        endif 
     endif                      ; savfits (clean search)


     ;; Loop to well-sample the parameter space (including the
     ;; wavelength coverage of the spectra)
     nbinz = ceil((dzspec[1]-dzspec[0])/config_strct.binz) 
     if (nbinz-1)*config_strct.binz+dzspec[0] ge dzspec[1] then $
        nbinz = nbinz - 1
     mxsep = round(config_strct.dvabs[1]-config_strct.dvabs[0] + $
                   3*config_strct.blim[1] + dvciv) > config_strct.binz*2.998e5
     nzsep = ceil(mxsep/(config_strct.binz*2.998e5)) + 1 ; force space
     nreal = nzsep*nbinncolm*nabs ; max number
     nsub = ceil(nbinz/double(nzsep))

     ;; Save names to pass to civ_complete
     trufits = strarr(nreal)
     recfits = strarr(nreal)
     nfits = 0
     
     zlimsub = [0,config_strct.binz] + dzspec[0]

     ;; Loop column density bins
     for incolm=0,nbinncolm-1 do begin
        nlimsub = config_strct.nlim[0]+[incolm,incolm+1]*$
                  config_strct.binncolm
        if nlimsub[1] gt config_strct.nlim[1] then $
           nlimsub[1] = config_strct.nlim[1]

        ;; _extra may be seed=, /calcewn, /debug
        civ_genprof,fullprof,nabs,zlimsub,nlimsub,$
                    siiv=siiv,config_fil=config_strct

        ;; Need to call sub to get components
        prof = civ_genprof_sub(fullprof,nabs,zlimsub,ncomp=ncomp) 

        ;; Loop absorbers
        for iabs=0,nabs-1 do begin

           ;; Loop redshift
           for iz=0,nzsep-1 do begin
              
              ;; Name
              num = 'z'+strtrim(iz+1,2)+'n'+strtrim(incolm+1,2)+$
                    'a'+strtrim(iabs+1,2)
              gd = where(prof[iabs].wrest gt 0.,ngd)
              unq = uniq(prof[iabs].wrest[gd],sort(prof[iabs].wrest[gd]))
              nionabs = n_elements(unq)
              ncompabs = ngd/nionabs

              ;; Space out feature across spectrum
              subprof = replicate(prof[iabs],nsub)
              for iion=0,nionabs-1 do begin
                 rng = lindgen(ncompabs)*nionabs + iion
                 for icomp=0,ncompabs-1 do begin
                    subprof.zabs[rng[icomp]] = $
                       subprof.zabs[rng[icomp]] + $
                       dindgen(nsub)*config_strct.binz*nzsep + $
                       iz*config_strct.binz
                    wcent = subprof.wrest[rng[icomp]]*(1+subprof.zabs[rng[icomp]])
                    wcent0 = prof[iabs].wrest[rng[icomp]]*$
                             (1+prof[iabs].zabs[rng[icomp]])
                    subprof.wv_lim[rng[icomp],0] = $ ; lower
                       wcent - (wcent0 - prof[iabs].wv_lim[rng[icomp],0])
                    subprof.wv_lim[rng[icomp],1] = $ ; upper
                       wcent + (prof[iabs].wv_lim[rng[icomp],1] - wcent0)
                 endfor         ; loop component ions 
              endfor            ; loop nionabs

              ;; Remove features if outside of dzspec (assuming 1550)
              bd = where(subprof.zabs[1] gt dzspec[1],nbd,complement=gd)
              if gd[0] eq -1 then begin
                 ;; Continuing should keep the density of lines
                 if nsub ne 1 then $
                    stop,'civ_synthspec: no good profiles to use'
                 continue 
              endif 
              if nbd ne 0 then subprof = subprof[gd]
              
              ;; Create vplin and synthesize noise
              synthspec_fil = outdir+'/synth'+num+'_'+instr_fil_nodir[ispec]
              civ_setvplin,subprof,vplin,dvabs=ncomp
              vpfx = civ_voigt(cleanspec.wave,vplin,fwhm=fwhm)
              civ_synthnoise,cleanspec,varflr,vpfx,synthspec,$
                             _extra=extra ;conti,seed,debug
              test = file_search(synthspec_fil,count=ntest)
              if keyword_set(resume) and ntest ne 0 then begin
                 if ispec le 6 then $
                    synthspec = xmrdfits(synthspec_fil,1,/silent) $
                 else begin
                    spos = strpos(synthspec_fil, 'f.fits')
                    sig_fil = strmid(synthspec_fil, 0, spos)+'e.fits'
                    test = x_readspec(synthspec_fil,wav=wv,sig=er,$
                                      fil_sig=sig_fil,npix=ntest)
                    synthspec = {wave:wv,flux:test,error:er}
                 endelse 
                 vpfx = xmrdfits(outdir+'/vp'+num+'_'+instr_fil_nodir[ispec],$
                                 /silent)
              endif else begin
                 if keyword_set(savfits) then begin
                    if ispec le 6 then $
                       mwrfits,synthspec,synthspec_fil,head,/create,/silent $
                    else begin
                       mwrfits,synthspec.flux,synthspec_fil,head,/create,/silent
                       spos = strpos(synthspec_fil, 'f.fits')
                       sig_fil = strmid(synthspec_fil, 0, spos)+'e.fits'
                       mwrfits,synthspec.error,sig_fil,head,/create,/silent
                    endelse
                    
                    mwrfits,vpfx,outdir+'/vp'+num+'_'+instr_fil_nodir[ispec],$
                            head,/create,/silent ; could append
                 endif 
                 civ_synthew,vpfx,subprof,/totew,wave=cleanspec.wave
              endelse 

              trufits[nfits] = outdir+'/true'+num+'_'+spec_root+$
                               strlowcase(dblt_name)+'.fits'
              ;; Must write out this file (will have EW)
              test = file_search(trufits[nfits],count=ntest)
              if keyword_set(resume) and ntest ne 0 then $
                 subprof = xmrdfits(trufits[nfits],1,/silent) $ ; use existing 
              else mwrfits,subprof,trufits[nfits],/create,/silent ; write new
              

              ;; Recover synthetic lines
              auto_fil = outdir+'/auto'+num+'_'+spec_root+'abslin.fits'
              civrec_root = outdir+'/rec'+num+'_'+spec_root+$
                            strlowcase(dblt_name)
              recfits[nfits] = civrec_root+'.fits'
              
              test = file_search(recfits[nfits],count=ntest)
              if keyword_set(resume) and ntest ne 0 then $
                 civstrct = xmrdfits(recfits[nfits],1,/silent) $
              else begin
                 ;; Have to search and save
                 searchspec,synthspec.wave,synthspec.flux,synthspec.error,3.,$
                            dopb,$
                            instr=2^ispec,dw0=inst_dw[ispec],w0=inst_w0[ispec],$
                            /shift,ostrct=auto,header=auto_hd,/silent
                 if keyword_set(savfits) then $
                    mwrfits,auto,auto_fil,auto_hd,/create,/silent
                 
                 if size(auto,/type) eq 8 then begin
                    civ_find,auto,ostrct=civstrct,zlim=dzspec,siiv=siiv,$
                             inspec=synthspec,instr_flg=2^ispec,$
                             _extra=extra ;/silent
                    if size(civstrct,/type) eq 8 then begin
                       ;; Other post-processing
                       civstrct.qso = qso
                       civstrct.zqso = zqso
                       civ_calcewn,civstrct,inspec=synthspec,instr_flg=2^ispec 
                       civ_aodm,civstrct,inspec=synthspec,/log,nbin=3
                       civ_flag,civstrct,dvlim=10.
                       ;; Must write out this file
                       mwrfits,civstrct,recfits[nfits],/create,/silent 
                    endif       ; CIV found
                 endif          ; features found
              endelse


              ;; Next file
              nfits = nfits + 1
           endfor               ; loop nzsep
        endfor                  ; loop nabs
     endfor                     ; loop ncolm

     ;; Completeness check (what is saved in trufits, recfits)
     ostrct_root = 'cmplt'+strlowcase(dblt_name)+'_'+spec_root
     civ_complete,trufits[0:nfits-1],recfits[0:nfits-1],$
                  ostrct_fil=ostrct_root+'.fits',$
                  binncolm=config_strct.binncolm, binz=config_strct.binz,$
                  siiv=siiv,$
                  spec_fil=instr_fil[ispec],instr_flg=2^ispec,$
                  _extra=extra  ;/silent,search_fil=
     ;; Two plots
     psfil = ostrct_root +'.ps' 
     prs = strsplit(psfil,'cmplt',/extract,count=nprs,/regex) 
     civ_complete_plot,ostrct_root+'.fits','cmpltn' + prs[nprs-1],0,$
                       yrng=[12.6,14.4],/pltnhist,siiv=siiv 
     civ_complete_plot,ostrct_root+'.fits','cmpltew' + prs[nprs-1],1,$
                       yrng=[30.,200.],/pltnhist,siiv=siiv 

  endfor                        ; loop instrument file

  ;; Return
  cd,curdir
end
