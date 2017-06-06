;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_genprof.pro             
; Author: Kathy Cooksey                      Date: 23 Aug 2011
; Project: 
; Description: Generate civcandstrct of randomly generated
;              values of z, [N, b, dv, ncomp]
; Input: 
;   nabs - number of absorption lines to make for each wrest
;   zlim - bounds of randomly generated redshifts
;   nlim - bounds of randomly generated column densities (logN)
;   blim - bounds of randomly generated Doppler parameters (km/s)
;   nlimobs - bounds of observed column densities that define f(N)
;   nabsobs - number of absorbers found in nlimobs (normalize f(N))
;   aobs - power-law exponent f(N) = b*N^(-aobs)
; Optional:
;   seed - IDL randomu random number generation seed
;   dvabs - bounds of randomly generated velocity offsets for 
;           components (km/s)
;   ndv - number of components, used for each input wrest, and
;         cannot exceed number civcandstrct holds
;   debug - stop at opportune times for debugging purposes
; Output: 
;   strct_fil - name of output civcandstrct, with nabs elements,
;               with n_elements(wrest) times ndv each
; Example:
;   sdss_genprof,'randciv.fits',50,[0,1],[12,15],[10,40],
;      [12.7,14.5],21,1.83,,dvabs=[-200,200],ndv=5
;   [Danfort & Shull 2007ph: 12.7 < ncolm < 14.5, alpha = 1.83+/-0.32,
;                            delta(z) = 2.22, nabsobs = 21 for CIV]
; History:
;   23 Mar 2008 - created by KLC
;   23 Aug 2011 - adapted from civ_genprof, KLC
;   28 Jul 2015 - change to use x_voigt() only, KLC
;   15 Mar 2016 - no longer need config_strct in sdss_genprof_voigt(),
;                 KLC
;   06 Sep 2016 - fix sdss_genprof_setnewz() (scrambled z_sys
;                 and wvlim_sys),
;                 sdss_genprof_reorder() (handle index_sys input and
;                 return consistent items), and
;                 sdss_genprof_mc() (clarify), DRM via KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@spec_ewlin

function sdss_genprof_config,config_fil,header=header
  ;; Read in and return structure from configuration file
  cfg = {nlim:fltarr(2), $
         ewlim:fltarr(2), $
         blim:fltarr(2), $
         dvlim:0., $
         ndv:0L, $
         flg_dvabs:0, $
         flg_capn:0, $
         dvabs:fltarr(2), $
         dndz:0., $
         frac_rm:0., $
         pixscale:0., $
         dvelo:0., $
         nabs:0L, $
         zlim:fltarr(2) $
        }
  if not keyword_set(config_fil) then return,cfg

  tags = tag_names(cfg)
  ntags = (size(tags,/dim))[0]

  readcol,config_fil,val,descript,format='a,a',delimiter='|',/silent
  descript = strtrim(descript,2)
  nval = n_elements(val)

  for ii=0,ntags-1 do begin
     mtch = where(tags[ii] eq strupcase(descript),nmtch)

     if nmtch eq 0 then begin
        print,'sdss_genprof_config(): value not in config file ',tags[ii]
        continue
     endif 
     
     if nmtch gt 1 then $
        stop,'sdss_genprof_config(): multiple values in config file ',tags[ii]
     
     mtch = mtch[0]

     if n_elements(cfg.(ii)) eq 2 then begin
        prs = strsplit(val[mtch],/extract)
        cfg.(ii) = [float(prs[0]),float(prs[1])]
        if keyword_set(header) then begin
           sxaddpar,header,tags[ii]+'MN',cfg.(ii)[0]
           sxaddpar,header,tags[ii]+'MX',cfg.(ii)[1]
        endif 
     endif else begin
        cfg.(ii) = float(val[mtch]) ; ndv should be long
        if keyword_set(header) then $
           sxaddpar,header,tags[ii],cfg.(ii)
     endelse 

  endfor                        ; loop ii=ntgas

  return,cfg
end                             ; sdss_genprof_config()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_setline, wrest,id_lin=id_lin
  if n_params() ne 1 then begin
     print,'sdss_genprof_setline(wrest,[id_lin=])'
     return,-1
  endif 
  
  strct = create_struct(x_setline(wrest),$ ; includes fval
                        'ID_SYS',0L,$      ; 0 to 2*config_strct.nabs
                        'NCOMP',0, $       ; in addition to primary
                        'ID_COMP',0, $     ; 0 = primary; sequential
                        'ID_LIN',0,$       ; 1 = wvI; 2 = wvII
                        'Z_SYS',0d,$       ; column-weighted
                        'N_SYS',0., $      ; column density sum
                        'WVLIM_SYS',[0.,0.],$ ; min/max empirically
                        'EW_SYS',0.,$         ; w/o noise
                        'EW_OBS',0.)          ; w/o noise but force to end at midpoint of doublet
  if keyword_set(id_lin) then strct.id_lin = id_lin

  return, strct 

end                             ; sdss_genprof_setline()



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_setnewz, vpstrct_fil, zlim, seed=seed, oseed=oseed, $
                               reorder=reorder,nabstot=nabstot,$
                               newindex=newindex, doppler=doppler
  ;; Change the redshifts for given redshifts limits and optionaly,
  ;; scramble the order of the systems to add more spice (and more
  ;; efficient than calling sdss_genprof_reorder() by itself)
  if n_params() ne 2 then begin
     print,'Syntax - sdss_genprof_setnewz(vpstrct_fil, zlim, [seed=, oseed=,'
     print,'                              /reorder,nabstot=, newindex=, /doppler])'
     return,-1
  endif 

  if size(vpstrct_fil,/type) eq 8 then $
     vpstrct = vpstrct_fil $
  else vpstrct = xmrdfits(vpstrct_fil,1,/silent) 

  if keyword_set(seed) then oseed = seed
  
  index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nsys)
  nlin_sys = 2*(vpstrct[index_sys].ncomp+1)

  if keyword_set(reorder) then begin
     ;; Change order of absorbers
     num0 = sort(randomu(oseed,nsys))
     if keyword_set(nabstot) then begin
        ;; Minimize looping by trunctating
        if nabstot lt nsys then nsys = nabstot
     endif                      ; nabstot=
  endif                         ; /reorder
  
  if keyword_set(doppler) then begin
     if keyword_set(reorder) then begin
        print,'sdss_genprof_setnewz(): WARNING! re-ordering is slow'
        for vv=0L,nsys-1 do begin
           ss = num0[vv] 
           sub = index_sys[ss] + lindgen(nlin_sys[ss])
           if vv eq 0 then newvpstr = vpstrct[sub] $
           else newvpstr = [newvpstr, vpstrct[sub]]
        endfor                  ; loop vv=nsys
     endif else newvpstr = vpstrct
     sub = where(newvpstr.id_lin eq 1,nsub) ; wvI
     num = randomu(oseed,nsub)
     newvpstr[sub].b = num*(zlim[1]-zlim[0]) + zlim[0] ; where zlim now blim
     newvpstr[sub+1].b = newvpstr[sub].b               ; match

  endif else begin
     ;; Scramble the redshifts (or Doppler parameters) for the given
     ;; redshift range Fix rest EW
     num = randomu(oseed,nsys)
     z_sys_new = num*(zlim[1]-zlim[0]) + zlim[0]

     ;; Loop and inject all profiles one-by-one to minimize overlap
     ;; (still will have overlap with real values)
     for vv=0L,nsys-1 do begin
        if keyword_set(reorder) then ss = num0[vv] $
        else ss = vv
        
        sub = index_sys[ss] + lindgen(nlin_sys[ss])
        subvpstr = vpstrct[sub]

        ;; DRM fix: ss reorders so vv and ss a bit scrambled
        z_new_over_old = (1. + z_sys_new[vv]) / (1. + vpstrct[index_sys[ss]].z_sys)
        
        ;; Set components first; preserve velocity structure
        subvpstr.zabs = (subvpstr.zabs - subvpstr[0].z_sys) * $
                        z_new_over_old + z_sys_new[vv]
        subvpstr.z_sys = z_sys_new[vv] ; by fiat
        
        subvpstr.wvlim_sys = subvpstr.wvlim_sys * z_new_over_old ; observed frame
        
        if vv eq 0 then newvpstr = subvpstr $
        else newvpstr = [newvpstr, subvpstr]
     endfor                     ; loop vv=nsys
  endelse
  
  newindex = where(newvpstr.id_comp eq 0 and newvpstr.id_lin eq 1,nsys)

  ;; DRM: Check:
  bd=where(((1.+newvpstr.zabs)*newvpstr.wrest lt newvpstr.wvlim_sys[0]) $
           or ((1.+newvpstr.zabs)*newvpstr.wrest gt newvpstr.wvlim_sys[1]),nbd)
  if nbd ne 0 then $
     stop,'sdss_genprof_setnewz() stop: components outside system wavelength bounds.'

  oseed = oseed[0]
  return, newvpstr
end                             ; sdss_genprof_setnewz()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_reorder, vpstrct, index_sys, seed=seed, oseed=oseed, $
                               newindex=newindex, nabstot=nabstot
  ;; Rearrange order of vpstrct array; better to do it while sorting
  if n_params() lt 1 then begin
     print,'Syntax - sdss_genprof_reorder(vpstrct,[index_sys, seed=, oseed=, newindex=,nabstot=])'
     return,-1
  endif

  if keyword_set(seed) then oseed = seed
  if not keyword_set(index_sys) then $
     index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nsys) $
  else nsys = (size(index_sys,/dim))[0] > 1
  nlin_sys = 2*(vpstrct[index_sys].ncomp+1)

  num = sort(randomu(oseed,nsys))
;  newindex = index_sys[num]     ; for return
  if keyword_set(nabstot) then begin
     if nabstot lt nsys then nsys = nabstot
;     newindex = newindex[0:nsys-1] ; for return??? 
  endif 
 
  for ss=0L,nsys-1 do begin
     nn = num[ss]
     rng = index_sys[nn]+lindgen(nlin_sys[nn])
     if ss gt 0 then $
        newvpstr = [newvpstr,vpstrct[rng]] $
     else newvpstr = vpstrct[rng]
  endfor                        ; loop ss=nsys

  newindex = where(newvpstr.id_comp eq 0 and newvpstr.id_lin eq 1) ; for return
  oseed = oseed[0]
  return, newvpstr
end                             ; sdss_genprof_reorder()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_voigt, wave, vpstrct, nosmooth=nosmooth
  ;; Take the input structure and make the actual profile
  if n_params() ne 2 then begin 
     print,'Syntax - sdss_genprof_voigt(wave, vpstrct, [\nosmooth])'
     return,-1
  endif 
  cinv = 3.3356410e-06          ; km^-1 s = 1./299792.4581d

  if keyword_set(nosmooth) then begin
     ;; Just want to call this (see sdss_genprof_calcew())
     vpfx = x_voigt(wave,vpstrct,fwhm=1,/nosmooth,_extra=extra)
     return, vpfx
  endif

;  npix = (size(wave,/dim))[0]
;
;  ;; need higher resolution wavelength array for x_voigt(); use
;  ;; it's version 
;  nsub = round((alog10(wave[npix-1]) - alog10(wave[0])) / 1.449E-6) + 1
;  subwv = 10^(alog10(wave[0]) + dindgen(nsub)*1.449E-6)
;
;  subvpfx0 = x_voigt(subwv,vpstrct,fwhm=1,/nosmooth,_extra=extra)
; 
;  ;; Better to smooth with convol, assuming constant resolution (not
;  ;; quite right) and knowing constant in velocity space; 162 to 136
;  ;; km/s 
;  ;; In spSpec*.fit 6th extension, there is a dispersion
;  ;; structure, which median resolution in pixels is 0.877945 for most
;  ;; spectra (so that fraction times the 69 km/s log binning)
;  nkpix = 31                    ; some odd number, 15*69 km/s = 1035 km/s
;  if config_strct.dvelo gt 0. then begin
;     velo = findgen(nkpix)*config_strct.dvelo  ; km/s
;     kernel = exp(-0.5*((velo - config_strct.dvelo*0.5*nkpix)/(20*69.))^2)
;  endif else begin
;     velo = findgen(nkpix)*152.5 ; ave 169 km/s and 136 km/s
;     kernel = exp(-0.5*((velo - 152.5*0.5*nkpix)/(20*69.))^2)
;  endelse 
;  kernel = kernel/total(kernel) ; normalize
;  subvpfx = convol(subvpfx0,kernel,/edge_truncate) 
;
;  ;; Put back to SDSS resolution; this looks like 
;  vpfx = rebin_spectrum(subvpfx,subwv,wave) ; could be slow !!!!
;  vpfx = vpfx / max(vpfx)                   ; "Normalize"
;  ;; Fix edges
;  vpfx[0] = vpfx[1]
;  vpfx[npix-1] = vpfx[npix-2]

  vpfx = x_voigt(wave,vpstrct,fwhm=2.2,_extra=extra)

  return,vpfx
end                             ; sdss_genprof_voigt()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_calcew, vpstrct_fil, config_fil, nosort=nosort
  ;; Purely generic plot of profiles
  ;; can also just instantiate EW
  ;; This is slow and almost can't help it
  if n_params() ne 2 then begin
     print,'Syntax - sdss_genprof_calcew(vpstrct_fil, config_fil)'
     return,-1
  endif 
  sdssdir = sdss_getsdssdir()

  if size(vpstrct_fil,/type) eq 7 then $
     vpstrct = xmrdfits(vpstrct_fil,1,/silent) $
  else vpstrct = vpstrct_fil
  if not keyword_set(nosort) then $
     vpstrct = sdss_srtvpstrct(vpstrct) ; will NOT be sorted

  if size(config_fil,/type) eq 7 then $
     config_strct = sdss_genprof_config(config_fil) $
  else config_strct = config_fil 

  ;; Generate high-resolution fiducial SDSS wavelength array
  wvmnx = alog10(sdss_getspecwave() + [-1000,1000]) ; extra width
  pixscale = 1.449E-6
  invpixscale = 690131.11
  npix = round((wvmnx[1] - wvmnx[0]) * invpixscale) + 1
  dpix = dindgen(npix)
  wave = 10.^(wvmnx[0] + dpix*pixscale)
  dwave = 10.^wvmnx[0]*(10.^((dpix+1)*pixscale) - 10.^(dpix*pixscale))

  index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nsys)

  for jj=0L,nsys-1 do begin
     nlin = vpstrct[index_sys[jj]].ncomp + 1
     sub = index_sys[jj] + lindgen(2*nlin)
     subvpstr = vpstrct[sub]

     rng = 2*lindgen(nlin)

     ;; Faster to calculate 
     iminI = round(( alog10(subvpstr[0].wvlim_sys[0])-wvmnx[0] ) * $
                   invpixscale ) ; wvI lower
     imaxI = round(( alog10(subvpstr[0].wvlim_sys[1])-wvmnx[0] ) * $
                   invpixscale ) ; wvI upper
     iminII = round(( alog10(subvpstr[1].wvlim_sys[0])-wvmnx[0] ) * $
                    invpixscale ) ; wvII lower
     imaxII = round(( alog10(subvpstr[1].wvlim_sys[1])-wvmnx[0] ) * $
                    invpixscale ) ; wvII upper
     wvmid = mean(subvpstr[0:1].wrest)*(1+subvpstr[0].z_sys)
     imid = round(( alog10(wvmid) - wvmnx[0] ) * invpixscale) ; midpoint

     ;; expand just a bit
     ;; what about rounding errors? should do partial pixels?
     pixmin = (iminI - 10) > 0
     pixmax = (imaxII + 10) < npix - 1
     
     ;; Previously had profiles generated for each line seperately but
     ;; physical world would never have one without the other;
     ;; so it's just a matter of where bounds are drawn.
     ;; Generate and calc EW for each line (unblended).
     vpfx = sdss_genprof_voigt(wave[pixmin:pixmax],subvpstr,/nosmooth)

     ;; WVI
     ;; assuming pixels defined on LHS
     ;; There will be overlap of EW because the profiles can blend
     ;; into each other
     subvpstr[rng].ew_sys = total((1. - vpfx[iminI-pixmin:imaxI-pixmin]) * $
                                  dwave[iminI:imaxI])
     ;; WVII
     subvpstr[rng+1].ew_sys = total((1. - vpfx[iminII-pixmin:imaxII-pixmin]) * $
                                    dwave[iminII:imaxII])

     ;; OBSERVED measurement, forcing no overlap so stop/start at
     ;; midpoint if wanting to avoid too much blending troubles
     imaxI = imaxI < imid
     iminII = iminII > imid
     subvpstr[rng].ew_obs = total((1. - vpfx[iminI-pixmin:imaxI-pixmin]) * $
                                  dwave[iminI:imaxI])
     subvpstr[rng+1].ew_obs = total((1. - vpfx[iminII-pixmin:imaxII-pixmin]) * $
                                    dwave[iminII:imaxII])

     ;; REST EW
     subvpstr.ew_sys = subvpstr.ew_sys / (1 + subvpstr.z_sys) 
     subvpstr.ew_obs = subvpstr.ew_obs / (1 + subvpstr.z_sys) 

     vpstrct[sub] = subvpstr
  endfor                        ; loop jj=nsys

  return, vpstrct               ; now sorted

end                             ; sdss_genprof_calcew()



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof_mc, dblt_name, config_strct, seed=seed, oseed=oseed,$
                          save_each=save_each, count=count, debug=debug
  ;; Generate all the structures necessary
  if n_params() ne 2 then begin
     print,'Syntax - sdss_genprof_mc(dblt_name, config_strct, [save_each=, seed=, oseed=,count=,/debug])'
     return,-1
  endif 

  ;; Parameters
  if keyword_set(seed) then oseed = seed
  if size(dblt_name,/type) eq 8 then $
     dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  c = 299792.4581d
  cinv = 3.3356410e-06          ; km^-1 s = 1./299792.4581d
  if config_strct.nlim[0] gt 20. then begin
     if keyword_set(debug) then $
        print,'sdss_genprof_mc() debug: config_strct.nlim already linear'
     linear_nlim = 1 
  endif else $
     config_strct.nlim = 10.^config_strct.nlim ; linear space

  rng = 2*lindgen(config_strct.nabs) ; to get both doublets
  count = config_strct.nabs          ; first

  vp_tmpltI = sdss_genprof_setline(dblt.wvI,id_lin=1) ; 0 to config_strct.ndv
  vp_tmpltII = sdss_genprof_setline(dblt.wvII,id_lin=2)
  vpstrct = replicate(vp_tmpltI,config_strct.nabs*2)  ; each line
  vpstrct[rng+1] = vp_tmpltII

  ;; Set System ID numbers
  vpstrct[rng].id_sys = lindgen(config_strct.nabs) ; "system" number
  vpstrct[rng+1].id_sys = vpstrct[rng].id_sys 

  ;; Generate and store random z, N, b for each doublet line
  num = randomu(oseed,config_strct.nabs*3,/double) ; want larger swath
  inum = 0L
  num2 = num[0:config_strct.nabs-1]
  inum = inum + config_strct.nabs
  vpstrct[rng].zabs = num2*(config_strct.zlim[1]-config_strct.zlim[0]) $
                      + config_strct.zlim[0]
  vpstrct[rng+1].zabs = vpstrct[rng].zabs

  num2 = num[inum:(inum+config_strct.nabs-1)]   
  inum = inum + config_strct.nabs
  vpstrct[rng].n = num2*(config_strct.nlim[1]-config_strct.nlim[0]) $
                   + config_strct.nlim[0]
  vpstrct[rng+1].n = vpstrct[rng].n

  num2 = num[inum:(inum+config_strct.nabs-1)] ; independent 
  inum = inum + config_strct.nabs ; scale with column if num2 not reset
  vpstrct[rng].b = num2*(config_strct.blim[1]-config_strct.blim[0]) $
                   + config_strct.blim[0]
  vpstrct[rng+1].b = vpstrct[rng].b


  if keyword_set(config_strct.ndv) then begin
     ;; ;;;;;;;
     ;; Add components
     if not keyword_set(config_strct.dvabs) then $
        stop,'sdss_genprof_mc(): must set config_strct.ndv and config_strct.dvabs[]'

     vpstrct[rng].ncomp = randomu(oseed,config_strct.nabs,$
                                  poisson=config_strct.ndv) 
     vpstrct[rng+1].ncomp = vpstrct[rng].ncomp
     
     sub = where(vpstrct[rng].ncomp gt 0,nsub) ; where to add?; rng limits to wvI
     if nsub gt 0 then begin
        ;; Create array
        sub = rng[sub]
        ncomptot = total(vpstrct[sub].ncomp)

        vpstrct_sub = replicate(vp_tmpltI,ncomptot*2) ; wvI, wvII
        rng = 2*lindgen(ncomptot)
        vpstrct_sub[rng+1] = vp_tmpltII ; instantiate

        num = randomu(oseed,ncomptot*3,/double) ; N, b, z
        if keyword_set(config_strct.flg_dvabs) then begin
           gnum = randomn(oseed,ncomptot*3,/double) ; N, b, z
           if keyword_set(debug) then $
              print,'sdss_genprof_mc() debug: using normally distributed dvabs'
        endif 
        inum = 0L

        ;; Loop over systems to add components too
        for ii=0L,nsub-1 do begin
           ss = sub[ii]
           id_comp = lindgen(vpstrct[ss].ncomp) 
           rng = 2*id_comp
           if ii gt 0 then $
              rng = rng + 2*total(vpstrct[sub[0:ii-1]].ncomp)
           vpstrct_sub[rng].id_sys = vpstrct[ss].id_sys ; duh!
           vpstrct_sub[rng+1].id_sys = vpstrct_sub[rng].id_sys
           vpstrct_sub[rng].id_comp = id_comp + 1
           vpstrct_sub[rng+1].id_comp = vpstrct_sub[rng].id_comp

           ;; Set dv offset
           ;; z = dv/c*(1+z0) + z0
           if keyword_set(gnum) then begin
              ;; Already determined in config_strct to be Gaussian
              ;; distribution around primary component
              num2 = gnum[inum:(inum+vpstrct[ss].ncomp-1)] ; Gaussian
              vpstrct_sub[rng].zabs = num2*config_strct.dvabs[1]*cinv*$
                                      (1 + vpstrct[ss].zabs) + vpstrct[ss].zabs
           endif else begin
              ;; Use uniformly distributed around primary component
              num2 = num[inum:(inum+vpstrct[ss].ncomp-1)]
              vpstrct_sub[rng].zabs = $
                 (num2*(config_strct.dvabs[1]-config_strct.dvabs[0]) $
                  + config_strct.dvabs[0])*cinv*$
                 (1 + vpstrct[ss].zabs) + vpstrct[ss].zabs
           endelse 
           vpstrct_sub[rng+1].zabs = vpstrct_sub[rng].zabs
           inum = inum+vpstrct[ss].ncomp ; access next chunk of random numbers
           
           ;; Total column originally selected is preserved; divided
           ;; between other components
           num2 = num[inum:(inum+vpstrct[ss].ncomp-1)] ; random number range for column
           inum = inum+vpstrct[ss].ncomp               ; make num2 "used" even if not really used 
           if keyword_set(config_strct.flg_capn) then begin
              jj = 0
              while jj lt vpstrct[ss].ncomp do begin
                 vpstrct_sub[rng[jj]].n = $
                    num2[jj]*(vpstrct[ss].n-config_strct.nlim[0]) $
                    + config_strct.nlim[0]
                 ;; Subtract it off to preserve total column
                 if vpstrct_sub[rng[jj]].n lt vpstrct[ss].n then begin
                    vpstrct[ss].n = vpstrct[ss].n - vpstrct_sub[rng[jj]].n
                    ;; But stay above min
                    if vpstrct[ss].n lt config_strct.nlim[0] then begin
                       vpstrct_sub[rng[jj]].n = 0.  ; exclude later
                       vpstrct[ss].n += vpstrct_sub[rng[jj]].n ; undo for realz (DRM)
                       bad = 1
                    endif else bad = 0
                 endif else bad = 1

                 if bad then begin 
                    if keyword_set(debug) then $
                       print,'sdss_genprof_mc() debug: truncating no. of components: ',$
                             string(vpstrct[ss].zabs,jj,vpstrct[ss].ncomp,$
                                    format='(f7.5,1x,i3,1x,i3)')
                    vpstrct[ss].ncomp = jj ; end while loop; one before
                 endif
                 jj++
              endwhile          ; loop jj < ncomp[ss]
           endif else begin     ; keyword_set(config_strct.flg_capn)
              vpstrct_sub[rng].n = num2*(config_strct.nlim[1]-config_strct.nlim[0]) $
                                   + config_strct.nlim[0]
           endelse              ; don't cap ncolm

           ;; Set Doppler; if locked (meaning no new num2), must come
           ;; after ncolm
           vpstrct_sub[rng].b = num2*(config_strct.blim[1]-config_strct.blim[0]) $
                                + config_strct.blim[0]
           
           ;; Save
           vpstrct[ss+1].n = vpstrct[ss].n ; DRM: sub instead of ss is overkill
           vpstrct_sub[rng+1].b = vpstrct_sub[rng].b
           vpstrct[ss+1].ncomp = vpstrct[ss].ncomp ; if not set earlier
           vpstrct_sub[rng+1].n = vpstrct_sub[rng].n
           vpstrct_sub[rng].ncomp = vpstrct[ss].ncomp ; possibly trimmed above
           vpstrct_sub[rng+1].ncomp = vpstrct_sub[rng].ncomp
        endfor                  ; loop ii=nsub

        ;; Append ones that have components instantiated
        gd = where(vpstrct_sub.n gt 0.,ngd)
        if gd[0] ne -1 then $
           vpstrct = [vpstrct,vpstrct_sub[gd]] ; don't sort!
        count = count + ngd / 2                ; update
     endif                      ; vpstrct.ncomp > 0

  endif                         ; config_strct.ndv > 0

  ;; Boundaries for little bits assuming +/-5 "sigma" (= b/c)
  wobs = vpstrct.wrest*(1+vpstrct.zabs)
  wvlimlo = wobs * (1 - 5*vpstrct.b*cinv)
  wvlimhi = wobs * (1 + 5*vpstrct.b*cinv)

  ;; Change to boundaries to be for system as a whole and since
  ;; components can be completely to blue or red of others, do min/max
  ;; on both limits of components
  index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nsys)
  for ii=0L,nsys-1 do begin
     ;; not sorted; and this will fail if I continue to have duplicate
     ;; id_sys problem
     rng = where(vpstrct[index_sys[ii]].id_sys eq vpstrct.id_sys and $
                 vpstrct.id_lin eq 1,nrng)
     if nrng ne vpstrct[index_sys[ii]].ncomp+1 then $
        stop,'sdss_genprof_mc() stop: WARNING! dupilcate id_sys'

     ;; WVI
     vpstrct[rng].wvlim_sys[0] = min(wvlimlo[rng])
     vpstrct[rng].wvlim_sys[1] = max(wvlimhi[rng])
     
     ;; redshift (rough flux weighting)
     wgt = vpstrct[rng].n ; DRM: b/c 1. - exp(-vpstrct[rng].n) ~ 1 (linear columns)
     vpstrct[rng].z_sys = total(vpstrct[rng].zabs*wgt)/total(wgt)
     vpstrct[rng+1].z_sys = vpstrct[rng].z_sys

     ;; WVII
     vpstrct[rng+1].wvlim_sys[0] = min(wvlimlo[rng+1])
     vpstrct[rng+1].wvlim_sys[1] = max(wvlimhi[rng+1])

     ;; Total column
     vpstrct[rng].n_sys = total(vpstrct[rng].n) ; already linear
     vpstrct[rng+1].n_sys = vpstrct[rng].n_sys
  endfor                        ; loop ii=maxid 


  ;; ;;;;;;;
  ;; EW: n needs to be log in sdss_genprof_calcew() and linear here
  ;; Comes back sorted
  config_strct.nlim = alog10(config_strct.nlim)
  vpstrct.n = alog10(vpstrct.n)
  vpstrct.n_sys = alog10(vpstrct.n_sys)
  vpstrct = sdss_genprof_calcew(vpstrct, config_strct); slow!
  config_strct.nlim = 10.^config_strct.nlim
  vpstrct.n = 10.^vpstrct.n
  vpstrct.n_sys = 10.^vpstrct.n_sys

  
  if keyword_set(config_strct.dvlim) then begin
     ;; ;;;;;;;
     ;; Efficient trimming of things too spread out
     cfg = config_strct
     cfg.dvlim = 0              ; make sure not to get here recursively
     
     ;; Remove things not within bounds and recursively call
     ;; sdss_genprof_mc() 
     dvlim = (vpstrct.zabs - vpstrct.z_sys) / (1 + vpstrct.z_sys) * c
     bd = where(abs(dvlim) gt config_strct.dvlim or $
                vpstrct.ew_obs lt config_strct.ewlim[0] or $
                vpstrct.ew_obs gt config_strct.ewlim[1],nbd)  
     if keyword_set(debug) then $
        print,'sdss_genprof_mc() debug: Components with width and/or EW out of bounds',$
              nbd,' of',2*count

     while nbd ne 0 do begin
        ;; New sample
        cfg.nabs = (10 * nbd) < 20000L ; just make more
        newvpstrct = sdss_genprof_mc(dblt_name,cfg,seed=oseed,oseed=oseed,$
                                     save_each=save_each,debug=debug)
        ;; Efficiently select from these
        newdvlim = (newvpstrct.zabs - newvpstrct.z_sys) / $
                   (1 + newvpstrct.z_sys) * c
        newgd = where(abs(newdvlim) le config_strct.dvlim and $
                      newvpstrct.ew_obs ge config_strct.ewlim[0] and $
                      newvpstrct.ew_obs le config_strct.ewlim[1])
        if newgd[0] eq -1 then begin
           if keyword_set(debug) then $
              stop,'sdss_genprof_mc() debug: no good new profiles'
           continue             ; fail; try drawing again
        endif 

        ;; Select the right newvpstrct systems by their ID
        ;; should be sorted
        unq = uniq(newvpstrct[newgd].id_sys) ;,sort(newvpstrct[newgd].id_sys))
        newgdid_sys = newvpstrct[newgd[unq]].id_sys
        nnewgd = (size(newgdid_sys,/dim))[0] > 1 ; foil singularity

        ;; Figure out stuff to excise and something good to put back
        ;; in with the right id_sys
        ;; preserve bdid_sys because vpstrct going to change, shuffle
        ;; have to use uniq() here b/c only one component could be
        ;; off!!! 
        unq = uniq(vpstrct[bd].id_sys,sort(vpstrct[bd].id_sys))
        bdid_sys = vpstrct[bd[unq]].id_sys
        nbdsys = (size(bdid_sys,/dim))[0] > 1 
        if keyword_set(debug) then $
           print,bdid_sys

        bb = 0L
        while (bb lt nbdsys and bb lt nnewgd) do begin
           sav = where(bdid_sys[bb] ne vpstrct.id_sys)
           newsub = where(newgdid_sys[bb] eq newvpstrct.id_sys)
           tmpstrct = newvpstrct[newsub]
           tmpstrct.id_sys = bdid_sys[bb] ; crucial to avoid confusion
           vpstrct = [vpstrct[sav],tmpstrct]
           bb++
        endwhile                 ; loop bb < nbdsys and bb < nnewgd

        ;; Set up new check
        dvlim = (vpstrct.zabs - vpstrct.z_sys) / (1 + vpstrct.z_sys) * c
        bd = where(abs(dvlim) gt config_strct.dvlim or $
                   vpstrct.ew_obs lt config_strct.ewlim[0] or $
                   vpstrct.ew_obs gt config_strct.ewlim[1],nbd)
        
        if keyword_set(debug) then $
           print,'                           ... new number out of bounds',$
                 nbd,' of',2*count

     endwhile                   ; nbd ne 0

     count = (size(vpstrct,/dim))[0] / 2 ; number of components
  endif                         ; config_strct.dvlim != 0 

  ;; Log space
  if not keyword_set(linear_nlim) then begin
     config_strct.nlim = alog10(config_strct.nlim)
     vpstrct.n = alog10(vpstrct.n)
     vpstrct.n_sys = alog10(vpstrct.n_sys)
  endif else $
     if keyword_set(debug) then $
        print,'sdss_genprof_mc() debug: returning linear column densities in vpstrct'
  ;; need to do this because of the recursive call that eliminates
  ;; profiles too wide

  oseed = oseed[0]              ; not the array

  ;; Now crucial to have this sorted
  srtstrct = sdss_srtvpstrct(vpstrct)

  if keyword_set(save_each) then begin
     if size(save_each,/type) eq 7 then ofil = save_each $
     else ofil = 'sdss_genprof_mc_intermediate.fit'
     mwrfits,srtstrct,ofil,/create,/silent
     mwrfits,config_strct,ofil,/silent
     print,'sdss_genprof_mc(): saved intermediate ',ofil
  endif 

;stop
  return, srtstrct
end                             ; sdss_genprof_mc()



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_genprof_plot, vpstrct, civ_fil=civ_fil, ewlim=ewlim
  ;; Purely generic plot of profiles
  ;; can also just instantiate EW
  if n_params() ne 1 then begin
     print,'Syntax - sdss_genprof_plot, vpstrct, [civ_fil=]'
     return
  endif 
  clr = getcolor(/load)
  sdssdir = sdss_getsdssdir()

  if keyword_set(civ_fil) then begin
     if size(civ_fil,/type) eq 7 then $
        civstr = xmrdfits(civ_fil,1,/silent) $
     else civstr = civ_fil
  endif 

  ;; Generate fiducial SDSS wavelength array
  wvmnx = sdss_getspecwave()
  pixscale = sdss_getspecpixscale(/loglam)
  npix = 4000
  wave = wvmnx[0]*10.^(dindgen(npix)*pixscale/alog(10))

  index_sys = where(vpstrct.id_comp eq 0 and vpstrct.id_lin eq 1,nsys)
  if keyword_set(ewlim) then begin
     case n_elements(ewlim) of
        1: gd = where(vpstrct[index_sys].ew_sys ge ewlim,ngd)
        2: gd = where(vpstrct[index_sys].ew_sys ge ewlim[0] and $
                      vpstrct[index_sys].ew_sys lt ewlim[1],ngd)
        else: stop,'sdss_genprof_plot stop: n_elements(ewlim) invalid'
     endcase
     if ngd eq 0 then $
        print,'sdss_genprof_plot: WARNING! no profiles meeting ewlim cut; using all' $
     else begin
        ;; Trim
        nsys = ngd
        index_sys = index_sys[gd]
     endelse
  endif 

  for jj=0L,nsys-1 do begin
     nlin = vpstrct[index_sys[jj]].ncomp + 1
     rng = 2 * lindgen(nlin)
     sub = index_sys[jj] + [rng, rng+1]

     mn = min(wave-vpstrct[sub[0]].wvlim_sys[0],pixmin,/abs)    ; I
     mn = min(wave-vpstrct[sub[nlin]].wvlim_sys[1],pixmax,/abs) ; II

     ;; expand
     pixmin = (pixmin - 50) > 0
     pixmax = (pixmax + 50) < npix - 1

     ;; Generate and plot for just the small region desired

     vpfx = sdss_genprof_voigt(wave[pixmin:pixmax],vpstrct[sub])

     wvr = wave[pixmin:pixmax] / (1. + vpstrct[sub[0]].z_sys)
     lgnd = 'MC VP'
     ttl = string(vpstrct[sub[0]].z_sys,vpstrct[sub[0]].ew_obs,$
                  vpstrct[sub[nlin]].ew_obs,$
                  format='("zvp=",f7.5,1x,"(EWI=",f5.2,"; EWII=",f5.2,")")')
     if keyword_set(civ_fil) then begin
        ;; Find closest on based on EW 
        mn = min(civstr.ew_orig[0]-vpstrct[sub[0]].ew_obs,imn,/abs)
        parse_sdss,sdssdir+civstr[imn].sdss_obs[0],fxtwo,wvobs,sig=sigthr
        wvtwo = wvobs / (1. + civstr[imn].zabs_orig[0])
        lgnd = [lgnd,civstr[imn].qso_name,'Error'] ; Info

        ;; Normalize
        cstrct = xmrdfits(sdssdir+civstr[imn].abslin_fil,1,/silent)
        cindx = fix(alog(cstrct.cflg)/alog(2))
        fxtwo = fxtwo / cstrct.conti[0:cstrct.npix-1,cindx]
        sigthr = sigthr / cstrct.conti[0:cstrct.npix-1,cindx]

        ;; Info to title
        ttl = ttl + ' Compared to ' + $
              string(civstr[imn].zabs_orig[0],$
                     civstr[imn].ew_orig[0:1],$
                     format='("zobs=",f7.5,1x,"(EWI=",f5.2,"; EWII=",f5.2,")")')
     endif else begin
        ;; don't plot anything but profile (less useful)
        wvtwo = 0
        fxtwo = 0
        sigthr = 0
     endelse 

     ;; Plot
     x_splot,wvr,vpfx,psym1=10,color1=clr.blue,$
             xtwo=wvtwo,ytwo=fxtwo,psym2=10,color2=clr.charcoal,$
             xthr=wvtwo,ythr=sigthr,psym3=10,color3=clr.red,$
             lgnd=lgnd,title=ttl,ymnx=[-0.05,1.2],/block
     
        
     ;; Option to exit
     if jj ne 0 and (jj mod 10) eq 0 then $
        stop,'sdss_genprof_plot stop: enter "return" now to exit; .c to continue'
  endfor                        ; loop ii=mxid

end                             ; sdss_genprof_plot



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_genprof_compare,list_fil,sdsssum,allciv_fil,dblt_name=dblt_name,$
                         clobber=clobber,debug=debug,istrt=istrt,$
                         in_fil=in_fil,plot=plot,_extra=extra
  ;; Process results of sdss_genprof_sanitycheck
  if n_params() ne 3 then begin
     print,'Syntax - sdss_genprof_compare, list_fil, sdssssum, allciv_fil,'
     print,'                               [dblt_name=, /clobber, /debug, '
     print,'                               istrt=,in_fil=,/plot,_extra=]'
     return
  endif 


  test = file_search(allciv_fil+'*',count=ntest)
  if ntest ne 0 and not keyword_set(clobber) then begin
     stop,'sdss_genprof_compare stop: file already exists but /clobber not set'
  endif 

  if not keyword_set(istrt) then istrt = 0L ; file
  if keyword_set(in_fil) then begin
     allciv = xmrdfits(in_fil,1,/silent)
     print,'sdss_genprof_compare: re-starting from ',in_fil
  endif 

  sdssdir = sdss_getsdssdir()

  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  
  ;; The original structure was suppose to be saved to the 2nd
  ;; extension of all the files but that didn't work so here's just
  ;; the structure
  ;; Need to correct all EW_ORIG because never needed to change them
  ;; in sdss_genprof_setnewz()
  vpstrct0 = xmrdfits(sdssdir+'civ/vpstrct50_sanitycheck.fit',1,/silent) 
  index_sys0 = where(vpstrct0.id_comp eq 0 and vpstrct0.id_lin eq 1,nvpstrct0)

  c = 299792.4581d
  cinv = 1./c                   ; km^-1 s 

  ;; Read in files
  readcol,list_fil,spec_fil,format='(a)',skip=1,/silent
  nfil = (size(spec_fil,/dim))[0] > 1 ; foils singularity
  sdsstab = xmrdfits(sdsssum,1,/silent)
  dblt_fil = sdss_getname(sdsstab,root=qso_name,mc=dblt.ion,dir=mcdir)
  dblt_fil = mcdir+dblt_fil

  ;; Want to know the results but would be exceedingly large to
  ;; concatenate the structures so just extract the values desired
  for ff=istrt,nfil-1 do begin
     test = file_search(sdssdir+dblt_fil[ff]+'*',count=ntest)
     if ntest eq 0 then continue ; SKIP b/c file DNE


     if keyword_set(in_fil) then begin
        ;; Check sightline not already in file before reading
        mtch = where(allciv.qso_name eq qso_name[ff])
        if mtch[0] ne -1 then begin
           print,'sdss_genprof_compare: QSO already concatenated ',qso_name[ff]
           continue             ; SKIP b/c QSO already appended
        endif                   ; mtch[0] ne -1
     endif                      ; in_fil= set

     ;; read in
     civstr = xmrdfits(sdssdir+dblt_fil[ff],1,/silent)

     ;; fix the problem with the old values not needing to be
     ;; corrected because rest EW is rest EW
     civstr.ew_orig[0] = vpstrct0[index_sys0].ew_obs
     civstr.ew_orig[1] = vpstrct0[index_sys0+1].ew_obs

     ;; Concatenate
     if keyword_set(allciv) then begin
        allciv = [allciv,civstr]
     endif else begin
        allciv = civstr
     endelse

     ;; Save intermediate
     if (ff mod 500) eq 0 and keyword_set(allciv) then begin
        mwrfits,allciv,allciv_fil,/create,/silent
        print,'sdss_genprof_compare: saved up to ff = ',ff
     endif 

  endfor                        ; loop ff=nfil
  
  ;; Write out (/clobber issue handled at the beginning)
  vpstrct0 = xmrdfits(sdssdir+dblt_fil[0],2,/silent) ; z and EW off
  mwrfits,allciv,allciv_fil,/create,/silent
  mwrfits,vpstrct0,allciv_fil,/silent
  spawn,'gzip -f '+allciv_fil
  print,'sdss_genprof_compare: created ',allciv_fil

  if keyword_set(plot) then begin
;     clr = getcolor(/load)
     x_splot,allciv.ew_orig[0],allciv.ew_final[0],psym1=3,$ ; period
             xtwo=[-10,10],ytwo=[-10,10],$
             xtitle='Original EW (Ang)',ytitle='Recovered EW (Ang)',$
             title='',/block
  endif 
  
end                             ;  sdss_genprof_compare


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_genprof_sanitycheck,list_fil,sdsssum,dblt_name=dblt_name,tcpu=tcpu,$
                             nabs=nabs,config_fil=config_fil,processor=processor,$
                             vpstrct_fil=vpstrct_fil,seed=seed,oseed=oseed,$
                             clobber=clobber,debug=debug,nosort=nosort,_extra=extra
  ;; Take the same systems and inject them into the given spectra
  ;; (list_fil) and see what the recovered parameters look like
  if n_params() ne 2 then begin
     print,'Syntax - sdss_genprof_sanitycheck,list_fil,sdsssum,[dblt_name=,'
     print,'                                  nabs=,config_fil=,processor=,'
     print,'                                  vpstrct_fil=,seed=,oseed=,'
     print,'                                  /clobber,/debug,_extra=]'
     return
  endif 

  sdssdir = sdss_getsdssdir()

  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  if not keyword_set(config_fil) then $
     config_fil = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civprof.config'
  if size(config_fil,/type) eq 8 then config_strct = config_fil $
  else config_strct = sdss_genprof_config(config_fil)
  if keyword_set(seed) then oseed = seed $
  else oseed = 1321.2528492     ; reproducibility; esp. for parallel jobs

  ;; Change params
  if keyword_set(nabs) then config_strct.nabs = nabs
  c = 299792.4581d
  cinv = 1./c                   ; km^-1 s 

  ;; Read in files
  readcol,list_fil,spec_fil,format='(a)',skip=1,/silent
  nfil = (size(spec_fil,/dim))[0] > 1 ; foils singularity
  sdsstab = xmrdfits(sdsssum,1,/silent)
  abslin_fil = sdss_getname(sdsstab,/abslin,root=qso_name,dir=absdir)
  abslin_fil = absdir + abslin_fil 
  dblt_fil = sdss_getname(sdsstab,mc=dblt.ion,dir=mcdir)
  if keyword_set(debug) then mcdir = 'test'+mcdir
  dblt_fil = mcdir+dblt_fil

  if keyword_set(processor) then begin
     sub = sdss_calcparalleljob(sdsstab,processor)
     istrt = sub[0]
     nfil = sub[1] + 1
     if keyword_set(debug) then $
        print,'sdss_genprof_sanitycheck debug: multi-processor run for just ',$
              istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif else istrt = 0L         ; default

  ;; Figure out wavelength bounds
  ;; _extra includes dvgal=, dvqso=
  tmp = sdss_measuresnr(spec_fil,wvlim_obs=wvlim_obs,/no_snr,$
                        dblt_name=dblt, zqso=sdsstab.z, _extra=extra)
  zlim_obs = wvlim_obs / dblt.wvI - 1

  ;; Generate profiles
  if keyword_set(vpstrct_fil) then begin
     if size(vpstrct_fil,/type) eq 7 then $
        vpstrct = xmrdfits(vpstrct_fil,1,/silent) $
     else vpstrct = vpstrct_fil
     if keyword_set(nosort) then $
        vpstrct = sdss_srtvpstrct(vpstrct)
  endif else $
     vpstrct = sdss_genprof_mc(dblt,config_strct,seed=oseed,oseed=oseed,$
                               debug=debug)

  index_sys = where(vpstrct.id_comp eq 0 and $
                    vpstrct.id_lin eq 1,nsys)
  nabs_orig = nabs_orig
  config_strct.nabs = 1         ; for loop

  ;; Where to save and instantiate all necessary
  civstrtmplt = replicate({sdsscivstrctslim},nsys) ; output; keep small
  civstrtmplt.wrest[0] = dblt.wvI
  civstrtmplt.wrest[1] = dblt.wvII
  civstrtmplt.cflg = sdss_getcflg(/hyb) ; would have to change if not true
  cindx = sdss_getcflg(/hyb,/index)
  civstrtmplt.rating[0] = sdss_getrating(/def)
  civstrtmplt.ncolm_orig[0] = vpstrct[index_sys].n_sys
  civstrtmplt.ncolm_orig[1] = vpstrct[index_sys+1].n_sys

  ;; Loop over all spectra and inject the same profiles and measure them
  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  for ff=istrt,nfil-1 do begin
     ;; Figure out if should do this loop
     test = file_search(sdssdir+dblt_fil[ff]+'*',count=ntest)
     if ntest eq 0 or keyword_set(clobber) then begin
        if ntest eq 0 then begin
           test2 = file_search(sdssdir+mcdir[ff],count=ntest2)
           if ntest2 eq 0 then $
              spawn,'mkdir -p '+sdssdir+mcdir[ff]
        endif 
     endif else begin
        print,'sdss_genprof_sanitycheck: skipping ',dblt_fil[ff]
        continue
     endelse 

     ;; Duplicate civstr template and modify for current LOS
     civstr = civstrtmplt
     civstr.qso_name = qso_name[ff]
     civstr.mjd = sdsstab[ff].smjd
     civstr.plate = sdsstab[ff].plate
     civstr.fiber = sdsstab[ff].fiber
     civstr.sdss_obs[0] = spec_fil[ff]
     civstr.z_qso = sdsstab[ff].z
     civstr.abslin_fil = abslin_fil[ff]

     ;; Scramble the redshfits for the given redshift range
     ;; Fix rest EW
     vpstrct = sdss_genprof_setnewz(vpstrct,zlim_obs[ff,*],seed=oseed,oseed=oseed)
     
     ;; Read in with desired format
     spec = xmrdfits(sdssdir+spec_fil[ff], 0, hdr, /silent)
     npix = (size(spec[*,0],/dim))[0] 
     wave = x_setwave(hdr,npix)
     cstrct = xmrdfits(sdssdir+abslin_fil[ff],1,/silent)

     ;; Loop and inject all profiles one-by-one to minimize overlap
     ;; (still will have overlap with real values)
     for vv=0L,nsys-1 do begin
        ;; This only works if vpstrct properly sorted
        sub = index_sys[vv] + lindgen(2*vpstrct[index_sys[vv]].ncomp) 
        subvpstr = vpstrct[sub]

        ;; Instantiate new values
        civstr[vv].ew_orig[0:1] = subvpstr[0:1].ew_obs  ; new rest frame
        civstr[vv].zabs_orig[0:1] = subvpstr[0].z_sys
        civstr[vv].wvlim_orig[0,*] = subvpstr[0].wvlim_sys
        civstr[vv].wvlim_orig[1,*] = subvpstr[1].wvlim_sys

        ;; Most efficient to measure EW here because having to
        ;; generate profiles anyway
        newspec = sdss_genprof(wave,dblt,config_strct,spec,$
                               vpstrct=subvpstr,debug=0,_extra=extra) 

        ;; Technically not the same without using snr_conv=
        ;; save to *_final tags
        tmp = sdss_ewciv(wave[cstrct.ipix0:cstrct.npix-1], $
                         newspec[cstrct.ipix0:cstrct.npix-1,0], $
                         newspec[cstrct.ipix0:cstrct.npix-1,2],$
                         cstrct.conti[cstrct.ipix0:cstrct.npix-1,cindx],$
                         dblt,civstr[vv].zabs_orig[0],$
                         /final,debug=0,istrct=civstr[vv])
        civstr[vv] = tmp

        vpstrct[sub] = subvpstr ; for writing out
     endfor                     ; loop vv=nsys


     ;; Write out
     civstr.cmplt_fil = dblt_fil[ff]

     mwrfits,civstr,sdssdir+dblt_fil[ff],/create,/silent
     mwrfits,vpstrct,sdssdir+dblt_fil[ff],/silent
     mwrfits,config_strct,sdssdir+dblt_fil[ff],/silent
     spawn,'gzip -f '+sdssdir+dblt_fil[ff]
     print,'sdss_genprof_sanitycheck: created ',dblt_fil[ff]

     if keyword_set(tcpu) and  ((ff+1) mod 500) eq 0 then begin
        tt = systime(/seconds)
        dt = (tt-tlast)         ; seconds
        print,'sdss_genprof_sanitycheck: Elapsed time (m) = ',dt/60.
        print,'sdss_genprof_sanitycheck: Average time per qq (s) = ',dt/1000L
        tlast = tt
     endif 

     if keyword_set(debug) then $
        stop,'sdss_genprof_sanitycheck debug stop: end of loop'
  endfor                        ; loop ff=nfil
  
  oseed = oseed[0]

  print, 'sdss_genprof_sanitycheck: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'sdss_genprof_sanitycheck: Elapsed time for '+strtrim(nfil-istrt,2)+$
        ' QSOs (m) = ',dt/60.
  print,'sdss_genprof_sanitycheck: Average time per QSO (s) = ',dt/(nfil-istrt)

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end                             ; sdss_genprof_sanitycheck



function sdss_genprof_getvpsubset, vpstrct, ewbinsize, n_per_bin, $
                                   min=min, max=max, config_fil=config_fil, $
                                   seed=seed, oseed=oseed, $
                                   index=index, count=count, $
                                   duplicate=duplicate, ewlim_dup=ewlim_dup, $
                                   debug=debug
  ;; vpstrct *has* to be sorted
  if n_params() ne 3 then begin
     print,'Syntax - '
     return,-1
  endif 

  sdssdir = sdss_getsdssdir()
  if not keyword_set(min) then begin
     if not keyword_set(config_fil) then $
        config_fil = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civprof.config'
     if size(config_fil,/type) eq 7 then cfg = sdss_genprof_config(config_fil) $
     else cfg = config_fil
     min = cfg.ewlim[0]
  endif 
  if not keyword_set(max) then max = max(vpstrct.ew_obs)
  if keyword_set(seed) then oseed = seed
  if not keyword_set(ewlim_dup) then ewlim_dup = 100. ; Ang; something large
  
  nloop = ceil((max - min) / ewbinsize)
  max_id_sys = max(vpstrct.id_sys)

  for ll=0L,nloop-1 do begin
     sub0 = where(vpstrct.ew_obs ge ll*ewbinsize + min and $
                  vpstrct.ew_obs lt (ll+1)*ewbinsize + min and $
                  vpstrct.id_lin eq 1) ; only wvI
     if sub0[0] eq -1 then begin
        print,'sdss_genprof_getvpsubset(): no profiles in ',$
              string(ll*ewbinsize+min,(ll+1)*ewbinsize+min,$
                     format='(f7.3," <= EW <",F7.3)')
        continue
     endif 
     sub0 = [sub0,sub0+1]       ; vpstrct *has* to be sorted
     sub0 = sub0[sort(sub0)]

     ;; Set up information
     tmpvpstr0 = vpstrct[sub0]
     indx_sys = where(tmpvpstr0.id_comp eq 0 and $
                      tmpvpstr0.id_lin eq 1,nsys0)
     
     if keyword_set(duplicate) and nsys0 lt n_per_bin and $
        min+ll*ewbinsize lt ewlim_dup then begin
        ;; Repeat profiles to make more for this bin
        nmore = ceil(n_per_bin / float(nsys0))
        if keyword_set(debug) then $
           print,'sdss_genprof_getvpsubset() debug: too few, duplicating ',$
                 nsys0, n_per_bin, nmore

        tmpvpstr = tmpvpstr0
        for mm=1,nmore-1 do begin
           nwtmpvpstr = tmpvpstr0
           ;; Have to keep systems spread out
           nwtmpvpstr.id_sys = nwtmpvpstr.id_sys + max_id_sys + 1L
           tmpvpstr = [tmpvpstr, nwtmpvpstr]
           max_id_sys = max(nwtmpvpstr.id_sys)
        endfor                  ; loop mm=nmore
        tmpvpstr0 = tmpvpstr
        indx_sys = where(tmpvpstr0.id_comp eq 0 and $
                         tmpvpstr0.id_lin eq 1,nsys0)
     endif                      ; /duplicate and nsys0 < n_per_bin


     if nsys0 le n_per_bin then begin
        ;; Just take everything
        nwtmpvpstr = tmpvpstr0  
        if keyword_set(debug) then $
           print,min+ll*ewbinsize,-1,nsys0,-1L ; match format below

     endif else begin
        ;; Take random assortment
        indx_sys = indx_sys[sort(randomu(oseed,nsys0))] ; random sorted order

        ;; Either going to have more to exclude or more to grab out; do
        ;; which ever is shortest
        if (nsys0-n_per_bin) gt n_per_bin then begin
           ;; Inclusion faster
           inclusive = 1 
           nsys = 0L            ; will add
        endif else begin
           ;; Exclusion faster
           inclusive = 0        
           nsys = nsys0         ; will subtract
        endelse 
        ss = 0L
        tmpvpstr = tmpvpstr0    ; may change
        while nsys ne n_per_bin do begin
           ;; This is why systems have to be spread out in number
           mtch = where(tmpvpstr.id_sys eq tmpvpstr0[indx_sys[ss]].id_sys,$
                        complement=unmtch)

           if inclusive then begin
              nsys++ 
              if ss eq 0 then sub = mtch $
              else sub = [sub, mtch]
;              print,'Inclusion: ss =',ss
;              printcol,tmpvpstr0[mtch].id_sys,tmpvpstr0[mtch].id_lin,$
;                       tmpvpstr0[mtch].id_comp,tmpvpstr0.ncomp
;              test = sdss_srtvpstrct(tmpvpstr0[sub])
           endif else begin
              nsys--
              tmpvpstr = tmpvpstr[unmtch] ; reduce
;              print,'Exclusion: ss =',ss
;              printcol,tmpvpstr.id_sys,tmpvpstr.id_lin,tmpvpstr.id_comp,$
;                       tmpvpstr.ncomp
;              test = sdss_srtvpstrct(tmpvpstr)
           endelse 
           
           ss++
        endwhile                ; loop nsys < n_per_bin

        if keyword_set(debug) then $
           print,min+ll*ewbinsize,inclusive,nsys0,ss

        if inclusive then nwtmpvpstr = tmpvpstr0[sub] $
        else nwtmpvpstr = tmpvpstr

;        print,'Finished bin'
;        test = sdss_srtvpstrct(nwtmpvpstr)
     endelse 

     if not keyword_set(subvpstr) then subvpstr = nwtmpvpstr $
     else subvpstr = [subvpstr, nwtmpvpstr]
  endfor                        ; loop ll=nloop
  
  oseed = oseed[0]
  srtvpstr = sdss_srtvpstrct(subvpstr) ; this is SLOW
  return, srtvpstr
end                             ; sdss_genprof_getvpsubset()




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof,wave, dblt_name, config_strct, cleanspec,$ 
                      seed=seed, oseed=oseed,$
                      conti=conti, vpstrct=vpstrct, $
                      debug=debug, _extra=extra
  if n_params() ne 4 then begin
     print,'Syntax - '+ $
           'sdss_genprof(wave, dblt_name, config_strct, cleanspec,'+$
           '            [seed=,oseed=,conti=,/debug]'
     return,-1
  endif                         ; param prompt

  ;; Default Values
  if keyword_set(seed) then oseed = seed
  if size(dblt_name,/type) eq 8 then $
     dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)

  ;; Generate profiles
  if not keyword_set(vpstrct) then $
     vpstrct = sdss_genprof_mc(dblt, config_strct, seed=oseed, oseed=oseed,$
                               debug=debug) $
  else if (size(vpstrct,/dim))[0] lt 2*config_strct.nabs then $
     stop,'sdss_genprof(): too few profiles in input vpstrct'

  ;; Generate Voigt profile and convolve with un-normalize spectrum
  vpfx = sdss_genprof_voigt(wave, vpstrct)

  ;; Save
  spectrum = cleanspec
  spectrum[*,0] = cleanspec[*,0]*vpfx

  ;; There's a floor to the error (some constant error) in SDSS
  ;; and will estimate that by taking the median of the lowest 5%
  ;; error > 0.
  ;; newvar = var0 - (var0 - varflr)*(1 - vpfx) 
  ;; Variance can change a lot over spectrum.
;  npix = round((size(cleanspec[*,0],/dim))[0] * 0.25)
;  var0 = cleanspec[*,2]^2
;  gdpix = where(cleanspec[*,2] gt 0)
;  srt = gdpix[sort(var0[gdpix])]
;  varflr = median(var0[srt[0:npix]],/even) 
;  spectrum[*,2] = sqrt(var0 - (var0 - varflr)*(1 - sqrt(vpfx)))


  if keyword_set(debug) then begin
     ;; Print
     print,'sdss_genprof() debug: synthetic line summary'
     print,'#','ID','Wobs','zabs','N','Ncomp','EWI','EWII',$
           format='(2(a3,1x),a7,1x,a7,1x,a5,1x,a5,2(1x,a5))'
     rng = where(abs(vpstrct.wrest-dblt.wvI) lt 1.e-4 and $
                 vpstrct.id_comp eq 0,nvpstrct)
     wobs = vpstrct.wrest*(1.+vpstrct.zabs)
     printcol,rng,vpstrct[rng].id_sys,wobs[rng],vpstrct[rng].zabs,$
              vpstrct[rng].n,vpstrct[rng].ncomp,$
              vpstrct[rng].ew_obs,vpstrct[rng+1].ew_obs,$
              format='(2(i3,1x),f7.2,1x,f7.5,1x,f5.2,1x,i3,2x,2(1x,f5.2))'
     print,''

     ;; Figure out good number to replicate as indication of centroids
     gd = where(vpfx lt 1.)
     fxabs = replicate(median(spectrum[gd,0]),nvpstrct)
     wobs = wobs[2*lindgen(nvpstrct)]
     ion = dblt.ion + ' '+strtrim(round(dblt.wvI),2)

     ;; Plot
     clr = getcolor(/load)
     subttl = string((1.+config_strct.zlim)*dblt.wvI,$
                     format='(2x,2(f7.2,1x))')
     wvbnds = [vpstrct[rng].wvlim_sys[0],vpstrct[rng].wvlim_sys[1]]
     fxbnds = [fxabs,fxabs]
     if keyword_set(conti) then $
        x_splot,wave,spectrum[*,0],psym1=10,$
                ytwo=spectrum[*,2],psym2=10,$
                ythr=vpfx*conti,psym3=10,color3=clr.blue,$
                xfou=wobs,yfou=fxabs,color4=clr.limegreen,psym4=1,$ ;cross
                xfiv=wvbnds,$
                yfiv=fxbnds,psym5=1,color5=clr.purple,$
                ymnx=[min(spectrum[*,2]),max(spectrum[*,0])],$
                title='Synthetic Spectrum and Voigt Profile'+subttl,/block,$
                lgnd=['Flux','Error','Voigt',ion] $
     else x_splot,wave,cleanspec[*,0],psym1=10,color1=clr.charcoal,$
                  ytwo=cleanspec[*,2],psym2=10,color2=clr.magenta,$
                  ythr=spectrum[*,0],psym3=10,color3=clr.black,$
                  yfou=spectrum[*,2],psym4=10,color4=clr.red,$
                  xfiv=wobs,yfiv=fxabs,color5=clr.limegreen,psym5=1,$ ;cross
                  xsix=wvbnds,ysix=fxbnds,psym6=1,color6=clr.purple,$
                  ymnx=[min(spectrum[*,2]),max(spectrum[*,0])],$
                  title='Cleaned and Synthetic Spectrum'+subttl,/block,$
                  lgnd=['Flux0','Error0','Flux','Error',ion]
;     stop
  endif

  if keyword_set(oseed) then oseed = oseed[0] ; not the array
  return,spectrum
end
