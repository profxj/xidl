;+ 
; NAME:
; sdss_gz.pro
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;    2-Nov-2011  created by KLC (removed from sdss_functions)
;-
;------------------------------------------------------------------------------

pro sdss_gz_readgzfil, gzcum_fil, gzglobal, zarr, ewarr, $
                       zbinsize=zbinsize, nzbin=nzbin, $
                       ewbinsize=ewbinsize, newbin=newbin, $
                       header=header, zlim=zlim, ewlim=ewlim
  if n_params() ne 4 then begin
     print,'Syntax - sdss_gz_readgzfil, gzcum_fil, gzglobal, zarr, ewarr, '
     print,'                            [zbinsize=, nzbin=, ewbinsize=, newbin=,'
     print,'                            header=, zlim=, ewlim=]'
     return
  endif 

  gzglobal = xmrdfits(gzcum_fil,0,header,/silent) 
  zarr = sdss_mkzarr(zlim,zbinsize,header=header,nzbin=nzbin)
  ewarr = sdss_mkewarr(ewlim,ewbinsize,header=header,$
                          newbin=newbin)
end                             ; sdss_gz_readgzfil


function sdss_gz_rebinz, gzcum_fil, zbinsize, zarr=zarr, $
                        nzbin=nzbin, zlim=zlim, $
                        header=header
  ;; Collapse gzglobal 2D grid to larger bins
  if n_params() ne 2 then begin
     print,'Syntax - '
     return, -1
  endif
  
  sdss_gz_readgzfil, gzcum_fil, gzglobal0, zarr0, ewarr0, $
                     zbinsize=zbinsize0, zlim=zlim0, nzbin=nzbin0, $
                     newbin=newbin0, header=header
  
  ;; Find nearest (lower) redshift bin
  if keyword_set(zlim) then begin
     izmin = floor((zlim[0] - zlim0[0]) / zbinsize0) 
  endif else begin
     izmin = 0
  endelse 
  sxaddpar,header,'ZMIN',zarr0[izmin]

  ;; Figure out number of z bins to sum over
  sumbin = floor(zbinsize / zbinsize0)
  zbinsize = sumbin * zbinsize0
  sxaddpar,header,'ZBINSZ',zbinsize

  ;; Figure out how many can fit this exactly
  nzbin = floor( (nzbin0 - izmin) / sumbin ) ; max out
  sxaddpar,header,'NAXIS1',nzbin
  sxaddpar,header,'ZMAX',(nzbin-1)*zbinsize

  ;; Instantiate (like in sdss_gz)
  gzglobal = fltarr(nzbin,newbin0,2) ; dz, dX
  for zz=0L,nzbin-1 do begin
     istart = izmin + sumbin*zz
     iend = istart + sumbin - 1
     gzglobal[zz,*,0] = total(gzglobal0[istart:iend,*,0],1)
     gzglobal[zz,*,1] = total(gzglobal0[istart:iend,*,1],1)
;     print,zz,istart,iend,gzglobal[zz,20,1],total(gzglobal0[istart:iend,20,1]),median(gzglobal0[istart:iend,20,1])
  endfor                        ; loop zz=nzbin
  zarr = sdss_mkzarr(zlim,zbinsize,header=header)

;  stop
  return, gzglobal
end                             ; sdss_gz_rebinz()


function sdss_gz_rebinew, gzcum_fil, ewbinsize, ewarr=ewarr, $
                        newbin=newbin, ewlim=ewlim, $
                        header=header
  ;; Collapse gzglobal 2D grid to larger bins
  if n_params() ne 2 then begin
     print,'Syntax - '
     return, -1
  endif
  
  sdss_gz_readgzfil, gzcum_fil, gzglobal0, zarr0, ewarr0, $
                     ewbinsize=ewbinsize0, ewlim=ewlim0, nzbin=nzbin0, $
                     newbin=newbin0, header=header
  
  ;; Find nearest (lower) redshift bin
  if keyword_set(ewlim) then begin
     iewmin = floor((ewlim[0] - ewlim0[0]) / ewbinsize0) 
  endif else begin
     iewmin = 0
  endelse 
  sxaddpar,header,'EWMIN',ewarr0[iewmin]

  ;; Figure out number of EW bins to sum over
  sumbin = floor(ewbinsize / ewbinsize0)
  ewbinsize = sumbin * ewbinsize0
  sxaddpar,header,'EWBINSZ',ewbinsize

  ;; Figure out how many can fit this exactly
  newbin = floor( (newbin0 - iewmin) / sumbin ) ; max out
  sxaddpar,header,'NAXIS2',newbin
  sxaddpar,header,'EWMAX',(newbin-1)*ewbinsize

  ;; Instantiate (like in sdss_gz)
  gzglobal = fltarr(nzbin0,newbin,2) ; dz, dX
  for ee=0L,newbin-1 do begin
     istart = iewmin + sumbin*ee
     iend = istart + sumbin - 1
     gzglobal[*,ee,0] = total(gzglobal0[*,istart:iend,0],2)
     gzglobal[*,ee,1] = total(gzglobal0[*,istart:iend,1],2)
     
;     print,ee,istart,iend,gzglobal[20,ee,1],total(gzglobal0[20,istart:iend,1]),median(gzglobal0[20,istart:iend,1])
  endfor                        ; loop ee=newbin
  ewarr = sdss_mkewarr(ewlim,ewbinsize,header=header)

;  stop
  return, gzglobal
end                             ; sdss_gz_rebinew()


function sdss_gz_dx, gzcum_fil, zval, ewval, dz=dz, $
                     zarr=zarr, ewarr=ewarr
  ;; Interpolate 2D grid to return dX (or dz) for given zval and ewval
  ;; (which may either both be arrays or one array and one single value)
  if n_params() ne 3 then begin
     print,'Syntax - sdss_gz_dx(gzcum_fil, zval, ewval, [/dz, zarr=, ewarr=])'
     return, -1
  endif 

  if keyword_set(dz) then idelta = 0 else idelta = 1
  if size(gzcum_fil,/type) eq 7 then begin
     gzglobal = xmrdfits(gzcum_fil, 0, header, /silent)
     zarr = sdss_mkzarr(zlim,zbinsize,header=header, nzbin=nzbin)
     ewarr = sdss_mkewarr(ewlim,ewbinsize,header=header, newbin=newbin)
  endif else begin
     if not keyword_set(zarr) or not keyword_set(ewarr) then begin
        print,'sdss_gz_dx(): must set zarr= and ewarr='
        return, -1              ; EXIT
     endif 
     gzglobal = gzcum_fil       ; assume vector
     nzbin = (size(zarr,/dim))[0]
     zlim = [zarr[0],zarr[nzbin-1]] ; assume sorted
     zbinsize = zarr[1] - zarr[0]   ; assume uniform [THIS MAY BE WRONG LATER]

     newbin = (size(ewarr,/dim))[0]
     ewlim = [ewarr[0],ewarr[newbin-1]] 
     ewbinsize = ewarr[1] - ewarr[0]
  endelse

  nzval = (size(zval,/dim))[0] > 1 ; handle singularity
  newval = (size(ewval,/dim))[0] > 1 

  if nzval ne newval then begin
     if nzval eq 1 or newval eq 1 then grid = 1 $
     else begin
        print,'sdss_gz_dx(): cannot handle mismatched zval, ewval sizes ',$
              nzval, newval
        return, -1
     endelse
  endif 

  ;; Interpolate in 2D
  ;; First must identify the fractional "indices" to use with
  ;; IDL's interpolate()
  ;; Linear interpolation for mid-point (xm, ym) between points
  ;; (x1,y2) and (x2,y2) 
  ;; ym =  (y2-y1)/(x2-x1) * (xm - x1) + y1

  iz1 = floor((zval - zlim[0])/zbinsize)
  iz2 = iz1 + 1
  izm = (iz2-iz1)/(zarr[iz2]-zarr[iz1])  * (zval - zarr[iz1]) + iz1
  
  test = where(izm lt 0,ntest)
  if ntest ne 0 then begin
     print,'sdss_gz_dx(): number of points below z bounds:',ntest
     izm[test] = 0              ; floor
  endif 
  test = where(izm ge nzbin,ntest)
  if ntest ne 0 then begin
     print,'sdss_gz_dx(): number of points above z bounds:',ntest
     izm[test] = nzbin-1              ; ceil     
  endif 

  iw1 = floor((ewval - ewlim[0])/ewbinsize)
  iw2 = iw1 + 1
  iwm = (iw2-iw1)/(ewarr[iw2]-ewarr[iw1])  * (ewval - ewarr[iw1]) + iw1

  test = where(iwm lt 0,ntest)
  if ntest ne 0 then begin
     print,'sdss_gz_dx(): number of points below EW bounds:',ntest
     iwm[test] = 0              ; floor
  endif 
  test = where(iwm ge newbin,ntest)
  if ntest ne 0 then begin
     print,'sdss_gz_dx(): number of points above EW bounds:',ntest
     iwm[test] = newbin-1              ; ceil     
  endif 

  
  ;; use indices to interpolate in gzglobal slice and set out-of-bound
  ;; values to -1
  dx = interpolate(gzglobal[*,*,idelta], izm, iwm, missing=-1, grid=grid) 

;  stop

  return, dx

end                             ; sdss_gz_dX()


function sdss_gz_sanitycheck, dblt_name=dblt_name, dz=dz, dztot=dztot, $
                              zlim=zlim, snrmin=snrmin, _extra=extra
  ;; Brute force way of making sure sdss_gz gets sensible results

  cosmology = sdss_setcosmology(_extra=extra) ; _extra = cosmology[3]

  sdssdir = sdss_getsdssdir()
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  if not keyword_set(snrmin) then snrmin = 4.
  snrstrct = xmrdfits(sdssdir+'inputs/dr7qso_noBAL_SNR.fit',1,/silent)
  tags = tag_names(snrstrct)
  isnr = (where(tags eq 'SNR_'+strupcase(dblt.ion)))[0]
  iwvobs = (where(tags eq 'WVOBS_'+strupcase(dblt.ion)))[0]
  if keyword_set(zlim) then begin
     if size(zlim,/dim) eq 2 then zrng = zlim $
     else zrng = [zlim,!values.d_infinity]
     gd = where(snrstrct.(isnr)[2] ge snrmin and snrstrct.z_qso ge zrng[0] $
                and snrstrct.z_qso lt zrng[1],ngd)
  endif else gd = where(snrstrct.(isnr)[2] ge snrmin,ngd)

  zmin = snrstrct[gd].(iwvobs)[0] / dblt.wvI - 1.
  zmax = snrstrct[gd].(iwvobs)[1] / dblt.wvI - 1.

  dxarr = cosm_xz(zmax,zmin=zmin,/silent,/exact,/noinit)
  dzarr = zmax - zmin

  dztot = total(dzarr)
  dxtot = total(dxarr)
  if keyword_set(dz) then begin
     tmp = dxtot
     dxtot = dztot
     dztot = tmp
  endif

  return,dxtot
     
end                             ; sdss_gz_sanitycheck()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_gz, list_fil, gzcum_fil, nsig=nsig, zbinsize=zbinsize, zlim=zlim,$
             ewbinsize=ewbinsize, ewlim=ewlim, $
             dblt_name=dblt_name, clobber=clobber, debug=debug
  ;; Rudimentary completeness test results
  if n_params() ne 2 then begin
     print,'Syntax - sdss_gz, list_fil, gzcum_fil, [nsig=, zbinsize=, zlim=, '
     print,'                  ewbinsize=, ewlim=,  dblt_name=, /clobber, /debug]'
     return
  endif

  cosmology = sdss_setcosmology(_extra=extra) ; _extra = cosmology[3]
  
  sdssdir = sdss_getsdssdir()
;  pixscale = sdss_getspecpixscale(/loglam) 

  ;; Defaults to make something sensible
  if not keyword_set(nsig) then nsig = 3. ; ~99% c.l.
  if not keyword_set(zbinsize) then zbinsize = 0.005
  if not keyword_set(zlim) then zlim = [1., 6.]
  if not keyword_set(ewbinsize) then ewbinsize = 0.05 ; Ang
  if not keyword_set(ewlim) then ewlim = [0.05, 5.]   ; Ang
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  dlambdar = dblt.wvII - dblt.wvI

  ;; For reference, on the c.l. of given nsig
  ;; gsigma = 1. + findgen(10)*0.5
  ;; gprob = [68.26895123013139,86.63856338891091,$
  ;;          95.4499772097615,98.75806921717441,$
  ;;          99.7300213908383,99.9534744959942,$
  ;;          99.99366582300915,99.9993204775126,$
  ;;          99.99994267123438,99.99999620223669]

  ;; Bins are defined on the left-hand side by convention
  zglobal = sdss_mkzarr(zlim, zbinsize, nzbin=nzbin)
  ewglobal = sdss_mkewarr(ewlim, ewbinsize, newbin=newbin)

  ;; Create place to store EVERYTHING
  gzglobal = fltarr(nzbin,newbin,2) ; dz and dX

  ;; Read in data and get names
  readcol,list_fil,spec_fil,format='a',/silent,skip=1
  nspec = (size(spec_fil,/dim))[0]

  ;;  gz array will be [npix, [zabs, dz, EWlim, mask, nsig]]
  gz_fil = sdss_getname(spec_fil,/spec,gz=dblt.ion,dir=gzdir)

  ;; Loop over spectra and sum up dz and dX
  for ss=0L,nspec-1 do begin 
     ;; Read it in
     gz_los = xmrdfits(sdssdir+gzdir[ss]+gz_fil[ss],0,/silent)
     npix = (size(gz_los[*,0],/dim))[0]

     ;; Scale the EW limit to desired Nsigma
     gz_los[*,2] = gz_los[*,2] * nsig / gz_los[*,4] 

     ;; Make partner pixel limit array
     ;; gz_los[*,1] = dz = dlambda / wvI
     indexII = lindgen(npix) + $
               round((1+gz_los[*,0])*dlambdar/(gz_los[*,1]*dblt.wvI))
     gd = where(indexII lt npix,ngd,complement=bd,ncomplement=nbd)
     ewlim_wvII = fltarr(npix,/nozero) 
     if ngd ne 0 then $
        ewlim_wvII[0:ngd-1] = gz_los[indexII[gd],2] ; shift value
     if nbd ne 0 then $
        ewlim_wvII[ngd:*] = !values.f_infinity ; default

     ;; Calculate the dX exactly for all pixels
     dxarr = cosm_xz(gz_los[*,0]+gz_los[*,1],zmin=gz_los[*,0],$
                     /silent,/exact,/noinit) ; calc
     
     ;; First and last redshift bins may be partials
     for zz=0L,nzbin-1 do begin
        ;; To minimize looping, just check redshift and goodness
        ;; (mask) 
        sub = where(gz_los[*,0] ge zglobal[zz] and $
                    gz_los[*,0] lt zglobal[zz]+zbinsize and $
                    gz_los[*,3] eq 1) ; dblt.wvI in range
        if sub[0] eq -1 then continue ; SKIP
        
        ;; Go through the EW bins
        for ww=0,newbin-1 do begin
           ;; Check that physical (> 0.) and less than the upper bin
           ;; boundary 
           ;; check both doublet lines conditions meet bin
           gd = where(gz_los[sub,2] lt ewglobal[ww]+ewbinsize and $
                      ewlim_wvII[sub] lt ewglobal[ww]+ewbinsize, $
                      npix_sub) ; mask
           
           if npix_sub eq 0 then $
              continue          ; SKIP
           sub = sub[gd]
           
           ;; Now can just add it up cumulatively 
           ;; What about fractional contributions?!!!
           gzglobal[zz,ww,0] = gzglobal[zz,ww,0] + $
                               total(gz_los[sub,1])
           gzglobal[zz,ww,1] = gzglobal[zz,ww,1] + $
                               total(dxarr[sub])
        endfor                  ; loop ww=ewbin
     endfor                     ; loop zz=zend
     
     if (ss mod 500L) eq 0L then begin
        maxdx = max(gzglobal[*,*,1])
        maxdz = max(gzglobal[*,*,0])
        print,'sdss_gz: ss = '+strtrim(ss,2),maxdz,maxdx
        if keyword_set(debug) then stop
     endif 
     
  endfor                        ; loop ss=nspec
  
  
  ;; Ready to write out
  fxhmake, header, gzglobal

  ;; Add keywords
  sxaddpar, header, 'NSPEC', nspec, 'Number of spectra'
  sxaddpar, header, 'NSIG', nsig, 'Assumed nsigma in EWlim est'
  sxaddpar, header, 'ZBINSZ', zbinsize, 'Redshift bin size'
  sxaddpar, header, 'ZMIN', zlim[0], 'Minimum redshift'
  sxaddpar, header, 'ZMAX', zlim[1], 'Maximum redshift'
  sxaddpar, header, 'EWBINSZ', ewbinsize, 'EWlim bin size (Ang)'
  sxaddpar, header, 'EWMIN', ewlim[0], 'Minimum EWlim (Ang)'
  sxaddpar, header, 'EWMAX', ewlim[1], 'Maximum EWlim (Ang)'
  sxaddpar, header, 'H0', cosmology[0], 'Adopted Hubble constant'
  sxaddpar, header, 'OmDM', cosmology[1], 'Cosmology: Matter density'
  sxaddpar, header, 'Lambda', cosmology[2], 'Cosmology: Vacuum energy density'

  test = file_search(sdssdir+gzcum_fil+'*',count=ntest)
  if ntest eq 0 or keyword_set(clobber) then begin
     mwrfits,gzglobal,sdssdir+gzcum_fil,header,/create,/silent
     spawn,'gzip -f '+sdssdir+gzcum_fil
     print,'sdss_gz: created ',gzcum_fil
  endif 

end                             ; sdss_gz


