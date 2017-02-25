;+ 
; NAME:
; fig_1darcfit
;     Version 1.1
;
; PURPOSE:
;-
;------------------------------------------------------------------------------

pro fig_coadd 

  if not keyword_set(SIDE) then side = 1
  if not keyword_set(PSFILE) then psfile='fig_coadd.ps'
  if not keyword_set(CSIZE) then csize = 1.5
  if not keyword_set(LSZ) then lsz = 1.2
  if not keyword_set(KEYINDX) then keyindx = 0L
  if not keyword_set(MEDINDX) then medindx = 100
  if not keyword_set(ENDTRM) then endtrm = 0
  if not keyword_set(ORDNM) then ordnm = 1
  if not keyword_set(SIGREJ) then sigrej=4.
  if not keyword_set(MINPIX1) then minpix1 = 100
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( CORDR ) then begin
      if side EQ 1 then CORDR = 85L else CORDR=55L
  endif
  if not keyword_set( SORDR ) then begin
      if side EQ 1 then SORDR = 70L else SORDR=30L
  endif
  if not keyword_set( SNRMIN ) then snrmin = 2.


  ;;
  readcol, 'Input/coadd.lst', files, FORMAT='A'
  nexp = n_elements(files)
  for q=0L,nexp-1 do begin
      tmp = xmrdfits(getenv('MIKE_PAP')+files[q], 1, STRUCTYP='mikeobjstrct', /silent)
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ 'a', ngd)
      ;; Add to total structure
      if q EQ 0 then mikeobj = tmp[gdobj] else mikeobj = [mikeobj, tmp[gdobj]]
  endfor

; Setup orders
  all_ordr = mikeobj.order
  all_ordr = all_ordr[uniq(all_ordr,sort(all_ordr))]

  ordrs = [min(all_ordr), max(all_ordr)]
  nspec = nexp

  echfspec = { mikefspecstrct }

  ;; Copy
  echfspec.nexp = nspec
  ;; Center ordr
  all_cen = where(mikeobj.order EQ CORDR, ncen)
  if ncen NE nspec then stop

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = mikeobj[all_cen[i]].exp

  ;; Other tags
  copy_struct, mikeobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  for i=0L,nspec-1 do echfspec.obj_fil[i] = files[i]

;  x_echcombspec, mikeobj, echfspec, ordrs, keyindx, SIGREJ=sigrej, $
;    ORDNM=ordnm, MCHK=mchk, NOFLUX=noflux, _EXTRA=extra

  allobj = mikeobj
  finspec = echfspec

;  Optional Keywords
  if size(ORDNM,/type) EQ 0 then ordnm = 1
  if not keyword_set(SIGREJ) then sigrej=4.
  if not keyword_set(MINPIX1) then minpix1 = 100
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( SNRMIN ) then snrmin = 2.
  if not keyword_set( OSTEP ) then ostep = 1

;  ordrs = [min(all_ordr), max(all_ordr)]

  nspec = finspec.nexp
  nexp = nspec

  svmedian = 1.
  svmedian = replicate(svmedian,nspec,ordrs[1]+1)
  svsnr = 1.
  svsnr = replicate(svsnr,nspec,ordrs[1]+1)
  ordr_typ = -1
  ordr_typ = replicate(ordr_typ,nspec,ordrs[1]+1)
  tmpfit = x_setfitstrct(NITER=3L, NORD=ordnm, $
                    HSIG=2.5, LSIG=2.5, /FLGREJ)
  fitstr = replicate(tmpfit,nspec,ordrs[1]+1)

  ;; Deal with FLUX vs not FLUX
  for qq=ordrs[0],ordrs[1],ostep do begin
      all_ordr = where(allobj.order EQ qq, nobj)
      for kk=0L,nobj-1 do begin
;          if nobj NE nspec then stop
          if keyword_set( NOFLUX ) then begin
              ;; Fill up flux array
;              print, 'hires_combspec: WARNING -- Filling flux ' + $
;                     'array with electrons'
              npix = allobj[all_ordr[kk]].npix
              if npix GT 0 then begin
                  allobj[all_ordr[kk]].flux[0:npix-1] = $
                    allobj[all_ordr[kk]].fx[0:npix-1] 
                  allobj[all_ordr[kk]].sig[0:npix-1] = $
                    sqrt(allobj[all_ordr[kk]].var[0:npix-1])
              endif
          endif
      endfor
  endfor

  ;; SHIFT and PAD
  if nexp GT 1 then begin
      for qq=ordrs[0],ordrs[1],ostep do begin
          all_ordr = where(allobj.order EQ qq, nobj)
          if nobj EQ 1 then continue
          ;; Allow for fewer exposures (i.e. non-overlapping order
          ;; limits)
          tmpk = keyindx < (nobj-1)
          fwv = allobj[all_ordr[tmpk]].wave[0]
          ;; Keyindx
          pad = (n_elements(allobj[all_ordr[tmpk]].flux)- $
                 allobj[all_ordr[tmpk]].npix)/2
          allobj[all_ordr[tmpk]].flux = $
            shift(allobj[all_ordr[tmpk]].flux, pad)
          allobj[all_ordr[tmpk]].wave = $
            shift(allobj[all_ordr[tmpk]].wave, pad)
          allobj[all_ordr[tmpk]].sig = $
            shift(allobj[all_ordr[tmpk]].sig, pad)
          ;; Other exposures
          for kk=0L,nobj-1 do begin
              if kk EQ tmpk then continue
              nbuff = n_elements(allobj[all_ordr[tmpk]].flux)- $
                allobj[all_ordr[tmpk]].npix
              ;; Min
              mn = min(abs(allobj[all_ordr[kk]].wave - fwv),imn)
              if imn NE 0 then begin
                  if (pad+imn) GT nbuff then stop
                  shft = pad-imn
              endif else begin
                  mn = min(abs(allobj[all_ordr[tmpk]].wave - $
                               allobj[all_ordr[kk]].wave[0]),shft)
              endelse
              
              ;; Shift
              allobj[all_ordr[kk]].flux = $
                shift(allobj[all_ordr[kk]].flux, shft)
              allobj[all_ordr[kk]].wave = $
                shift(allobj[all_ordr[kk]].wave, shft)
              allobj[all_ordr[kk]].sig = $
                shift(allobj[all_ordr[kk]].sig, shft)
;              if qq EQ 79 then begin
;                  x_splot, allobj[all_ordr[tmpk]].wave, $
;                    allobj[all_ordr[tmpk]].flux, $
;                    xtwo=allobj[all_ordr[kk]].wave, $
;                    ytwo=allobj[all_ordr[kk]].flux, /bloc
;              endif
          endfor
      endfor
  endif
  
  ;; Loop on Orders
  for qq=ordrs[0],ordrs[1],ostep do begin
      all_ordr = where(allobj.order EQ qq, nobj)
      finspec.phys_ordr[qq] = qq
      if nexp EQ 1 then begin
          npix = allobj[all_ordr].npix
          finspec.wave[0:npix-1,qq] = allobj[all_ordr].wave[0:npix-1]
          finspec.fx[0:npix-1,qq] = allobj[all_ordr].flux[0:npix-1]
          finspec.var[0:npix-1,qq] = (allobj[all_ordr].sig[0:npix-1])^2.
          ordr_typ[qq] = -1 
          continue
      endif 
      ;; Loop on exposures
      for kk=0L,nobj-1 do begin
          if kk EQ keyindx then continue
        
;        gdsig = dindgen(5000L)*0.-1.
;        gdnum = findgen(5000L)*0.
        
          ;; Find median ratio of flux and SNR
;          tst = where(allobj[all_ordr[keyindx]].wave NE 0.)
;          tstlen = size(tst,/dimensions)
;          tstwv = allobj[all_ordr[keyindx]].wave
;          tstwv[0:endtrm] = 0.
;          tstwv[tstlen[0]-endtrm-1:tstlen[0]-1] = 0.
          highsnr = where(allobj[all_ordr[keyindx]].sig GT 0. AND $
                          allobj[all_ordr[keyindx]].flux / $
                          allobj[all_ordr[keyindx]].sig GT SNRMIN  AND $
                          allobj[all_ordr[kk]].sig GT 0. AND $
                          allobj[all_ordr[kk]].flux GT 0. AND $
                          allobj[all_ordr[kk]].flux / $
                          allobj[all_ordr[kk]].sig GT SNRMIN AND $
                          allobj[all_ordr[keyindx]].wave GT 0., nhigh)
          if nhigh GT 100 then $
            med0 = median(allobj[all_ordr[keyindx]].flux[highsnr]$
                          /allobj[all_ordr[keyindx]].sig[highsnr]) $
          else med0 = 1.

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; CALCULATE RATIO AND SNR
          if not keyword_set(SILENT) then begin
              print, 'hires_combspec: order = ', strtrim(qq,2), $
                ' nhigh = ', strtrim(nhigh,2), ' med0 = ', med0
          endif
          if nhigh LT minpix1 then begin
              ;; Fill with exp. time ratio, may replace after other orders fit
              med_rto = allobj[all_ordr[keyindx]].exp/ $
                allobj[all_ordr[kk]].exp
              snr_rto = sqrt(allobj[all_ordr[keyindx]].exp/ $
                             allobj[all_ordr[kk]].exp)
              svmedian[kk,qq] = med_rto
              svsnr[kk,qq] = snr_rto
              
              ordr_typ[kk,qq] = 0 
          endif else begin
              svfx = median(allobj[all_ordr[keyindx]].flux[highsnr],medindx)
              svsig = median(allobj[all_ordr[keyindx]].sig[highsnr],medindx)
              if nhigh LT 300 then begin
                  ;; Take one value
                  ;; Ratio
                  rtio = (svfx/$
                          median(allobj[all_ordr[kk]].flux[highsnr],medindx))
                  djs_iterstat, rtio, maxiter=7L, median=medv, sigrej=3.
                  med_rto = medv
                  ;; SNR
                  snr = allobj[all_ordr[kk]].flux[highsnr]/ $
                    allobj[all_ordr[kk]].sig[highsnr]/  (svfx/svsig)
                  djs_iterstat, snr, maxiter=7L, median=medv, sigrej=3.
                  snr_rto = medv
                  ;; Save
                  svmedian[kk,qq] = med_rto
                  svsnr[kk,qq] = snr_rto
                  ordr_typ[kk,qq] = 1 
              endif else begin
                  ;; Fit a 2nd order POLY to the ratio 
                  ;; Ratio
                  rtio = (svfx /  $
                          median(allobj[all_ordr[kk]].flux[highsnr],medindx))
                  fitwv = allobj[all_ordr[kk]].wave[highsnr]
                  fit = x_fitrej(allobj[all_ordr[kk]].wave[highsnr], $
                                 rtio, FITSTR=fitstr[kk,qq])
                  med_rto = median(fit)
                  ;;if keyword_set(CHK) then begin
                  ;;  x_splot, rtio, ytwo=fit, /block, psym2=10, TITLE='rtio,fit'
                  ;;  stop
                  ;;endif
                  ;; SNR
                  snr = allobj[all_ordr[kk]].flux[highsnr]/ $
                    allobj[all_ordr[kk]].sig[highsnr]/ $
                    svfx * $
                    svsig
                  djs_iterstat, snr, maxiter=7L, median=medv, sigrej=3.
                  snr_rto = medv
                  ;; Save
                  svmedian[kk,qq] = med_rto
                  svsnr[kk,qq] = snr_rto
                  ordr_typ[kk,qq] = 2 

                  ;; 
                  if qq EQ 91 then begin
                      x_psopen, psfile, /maxs
                      clr = getcolor(/load)
                      wv = allobj[all_ordr[kk]].wave[highsnr]
                      yrng = [0.6, 1.0]
                      plot, [wv], [rtio], psym=1, $
                            background=clr.white, color=clr.black, $
                            xtitle='Wavelength (Ang)', ytitle='Median Smoothed Ratio', $
                            xmargin=[8,1], $
                            ymargin=[4,0], yrange=yrng, xrange=[3880., 3950.], $
                            xstyle=1, ystyle=1, charsize=csize
                      oplot, wv, fit, color=clr.red, linesty=2
                      x_psclose
                      !p.multi=[0,1,1]
                      ;stop
                      ;x_splot, allobj[all_ordr[0]].flux, ytwo=allobj[all_ordr[1]].flux*0.75, /bloc
                  endif
              endelse
          endelse
      endfor
  endfor


  return
end

