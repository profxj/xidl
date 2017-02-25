;+ 
; NAME:
; x_echcombspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run even on an object with a single exposure
;
; CALLING SEQUENCE:
;  x_echcombspec, allobj, finspec, ordrs, keyindx
;
; INPUTS:
;   allobj   -- Array of echspec structures
;   ordrs    -- Orders to coadd 
;   keyindx  -- Exposure to serve as the fiducial for stats
;
; RETURNS:
;
; OUTPUTS:
;   finspec  --  echfspec structure containing the combined spectrum
;
; OPTIONAL KEYWORDS:
;    /SILENT   - No text output
;    /NOFLUX   - Use the fluxed array
;    MINPIX1=  - Minimum 'good' pixels to calculate fitting profile
;    ORDNM=    - Fitting profile order [default: 1]
;    SNRMIN=   - Minimum S/N per pixel for stats [Default: 2]
;    REJSIG=   - Parameter passsed to x_combine [default = 4.]
;    MEDINDX=  - Value for median smoothing [default: 100]
;    /NOREJ    - Do not search for and reject bad pixels (e.g. cosmic
;    /ACHK     - Force check of all orders
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hires_combspec, hires, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;   10-Jun-2004 Revisited by GEP
;   Sep-2005 Revisited by JXP
;   17-Jun-2008 Modified for MAGE by CLW
;-
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro mage_echcombspec, allobj, finspec, ordrs, keyindx, SILENT=silent, $
                    SNRMIN=snrmin, MEDINDX=medindx, $
                    NOFLUX=noflux, FCHK=fchk, $
                    MINPIX1=minpix1, ORDNM=ordnm, $
                    CHK=chk, SIGREJ=sigrej, MCHK=mchk, ACHK=achk, $
                    OSTEP=ostep, REJSIG=rejsig

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'x_echcombspec, allobj, finspec, ordrs, keyindx, /SILENT, ' + $
      'OBJ_NM=, /PATCH, /CHK, /MCHK  [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(MEDINDX) then medindx = 100
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
                  allobj[all_ordr[kk]].nosig[0:npix-1] = $
                    sqrt(allobj[all_ordr[kk]].novar[0:npix-1])
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
          cc = getcolor(/load)

          tmpk = keyindx < (nobj-1)
          fwv = min(allobj[all_ordr[tmpk]].wave)
          ;; Keyindx
          pad = (n_elements(allobj[all_ordr[tmpk]].flux)- $
                 allobj[all_ordr[tmpk]].npix)/2
          allobj[all_ordr[tmpk]].flux = $
            shift(allobj[all_ordr[tmpk]].flux, pad)
          allobj[all_ordr[tmpk]].wave = $
            shift(allobj[all_ordr[tmpk]].wave, pad)
          allobj[all_ordr[tmpk]].sig = $
            shift(allobj[all_ordr[tmpk]].sig, pad)
          allobj[all_ordr[tmpk]].nosig = $
            shift(allobj[all_ordr[tmpk]].nosig, pad)
;          plot, allobj[all_ordr[tmpk]].wave, allobj[all_ordr[tmpk]].flux

          ;; Other exposures
          for kk=0L,nobj-1 do begin
              if kk EQ tmpk then continue
              nbuff = n_elements(allobj[all_ordr[tmpk]].flux)- $
                allobj[all_ordr[tmpk]].npix
              ;; Min
              mn = min(abs(allobj[all_ordr[kk]].wave - fwv),imn)
              if imn NE 0 then begin
;                  if (pad+imn) GT nbuff then stop
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
              allobj[all_ordr[kk]].nosig = $
                shift(allobj[all_ordr[kk]].nosig, shft)
;              oplot, allobj[all_ordr[kk]].wave, allobj[all_ordr[kk]].flux, color=cc.red
;          stop

          endfor
      endfor
  endif
  
  ;; Loop on Orders
  for qq=ordrs[0],ordrs[1],ostep do begin
      all_ordr = where(allobj.order EQ qq, nobj)
      finspec.phys_ordr[qq] = qq
      if nexp EQ 1 then begin
          npix = allobj[all_ordr].npix
          if npix GT 0 then begin
              finspec.wave[0:npix-1,qq] = allobj[all_ordr].wave[0:npix-1]
              finspec.fx[0:npix-1,qq] = allobj[all_ordr].flux[0:npix-1]
              finspec.var[0:npix-1,qq] = (allobj[all_ordr].sig[0:npix-1])^2.
              finspec.novar[0:npix-1,qq] = (allobj[all_ordr].nosig[0:npix-1])^2.
          endif
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
              print, 'mage_combspec: order = ', strtrim(qq,2), $
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
                    svfx * svsig
                  djs_iterstat, snr, maxiter=7L, median=medv, sigrej=3.
                  snr_rto = medv
                  ;; Save
                  svmedian[kk,qq] = med_rto
                  svsnr[kk,qq] = snr_rto
                  ordr_typ[kk,qq] = 2 
              endelse
          endelse
      endfor
  endfor

  medianave = replicate(1.,nexp)
  mediansig = replicate(0.,nexp)
  ;; Have probably put a bug in here relating to keyindx
  for kk=0,nexp-1 do begin
      ;; fix BAD orders
      for qq=ordrs[0],ordrs[1],ostep do begin
          all_ordr = where(allobj.order EQ qq, nobj)
          if nobj EQ 1 then continue
          if ordr_typ[kk,qq] EQ 0 then begin
              useordr = where(abs(lindgen(ordrs[1]+1)-qq) LE 5 $
                              AND ordr_typ[kk,*] GE 1, nuse)
              ord_diff = 10
              while nuse LE 5 AND ord_diff LE 35 do begin
                  useordr = where(abs(lindgen(ordrs[1]+1)-qq) LE ord_diff $
                                  AND ordr_typ[kk,*] GE 1, nuse)
                  ord_diff = ord_diff + 5
              endwhile
              if nuse LT 5 then begin
                  med_rto = allobj[all_ordr[keyindx]].exp/ $
                    allobj[all_ordr[kk]].exp
                  snr_rto = sqrt(allobj[all_ordr[keyindx]].exp/ $
                                 allobj[all_ordr[kk]].exp)
                  svmedian[kk,qq] = med_rto
                  svsnr[kk,qq] = snr_rto
              endif else begin
                  svsnr[kk,qq] = median(svsnr[kk,useordr])
                  svmedian[kk,qq] = median(svmedian[kk,useordr])
              endelse
          endif
      endfor 
      ;; Get Median, sig
      medianave[kk] = median(svmedian[kk,ordrs[0]:ordrs[1]])
      for qq=ordrs[0],ordrs[1],ostep do $
        mediansig[kk] = mediansig[kk] + (medianave[kk]-svmedian[kk,qq])^2.
      mediansig[kk] = sqrt(mediansig[kk]/(ordrs[1]-ordrs[0]+1))
      ;; Check for outliers
;      if keyword_set(MCHK) then begin
;          xcol = lindgen(ordrs[1]-ordrs[0]+1)+ordrs[0]
;          ptit = 'Median Ratios for Exposure ' + strtrim(string(kk),2) + $
;            ' median = ' + strtrim(string(medianave[kk]),2) + $
;            ' sigma = ' + strtrim(string(mediansig[kk]),2)
;          x_splot, xcol, svmedian[kk,ordrs[0]:ordrs[1]], /block, TITLE=ptit
;      endif
      ;; Replace
      for qq=ordrs[0],ordrs[1],ostep do begin
          all_ordr = where(allobj.order EQ qq, nobj)
          if nobj NE nexp then continue
          if abs(svmedian[kk,qq] - medianave[kk]) $
            GT 2.*mediansig[kk] then begin
              ordr_typ[kk,qq] = 0
              useordr = where(abs(lindgen(ordrs[1]+1)-qq) LE 5 $
                              AND ordr_typ[kk,*] GE 1, nuse)
              ord_diff = 10
              while nuse LE 5 AND ord_diff LE 35 do begin
                  useordr = where(abs(lindgen(ordrs[1]+1)-qq) LE ord_diff $
                                  AND ordr_typ[kk,*] GE 1, nuse)
                  ord_diff = ord_diff + 5
              endwhile
              if nuse LT 5 then begin
                  med_rto = allobj[all_ordr[keyindx]].exp/ $
                    allobj[all_ordr[kk]].exp
                  snr_rto = sqrt(allobj[all_ordr[keyindx]].exp/ $
                                 allobj[all_ordr[kk]].exp)
                  svmedian[kk,qq] = med_rto
                  svsnr[kk,qq] = snr_rto
              endif else begin
                  svsnr[kk,qq] = median(svsnr[kk,useordr])
                  svmedian[kk,qq] = median(svmedian[kk,useordr])
              endelse
          endif
      endfor
  endfor
  if keyword_set(MCHK) then begin
      xcol = lindgen(ordrs[1]-ordrs[0]+1)+ordrs[0]
      for qq=0,nexp-1 do begin
          if qq EQ keyindx then continue
          ptit = 'Final Median Ratios for Exposure ' + strtrim(string(qq),2) + $
            ' median = ' + strtrim(string(medianave[qq]),2) + $
            ' sigma = ' + strtrim(string(mediansig[qq]),2)
         x_splot, xcol, svmedian[qq,ordrs[0]:ordrs[1]], /block, TITLE=ptit
      endfor
  endif
    
  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Scale
  for kk=0,nexp-1 do begin
      for qq=ordrs[0],ordrs[1],ostep do begin
          all_ordr = where(allobj.order EQ qq, nobj)
;          if nobj NE nspec then continue
          finspec.phys_ordr[qq] = qq
          ;; Rescale by the median (one number)
          if ordr_typ[kk,qq] EQ 1 OR ordr_typ[kk,qq] EQ 0 then begin
              allobj[all_ordr[kk]].flux = $
                allobj[all_ordr[kk]].flux * svmedian[kk,qq] 
              allobj[all_ordr[kk]].sig = $
                allobj[all_ordr[kk]].sig * svmedian[kk,qq]
              allobj[all_ordr[kk]].nosig = $
                allobj[all_ordr[kk]].nosig * svmedian[kk, qq]

              if keyword_set(CHK) and kk EQ (nexp-1) then begin
                wuse = where(allobj[all_ordr[keyindx]].wave NE -1000.)
                ptit = 'Rescaled 0, 1, Order: ' + strtrim(string(qq),2) + $
                       ' and Exposure (red): ' + strtrim(string(kk),2)
                
                x_splot, allobj[all_ordr[keyindx]].wave[wuse],allobj[all_ordr[keyindx]].flux[wuse], TITLE=ptit,$
                   ytwo=allobj[all_ordr[kk]].flux[wuse], /block
              endif
          endif 
          if ordr_typ[kk,qq] EQ 2 then begin
              ;; Rescale by a linear correction
;              tstwv = allobj[all_ordr[keyindx]].wave
;              tstwv[0:endtrm] = 0.
;              tstwv[tstlen[0]-endtrm-1:tstlen[0]-1] = 0.
              highsnr = where(allobj[all_ordr[keyindx]].sig GT 0. AND $
                              allobj[all_ordr[keyindx]].flux / $
                              allobj[all_ordr[keyindx]].sig GT SNRMIN  AND $
                              allobj[all_ordr[kk]].sig GT 0. AND $
                              allobj[all_ordr[kk]].flux / $
                              allobj[all_ordr[kk]].sig GT SNRMIN AND $
                              allobj[all_ordr[keyindx]].wave GT 0., nhigh)
              svfx = median(allobj[all_ordr[keyindx]].flux[highsnr],medindx)
              svsig = median(allobj[all_ordr[keyindx]].sig[highsnr],medindx)
              rtio = (svfx /  $
                        median(allobj[all_ordr[kk]].flux[highsnr],medindx))
              gdwv = where(allobj[all_ordr[kk]].wave GT 0.)
              fitwv = allobj[all_ordr[kk]].wave[highsnr]
              tmpfit = x_setfitstrct(NITER=3L, NORD=ordnm, $
                    HSIG=2.5, LSIG=2.5, /FLGREJ)
              fit = x_fitrej(allobj[all_ordr[kk]].wave[highsnr], $
                             rtio, FITSTR=tmpfit)
              if (keyword_set(CHK) and kk EQ (nexp-1)) OR $
                keyword_set(ACHK) then begin
                ptit = 'RTIO,FIT for Order: ' + strtrim(string(qq),2) + $
                       ' and Exposure: ' + strtrim(string(kk),2)
                x_splot, highsnr, rtio, ytwo=fit, /block, psym2=10, TITLE=ptit
              ;  stop
              endif
              allfit = x_calcfit(allobj[all_ordr[kk]].wave[gdwv], $
                                 FITSTR=tmpfit)
                                 ;FITSTR=svfit[qq])
              allobj[all_ordr[kk]].flux[gdwv] = $
                allobj[all_ordr[kk]].flux[gdwv] * allfit
              allobj[all_ordr[kk]].sig[gdwv] = $
                allobj[all_ordr[kk]].sig[gdwv] * allfit
              allobj[all_ordr[kk]].nosig[gdwv] = $
                allobj[all_ordr[kk]].nosig[gdwv] * allfit

              if (keyword_set(CHK) and kk EQ (nexp-1)) OR $
                keyword_set(ACHK) then begin
                wuse = where(allobj[all_ordr[keyindx]].wave NE -200.)
;                w2use = where(allobj[all_ordr[kk]].wave NE 0.)
                ptit = 'Rescaled 2, Order: ' + strtrim(string(qq),2) + $
                       ' and Exposure (red): ' + strtrim(string(kk),2)
                plot, allobj[all_ordr[keyindx]].wave[wuse], allobj[all_ordr[keyindx]].flux[wuse]
                oplot, allobj[all_ordr[kk]].wave[wuse],allobj[all_ordr[kk]].flux[wuse], color=cc.red
                stop
             endif
          endif 
      endfor
    endfor

    ;;;;;;;;;;;;;;;;;;;;;;
    ;; Combine
    for qq=ordrs[0],ordrs[1],ostep do begin
      all_ordr = where(allobj.order EQ qq, nobj)
      if nobj EQ 1 then begin
          gd = where(allobj[all_ordr[0]].var GT 0., ngd)
          if ngd NE 0 then begin
              fpx = gd[0]
              epx = gd[ngd-1]
              npix = epx-fpx+1
              ;; Fill it up
              finspec.wave[0:npix-1,qq] = allobj[all_ordr[0]].wave[fpx:epx]
              finspec.fx[0:npix-1,qq] = allobj[all_ordr[0]].flux[fpx:epx]
              finspec.var[0:npix-1,qq] = allobj[all_ordr[0]].sig[fpx:epx]^2
              bdpix=where(allobj[all_ordr[0]].sig[fpx:epx] lt 0,nbd)
              if nbd gt 0 then finspec.var[bdpix,qq]=-1
              finspec.novar[0:npix-1,qq] = allobj[all_ordr[0]].nosig[fpx:epx]^2
          endif
          continue
      endif
          
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Set bad pix to -1
      sig = double(allobj[all_ordr].sig)
      nosig = double(allobj[all_ordr].nosig)
      a = where(allobj[all_ordr].sig LE 0., na)
      if na GT 0 then begin
          sig[a] = -1.
          nosig[a] = -1.
      endif
      if not keyword_set(SILENT) then print, 'hires_combspec: Order ', qq

      ;; FIX SIGMA FOR 2 EXPOSURE CASE
      if nobj EQ 2 then begin
         sigma = 0.
         gdvar = where(allobj[all_ordr[0]].sig GT 0. AND $
                       allobj[all_ordr[1]].sig GT 0., ngdvar)
         if ngdvar NE 0 then begin
             d = (allobj[all_ordr[0]].flux[gdvar]- $
                  allobj[all_ordr[1]].flux[gdvar])/ $
                 sqrt(allobj[all_ordr[0]].sig[gdvar]^2 + $
                      allobj[all_ordr[1]].sig[gdvar]^2)
	     djs_iterstat, d, sigma=sigma, sigrej=sigrej, maxiter=5
             sig = sig*sigma
             nosig = nosig*sigma
             d = (allobj[all_ordr[0]].flux[gdvar]- $
                  allobj[all_ordr[1]].flux[gdvar])/ $
                 sqrt(sigma^2*allobj[all_ordr[0]].sig[gdvar]^2 $
                      + sigma^2*allobj[all_ordr[1]].sig[gdvar]^2)
	     djs_iterstat, d, sigma=sigma2, sigrej=sigrej, maxiter=5
         endif 
      endif 

      ;; COMBINE
      sz = size(allobj[all_ordr].flux,/dimensions)
      mask = replicate(0,sz[0],sz[1])
      svmedian[*,qq] = 1.

;      plot, allobj[all_ordr[0]].wave, allobj[all_ordr[0]].flux
;      oplot, allobj[all_ordr[1]].wave, allobj[all_ordr[1]].flux, color=cc.red

;      stop
      x_combspec, allobj[all_ordr].flux, (sig>0.)^2, $
                  fflux, fvar, WAVE=allobj[all_ordr[0]].wave, $
                  SNR=svsnr[*,qq], SCALE=svmedian[*,qq], mask=mask, $
                  NSIG=rejsig, novar = (nosig > 0.)^2, fnovar = fnovar
      plot, allobj[all_ordr[0]].wave

      ;; Find first pix
      gd = where(fvar GT 0., ngd)
      if ngd NE 0 then begin
          fpx = gd[0]
          epx = gd[ngd-1]
          npix = epx-fpx+1
          ;; Fill it up
          finspec.wave[0:npix-1,qq] = allobj[all_ordr[0]].wave[fpx:epx]
          finspec.fx[0:npix-1,qq] = temporary(fflux[fpx:epx])
          finspec.var[0:npix-1,qq] = temporary(fvar[fpx:epx])
          finspec.novar[0:npix-1,qq] = temporary(fnovar[fpx:epx])
          bd = where(finspec.wave[*,qq] LE 0,nbd)
          finspec.var[bd,qq] = 0.
          finspec.novar[bd, qq] = 0.
      endif

;      gdp = where(allobj[all_ordr[0]].wave GT 0.,ngdp)
;      fpx = gdp[0]
;      epx = gdp[ngdp-1]
;      npix = epx-fpx+1
;      x_splot, finspec.wave[gdp,qq], finspec.fx[gdp,qq], /block, TITLE='sigma'

      ;; FIX SIGMA FOR > 2 EXPOSURE CASE
      if nobj GT 2 then begin
          usig = allobj[all_ordr[0]].sig
          bad = where(mask[*,0] EQ 0,nbad)
          if nbad GT 0 then usig[bad] = -1.
          gd = where(usig GT 0., ngd)
          ntot = ngd
          if ngd GT 0 then $
            chi = (allobj[all_ordr[0]].flux[gd]-fflux[gd])/usig[gd]
          for i=1, nobj-1 do begin
              usig = allobj[all_ordr[i]].sig
              bad = where(mask[*,i] EQ 0,nbad)
              if nbad GT 0 then usig[bad] = -1.
              gd = where(usig GT 0.,ngd)
              if ngd GT 0 AND ntot GT 0 then $
                chi = [chi,(allobj[all_ordr[i]].flux[gd]-fflux[gd])/usig[gd]]
              if ngd GT 0 AND ntot EQ 0 then $
                chi = (allobj[all_ordr[i]].flux[gd]-fflux[gd])/usig[gd]
              ntot = ntot+ngd
          endfor
          djs_iterstat, chi, sigma=sigma, sigrej=sigrej, maxiter=5
;;;TEST
          if ntot GT 0 then ochihist = histogram(chi,min=-5,max=5,binsize=0.1)
          usig = sigma*allobj[all_ordr[0]].sig
          bad = where(mask[*,0] EQ 0,nbad)
          if nbad GT 0 then usig[bad] = -1.
          ntot = 0
          gd = where(usig GT 0., ngd)
          ntot = ntot+ngd
          if ngd GT 0 then chi = (allobj[all_ordr[0]].flux[gd]-fflux[gd])/usig[gd]
          for i=1, nobj-1 do begin
              usig = sigma*allobj[all_ordr[i]].sig
              bad = where(mask[*,i] EQ 0,nbad)
              if nbad GT 0 then usig[bad] = -1.
              gd = where(usig GT 0.,ngd)
              ntot = ntot+ngd
              if ngd GT 0 then $
                chi = [chi,(allobj[all_ordr[i]].flux[gd]-fflux[gd])/usig[gd]]
          endfor
          if ntot GT 0 then begin
              chihist = histogram(chi,min=-5,max=5,binsize=0.1)
              djs_iterstat, chi, sigma=sigma2, sigrej=sigrej, maxiter=5
              title = 'Sigma Original: ' + strtrim(string(sigma),2) + $
                      ' New Sigma: ' + strtrim(string(sigma2),2)
              ;;   x_splot, findgen(101)*0.1-5, ochihist, ytwo=chihist, /block, TITLE=title

              x_combspec, allobj[all_ordr].flux, ((sig>0.)*sigma)^2, $
                          nflux, nvar, WAVE=allobj[all_ordr[0]].wave, $
                          SNR=svsnr[*,qq], SCALE=svmedian[*,qq], mask=mask, $
                          NSIG=REJSIG, novar = ((nosig>0.)*sigma)^2, $
                          fnovar = novar
;          gdp = where(nvar GT 0. and finspec.wave GT 0.,ngdp)
              gdp = where(nvar GT 0. and allobj[all_ordr[0]].wave GT 0.,ngdp)
              fpx = gdp[0]
              epx = gdp[ngdp-1]
              npix = epx-fpx+1
;              x_splot, finspec.wave[gdp,qq], fflux[gdp], ytwo=nflux[gdp], /block, TITLE='Flux'
;              x_splot, finspec.wave[gdp,qq], sqrt(fvar[gdp]), ytwo=sqrt(nvar[gdp]), /block, TITLE='sigma'
              finspec.wave[0:npix-1,qq] = allobj[all_ordr[0]].wave[fpx:epx]
              finspec.fx[0:npix-1,qq] = temporary(nflux[fpx:epx])
              finspec.var[0:npix-1,qq] = temporary(nvar[fpx:epx])
              finspec.novar[0:npix-1, qq] = temporary(novar[fpx:epx])
          endif
      endif 

      ;;
      if keyword_set( FCHK ) then begin
          wuse = where(finspec.wave[*,qq] GT 1000.)
          x_splot, finspec.wave[wuse,qq], finspec.fx[wuse,qq], $
                   TITLE='Combined -- Order: '+strtrim(qq,2), $
                   ytwo=sqrt(finspec.var[wuse,qq]), /block
;                   ythr=allobj[all_ordr[1]].flux[wuse], psym3=-3
      endif
  endfor

  return
end
