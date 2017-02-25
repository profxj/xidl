;+ 
; NAME:
; sdss_omgdla
;
; PURPOSE:
;    Calculates the first moment of f(N,X) [i.e. Omega_DLA]
;    Includes He as the default
;
; CALLING SEQUENCE:
;  sdss_omgdla, GZFIL=, PEROUX=, MODELS=, STRCT=, BINS=, VPROX=, DR3=, ALL=
;
; INPUTS:
;  GZFIL= -- Filename of the g(z) file for the quasars
;
; RETURNS:
;
; OUTPUTS:
;  STRCT=  -- Structure summarizing the calculation
;
; OPTIONAL KEYWORDS:
;  BINS=  -- Redshift intervals, e.g. [ [1,2], [2,3]]
;  /PEROUX -- Include DLAs in Peroux03 in the stats
;  /DR3  -- Restrict the DLAs to DR3
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Nov-2006 Written by JXP  (based on fig_omega)
;-
;------------------------------------------------------------------------------
pro sdss_omgdla_boot, NHI, N1, N2, NTRIAL=ntrial, POISS=poiss
  ;; Ntrial
  if not keyword_set( NTRIAL ) then ntrial = 1000L
  if not keyword_set( SEED ) then seed = -1997

  nval = n_elements(NHI)

  rnd = fix(randomu(seed, ntrial, 10*nval)*float(nval))

  sumNH = dblarr(ntrial)
  if not keyword_set( POISS ) then begin
      for qq=0L,ntrial-1 do sumNH[qq] = total( 10^NHI[ rnd[qq,*] ], /double)
  endif else begin
      rndn = randomn(seed, ntrial)
      for qq=0L,ntrial-1 do begin
          nsys = (nval + round(sqrt(nval)*rndn[qq])) > 1
          sumNH[qq] = total( 10^NHI[ rnd[qq,0:nsys-1] ], /double)
      endfor
  endelse
  
  ;; Stats
  srt = sort(sumNH)
  
  N1 = sumNH[ srt[ round(ntrial*0.165) ] ]
  N2 = sumNH[ srt[ round(ntrial*(1-0.165)) ] ]

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_omgdla, GZFIL=gzfil, PEROUX=peroux, GZSTR=gzstr, $
                 MODELS=models, STRCT=strct, BINS=bins, HUB=hub, $
                 VPROX=vprox, DR3=dr3, ALL=all, CST=cst, XZFIL=xzfil,$
                 VMIN=vmin, ADR=alldr

  ;; Get structure if necessary
  if not keyword_set( GZFIL ) then $
    gzfil = '/u/xavier/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'

  ;; read gzfil
  if not keyword_set(GZSTR) then gzstr = xmrdfits(gzfil, 1, /silent)

  ;; SDSS
  if not keyword_set(ALLDR) then $
     sdss_dlastrct, alldr, /stat, GZSTR=gzstr, VPROX=vprox, ADR3=DR3, ALL=all,$
                    GZFIL=gzfil, VMIN=vmin

  ;; Read DLA
  ndr = n_elements(alldr)
  print, 'NDLA = ', ndr

  ;; Bins
  if keyword_set(PEROUX) then prxdla = peroux_dla(ALLDR=alldr)

  ;; 
  c = x_constants()
  AMU=1.3
  if not keyword_set(HUB) then HUB = 70.
  Cst=(amu*hub*1e5/c.Mpc)*c.mp/c.c/(c.rhoc*hub/100.*hub/100.)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; TOTAL
  if not keyword_set( BINS ) then $
    bins = [ [1.7,2.2], [2.2, 2.5], [2.5, 3.0], [3.0, 3.5], $
             [3.5,4.0], [4.,5.3]]
  if n_elements(bins) EQ 2 then sz = [2,1] else sz = size(bins,/dimensions)
  nbin = sz[1]

  ;; Structrue
  tmp = { $
          dX: 0., $
          pdX: 0., $
          sall: 0., $
          pall: 0., $
          medz: 0., $
          bins: fltarr(2), $
          omega: 0.d, $
          somg: dblarr(2) $
        }
  strct = replicate(tmp, nbin)

  ;; LOOOP
  for ii=0L,nbin - 1 do begin
      dX = sdss_dladx( bins[0,ii], bins[1,ii], vprox, $
                       GZSTR=gzstr, XZFIL=xzfil, VMIN=vmin)
      if keyword_set(PEROUX) then $
        pdX = peroux_dx( bins[0,ii], bins[1,ii], /uniq, GZSTR=gzstr) $ 
      else pdX = 0.

      ;; SDSS
      gd = where(alldr.zabs GE bins[0,ii] AND alldr.zabs LT bins[1,ii] AND $
                 alldr.NHI GE 20.3, nsdss)
      
      ;; Peroux
      if keyword_set(PEROUX) then $
        gdp = where(prxdla.zabs GE bins[0,ii] and prxdla.zabs LT bins[1,ii] AND $
                    prxdla.NHI GE 20.3, nprx)

      ;; Sum SDSS
      if gd[0] NE -1 then begin
          allNH = total(10^alldr[gd].NHI, /double) 
          avsdss = allNH / float(nsdss)
      endif else begin
          dX = 0.
          allNH = 0.
          avsdss = 0.
          ngd = 0
      endelse


      ;; Lanzetta et al 1991
;      varsdss = total( (10.d^dr3[gd].NHI-avsdss)^2, /double)  / $
;        (1-1/float(ngd)) / dX^2 / float(ngd)

      if keyword_set(PEROUX) then begin
          prxNH = total(10^prxdla[gdp].NHI,/double)
          avprx = prxNH / float(nprx)
          varprx = total( (10.d^prxdla[gdp].NHI-avprx)^2 ,/double) / pdX^2 / $
                   (1-1/float(nprx)) / float(nprx)
      endif else prxNH = 0.
      ;; Structure
      strct[ii].dX = dX
      strct[ii].pdX = pdX
      strct[ii].sall = allNH
      strct[ii].pall = prxNH

      ;; Total
      allNH= allNH + prxNH
      tdx = dX + pdX

      
      ;; Median z
      if keyword_set(PEROUX) then begin
          if gd[0] NE -1 then $
            medz = (total(alldr[gd].zabs*10^alldr[gd].NHI) $
                    + total(prxdla[gdp].zabs*10^prxdla[gdp].NHI)) / allNH $
          else medz = total(prxdla[gdp].zabs*10^prxdla[gdp].NHI) / allNH
      endif else begin
          medz = total(alldr[gd].zabs*10^alldr[gd].NHI) / allNH 
      endelse
          
      ;; Omega
      omega = allNH*Cst/tdX
;      print, 'Const: cst', cst

      ;; Sig
      if keyword_set( PEROUX ) then begin
          if gd[0] NE -1 then $
            sdss_omgdla_boot, [alldr[gd].NHI, prxdla[gdp].nhi], sig1, sig2, /POISS $
          else sdss_omgdla_boot, [prxdla[gdp].nhi], sig1, sig2, /POISS
      endif else begin
          sdss_omgdla_boot, [alldr[gd].NHI], sig1, sig2, /POISS 
      endelse
      sig1 = sig1*Cst/tdX
      sig2 = sig2*Cst/tdX

      strct[ii].medz = medz
      strct[ii].bins = bins[*,ii]
      strct[ii].omega = omega
      strct[ii].somg = [sig2,sig1]

      midz = (bins[0,ii]+bins[1,ii])/2.
  endfor

  print, 'sdss_omgdla: All done!'

  return
end
      
      
