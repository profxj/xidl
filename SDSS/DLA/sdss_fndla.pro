;+ 
; NAME:
; sdss_fndla
;
; PURPOSE:
;    Calculates f(N,X) for a set of DLA
;
; CALLING SEQUENCE:
;
; INPUTS:
;  GZSTR= -- Structure summarizing g(z) for the quasars
;  GZFIL= -- Filename of the g(z) file for the quasars
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /DR3  -- Restrict to DR3
;  /NODBL -- Do not solve for double power-law
;  CL=   -- Confidence limits for error calculations
;  ZBIN= -- redshift interval [default: 2.2 to 5.5]
;  /SALL  -- Use all of the DR releases
;  BINS  -- NHI bins
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
;   20-Aug-2006 Written by JXP  (based on fig_fndla)
;-
;------------------------------------------------------------------------------
pro sdss_fndla, GZFIL=gzfil, DLALST=dlalst, ZBIN=zbin, VPROX=vprox, $
                STRCT=strct, GZSTR=gzstr, PRX=prx, SALL=sall, BINS=bins, $
                stp=stp, NODBL=nodbl, OUTFIL=outfil, NMAX=nmax, NPLT=nplt, $
                DR3=dr3, ALLDR=alldr, CL=cl, GMMNHI=gmmnhi, NMIN=nmin, DX=dx,$
                XZFIL=xzfil

  ;; Get structure if necessary
;  if not keyword_set( GZFIL ) then $
;    gzfil = getenv('SDSSPATH')+'/DR5_QSO/dr5_dlagz_s2n4.fits'
  if not keyword_set( NPLT ) then nplt = 18L
  if not keyword_set( STP ) then stp = 0.1
  if not keyword_set( NMIN ) then nmin = 20.3
  if not keyword_set( ZBIN ) then zbin = [2.2, 5.5]
  if not keyword_set(GMMNHI) then GMMNHI = [20.6, 21.9]
  
  ;; GZ
;  if not keyword_set(DR3) then ALL=1 else all=0

  ;; SDSS
  if not keyword_set(ALLDR) then $
    sdss_dlastrct, alldr, /stat, GZSTR=gzstr, VPROX=vprox, ADR3=DR3, $
                   ALL=sall, GZDR3=DR3 $
  else print, 'sdss_fndla: WARNING -- Using input ALLDR'
  
  i2 = where(alldr.zabs GE zbin[0] AND $
             alldr.zabs LE zbin[1], nsdss)
  if nsdss EQ 0 then return
  alldr = alldr[i2]
  
  ;; Peroux DLA
  if keyword_set(PRX) then begin
      stop
      allprx = peroux_dla(zbin, ALLDR=alldr)
      if nsdss EQ 0 then alldr = allprx else $
        alldr = [alldr, allprx]
  endif
  
  ;; dX
  if not keyword_set(dX) then $
    dX = sdss_dladx( zbin[0], zbin[1], vprox, XZFIL=xzfil, $
                     GZSTR=gzstr, GZFIL=gzfil, DZ=dz)
  if keyword_set(PRX) then begin
      pdX = peroux_dx( zbin[0], zbin[1], /uniq, GZSTR=gzstr, GZFIL=gzfil)
      dX = dX + pdX
  endif
  ndr = n_elements(alldr)

  ;; f(N) plot
  if not keyword_set(BINS) then begin
      bins = 20.3 + findgen(nplt)*stp
  endif else begin
      szbb = size(bins,/dimens)
      nplt = szbb[1]
  endelse
  yplt = dblarr(nplt)
  yerr1 = dblarr(nplt)
  yerr2 = dblarr(nplt)
  xplt = dblarr(nplt)
  o3xplt = dblarr(nplt)
  xval = dblarr(nplt)
  x1 = dblarr(nplt)
  x2 = dblarr(nplt)
  svdn = dblarr(nplt)
  svm = lonarr(nplt)

  for qq=0L,nplt-1 do begin
      xval[qq] = bins[qq] + stp/2.
      gd = where(alldr.NHI GE (bins[qq] - 1e-5) AND $
                 alldr.NHI LT (bins[qq]+stp - 1e-5),ngd)
      dN = 10^(bins[qq]+stp) - 10^bins[qq]
      svdN[qq] = dN
      case ngd of
          0: xplt[qq] = xval[qq]
          1: xplt[qq] = alldr[gd].NHI
          else: begin
              xplt[qq] = alog10(mean(10^alldr[gd].NHI))
              o3xplt[qq] = mean(alldr[gd].NHI)
          end
      endcase
      svm[qq] = ngd
      if ngd NE 0 then begin
          yplt[qq] = alog10(ngd/dN/dX) 
          ;; Poisson 68.3% (1sigma)
          val = x_poisscl(double(ngd), sigma=1)
          ;; Fill
          yerr1[qq] = alog10(val[0]/dN/dX) 
          yerr2[qq] = alog10(val[1]/dN/dX) 
      endif else begin
          yplt[qq] = -99.
          val = x_poisscl(0., sigma=2)
          yerr1[qq] = alog10(val[0]/dN/dX) 
      endelse
  endfor
  if abs(total(svm)-n_elements(alldr)) GT 1e-5 then stop

  ;; Output
  if keyword_set(OUTFIL) then $
    writecol, outfil, xplt, yplt, yerr1, yerr2

  ;; Function
  fplt = 20.2 + (22-20.2)*findgen(1000L)/1000.


  ;; Single power law
  maxb = x_maxsngpow(10^alldr.NHI, NMIN=10.^NMIN, SMM=smmnh, $
                    KSPROB=ksprob, NOISE=0.05, ERR=errs, CL=0.683, $
                    ARNG=[-7., -1.1], NMAX=nmax)
  if not keyword_set(NMAX) then $
    k1 = -1.*float(ndr/dX) * (maxb+1.) / (10.^Nmin)^(1.+maxb) $
  else $
    k1 = -1.*float(ndr/dX) * (maxb+1.) / $
         ((10.^Nmin)^(1.+maxb) - (10.^Nmax)^(1.+maxb))

  sigk1 = x_poisscl(double(ndr), sigma=1)*k1/float(ndr)
  thy = k1 * (10.^fplt)^(maxb)

  print, 'Maxb = ', maxb
  if keyword_set(KSPROB) then print, 'Single KS = ', ksprob

  ;; Gamma function
  x_maxschecht, 10^alldr.NHI, GMMNHI, $
    [-2.3, -1.35], val, Nmin=2e20, LIK=lik, KSPROB=ksprob, NOISE=0.05, $
    ERR=err, /log, CL=CL
  print, 'Gamma val = ', val
  print, 'Gamma KS = ', ksprob
  k2 = float(ndr/dX) / val[0] / (gamma(val[1]+1.)* $
                              (1.-igamma(val[1]+1.,2e20/val[0])))
  errk2 = x_poisscl(double(ndr), sigma=1)*k2/float(ndr)
  gerr = dblarr(2,3)
  gerr[*,0:1] = err
  gerr[*,2] = errk2
  thy = k2 * (10.^fplt/val[0])^(val[1])*exp(-1.*(10.^fplt/val[0]))

  ;; Double Power-Law
  if not keyword_set( NODBL ) then begin
      x_maxdblpow, 10^alldr.NHI, [21., 21.9], [-3.,-1.1], [-10.,-1.8], vald, $
        STN=0.005d, NB1=500L, NB2=500L, NMIN=(10^20.3), KSPROB=ksprob, $
        NOISE=0.05, ERR=err
      print, 'Double power law: ', vald
      print, 'Double power KS = ', ksprob
      ns10 = 10.^vald[0]
      if not keyword_set(NMAX) then nmax = 1e99
      k3 = float(ndr/ns10/dX) / $
        ( (1. - (10^Nmin/Ns10)^(vald[1]+1))/(1.+vald[1]) + $
          ((Nmax/Ns10)^(vald[2]+1) - 1.)/(vald[2]+1)) 
      errk3 = x_poisscl(double(ndr), sigma=1)*k3/float(ndr)
      derr = dblarr(2,4)
      derr[*,0:2] = err
      derr[*,3] = errk3
      thy = fltarr(n_elements(fplt))
      lw = where(fplt LT vald[0], complement=hi, ncomplement=numhi)
      thy[lw] = k3 * (10^fplt[lw]/ns10)^vald[1]
      if numhi NE 0 then thy[hi] = k3 * (10^fplt[hi]/ns10)^vald[2]
  endif else begin
      k3 = 0.
      vald = fltarr(3)
      derr = fltarr(3,3)
  endelse

  ;; Structure
  strct = { $
          zmean: mean(alldr.zabs), $ ; Mean redshift
          xval: xval, $      ; Center of the bin
          xplt: xplt, $      ; Log of the mean
          o3xplt: o3xplt, $  ; Mean of the log values
          yplt: yplt, $
          mdla: svm, $
          bins: bins, $
          yerr1: yerr1, $
          yerr2: yerr2, $
          dX: dX, $
          dz: dz, $
          dN: svdN, $
          k1: k1, $
          sigk1: sigk1, $
          a1: maxb, $
          siga1: errs, $
          k2: k2, $
          a2: val[1], $
          Ng: val[0], $
          gerr:gerr, $
          k3: k3, $
          Nd: vald[0], $
          a3: vald[1], $
          a4: vald[2], $
          derr:derr $
          }
            
    
  print, 'sdss_fndla: All done!'

  return
end
      
      
