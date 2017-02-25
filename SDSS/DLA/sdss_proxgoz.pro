;+ 
; NAME:
; sdss_proxgoz
;
; PURPOSE:
;    Determines the g(z) path for each QSO in the inputted list.
;    Restrict to the proximate DLAs
;
; CALLING SEQUENCE:
;   sdss_proxgoz
;
; INPUTS:
;  qalfil -- Filename of QSO structure for g(z) calculation
;  XZFIL= -- Filename of the dX calculation as a function of z.
;
; RETURNS:
;
; OUTPUTS:
;  OUTFIL= -- Filename for g(z) FITS file
;
; OPTIONAL KEYWORDS:
;  ORIG_ZEM=  -- Use the original redshifs from SDSS (in this file)
;  SNRLMT=    -- SNR limit for DLA search [default: 4]
;
; OPTIONAL OUTPUTS:
;  /CHK  -- Plot g(z) to the screen
;
; COMMENTS:
;
; EXAMPLES:
;  sdss_proxgoz, 'sdss_dr3_QAL.fits', 'dr3_proxdlagz_s2n4.fits'
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by SHF
;   14-Oct-2004 Modified by JXP
;-
;------------------------------------------------------------------------------
pro sdss_proxgoz, qalfil, outfil, ntrials=ntrials, XZFIL=xzfil, $
                  SNRLMT=snrlmt, CHK=chk, ISTRT=istrt, ORIG_ZEM=orig_zem

  if not keyword_set( XZFIL ) then $
    xzfil = '/u/xavier/SDSS/DR2_QSO/xz_val.fits'
  if not keyword_set( OFFL ) then offl = 0.  ;; Offset to start wave
  if not keyword_set( SNRLMT ) then snrlmt = 4.
  if not keyword_set( ISTRT ) then istrt = 0L

 ;;read in spectra
  x=xmrdfits(qalfil,1,/silent)
  nspectra=n_elements(x)
  f=findgen(5500L)*.001
  g=lonarr(5500L)

  ;; XZFIL
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)

  ;; RA
  ra=dblarr(nspectra)
  dec=dblarr(nspectra)
  zem=dblarr(nspectra)
  if not keyword_set(ntrials) then ntrials=nspectra

  ;; zem
  if keyword_set(ORIG_ZEM) then qsos = xmrdfits(orig_zem,1)

  ;; 8000 km/s deal
  c = x_constants()
  vc = 8000. / (c.c / 1e5)
  ovc = 1. - vc


  ;;loop over spectra 
  gdqso = [-1L]
  allz1 = [0.d]
  allz2 = [0.d]
  dX = [0.d]
  svplt = [0L]
  svfib = [0L]
  svra = [0.d]
  svbal = [0]
  svdec = [0.d]
  svzem = [0.d]
  svmag = [0.d]
  svs2n = [0.]
  for i=istrt,ntrials-1 do begin
      if i MOD 100 EQ 0 then print, 'i = ', i
      if x[i].z_qso LE 2.15 then continue

      name=x[i].file_name
      
      ra[i]=x[i].ra
      dec[i]=x[i].dec
      zem[i]=x[i].z_qso
      
      if keyword_set(ORIG_ZEM) then begin
          mt = where(x[i].plate EQ qsos.plate and x[i].fiberid EQ qsos.fiberid, nmt)
          if nmt EQ 1 then zem[i] = qsos[mt].z 
      endif

      ;; Check for double
      if i gt 0 then begin   
          dbl=where(abs(ra[i]-ra[0:i-1]) le 1E-3 $
                    and abs(dec[i]-dec[0:i-1]) le 1E-3 and $
                    abs(zem[i]-zem[0:i-1]) le .1,ndbl) 
          if ndbl NE 0 then begin
              print, 'sdss_dr2_goz: Repeat quasar!'
          endif
      endif else ndbl=0

      ;; Recalculate start wave
      istrt = x[i].start_wave
      parse_sdss, getenv('SDSSPATH')+x[i].file_name, flux, wave, SIG=sig
      osig = sig
      s=where(flux GT 0. AND sig GT 0., complement=bad, ncomplement=nbad)
      if nbad NE 0 then begin 
          flux[bad] = 0.
          sig[bad] = -1.
      endif
      
      ;; Lyb emission (and OVI)
      msk = replicate(1B, n_elements(flux))
      msk[0:10] = 0B  ;; Deals with median issue
      msk2 = msk
      vel = (wave/(1025.7*(1+x[i].z_qso)) - 1.)*2.997924d5


      pixlyb = where(abs(vel) LT 5000., npix)
      if npix NE 0 then begin
          ;; Old way and was previously bugged!!  
          ;; Thanks to MM for cathching this (after the paper was
            ; accepted!)
;          bad = where(median(flux[pixlyb]/sig[pixlyb],20) LT SNRLMT*3.,nbad)
;          if nbad NE 0 then msk[pixlyb[bad]] = 0B
          msk[pixlyb] = 0B
      endif

      ;; Start point
      start=where(median(flux/sig, 20) GT SNRLMT $
                  AND sig GT 0. AND msk, nstart)

      if nstart EQ 0 then x[i].start_wave = -1. $
      else x[i].start_wave = wave[start[0]]
      if keyword_set( CHK ) then begin
          if abs(istrt-x[i].start_wave) GT 5. and istrt NE 0. $
            and start[0] GT 15 and istrt GT 3833. then begin
              print, i, istrt, x[i].start_wave, x[i].z_qso
              stop
          endif
      endif

      ;; S2N
      if x[i].z_qso LT 5. then begin
          if x[i].flg_bal EQ 0 then lim = [1440., 1490.] $
          else lim = [1440., 1470.]
      endif else begin ;; High z
          lim = [1287., 1366.]
      endelse
      lim = lim*(1+x[i].z_qso)
      gd = where(wave GT lim[0] AND wave LT lim[1] AND sig GT 0.,ngd)
      if ngd EQ 0 then x[i].snr = 0. else begin
          djs_iterstat, flux[gd]/sig[gd], mean=mean, maxiter=5, sigrej=2.
          x[i].snr = mean
      endelse

      ;; Look for stretch of bad data
      bmed = median(osig, 50L)
      bad = where(bmed LE 0.,  nbad)
      flg_bad = 0L
      if nbad NE 0 then begin
          bb = where(bad GT 50L and bad LT (n_elements(sig)-50L),nbb)
          if nbb NE 0 then begin
              flg_bad = 1L
              l2b = wave[bad[bb[0]]]
          endif
      endif 

      ;; Calculate z1, z2
      if ndbl le 0 and x[i].start_wave GT 0. then begin
          flg = 0
          zstrt = (x[i].start_wave)/1215.6701 -1
          case x[i].flg_bal of
              0: begin
                  z2 = x[i].z_qso   ; 
                  z1 = zstrt*1.005 + 0.005  ; Add 1500 km/s to start wave
              end
              1: begin  ;; Modest BAL
                  z1=(((1060.*(1+x[i].z_qso))/1215.6701)-1) $
                    > (zstrt*1.005 + 0.005)
                  z2=x[i].z_qso 
              end
              2: begin  ;; Rejected BAL
                  z1 = -1.
                  z2 = -1.
                  flg = -1
              end
              else: flg = -1
          endcase

          ;; Min on z1 of 2.2 for SDSS
          z1 = z1 > 2.2

          ;; Require z1 be 8000 km/s from z2
          if z1 GT (z2*ovc - vc) then flg = -1

          ;; Deal with 'Null data'
          if flg_bad EQ 1 then z2 = z2 < (((l2b/1215.6701) -1.)*0.995-0.005)

          ;; gz
          gdi=where(f ge z1 and f le z2, ngd)
          if ngd ne 0 then g[gdi]=g[gdi]+1 else flg=-1
          if flg NE -1 then begin
              ;;
              gdqso = [gdqso,i]
              allz1 = [allz1, z1]
              allz2 = [allz2, z2]
              svbal = [svbal, x[i].flg_bal]
              svra = [svra, ra[i]]
              svdec = [svdec, dec[i]]
              svzem = [svzem,zem[i]]
              svplt = [svplt, x[i].plate]
              svfib = [svfib, x[i].fiberid]

              ;; Magnitude and s2n
              svmag = [svmag,x[i].qso_mag]
              svs2n = [svs2n, x[i].snr]

              ;; Calculate dX
              mn = min(abs(xzz-z1),imn1)
              mn = min(abs(xzz-z2),imn2)
              dX = [dx, (xzv[imn2]-xzv[imn1])>0.] 
          endif
      endif
  endfor

  ;;plot
  x_splot,f,g,/block
  maxg=max(g)
;stop
  old_device=!d.name
  set_plot,'ps'
  device,filename='prox_goz.ps'
;device,xsize=7,ysize=7,/inches
  plot,f,g,yrange=[0,maxg+1]
  device,/close_file
  set_plot,old_device

  ;; OUT
  if keyword_set( OUTFIL ) then begin
      strct = { $
                gzz: f, $
                gzv: g, $
                plate: svplt[1:*], $
                fib: svfib[1:*], $
                ra: svra[1:*], $
                dec: svdec[1:*], $
                flg_bal: svbal[1:*], $
                iqso: gdqso[1:*], $
                mag: svmag[1:*], $
                s2n: svs2n[1:*], $
                z1: allz1[1:*], $
                z2: allz2[1:*], $
                zem: svzem[1:*], $
                dX: dX[1:*] $
              }
      mwrfits,strct,outfil,/create
      spawn, 'gzip -f '+outfil
  endif
  print, 'sdss_goz: All done!'

end

