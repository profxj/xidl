;+ 
; NAME:
; sdss_llsgoz
;
; PURPOSE:
;    Determines the g(z) path for a LLS search
;
; CALLING SEQUENCE:
;   sdss_llsgoz
;
; INPUTS:
;  qalfil -- Filename of QSO structure for g(z) calculation
;
; RETURNS:
;
; OUTPUTS:
;  OUTFIL= -- Filename for g(z) FITS file
;
; OPTIONAL KEYWORDS:
;  ORIG_ZEM=  -- Use the original redshifs from SDSS (in this file)
;  SNRLMT=    -- SNR limit for DLA search [default: 2]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  sdss_llsgoz, '../dr5_qso.fits', 'dr5_lls_goz_sn2.fits'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by SHF
;   14-Oct-2004 Modified by JXP
;-
;------------------------------------------------------------------------------
pro sdss_llsgoz, qso_fil, outfil, ntrials=ntrials, XZFIL=xzfil, BALFIL=balfil, $
                  SNRLMT=snrlmt, CHK=chk, ISTRT=istrt, ORIG_ZEM=orig_zem

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_llsgoz, qso_fil, outfil, ORIG_ZEM='
    return
  endif 

  if not keyword_set(MINZ1) then minz1 = 3.28
  if not keyword_set(SNRLIM) then snrlim = 2.
  if not keyword_set(DRPATH) then drpath = getenv('SDSSPATH')+'/DR7_QSO/'
  if not keyword_set(QSOtempwav0) then QSOtempwav0 = 2.8
  if not keyword_set( XZFIL ) then $
    xzfil = getenv('SDSSPATH')+'/DR7_QSO/xz_val.fits'

  if not keyword_set(BALFIL) then $
    balfil = getenv('SDSSPATH')+'/DR5_QSO/LLS/lls_bal.lst'

 ;;read in spectra
  x=xmrdfits(qso_fil,1,/silent)
  nspectra=n_elements(x)
  f=findgen(5500L)*.001
  g=lonarr(5500L)

  ;; XZFIL
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)


  qsodr=xmrdfits(qso_fil,1,/silent)

  szstart = [3.2, 3.4, 3.7, 4.0, 4.4]
  szend =   [3.4, 3.7, 4.0, 4.4, 5.0]
  tfil = ['LLS/hiztemplate_32_34.fits', $
          'LLS/hiztemplate_34_37.fits', $
          'LLS/hiztemplate_37_40.fits', $
          'LLS/hiztemplate_40_44.fits', $
          'LLS/hiztemplate_44_50.fits']

  nzcut = n_elements(szstart)

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

  ;; REMOVE BAL + BAD QSOs
  readcol, BALFIL, bad_plate, bad_fib, format='L,L', delim=','
  nbad = n_elements(bad_plate)

  ;; Big Loop
  for kk=0L,nzcut-1 do begin
      
      good = where(qsodr.z GE szstart[kk] AND qsodr.z LT szend[kk],nqsos)

      template=xmrdfits(drpath+tfil[kk],/silent)
      
      clr = getcolor(/load)
      
      ;; LOOP
      for qq=0L,nqsos-1 do begin

          ;; Bad?
          bada = where(qsodr[good[qq]].plate EQ bad_plate AND $
                       qsodr[good[qq]].fiberid EQ bad_fib, nbada)
          if nbada NE 0 then continue

          ;;
          if qq MOD 100 EQ 0 then print, qq
          IF(qsodr[good[qq]].plate LT 1000) THEN BEGIN
              splate=strcompress('0'+string(qsodr[good[qq]].plate),/remove_all)
          ENDIF ELSE BEGIN
              splate=strcompress(string(qsodr[good[qq]].plate))
          ENDELSE
          
          smjd=string(qsodr[good[qq]].mjd, format='(i5.5)')
          
          sfiber=string(qsodr[good[qq]].fiberid, format='(i3.3)')
          
          filename = strcompress('spSpec-'+smjd+'-'+splate+'-'+sfiber+'.fit.gz', $
                                 /remove_all)
          filedir = strcompress(drpath+'spectro/1d_26/'+splate+'/1d/', $
                                /remove_all)
          file=filedir+filename
          
          ;; Data
          spspec=xmrdfits(file,0,hdr,/silent)
          spec = spspec[*,0]
          specerr = spspec[*,2]
          speczem = qsodr[good[qq]].z
          specwav0 = sxpar(hdr,"COEFF0")
          logwav = findgen(n_elements(spec))*0.0001 + specwav0
          
          ;; Template
          templ_spec = sdss_qsotempl(template, QSOtempwav0, logwav, $
                                     spec, specerr, speczem, $
                                     EPX=epx, FPX=fpx)
          
          wv = 10^logwav[0:epx-fpx]

          ;; Smooth sigma
          sig = median(specerr,15)
          sig[0:10] = specerr[0:10]
          sig = sig[0:epx-fpx]

          ;; Cut
          a = where(templ_spec GT sig*SNRLIM and specerr[0:epx-fpx] GT 0., na)

          ;; Plot
          if keyword_set(PLOT) then begin
              plot, wv, spec[0:epx-fpx], color=clr.black, background=clr.white, $
                    xrange=[3800.,6000.]
              oplot, wv, sig, color=clr.blue, linestyle=10
              oplot, wv, templ_spec, color=clr.red
              if na NE 0 then $
                oplot, replicate(wv[a[0]],2), [-1e9,1e9], color=clr.green, linesty=2
              wait, 1.5
          endif

          if na ne 0 then begin
              z1 = (wv[a[0]] / 911.764 - 1) > MINZ1
              z2 = .99*speczem - 0.01 ; 3000 km/s blue of Lya
              ;; Need to deal with BAL
              gdi=where(f ge z1 and f le z2, ngd)
              if ngd ne 0 then begin
                  ;; Fill up
                  g[gdi]=g[gdi]+1 

;                  gdqso = [gdqso,i]
                  allz1 = [allz1, z1]
                  allz2 = [allz2, z2]
;                  svbal = [svbal, x[i].flg_bal]
                  svra = [svra, qsodr[good[qq]].raobj]
                  svdec = [svdec, qsodr[good[qq]].decobj]
                  svzem = [svzem,speczem]
;                  if qsodr[good[qq]].plate EQ 1958 AND $
;                    qsodr[good[qq]].fiberid EQ 449 then stop
                  svplt = [svplt, qsodr[good[qq]].plate]
                  svfib = [svfib, qsodr[good[qq]].fiberid]
                  
                  ;; Magnitude and s2n
                  svmag = [svmag,qsodr[good[qq]].psf_r]
;                  svs2n = [svs2n, x[i].snr]

                  ;; Calculate dX
                  mn = min(abs(xzz-z1),imn1)
                  mn = min(abs(xzz-z2),imn2)
                  dX = [dx, (xzv[imn2]-xzv[imn1])>0.] 
              endif
          endif
      endfor
  endfor

  if not keyword_set(NOPLOT) then x_splot,f,g,/block
  maxg=max(g)

  ;; OUT
  nqso = n_elements(svzem[1:*])
  strct = { $
          gzz: f, $
          gzv: g, $
          plate: svplt[1:*], $
          fib: svfib[1:*], $
          ra: svra[1:*], $
          dec: svdec[1:*], $
          mag: svmag[1:*], $
          z1: allz1[1:*], $
          z2: allz2[1:*], $
          zem: svzem[1:*], $
          flg_bal: replicate(0,nqso), $
          dX: dX[1:*] $
          }
;          flg_bal: svbal[1:*], $
;          iqso: gdqso[1:*], $
;          s2n: svs2n[1:*], $

  ;; Out
  mwrfits,strct,outfil,/create
  spawn, 'gzip -f '+outfil

end

