;+ 
; NAME:
; sdss_qsolin   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to SDSS QSO spectroscopic data interactively
;    using a PCA as a guess for the Lya and CIV emission lines.  The
;    code then automoatically searches for absorption line features
;    redward of Lya.
;
; CALLING SEQUENCE:
;   sdss_qsolin, sdss_list
;
; INPUTS:
;  sdss_list -- List of QSOs to analyze
;
; RETURNS:
;
; OUTPUTS:
;  One FITS file per QSO which contains the continuum and a list of
;  absorption lines detected
;
; OPTIONAL KEYWORDS:
;  LSNR= -- Statistical significance for an absorption line feature
;           [default: 3.5]
; /NOPLOT -- Do not plot to screen
; /REDO  -- Redo the search even if the file exists
; /NOCOMPRESS -- Do not gzip suprress the output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_sdss, fil
;       sdss_qsolin, 'dr1_fit.lst', OUTDIR='ABSLIN/', /noplt
;       sdss_qsolin, 'dr3_fit.list', OUTDIR='ABSLIN/', /noplt
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   09-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_qsolin, sdss_list, LSNR=lsnr, CCHK=cchk, PSFIL=psfil, $
                 OUTDIR=outdir, NOPLT=noplt, REDO=REDO, ISTRT=istrt, $
                 NOCOMPRESS=nocompress

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'parse_sdssqso, sdss_list, PSFIL=, OUTDIR=, /NOCOMPR, ISTRT= [v1.0]'
    return
  endif 

  if not keyword_set( LSNR ) then lsnr = 3.5
  if not keyword_set( WV_PCA ) then wv_pca = alog10(1230.)
  if not keyword_set( ISTRT ) then istrt = 0L

  if not keyword_set( PCAFIL ) then pca_fil = $
          getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'

; Open pca
  pca = xmrdfits(pca_fil, 0, pca_head, /silent)
  sz = size(pca,/dimensions)
  if sz[0] NE 2000L then stop
  c0 = sxpar(pca_head, 'COEFF0')
  if c0 NE 3.0 then stop

  xmed = dindgen(sz[0])*0.0001 + 3.0d 
  ymed = pca[*,0]

  ;; Read SDSS files
  readcol, sdss_list, fil_sdss, format='A'

  ;; Parse out first line
  print, 'Note:  Ignoring first line!'
  fil_sdss = fil_sdss[1:*]

  ;; Create PCA arrays
  delta = fltarr(sz[0])
  d_ivar = fltarr(sz[0]) 

  if keyword_set( PSFIL ) then begin
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 5
      device, decompose=0
      ps_open, filename=psfil, /color, xs=10.8, ys=8.2
  endif ;else wset,2


  clr = getcolor(/load)
  if not keyword_set( NOPLT ) then !p.multi=[0,1,8]

; LOOP
  nfil = n_elements(fil_sdss)

  for qq=istrt,nfil-1 do begin

      ;; Get header
      sdss_head = xheadfits(fil_sdss[qq])
      sdss_name = strtrim( sxpar(sdss_head, 'NAME'), 2)+'-'+ $
        strtrim( sxpar(sdss_head, 'FIBERID'), 2)
      outfil = OUTDIR+sdss_name+'.fits'

      ;; REDO?
      if x_chkfil(outfil+'*', /silent) NE 0 AND not keyword_set( REDO ) then begin
          print, 'Skipping ', qq, outfil
          continue
      endif

      ;; Parse SDSS
      parse_sdss, fil_sdss[qq], flux, wave, ZQSO=zqso, IVAR=ivar, NPIX=npix, $
        HEAD=sdss_head, SIG=sig

      ;; 


      print, 'sdss_qsolin: qq ', qq, ' ', fil_sdss[qq]
      ;; rlam
      rlam = alog10( wave / (1.+zqso) )

      ;; Call qso_tilt
      pca_qsotilt, rlam, flux, ivar, xmed, ymed, tilt, tiltivar, $
        ANSWER=tilt_ans, FLG_TILT=flg_tilt

      ;; Find wave length edges for PCA analysis
      mn = min( abs(rlam - WV_PCA), ipstrt)
      mn = min( abs(rlam - xmed[sz[0]-1]), ipend)

      ;; Check for bad flux at beginning
      gdiv = where(ivar GT 0., ngdiv)
      if gdiv[0] GT ipend then ipend = 0L

      ;; Fit when 1584 is at lambda > 3800A
      if ipend GT 200 then begin

          ;; Find corresponding spots in xmed
          mn = min(abs(rlam[ipstrt]-xmed), ipc0)
          mn = min(abs(rlam[ipend]-xmed), ipc1)
          
          if (ipc1-ipc0) NE (ipend-ipstrt) then stop

          ;; Iterate
          for mm=0,2 do begin
          
              delta[*] = 0.
              delta[ipc0:ipc1] = tilt[ipstrt:ipend]/ymed[ipc0:ipc1] - 1.
              d_ivar[*] = 0.
              d_ivar[ipc0:ipc1] = tiltivar[ipstrt:ipend] * ymed[ipc0:ipc1]^2
              
              eigenspec = double(pca[*,1:*])
              
              ;; pca_weight
              pca_weight = transpose( eigenspec ) * $
                           (d_ivar ## replicate(1,20) )
              alpha = pca_weight # eigenspec
              beta = (pca_weight # delta)[*]
              
              ;; Linear algebra
              choldc, alpha, p
              result = cholsol(alpha, p, beta)
              pca_fit = eigenspec # result
              
              
              ;; PCA Continuum
              conti = (1.+pca_fit[ipc0:ipc1]) * $
                      ymed[ipc0:ipc1] * tilt_ans[ipstrt:ipend]
              
              ;; Smooth
              conti = smooth(conti,7)

              ;; Reject 
              bad = where((conti-flux[ipstrt:ipend]) GT 2.*sig[ipstrt:ipend],nb)
              if nb NE 0 then begin
                  nip = ipend-ipstrt+1
                  grow = [bad, [(bad+1) < (nip-1)], [(bad-1) > 0]]
                  ugrow = grow[uniq(grow,sort(grow))]
                  tiltivar[ipstrt+ugrow] = 0.
;                  x_splot, flux[ipstrt:ipend], ytwo=conti, /block, $
;                           ythr=tiltivar[ipstrt:ipend]
              endif
          endfor

          ;; Save
          fconti = tilt_ans
          fconti[ipstrt:ipend] = conti
          
          if keyword_set( CCHK) then x_splot, $
            flux[ipstrt:ipend], ytwo=conti, /block, ythr=tiltivar

          ;; Normalize
          norm = tilt[*,0]
          norm[ipstrt:ipend] = flux[ipstrt:ipend]/conti
          
          nrm_ivar = tiltivar[*,0]
          nrm_ivar[ipstrt:ipend] = ivar[ipstrt:ipend]*conti^2

      endif else begin  ; Low z QSO
          norm = flux[*]
          nrm_ivar = ivar[*]
          fconti = replicate(1.,  n_elements(flux))
      endelse

      ;; BSPLINE the red stuff

      cut = norm[ipstrt:*]
      cutwv = rlam[ipstrt:*]
      cutiv = nrm_ivar[ipstrt:*]
      tt = bspline_iterfit(cutwv, cut, invvar=cutiv, everyn=25, yfit=yfit, $
                           /groupbadpix, maxrej=20, /silent, /sticky, $
                           lower=2.0, upper=3., outmask=omsk, nord=3)
  
      ;; CHK
      fconti[ipstrt:*] = fconti[ipstrt:*]*yfit


      fx = cut/yfit
      ;; Sig
      sig = replicate(9.e9, n_elements(fx))
      gd = where(cutiv GT 0., ngd)
      sig[gd] = 1. / sqrt(cutiv[gd]*yfit^2)

      if keyword_set(CCHK) then $
        x_splot, wave, flux, ytwo=fconti, /bloc

      ;; Find lines
      x_findgauss, 1.-fx, 1./sig^2, xpeak=xpeak, ypeak=ypeak, $
        xsn=xsn, sn=sn, nfind=300L
      good = where(sn GT lsnr,ngood)
      if ngood NE 0 then begin
          gpk = xpeak[good]
          
          ;; Cut out bad sky regions (identify by too many lines at
        ;;; lambda>8000A)
          high = where(wave[ipstrt+gpk] GT 8000., complement=low, nhigh)
          if nhigh GT 10 then begin
;              gpk = gpk[low]
              flg_sky = 1
          endif else flg_sky = 0
          npk = n_elements(gpk)
          
          ;; Also 5575 and 6300
          if npk GT 0 then begin
              mn = min(abs(wave[ipstrt:*]-5579.), i55)
              mn = min(abs(wave[ipstrt:*]-6300.), i63)
              gd = where( abs(gpk-i55) GT 5 AND abs(gpk-i63) GT 2, npk)
              if npk NE 0 then gpk = gpk[gd] else gpk = 0L
          endif else gpk = 0L
      endif else gpk = 0L
          

;  x_splot, flux[l273:*,0], xtwo=xpeak[good], ytwo=fltarr(ngood), /block, $
;    psym2=1
;  x_splot, fx, xtwo=xpeak[good], ytwo=fltarr(ngood)+0.5, /block, $
;    psym2=1

;  x_splot, osnr[1000:*], ytwo=fsnr, /block

      if not keyword_set( NOPLT ) then begin
          plot, wave[ipstrt:*], [0.1, 1.5], /nodata, background=clr.white, $
            color=clr.black, xstyle=1, ystyle=1, charsize=1.3, ymargin=[3,0], $
            xmargin=[3,0], yticks=3
          oplot, wave[ipstrt:*], norm[ipstrt:*], color=clr.black
          oplot, wave[ipstrt:*], yfit, color=clr.red
          if ngood NE 0 then $
            oplot, [wave[ipstrt+gpk]], (norm[ipstrt+gpk]/2. > 0.2), $
            psym=1, color=clr.blue

          ;; Label
          xyouts, wave[ipstrt]+80., -0.1, sdss_name, charsize=1.5, color=clr.green
      endif
      
;      if qq EQ 12 then begin
;          !p.multi=[0,1,1]
;          stop
;      endif
      if not keyword_set( PSFIL ) and not keyword_set( NOPLT ) then begin
          if (qq+1) MOD 8 EQ 0 then stop
      endif

      ;; Output!
      if keyword_set( OUTDIR ) then begin
          outfil = OUTDIR+sdss_name+'.fits'
          mwrfits, fconti, outfil, /create
          ;; Just the unique ones!
          fwv = wave[ipstrt+gpk]
          if ngood NE 0 then $
            uwv = fwv[ uniq( fwv, sort(fwv) ) ] else uwv = [0.]
          if flg_sky EQ 1 then uwv = [uwv, -9999]
          mwrfits, uwv, outfil
          if not keyword_set(NOCOMPRESS) then spawn, 'gzip -f '+outfil
      endif

  endfor
  !p.multi=[0,1,1]
  if keyword_set( PSFIL ) then begin
      ps_close, /noprint, /noid
      device, decomposed=1
      !x.thick = 1
      !y.thick = 1
      !p.charthick = 1
      spawn, 'gzip -f '+psfil
  endif

  print, 'sdss_qsolin: All done!'
  return
end
