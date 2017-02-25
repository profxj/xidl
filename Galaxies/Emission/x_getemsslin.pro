;+ 
; NAME:
; x_getemsslin
;   Version 1.0
;
; PURPOSE:
;    Measure Galaxy emission lines given an emission line structure
;
; CALLING SEQUENCE:
;   x_getemsslin, infil, outfil, emstrct
;
; INPUTS:
;    INFIL -- 
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL -- File of emission line fluxes and EW values
;
; OPTIONAL KEYWORDS:
;  LINLIST= -- Line list of emission line features [default:
;              gal_vac.lst]
;  FNORM=   -- Normalization value for the data [default: 1e-17]
;  HELIO=   -- Perform a heliocentric correction [Requires input
;              values]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_getemsslin, infil, outfil, emstrct, LINLIST=linlist, SETEW=setew, $
                  PLOT=plot, ZSYS=zsys, HELIO=helio, CHK=chk, $
                  DATFLG=datflg, FNORM=fnorm

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_getemsslin, infil, outfil, LINLIST= (v1.0)'
    return
  endif 

;  Optional Keywords
  if not keyword_set( FNORM ) then fnorm = 1e17
  if not keyword_set( DATFLG ) then datflg = 0L
  if not keyword_set( ZSYS ) then zsys = 1.
  if not keyword_set( LINLIST ) then linlist = $
    getenv('XIDL_DIR')+'/Spec/Lines/Lists/gal_vac.lst'

  ;; Read in the line list
  lines = x_setllst(linlist, 1)


  ;; Read input file
  close, /all
  openr, 1, infil
  datfil = ''
  readf, 1, datfil

  ;; Open data file
  case datflg of
      0: begin
          fx = x_readspec(datfil, wav=wv, sig=sig, /autofsig)
          dat = { $
                  fx: fx, $
                  wv: wv, $
                  sig: sig }
      end
      1: dat = xmrdfits(datfil, 1, /silent)
      else: stop
  endcase

  ;; Helio correction
  if keyword_set( HELIO ) then begin
      helio = x_keckhelio(helio[0], helio[1], helio[2], jd=helio[3], $
                          altitude=2282., $
                          longitude=360.-70.70, $
                          latitude=-29.08333)
      print, 'x_getemsslin: Helio correction of ', helio, ' km/s'
      hel_corr = sqrt( (1.d + helio/299792.458d) / (1.d - helio/299792.458d) )
      dat.wv = dat.wv*hel_corr
  endif
      
  if keyword_set( PLOT ) then begin
      x_specplot, dat.fx, dat.sig, wave=dat.wv, inflg=4, /block, /gal
  endif

  ;; Create emiss structure
  nlin = 0L
  readf, 1, nlin

  tmp = {emissstrct}
  emstrct = replicate(tmp, nlin)

  ;; Continuum
  cnti = dblarr(2,2)
  wvlim = dblarr(2)
  dumt = dblarr(2)
  contistr = x_setfitstrct(flgrej=1, nord=1)

  ;; Loop on emission lines
  for qq=0L, nlin-1 do begin
      readf, 1, wval, format='(f)'
      mn = min(abs(wval-lines.wave), imn)
      ;; 
      emstrct[qq].ion = lines[imn].name
      emstrct[qq].wrest = lines[imn].wave
      ;; Input
      readf, 1, dumt 
      cnti[0,*] = dumt
      readf, 1, dumt 
      wvlim[*] = dumt
      readf, 1, dumt 
      cnti[1,*] = dumt

      ;; Measure continuum
      ;; Grab pixels
      cpix = where( (dat.wv GT cnti[0,0] AND dat.wv LT cnti[0,1]) OR $
                    (dat.wv GT cnti[1,0] AND dat.wv LT cnti[1,1]), npix)
      if npix EQ 0 then stop
      fit = x_fitrej(dat.wv[cpix], dat.fx[cpix], FITSTR=contistr)

      ;; Measure the line
      lpix = where( dat.wv GT wvlim[0] AND dat.wv LT wvlim[1], nem)
      if nem EQ 0 then stop

      ;; Total
      dwv = (shift(dat.wv[lpix],-1)-shift(dat.wv[lpix], 1))/2.
      dwv[0] = dat.wv[lpix[1]]-dat.wv[lpix[0]]
      dwv[nem-1] = dat.wv[lpix[nem-1]]-dat.wv[lpix[nem-2]]
      ;; Calc conti
      cval = x_calcfit(dat.wv[lpix], fitstr=contistr)

      ;; Flux
      emstrct[qq].flux = total( (dat.fx[lpix]-cval)*dwv ) / (1.+zsys)
      emstrct[qq].fsig = sqrt( total( (dat.sig[lpix]*dwv)^2 ) ) / (1.+zsys)
      ;; EW
      emstrct[qq].ew = total( (dat.fx[lpix]/cval)*dwv ) / (1.+zsys)
      emstrct[qq].ewsig = sqrt( total( (dat.sig[lpix]*dwv/cval)^2 ) ) / (1.+zsys)
      ;; Wave centroid
      wcen = total(dat.fx[lpix]*dat.wv[lpix])/total(dat.fx[lpix])

      ;; Gaussian
      yfit = gaussfit(dat.wv[lpix], dat.fx[lpix]*FNORM, acoeff, $
                      estimates=[max(dat.fx[lpix])*FNORM, wcen, $
                                 1., median(cval)*FNORM], $
                      sigma=sigma, nterms=4)
      emstrct[qq].wcen = acoeff[1]
      emstrct[qq].sig_wc = sigma[1]

      emstrct[qq].vsigma = acoeff[2]/acoeff[1]*3e5
      emstrct[qq].sig_vsig = sigma[2]/acoeff[1]*3e5

      ;; Wavelength
      if keyword_set( CHK ) then $
        x_splot, dat.wv[lpix], dat.fx[lpix]*FNORM, psym1=1, ytwo=yfit, /block

      ;; Redshift
      emstrct[qq].zem = emstrct[qq].wcen/emstrct[qq].wrest - 1.
      printcol, emstrct[qq].wrest, emstrct[qq].zem, acoeff[2]/acoeff[1]*3e5
      
      if keyword_set( PLOT ) then begin
          ymin = min([dat.fx[cpix],dat.fx[lpix]], max=ymax)
          clr = getcolor(/load)
          plot, dat.wv[cpix], dat.fx[cpix], /nodata, xstyle=1, ystyle=1, $
            background=clr.white, color=clr.black, xtitle='!17Wavelength', $
            ytitle='!17Flux', xmargin=[7,1], ymargin=[3,1], $
            yrange=[-ymax*0.1,ymax*1.1]
          oplot, dat.wv[cpix], dat.fx[cpix], color=clr.green, psym=1
          oplot, dat.wv[lpix], dat.fx[lpix], color=clr.blue, psym=10
          oplot, dat.wv[lpix], cval, color=clr.red, psym=1
          xyouts, 0.2, 0.1, emstrct[qq].ion, color=clr.black, charsize=1.5, $
            /normal
          stop
      endif

      ;; Gaussian fluxes and EW
      fgauss = acoeff[0] * acoeff[2] * sqrt(2.*!pi)
      
      
  endfor

  print, 'x_getemsslin:  All done'
  mwrfits, emstrct, outfil, /create

  return

end
