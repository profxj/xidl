;+ 
; NAME:
; fig_fittrcarc
;     Version 1.1
;
; PURPOSE:
;   To fit the slope of the arc lines as a function of order number
;   and y position on the CCD.  This information is then used to
;   construct a 2D wavelength image.  The fitting routine is the usual
;   least-squares algorithm with two rounds of rejection.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  Fits file with the coefficients of the 2D fit.  Filename like
;  'Arcs/TRC/Arc_mb0539_F.fits' 
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Plots residuals
;  /CLOBBER -- Overwrite previous solution
;  /ORDRCLOB -- Overwrite arc_m in the order structure
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
;   28-Apr-2003 Written by SB
;   Feb-2005 Ported to XIDL by JXP
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fig_fittrcarc, arc_fil, trc_fil, ordr_str, $
                      CHK=chk, CLOBBER=clobber, ORDR_FIL=ordr_fil, $
                      ORDRCLOB=ordrclob, NYCOEFF=nycoeff, NOCOEFF=nocoeff, $
                      _EXTRA=extra
;
;  if  N_params() LT 4  then begin 
;      print,'Syntax - ' + $
;        'rslt = x_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, [qafil], ' + $
;        '/CHK, /CLOBBER /ORDRCLOB, NYCOEFF=, NOCOEFF=) [v1.1]'
;      return
;  endif 

;  Optional Keywords

  if keyword_set( CLOBBER ) then ordrclob = 1
  if NOT keyword_set(nycoeff) then nycoeff = 2
  if NOT keyword_set(nocoeff) then nocoeff = 3

  if not keyword_set(PSFIL) then psfil='fig_fittrcarc.ps'
  if not keyword_set(ARC_FIL) then arc_fil = $
    getenv('MIKE_PAP')+'Arcs/Arc_mb0005.fits'
  if not keyword_set(TRC_FIL) then trc_fil = $
    getenv('MIKE_PAP')+'Arcs/TRC/Arc_mb0005_T.fits'
  if not keyword_set(OSTRFIL) then ostrfil = $
    getenv('MIKE_PAP')+'Flats/OStr_B_01.fits'
  if not keyword_set(ARC_INFO) then arc_info = $
    getenv('MIKE_PAP')+'Arcs/Fits/mb0005_fit.idl'
  lsz = 1.3
  csize = 2.1

  
  ;; 
  ordr_str = xmrdfits(ostrfil, 1, /silent)
  nordr = n_elements(ordr_str)

  ;; Grab Arc Trace
  trc_arc = xmrdfits(trc_fil, 1, /silent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SETUP THE DATA

  print, 'x_fittrcarc: Setting up the values (normalizing)'
  npix = round(total(trc_arc.ngood))
  t = dblarr(npix)
  ;; t, PIX and slope
  cnt = 0L
  for j=0L,nordr-1 do begin
      ngd = trc_arc[j].ngood
      ;; PIX
      if ngd GT 0 then begin
        if keyword_set(all_pix) EQ 0 then $
             all_pix = [trc_arc[j].xguess[0:ngd-1]] $
        else all_pix = [all_pix,trc_arc[j].xguess[0:ngd-1]]
        ;; SLOPE
        if keyword_set(all_slope) EQ 0 then $
          all_slope = [(2.*trc_arc[j].coeff[1,0:ngd-1]/ $
                      (trc_arc[j].xmax - trc_arc[j].xmin))[*]] $
        else $
          all_slope = [all_slope, $
                     (2.*trc_arc[j].coeff[1,0:ngd-1]/ $
                      (trc_arc[j].xmax - trc_arc[j].xmin))[*]]
        ;; Order #
        t[cnt:cnt+ngd-1] = ordr_str[j].order
        cnt = cnt + ngd
     endif
  endfor
          
  nrm = dblarr(2)
  ;; NORMALIZE PIX
  mnx = min(all_pix, MAX=mxx)
  nrm[0] = 0.5 * (mnx + mxx)
  nrm[1] = mxx - mnx
  pix_nrm = 2. * (all_pix - nrm[0])/nrm[1]
  
  ;; NORMALIZE ORDER
  nrmt = dblarr(2)
  mnx = min(t, MAX=mxx)
  nrmt[0] = 0.5 * (mnx + mxx)
  nrmt[1] = mxx - mnx
  t_nrm = 2. * (t - nrmt[0])/nrmt[1]
  
  invvar = replicate(1., npix)
  
  ;;  Setup the Functions
  work2d = dblarr(npix,nycoeff*nocoeff)
  worky = flegendre(pix_nrm[*], nycoeff)
  workt = flegendre(t_nrm[*], nocoeff)
  
  for i=0,nocoeff-1 do begin
      for j=0,nycoeff-1 do begin
          work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
      endfor
  endfor
  
  ;; Do the matrix algebra
  work2di = transpose(work2d * sqrt(invvar[*] # replicate(1,nocoeff*nycoeff)))
  alpha = work2di # transpose(work2di)
  beta = work2di # (all_slope * sqrt(invvar[*]))
  choldc, alpha, p, /double
  res = cholsol(alpha,p,beta, /double)
  slop_mod = dblarr(npix)
  slop_mod[*] = work2d # res
  
  ;; Get RMS
  gd_wv = where(invvar GT 0.0, ngd)
  ;; REJECT
  diff = (slop_mod-all_slope)
  djs_iterstat, diff[gd_wv], sigrej=4.0, sigma=rms
  msk = (abs(diff) LT 4.0*rms)*(invvar GT 0.0)
  gd = where(msk EQ 1B, complement=bad)
  print, 'x_fittrcarc:  RMS1 = ', rms
  ;; RESET invvar
  if bad[0] NE -1 then invvar[bad] = 0.
          
  ;; Do the matrix algebra
  work2di = transpose(work2d * sqrt(invvar[*] # replicate(1,nocoeff*nycoeff)))
  alpha = work2di # transpose(work2di)
  beta = work2di # (all_slope * sqrt(invvar[*]))
  choldc, alpha, p, /double
  res = cholsol(alpha,p,beta, /double)
  slop_mod = dblarr(npix)
  slop_mod[*] = work2d # res
      
  ;; MSK
  gd_wv = where(invvar GT 0.0, ngd)
  diff = (slop_mod-all_slope)
  djs_iterstat, diff[gd_wv], sigrej=2.5, sigma=rms

  msk = (abs(diff) LT 2.5*rms) * (invvar GT 0.0)
  gd = where(msk EQ 1B, complement=bad)
  print, 'x_fittrcarc:  RMS2 = ', rms
          
  ;; One more time
  
  ;; RESET invvar
  if bad[0] NE -1 then invvar[bad] = 0.
  
  ;; Do the matrix algebra
  work2di = transpose(work2d * sqrt(invvar[*] # $
                                replicate(1,nocoeff*nycoeff)))
  alpha = work2di # transpose(work2di)
  beta = work2di # (all_slope * sqrt(invvar[*]))
  choldc, alpha, p, /double
  res = cholsol(alpha,p,beta, /double)
  slop_mod = dblarr(npix)
  slop_mod[*] = work2d # res
  
  gd_wv = where(invvar GT 0.0, ngd)
  diff = (slop_mod-all_slope)
  djs_iterstat, diff[gd_wv], sigrej=2.5, sigma=rms
  msk = (abs(diff) LT 2.5*rms) * (invvar GT 0.0)
  gd = where(msk EQ 1B, complement=bad)

  print, 'x_fittrcarc:  RMS3 = ', rms
  if keyword_set( CHK ) then $
    x_splot, t, diff, psym1=3, xthr=[min(t),max(t)], ythr=[0,0], psym3=-1, $
             xtwo=t[bad], ytwo=diff[bad],psym2=4, /block
;             xtitle='Order Number', ytitle='Slope errors'
          
          
  ;; Calculate slope at all y values for each order
  ny = n_elements(ordr_str[0].lhedg)
  pix_nrm = 2. * (dindgen(ny) - nrm[0])/nrm[1] # replicate(1., nordr)
  t_nrm = 2.d * (ordr_str.order - nrmt[0])/nrmt[1] ## replicate(1., ny)
  
;      invvar = replicate(1., ny*nordr)
  nwork2d = dblarr(ny*nordr,nycoeff*nocoeff)
  worky = flegendre(pix_nrm[*], nycoeff)
  workt = flegendre(t_nrm[*], nocoeff)
  for i=0,nocoeff-1 do begin
      for j=0,nycoeff-1 do begin
          nwork2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
      endfor
  endfor
  nslop_mod = dblarr(ny,nordr)
  nslop_mod[*] = nwork2d # res

  ;; Fill up the order structure (if not done before)
  if ordr_str[0].arc_m[0] EQ 0. OR keyword_set( ORDRCLOB ) then begin
      print, 'x_fittrcarc:  Writing arc_m in order structure..'
      ordr_str.arc_m = nslop_mod
      flg_owrite = 1
  endif else flg_owrite = 0

  ;; Output
;  print, 'x_fittrcarc:  Saving fit info to ', out_fil
  slope_str = { $
                res: res, $
                nrm: nrm, $
                nrmt: nrmt, $
                nycoeff: nycoeff, $
                nocoeff: nocoeff $
              }
;  mwrfits, slope_str, out_fil, /create
  if flg_owrite EQ 1 then begin
      print, 'x_fittrcarc:  Overwriting order structure..'
      if not keyword_set(ORDR_FIL) then stop
      mwrfits, ordr_str, ordr_fil, /create
  endif

  ;; QA
  head = xheadfits(arc_fil, /silent)
  sz = lonarr(2)
  sz[0] = sxpar(head,'NAXIS1')
  sz[1] = sxpar(head,'NAXIS2')

  x_psopen, psfil, /maxs
  clr = getcolor(/load)
;  !p.multi=[0,3,2]
      
  ;; NORMALIZE PIX
  plt_pix = dindgen(round(sz[1]/8.))*8.
  npix = n_elements(plt_pix)
  pix_nrm = 2. * (plt_pix - nrm[0])/nrm[1]
  worky = flegendre(pix_nrm[*], nycoeff)
  
  ;; Main plot
  mn = min(slop_mod*max(trc_arc.xmax), max=mx) < 0.
;  plot, [0.], [0.], color=clr.black, thick=5, $
;        background=clr.white, charsize=1.5, yrange=[mn, mx], $
;        xrange=[0.,sz[1]], xstyle=1, ystyle=1, xtitle='Row', $
;        ytitle='Arc Line Tilt', xmargin=[11,2], ymargin=[5,1], /nodata
;  
;  xyouts, 0.5, 0.96, 'RMS = '+string(rms,format='(f7.5)'),$
;          /normal, alignment=0.5, color=clr.black, charsize=2.5
;  
;  !p.thick = 1
  nordr = n_elements(ordr_str)
  for jj=0L,nordr-1 do begin
      ;; NORMALIZE ORDER
      ii = ordr_str[jj].order
      tsub = replicate(float(ii), npix)
      t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
          
      ;; work2d and wv
      work2d = dblarr(npix,nycoeff*nocoeff)
      workt = flegendre(t_nrm[*], nocoeff)
      
      for i=0,nocoeff-1 do begin
          for j=0,nycoeff-1 do begin
              work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
          endfor
      endfor
      
      ;; Model
      slope = dblarr(npix)
      slope[*] = work2d # slope_str.res
;      oplot, plt_pix, slope*trc_arc[jj].xmax, color=clr.black
      
      ;; Resid
      pts = where(t EQ ii, npts)
      if npts NE 0 then begin
          sres = (slop_mod[pts] - all_slope[pts])*trc_arc[jj].xmax
;          oplot, all_pix[pts],  slop_mod[pts]*trc_arc[jj].xmax - sres, $
;                 psym=1, color=clr.blue, symsize=0.5
      endif
  endfor
      
  ;;;;;;;;;
  ;; Individual plots
  !p.multi=[0,3,2]
  for jj=12L,18-1 do begin
      ;; NORMALIZE ORDER
      ii = ordr_str[jj].order
      tsub = replicate(float(ii), npix)
      t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
      
      ;; work2d and wv
      work2d = dblarr(npix,nycoeff*nocoeff)
      workt = flegendre(t_nrm[*], nocoeff)
      
      for i=0,nocoeff-1 do begin
          for j=0,nycoeff-1 do begin
              work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
          endfor
      endfor
      
      ;; Model
      slope = dblarr(npix)
      slope[*] = work2d # slope_str.res
      plot, plt_pix, slope*trc_arc[jj].xmax, color=clr.black, $
            background=clr.white, charsize=csize, yrange=[-2.4, -1.1], $
            xrange=[0.,sz[1]], xstyle=1, ystyle=1, thick=1, xtitle='CCD Row', $
            ytitle='Arc Line Tilt (pixels per half-slit)', $
            xmargin=[8,1], ymargin=[4,0.5]
          
      ;; Resid
      pts = where(t EQ ii, npts)
      if npts NE 0 then begin
          sres = (slop_mod[pts] - all_slope[pts])*trc_arc[jj].xmax
          oplot, all_pix[pts],  slop_mod[pts]*trc_arc[jj].xmax - sres, $
                 psym=1, color=clr.blue, symsize=0.5
          ;; RMS
          rms = sqrt( total( sres^2 ) / float(npts))
          
          ;; Rejected
          bad = where(invvar[pts] LE 0., nbad)
          if nbad NE 0 then begin
              oplot, [all_pix[pts[bad]]],  $
                     [slop_mod[pts[bad]]*trc_arc[jj].xmax - sres[bad]], $
                     psym=2, color=clr.red, symsize=0.5
          endif
      endif

      ;; Label
      xyouts, 100, -1.25, 'Order '+strtrim(ordr_str[jj].order,2), $
              color=clr.black, charsiz=1.6
      
  endfor
  
  x_psclose
  !p.multi=[0,1,1]
      
  return
  end
