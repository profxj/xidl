;+ 
; NAME:
; x_extobjbox   
;    Version 1.0
;
; PURPOSE:
;    Given the flux and wave images and the starting position of the
;       object to extract (x,y),  do a boxcar extraction
;
; CALLING SEQUENCE:
;   
;   x_extobjbox, slitstr, map
;
; INPUTS:
;   sub_fil      - Sky subtracted image file
;   slit_fil     - Slit file
;   obj_fil      - Object structure file
;   [fin_fil]    - Flux,sig,wave file
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;   WVMNX - Endpoints for extraction (default: [3200., 11000.])
;   COLLMNX - Endpoints for collapsing the spectrum prior to extraction
;              (default: [3400., 8000.])
;   REDBLUE - Spectrum runs from red to blue
;   SKYSUB  - Data has been sky subtracted (helps with aperture)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_extobjbox, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_extobjbox, fx, wv, xyguess, fin_spec, MSK=msk, VAR=var, WVMNX=wvmnx, $
                 NMED=nmed, PIX=pix, FRAC=frac, RADIUS=radius, $
                 APER=aper, TOT_TRC=tot_trc, CRVAL1=crval1, CDELT=cdelt, $
                 NPIX=npix, COLLMNX=collmnx, DEBUG=debug, REDBLUE=redblue, $
                 SKYSUB=skysub, REJ_CR=rej_cr, MXSHIFT=MXSHIFT, TRADIUS=tradius,$
                 NEWWV=newwv, TRC_ORD=trc_ord, NOCRREJ=nocrrej, REJSIG=rejsig,$
                 REBINC=rebinc, BKAPER=bkaper, SIG_COLL=sig_coll, SILENT=silent



;  Error catching
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_extobjbox, fx, wv, xyguess, fin_spec, MSK=, VAR=, '
    print, '       WVMNX=, NMED=, PIX=, FRAC=j, RADIUS=, APER=, BKAPER='
    print, '       TOT_TRC=, CRVAL1=, CDELT=, /REDBLUE, /CREBIN [v1.0]'
    return
  endif 


;  Optional Keywords

  sz = size(fx, /dimensions)
  if not keyword_set(WVMNX) then wvmnx = [3200., 11000.]
  if not keyword_set(COLLMNX) then collmnx = [3400., 8000.]
  if not keyword_set(NMED) then nmed = 5L
  if not keyword_set(NAVE) then nave = 5L
  if not keyword_set(RADIUS) then radius = 10L
  if not keyword_set(FRAC) then frac = 0.025
  if not keyword_set(REJSIG) then rejsig = 7.
  if not keyword_set(TRADIUS) then tradius = 3.
  if not keyword_set(TRC_ORD) then trc_ord = 5L
  if not keyword_set(SIG_COLL) then sig_coll = 2.5
  if not keyword_set( VAR ) then begin
      var = fltarr(sz[0],sz[1])
      var = (fx>0.) + 25.  ; Assumes readnoise = 25
  endif

;  Find all pixels in that slit
  if keyword_set( MSK ) then nwmsk = msk else nwmsk = bytarr(sz[0],sz[1])+1 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; TRACE ;;;;;;;;;
; Transpose for tracecrude

  if not keyword_set( TOT_TRC ) then begin
      nwmsk = transpose(nwmsk)
      tfx = transpose(fx)
      gdpix = where(nwmsk NE 0)

  ; Calculate subtracted ivar for the trace [dont use input var]
      ivar = fltarr(sz[1],sz[0]) 
      ivar[gdpix] = 1./(abs(tfx[gdpix]))

  ; Tracecrude
      newx = trace_fweight(tfx, xyguess[1], long(xyguess[0]), invvar=ivar, $
                           xerr=sigx)
      newx = trace_fweight(tfx, newx, long(xyguess[0]), invvar=ivar, $
                           xerr=sigx)
      xcen_pos = trace_crude(tfx, yset=ycen_pos, xstart=newx, nmed=nmed, $
                             nave=nave, radius=tradius, ystart=xyguess[0], $
                             xerr=xerr_pos, MAXSHIFTE=mxshift)

  ; Fit to the trace 
      gdtrc = where(xerr_pos LT 0.1, ntrc)

      if ntrc LE 5 then begin
          print, 'x_extobjbox: No object to extract!'
          fin_spec = { npix: 0 }
          return
      endif
      
      fitstr = { fitstrct }
      fitstr.func = 'POLY'
      fitstr.nord = TRC_ORD
      fitstr.hsig = 3.
      fitstr.lsig = 3.
      fitstr.maxrej = ntrc/5
      fitstr.niter = 2
      fitstr.minpt = 20
      
      if ntrc LT 40 then fitstr.nord = 3L

      trc = x_fitrej(float(gdtrc), xcen_pos[gdtrc], FITSTR=fitstr)
      tot_trc = x_calcfit(findgen(sz[0]), FITSTR=fitstr)
  endif 

  if keyword_set(DEBUG) then begin
      tmpfx = fx
      tmpfx[lindgen(sz[0]),tot_trc[lindgen(sz[0])]] = -100
      xatv, tmpfx, /block, min=-50, max=50, wvimg=wv, sigimg=var
      stop
  endif

;;;;; APERTURE ;;;;;;;;;
; Collapse along the trace (restrict to 3200 - 8000A by default)
  rnd_trc = round(tot_trc)
  collpix = where(wv[lindgen(sz[0]),rnd_trc] GT collmnx[0] AND $
                  wv[lindgen(sz[0]),rnd_trc] LT collmnx[1], ncoll)
  ; Pad with zeros before shifting (key for the edges)
  pad_img = fltarr(ncoll, sz[1]+200L)
  pad_var = fltarr(ncoll, sz[1]+200L)
  pad_img[*,100:sz[1]+99] = fx[collpix,*]
  pad_var[*,100:sz[1]+99] = var[collpix,*] > 0. ; No neg values
  ;; Create straightened image
  coll_img = fltarr(ncoll,sz[1])
  coll_var = fltarr(ncoll,sz[1])
  for i=0L,ncoll-1 do begin
      yshft = -round(tot_trc[collpix[i]]-tot_trc[collpix[0]])
      coll_img[i,*] = (shift(pad_img[i,*], 0, yshft))[*,100:sz[1]+99]
      coll_var[i,*] = (shift(pad_var[i,*], 0, yshft))[*,100:sz[1]+99]
  endfor
  ;; Total along columns
  ymn = (round(tot_trc[collpix[0]]) - radius) > 0L
  ymx = (round(tot_trc[collpix[0]]) + radius) < (sz[1]-1)
  sum_cimg = total(coll_img[*,ymn:ymx], 2)
  sum_cvar = total(coll_var[*,ymn:ymx], 2)
  a = where(sum_cvar GT 0.)
  ;; Calculate S/N
  sn_coll = fltarr(ncoll)
  sn_coll[a] = sum_cimg[a]/sqrt(sum_cvar[a])
  gd_coll = where(sn_coll GT sig_coll, ngd)

  if ngd LT ncoll/20 then flg_smsh = 0 else begin
      flg_smsh = 1
      smsh = djs_median(coll_img[gd_coll,*], 1)
  endelse
      
  if keyword_set(DEBUG) then stop
; Find or set aperture 
  if not keyword_set( APER ) then begin
      if flg_smsh NE 1 then begin
          if keyword_set( BKAPER ) then aper = bkaper else aper = [5.,5.]
          print, 'x_extobjbox: No flux to set aperture for! Using backup!',$
            aper
      endif else begin
          aper = x_setaper(smsh, tot_trc[collpix[0]], frac, RADIUS=radius, $
                           SKYSUB=skysub)
      endelse
  endif

; Spline the aperture
  if flg_smsh EQ 1 then begin
      sv_cen = tot_trc[collpix[0]]
      amin = round(sv_cen - aper[0] - 1) > 0
      amax = round(sv_cen + aper[1] + 1) < (n_elements(smsh)-1)
      ax = amin + findgen(amax-amin+1)
      ay = smsh[amin:amax]
      splin = spl_init(ax, ay, /double)
      ;; Find the area of the spline
      dx = (amax-amin)/500.
      tx = (amin + findgen(500)*dx)
      fspln = spl_interp(ax,ay,splin,tx,/double)
      spln_area = total( fspln*dx )
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;  EXTRACT ;;;;;;;

  fx_img = fx
  msk = bytarr(sz[0],sz[1])

  ; LOOP
  for i=0L,sz[0]-1 do begin
      ; Check wavelength
      a = where(wv[i,*] GT wvmnx[0] AND wv[i,*] LT wvmnx[1], nwv)
      if nwv EQ 0 then continue
      ; Find aperture edges
      lmin = (sz[1]-1) < (long(tot_trc[i]-aper[0]-0.5)+1) > 0L
      lmax = 0L > (long(tot_trc[i]+aper[1]+0.5) < (sz[1]-1)) 
      msk[i,lmin:lmax] = 1
      ; Deal with endpoints
      if (lmin-0.5) LT tot_trc[i]-aper[0] then $
        fx_img[i,lmin] = fx[i,lmin]*(1.-(tot_trc[i]-aper[0]-lmin+0.5))
      if lmax GT tot_trc[i]+aper[1] then $
        fx_img[i,lmax] = fx[i,lmax]*(tot_trc[i]+aper[1]+0.5-lmax)
      ; Find total 
      tot = total(fx_img[i,lmin:lmax])
      ; Ignore tot LT 0 case
      if tot GT 0. AND not keyword_set( NOCRREJ ) AND flg_smsh EQ 1 then begin
          xpix = lmin + findgen(lmax-lmin+1) - tot_trc[i] + sv_cen
          ; Get spline values
          spval = spl_interp(ax,ay,splin,xpix,/double)/spln_area
          ; Calc number of sigma [need to do correct variance]
          nsig = abs(fx_img[i,lmin:lmax]/tot - spval) / $
            (sqrt(var[i,lmin:lmax]) / tot)
          ; Reject (require >2 factor difference)
          rej = where(nsig GT rejsig AND $
                      abs(fx_img[i,lmin:lmax]/tot/spval) GT 2, nrej)
          if nrej NE 0 then var[i,lmin+rej] = -1.
;          if wv[i,lmin] GT 7800. and nrej NE 0 then begin
;              plot, fx_img[i,lmin:lmax]/tot
;              oplot, spval
;              print, i
;              stop
;          endif
      endif
  endfor
  gdpix = where(msk EQ 1) 


;;;;;;;; REBIN ;;;;;;;;

  ; Set wavelength scale (constant velocity (resolution) as default)
  if not keyword_set( NEWWV ) then begin
      if not keyword_set( CRVAL1 ) then $
        crval1 = double(alog10( wv[0,tot_trc[0]] < wv[sz[0]-1,tot_trc[sz[0]-1]] ))
      if not keyword_set( CDELT ) then $
        cdelt = (alog10( wv[0,tot_trc[0]] > wv[sz[0]-1,tot_trc[sz[0]-1]] ) - $
                 crval1)/double(sz[0])
      if not keyword_set( NPIX ) then npix = 2000L
      newwv = 10^( crval1 + dindgen(npix)*cdelt)
  endif else npix = n_elements(newwv)

  ; Check wavelength direction
  if not keyword_set( REDBLUE ) then begin
      if wv[sz[0]/2,tot_trc[sz[0]/2]] GT wv[sz[0]/2 + 1, tot_trc[sz[0]/2 + 1]] $
        then redblue = 1 else redblue = 0
  endif
  
  ; Rebin
  if not keyword_set( SILENT ) then $
    print, 'x_extobjbox: Rebinning...'
  x_rebin2dspec, wv, fx, newwv, newfx, GDPIX=gdpix, $
    VAR=var, NWVAR=newvar, SILENT=silent, REDBLUE=redblue, CR=rej_cr, $
    REBINC=rebinc

;;;;;; PASS IT BACK ;;;;;;;
  
  fin_spec = { $
               npix: npix, $
               wv: newwv, $
               fx: newfx, $
               var: newvar, $
               trc: tot_trc, $
               aper: aper $
             }

  return
end
