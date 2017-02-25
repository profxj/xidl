;+ 
; NAME:
; x_extapprof
;    Version 1.0
;
; PURPOSE:
;    Given the flux and wave images and the starting position of the
;       object to extract (x,y),  do a trace (if necessary) and then
;       calcualte the aperture and a Gaussian fit to the aperture.
;
; CALLING SEQUENCE:
;   
;   x_extapprof, fx, wv, xyguess
;
; INPUTS:
;   fx     -- Flux image
;   wv     -- Wavelength image
;  xyguess -- Guess at the trace
;
; RETURNS:
;
; OUTPUTS:
;   apstrct -- Structure containing the key info for the aperture and
;              profile
;
; OPTIONAL KEYWORDS:
;   WVMNX - Endpoints for extraction (default: [3200., 11000.])
;   COLLMNX - Endpoints for collapsing the spectrum prior to extraction
;              (default: [3400., 8000.])
;   VAR=  -- Variance array
;   FRAC=  -- Fraction of flux NOT to include in aperture [default:
;             0.025, which corresponds to 5% total]
;   TRC_ORD= -- Polynomial order for fitting the trace [default: 5]
;   /SKYSUB  -- Perform simple sky subtraction before calculating the
;               profile
;  PARAMETERS to TRACE_CRUDE:  NMED, NAVE, TRADIUS, MXSHIFT
;
; OPTIONAL OUTPUTS:
;  TOT_TRC -- Final trace (can be an input too)
;  FLG_APER= -- Flag describing the success (1= good)
;
; COMMENTS:
;
; EXAMPLES:
;   x_extapprof, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_extapprof, fx, wv, xyguess, MSK=msk, VAR=var, $
                 NMED=nmed, FRAC=frac, RADIUS=radius, $
                 APER=aper, TOT_TRC=tot_trc, NAVE=nave, $
                 COLLMNX=collmnx, DEBUG=debug, $
                 SKYSUB=skysub, MXSHIFT=MXSHIFT, TRADIUS=tradius,$
                 TRC_ORD=trc_ord, APSTRCT=apstrct,$
                 SIG_COLL=sig_coll, SILENT=silent, CHK=chk, $
                 FLG_APER = flg_aper, USEGAU = USEGAU, FIT_APER = FIT_APER

;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_extapprof, fx, wv, xyguess, MSK=, VAR=, '
    print, '       WVMNX=, NMED=, PIX=, FRAC=j, RADIUS=, APER=, '
    print, '       TOT_TRC=   [v1.0]'
    return
  endif 


;  Optional Keywords

  sz = size(fx, /dimensions)
  if not keyword_set(COLLMNX) then collmnx = [3400., 8000.]
  if not keyword_set(NMED) then nmed = 5L
  if not keyword_set(NAVE) then nave = 5L
  if not keyword_set(RADIUS) then radius = 10L
  if not keyword_set(FRAC) then frac = 0.025
  if not keyword_set(TRADIUS) then tradius = 2.
  if not keyword_set(TRC_ORD) then trc_ord = 5L
  if not keyword_set(SIG_COLL) then sig_coll = 2.5
  if not keyword_set( VAR ) then begin
      var = fltarr(sz[0],sz[1])
      var = (fx>0.) + 25.  ; Assumes readnoise = 25
  endif

;  Find all pixels in that slit
  if keyword_set( MSK ) then nwmsk = msk else nwmsk = bytarr(sz[0],sz[1])+1 

  flg_aper = 1L

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
          flg_aper = -1L
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
      
; RMB -> really bad traces ...
      if ntrc LT 40 then fitstr.nord = 3L
      if ntrc LT 10 then fitstr.nord = 2L

      trc = x_fitrej(float(gdtrc), xcen_pos[gdtrc], FITSTR=fitstr)
      tot_trc = x_calcfit(findgen(sz[0]), FITSTR=fitstr)
  endif

  if keyword_set(DEBUG) then begin
      tmpfx = fx
      tmpfx[lindgen(sz[0]),tot_trc[lindgen(sz[0])]] = -100
      xatv, tmpfx, /block, min=-50, max=50, wvimg=wv, sigimg=var
  endif

;;;;; APERTURE ;;;;;;;;;
; Collapse along the trace (restrict to 3200 - 8000A by default)
  rnd_trc = round(tot_trc)
  collpix = where(wv[lindgen(sz[0]),rnd_trc] GT collmnx[0] AND $
                  wv[lindgen(sz[0]),rnd_trc] LT collmnx[1] AND $
                  rnd_trc LT sz[1], ncoll)
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
  sv_ngcoll = ngd
      
  if keyword_set(DEBUG) then stop
; Find or set aperture 
  if not keyword_set( APER ) then begin
      if flg_smsh NE 1 then begin
          if keyword_set( BKAPER ) then aper = bkaper else aper = [5.,5.]
          print, 'x_extobjbox: No flux to set aperture for! Using backup!',$
            aper
      endif else begin
            aper = x_setaper(smsh, tot_trc[collpix[0]], frac, RADIUS = radius, $
                             /SKYSUB)
            IF KEYWORD_SET(fit_aper) THEN BEGIN
                aper[0] = aper[0] >  fit_aper[0]
                aper[1] = aper[1] >  fit_aper[1]
            ENDIF
            ;; JFH 05/07 force aper to be at least as large as that requested
            ;; for the fitting. 
        endelse
    endif

; Spline the aperture
  spln_area  = 0.
  splin = 0.
  gcoeff = fltarr(3)
  sv_cen = 0.
  ax = 0.
  ay= 0.
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
      ;; Gaussian
      ;stop
      ;pk = max(ay,ipk)
      pk = interpol(ay, ax, sv_cen)
      estim = [pk, sv_cen, 2.5]
      IF KEYWORD_SET(FIT_APER) THEN BEGIN 
          afitmin = round(sv_cen - fit_aper[0] - 1) > 0
          afitmax = round(sv_cen + fit_aper[1] + 1) < (n_elements(smsh)-1)
          indfit = WHERE(ax GE afitmin AND ax LE afitmax)
      ENDIF ELSE indfit = lindgen(n_elements(ax))
      gfit = gaussfit(ax[indfit], ay[indfit], gcoeff, estim = estim $
                      , nterms = 3) ;, sigma=gsig)
      Z = (ax-gcoeff[1])/gcoeff[2]  ;GET Z
      geval = gcoeff[0]*EXP(-Z^2/2.)          ;GAUSSIAN PART
      if gcoeff[2] LT 0.2 then stop
;      if gcoeff[2]/gsig[2] LT 5. then $
;        print, 'x_extapprof: Warning -- Gaussian not well defined!', gcoeff[2],$
;        gsig[2]

      ;; Plot

      if keyword_set( CHK ) then begin
          clr = getcolor(/load)
          plot, ax, ay, color=clr.black, background=clr.white, psym=10
          oplot, ax[indfit], gfit, color = clr.red
          IF n_elements(indfit) NE n_elements(ax) THEN $
            oplot, ax, geval, color = clr.red, linestyle = 2
          IF KEYWORD_SET(USEGAU) THEN BEGIN
              usefit = gauss1(ax, [sv_cen, usegau, 1.0d])
              usefit = usefit/max(usefit)*pk
              oplot, ax, usefit, color = clr.blue
          ENDIF
      endif
  endif

  if not keyword_set( GFIT ) then gfit = -1

  apstrct = { $
              aper: aper, $
              gcoeff: gcoeff, $
              flg_smsh: flg_smsh, $
              ngcoll: sv_ngcoll, $
              sv_cen: sv_cen, $
              spln_area: spln_area, $
              ax: ax, $
              ay: ay, $
              gfit: gfit, $
              splin: splin $
            }
  
  

;              gsig: gsig, $

              

  return
end
