;+ 
; NAME:
; x_echfitstd   
;   Version 1.3
;
; PURPOSE:
;    GUI for tweaking a sensitivity function for an Echelle spectrum
;
; CALLING SEQUENCE:
;   
;   fit = x_echfitstd([xdat],ydat,func=,nord=, /inter, xsize=, ysize=)
;
; INPUTS:
;   xdat       - Values along one dimension [optional]
;   ydat       - Values along the other
;
; RETURNS:
;   fit        - Values at each xdat
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   func       - String for Fitting function (POLY, LEGEND, BSPLINE,
;                GAUSS, CHEBY)
;   nord       - Order of the fit
;   /inter      - Interactive fitting
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   sig        - Errors in the points
;   reg        - Regions of data to fit
;   LSIG       - Lower SIGMA
;   HSIG       - High SIGMA
;   /REJ        - Turn rejection on
;   DELPTS     - Array of user-deleted points
;   MSK        - Array of 0,1 values [0=do not include]
;
; OPTIONAL OUTPUTS:
;   bset       - Bspline info
;   rms        - RMS of the fit (only valid for fits with rejection)
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_echfitstd(x, y, 'POLY', nord=5, /inter)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  POLY_FIT
;  POLY
;  SVDFIT
;  SVLEG
;  BSPLINE
;  GAUSS
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;  GETCOLOR (coyote package)
;
; REVISION HISTORY:
;   22-June-2001 Written by JXP
;   28-June-2001 Added LEGENDRE fits (JXP)
;   23-Nov-2001 Added BSPLIN (SB), GAUSS, REJECTION (JXP)
;   31-Jan-2002 Added CHEBYSHEV
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_echfitstd_initcommon, inord, fitstr

;

common x_echfitstd_fit, nord, tot_fit, sens_func

  nord = inord
  if not keyword_set( FITSTR ) then $
    fitstr = x_setfitstrct(FUNC='LEGEND', NORD=9L, LSIG=3., HSIG=3., $
                           FLGREJ=1)
  
  tot_fit = replicate(fitstr, nord)

end

;;;;
; Events
;;;;

;function x_echfitstd_ev, ev
pro x_echfitstd_event, ev

common x_echfitstd_fit 

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'FUNCLIST' : begin
          nfunc = ev.index
          case nfunc of 
              0 : fitprm.func = 'POLY'
              1 : fitprm.func = 'LEGEND'
              2 : fitprm.func = 'CHEBY'
              3 : fitprm.func = 'BSPLIN'
              4 : fitprm.func = 'GAUSS'
              else:
          endcase
          ; Reset rejected
          rej = where(state.gdpix[0:state.norg-1] MOD 8 GE 4, nrej)
          if nrej NE 0 then state.gdpix[rej] = state.gdpix[rej] - 4
          ; Update fit
          x_echfitstd_UpdateFit, state
      end
      'UP' : begin
          fitprm.nord = fitprm.nord+1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm.nord, format='(i4)')
      end
      'DOWN' : begin
          fitprm.nord = fitprm.nord-1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm.nord, format='(i4)')
      end
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Delete a region
                          x_echfitstd_zeroreg, state
                          if state.flg_zero EQ 1 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state, /no_copy
                              return
                          endif 
                      end 
                      4 : x_echfitstd_setordr, state
                      else :
                  endcase
              end
              1 : begin ; Button release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          if (state.flg_zero EQ 1 AND eventch NE 's') then begin
              print, 'x_echfitstd: Expecting another s or LMB click !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
           endif
           if (state.flg_unzero EQ 1 AND eventch NE 'x') then begin
              print, 'x_echfitstd: Expecting another x'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
           endif
          case eventch of
              'u': begin
                  if state.currordr LT 1 then return
                  gdo = where(state.objstr.order EQ state.currordr)
                  tot_fit[gdo].nord = tot_fit[gdo].nord + 1
                  print, 'x_echfitstd: Fit order = ', tot_fit[gdo].nord
;                  widget_control, state.lblordr_id, $
;                    set_value=string(fitprm.nord, format='(i4)')
                  x_echfitstd_UpdateFit, state, ordr=state.currordr
              end
              'd': begin
                  if state.currordr LT 1 then return
                  gdo = where(state.objstr.order EQ state.currordr)
                  tot_fit[gdo].nord = (tot_fit[gdo].nord - 1) > 0
                  print, 'x_echfitstd: Fit order = ', tot_fit[gdo].nord
;                  widget_control, state.lblordr_id, $
;                    set_value=string(fitprm.nord, format='(i4)')
                  x_echfitstd_UpdateFit, state, ordr=state.currordr
              end
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
;              ']': x_specpan, state
;              '[': x_specpan, state, /left
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
;              'z': begin  ; Zoom
;                  ximgd_setzoom, state, /plot
;                  if state.flg_zero EQ 1 then begin
;                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
;                      return
;                  endif
;              end
              's': begin  ; Region
                  x_echfitstd_zeroreg, state
                  if state.flg_zero EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, $
                        set_uvalue = state, /no_copy
                      return
                  endif 
               end
              'x': begin        ; Region
                  x_echfitstd_unzeroreg, state
                  if state.flg_unzero EQ 1 then begin
                     WIDGET_CONTROL, state.base_id, $
                        set_uvalue = state, /no_copy
                      return
                  endif 
               end
              'f': begin        ; Region
                  x_echfitstd_forcereg, state
                  if state.flg_force EQ 1 then begin
                     WIDGET_CONTROL, state.base_id, $
                        set_uvalue = state, /no_copy
                      return
                  endif 
               end
              'k': state.res1ymnx = state.res1ymnx/2.
              'K': state.res1ymnx = state.res1ymnx*2.
              'F': begin ;; Reset y values on Flux window
                  gd = where(state.std_wv GT state.xymnx[0] AND $
                             state.std_wv LT state.xymnx[2], ngd)
                  if ngd NE 0 then begin
                      state.stdymnx = [0.9*min(state.std_fx[gd], max=mx), $
                                       mx*1.1] 
                  endif
              end
              'W': state.xymnx = state.svxymnx   ; Reset the screen
;              'D': x_echfitstd_Delete, state ; Delete/Undelete a data point
;              'x': x_echfitstd_DelReg, state        ; Delete one region
;              'S': x_echfitstd_DelReg, state, /all  ; Delete all regions and reset
              'R': begin  ;; Reset sigma of current order
                  if state.currordr LT 1 then return
                  gdo = where(state.objstr.order EQ state.currordr)
                  state.objstr[gdo].sig = (state.svsig)[*,gdo]
                  x_echfitstd_UpdateFit, state, ordr=state.currordr
              end
;              'V': x_echfitstd_SetValues, state     ; Plot values
;              'C': x_echfitstd_ClearRej, state ; Clear all rejected points
              'H': x_helpwidg, state.help
              'q': begin
;                  x_echfitstd_setpnt, state
                  widget_control, ev.top, /destroy
                  return
              end
              else:  ; Nothing
          endcase
      end
;            REJECTION
      'HSIGVAL': begin
          fitprm.hsig = ev.value
          if fitprm.niter EQ 0 then fitprm.niter=3
          x_echfitstd_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.lsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'LSIGVAL': begin
          fitprm.lsig = ev.value
          if fitprm.niter EQ 0 then fitprm.niter=3
          x_echfitstd_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.hsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'DONE' : begin
          x_echfitstd_setpnt, state
          widget_control, ev.top, /destroy
          return
      end
      else:
  endcase

; Update Plot
  x_echfitstd_UpdatePlot, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_echfitstd_UpdatePlot, state, ALL=all
  
common x_echfitstd_fit


  clr = getcolor(/load)

  ;; Plot Data
  widget_control, state.draw_id, get_value=wind
  wset, wind

  gdw = where( state.objstr.box_wv GT state.xymnx[0] AND $
               state.objstr.box_wv LT state.xymnx[2], ngdw )
  if ngdw LE 1 then return

  ;; Plot good and bad data
  gd1 = where( (state.objstr.sig)[gdw] GT 0. $
               AND ((state.ordnum)[gdw] MOD 2 EQ 0), ngd1)
  IF ngd1 NE 0 THEN gd1 = gdw[gd1] ; added by JFH to fix bug
  gd2 = where( (state.objstr.box_wv)[gdw] GT 0. $
               and (state.objstr.sig)[gdw] GT 0. $
               AND ((state.ordnum)[gdw] MOD 2 EQ 1), ngd2)
  IF ngd2 NE 0 THEN gd2 = gdw[gd2] ; added by JFH to fix bug

  if state.currordr GT 0 then begin
      gdo = where( (state.objstr.box_wv)[gdw] GT 0. $
                   and (state.objstr.sig)[gdw] GT 0. $
                   AND ((state.ordnum)[gdw] EQ state.currordr) , ngdw)
      if ngdw NE 0 then gdo = gdw[gdo] else gdo = [0L,0L]
  endif
  bad = where((state.objstr.sig)[gdw] LE 0., nbad)
  force = where((state.objstr.sig)[gdw] GT 1.0d8, nforce)
  IF nbad NE 0 THEN bad = gdw[bad]
  IF nforce NE 0 THEN force = gdw[force]
  ;; Data
  plot, [0.], [0.], $
     xrange=[state.xymnx[0],state.xymnx[2]], $
     yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
     charsize=1.2, psym=1, background=clr.white, color=clr.black, $
     xtitle='!17Wavelength', xmargin=[0,0], ymargin=[0,0], xstyle=1,$
     ystyle=1, /nodata

 ;;  plot, (state.objstr.box_wv)[gd1], (state.objstr.flux)[gd1], $
;;     xrange=[state.xymnx[0],state.xymnx[2]], $
;;     yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
;;     charsize=1.2, psym=1, background=clr.white, color=clr.black, $
;;     xtitle='!17Wavelength', xmargin=[0,0], ymargin=[0,0], xstyle=1,$
;;     ystyle=1, /nodata

  ;; Plot the good
  IF ngd1 GT 0 THEN $
     oplot, (state.objstr.box_wv)[gd1], (state.objstr.flux)[gd1], psym = 1 $
            , color = clr.red
  
  IF ngd2 GT 0 THEN $
     oplot, (state.objstr.box_wv)[gd2], (state.objstr.flux)[gd2], psym = 1 $
            , color = clr.blue

  IF state.currordr GT 0 then BEGIN
      oplot, (state.objstr.box_wv)[gdo], (state.objstr.flux)[gdo], $
             psym = 1, color = clr.black
      oplot, (state.objstr.box_wv)[gdo], (state.fit)[gdo], $
             color = clr.magenta, thick = 3.0
  END

  ;; Plot the bad
  IF nbad GT 0 THEN $
     oplot, (state.objstr.box_wv)[bad], (state.objstr.flux)[bad], psym = 1 $
            , color = clr.green
  ;; Plot the force
  IF nforce GT 0 THEN $
     oplot, (state.objstr.box_wv)[force], (state.objstr.flux)[force], psym = 1 $
            , color = clr.yellow

  
  ;; Label
  if state.currordr NE 0 then $
    xyouts, 0.82, 0.9, 'Order# '+string(state.currordr,format='(i3)'), $
    /normal, color=clr.black, charsize=2.

  ;; Plot the residuals
  widget_control, state.resid1_id, get_value=wind
  wset, wind

  plot, [0.], [0.], $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=state.res1ymnx, $
    charsize=1.2, background=clr.white, color=clr.black, $
    xmargin=[8,0], ymargin=[0,0], xstyle=1, ystyle=1, /nodata
  IF ngd1 GT 0 THEN oplot, (state.objstr.box_wv)[gd1], (state.resid1)[gd1] $
                           , psym = 1, color = clr.red
  IF ngd2 GT 0 THEN oplot, (state.objstr.box_wv)[gd2], (state.resid1)[gd2] $
                           , psym = 1, color = clr.blue
  if state.currordr GT 0 then $
     oplot, (state.objstr.box_wv)[gdo], (state.resid1)[gdo] $
            , psym = 1, color = clr.black
  
  ;; Plot the standard
  widget_control, state.stdstr_id, get_value=wind
  wset, wind

  plot, [0.], [0.], xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=state.stdymnx, $
    charsize=1.2, background=clr.white, color=clr.black, $
    xmargin=[8,0], ymargin=[0,0], xstyle=1, ystyle=1, /nodata
  IF ngd1 GT 0 THEN oplot, (state.objstr.box_wv)[gd1], (state.newflux)[gd1] $
                           , psym = 1, color = clr.red
  IF ngd2 GT 0 THEN oplot, (state.objstr.box_wv)[gd2], (state.newflux)[gd2] $
                           , psym = 1, color = clr.blue
;  if state.currordr GT 0 then $
;    oplot, (state.objstr.box_wv)[gdo], alog10((state.newflux)[gdo]), $
;    psym=1, color=clr.black
  oplot, state.std_wv, state.std_fx, color=clr.black

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Fit
;;;;;;;;;;;;;;;;;;;;

pro x_echfitstd_UpdateFit, state, ALL=all, ORDR=ordr

common x_echfitstd_fit

  if keyword_set(ALL) then begin
      iord = 0L
      ford = nord-1
  endif else begin
      ii = where(state.objstr.order EQ ordr)
      iord = ii[0]
      ford = ii[0]
  endelse

  ;; Loop
  for qq=iord,ford do begin
      gd = where(lindgen(state.objstr[qq].npix) LT state.objstr[qq].npix $
                 AND state.objstr[qq].sig GT 0., complement = bad)
      tmp = tot_fit[qq]
      ;; Fit
      sig_fit  = state.objstr[qq].sig[gd]
      force_pix = WHERE(sig_fit GE 1.0d8, nforce)
      ;; Give high weight to force pix
      gd_fit =  where(lindgen(state.objstr[qq].npix) LT state.objstr[qq].npix $
                      AND state.objstr[qq].sig GT 0. AND $
                      state.objstr[qq].sig LE 1.0d6)
      IF nforce NE 0 THEN sig_fit[force_pix] = $
         1.0d-5*djs_median(sig_fit[gd_fit])
      fit = x_fitrej(state.objstr[qq].box_wv[gd], $
                     state.objstr[qq].flux[gd], $
                     sig = sig_fit[gd], $
                     fitstr=tmp, rejpt=rejpt)
;                     state.objstr[qq].flux[gd], $
;                     sig=sqrt(state.objstr[qq].box_var[gd]>0), $

      if rejpt[0] NE -1 then state.objstr[qq].sig[gd[rejpt]] = 0.
      tot_fit[qq] = tmp

      ;; Residuals
      allf = x_calcfit(state.objstr[qq].box_wv, fitstr=tmp)
      state.resid1[*,qq] = state.objstr[qq].flux - allf

      ;; Std flux
      state.newflux[*,qq] = state.objstr[qq].box_fx / allf
      state.fit[*, qq] = allf
      ;; Sensitivity function
;      gd = where(state.objstr[qq].wave GT 0.)
;      linterp, state.std_wv, state.std_fx, $
;        state.objstr[qq].wave[gd], standard
;      sens_func[gd,qq] = standard / $
;        x_calcfit(state.objstr[qq].wave[gd], FITSTR=tmp)
  endfor
  gd = where(state.objstr.flux NE 0)
  djs_iterstat, (state.resid1)[gd], sigma=sigr, sigrej=2.5
  state.res1ymnx = 3*sigr*[-1.,1]

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Set order
;;;;;;;;;;;;;;;;;;;;
pro x_echfitstd_setordr, state

  ;; Grab x position
  xv = xgetx_plt(state,/strct)

  ;; Grab closest pixel
  mn = min(abs(state.objstr.box_wv-xv),imn)
  dwv = abs((state.objstr.box_wv)[imn] - (state.objstr.box_wv)[imn+1])

  ;; Grab closest within 10pixels
  gd = where(abs(state.objstr.box_wv-xv) LT 10*dwv)
  mx = max( (state.objstr.flux)[gd], imx)

  ;; Order
  state.currordr =  (state.ordnum)[gd[imx]]
  return
end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_echfitstd_Reset, state

  ;; Set the flux
  for qq=0L,n_elements(state.objstr)-1 do begin
      ;; gd
      gd = where(state.objstr[qq].box_wv GT 0.)
      linterp, state.std_wv, state.std_fx, state.objstr[qq].box_wv[gd], fx
      state.objstr[qq].flux[gd] = state.objstr[qq].box_fx[gd] / fx
      state.objstr[qq].sig[gd] = sqrt(state.objstr[qq].box_var[gd]) / fx
  endfor
  state.svsig = state.objstr.sig

  ;; Good ones
  gd = where(state.objstr.box_wv GT 0. and state.svsig GT 0. )

  ;; Set xdat, ydat
  xdat = (state.objstr.box_wv)[gd]
  ydat = (state.objstr.flux)[gd]

  ;; Now svxymnx
  state.svxymnx =  [min(xdat)-0.01*abs(max(xdat)-min(xdat)), $
                        min(ydat)-0.01*abs(max(ydat)-min(ydat)), $
                        max(xdat)+0.01*abs(max(xdat)-min(xdat)), $
                        max(ydat)+0.01*abs(max(ydat)-min(ydat))]
  state.xymnx = state.svxymnx


  ;; Standard star
  gd = where(state.std_wv GT state.xymnx[0] AND $
             state.std_wv LT state.xymnx[2])
  state.stdymnx = [min(state.std_fx[gd], max=mx), mx] 
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;  Zero Regions
;;;;

pro x_echfitstd_zeroreg, state

  ;; Check current order
  if state.currordr LT 1 then begin
      print, 'x_echfitstd:  Need to set the order with RMB first!'
      return
  endif

  ;; Set the flag
  if state.flg_zero MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_zero = 1 
      state.zero[0]= xgetx_plt(state, /strct)
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.zero[0] then begin
          state.zero[1] = tmp 
      endif else begin
          state.zero[1] = state.zero[0]
          state.zero[0] = tmp
      endelse
      ;; Find the pixels
      gdo = where(state.objstr.order EQ state.currordr)
      zpix = where(state.objstr[gdo].box_wv GT state.zero[0] AND $
                   state.objstr[gdo].box_wv LT state.zero[1], nzp)
      if nzp NE 0 then $
        state.objstr[gdo].sig[zpix] = 0.
                   
      ;; Update fit
      x_echfitstd_UpdateFit, state, ordr=state.currordr
      state.flg_zero = 0
  endelse

end

; This routine added by JFH to unzero regions  
pro x_echfitstd_unzeroreg, state

  ;; Check current order
  if state.currordr LT 1 then begin
      print, 'x_echfitstd:  Need to set the order with RMB first!'
      return
  endif

  ;; Set the flag
  if state.flg_unzero MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_unzero = 1 
      state.unzero[0]= xgetx_plt(state, /strct)
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.unzero[0] then begin
          state.unzero[1] = tmp 
      endif else begin
          state.unzero[1] = state.unzero[0]
          state.unzero[0] = tmp
      endelse
      ;; Find the pixels
      gdo = where(state.objstr.order EQ state.currordr)
      zpix = where(state.objstr[gdo].box_wv GT state.unzero[0] AND $
                   state.objstr[gdo].box_wv LT state.unzero[1], nzp)
      if nzp NE 0 then $
        state.objstr[gdo].sig[zpix] = state.svsig[zpix, gdo[0]]
                   
      ;; Update fit
      x_echfitstd_UpdateFit, state, ordr=state.currordr
      state.flg_unzero = 0
  endelse

end

; This routine added by JFH to unzero regions  
pro x_echfitstd_forcereg, state

  ;; Check current order
  if state.currordr LT 1 then begin
      print, 'x_echfitstd:  Need to set the order with RMB first!'
      return
  endif

  ;; Set the flag
  if state.flg_force MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_force = 1 
      state.force[0]= xgetx_plt(state, /strct)
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.force[0] then begin
          state.force[1] = tmp 
      endif else begin
          state.force[1] = state.force[0]
          state.force[0] = tmp
      endelse
      ;; Find the pixels
      gdo = where(state.objstr.order EQ state.currordr)
      zpix = where(state.objstr[gdo].box_wv GT state.force[0] AND $
                   state.objstr[gdo].box_wv LT state.force[1], nzp)
      if nzp NE 0 then $
         state.objstr[gdo].sig[zpix] = 1.0d10 ;*djs_median(state.svsig[*, gdo])
      ;; choose the weight to be 1d6 
      ;; Update fit
      x_echfitstd_UpdateFit, state, ordr=state.currordr
      state.flg_force = 0
   endelse

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Reset SVXY
pro x_echfitstd_UpdateSVXY, state
  ; Set xdat, ydat
  xdat = state.xtot[0:state.ntot-1]
  ydat = state.ytot[0:state.ntot-1]
  ; Now svxymnx
  state.dat_svxymnx =  [min(xdat)-0.01*abs(max(xdat)-min(xdat)), $
                        min(ydat)-0.01*abs(max(ydat)-min(ydat)), $
                        max(xdat)+0.01*abs(max(xdat)-min(xdat)), $
                        max(ydat)+0.01*abs(max(ydat)-min(ydat))]
  state.svxymnx = state.dat_svxymnx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_echfitstd, objstr, std_wv, std_fx, outfil, XSIZE=xsize, FITSTR=fitstr,$
                 INFIL = infil, savefil = savefil, insave = insave

common x_echfitstd_fit

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'fit = x_echfitstd(objstr, std_wv, std_fx, outfil, ' + $
      'FITSTR=,XSIZE=, INFIL=) [v1.1]'
    return
  endif 
  
  IF KEYWORD_SET(INSAVE) THEN restore, insave

  ;; STDSTRFIL
;  readcol, stdfil, swv, sfx
  ;; Initialize the common block
  x_echfitstd_initcommon, n_elements(objstr), fitstr
  sens_func = objstr.fx
  sens_func[*] = 0.

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200

  
  ;;  STATE
  state = { $
            objstr: objstr, $
            svsig: objstr.sig, $
            resid1: objstr.box_var, $
            newflux: objstr.box_var, $
            fit: objstr.box_var, $
            std_wv: std_wv, $
            std_fx: std_fx, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_zero: 0, $
            zero: fltarr(2), $
            flg_unzero: 0, $
            unzero: fltarr(2), $
            flg_force: 0, $
            force: fltarr(2), $
            ordnum: replicate(1., n_elements(objstr[0].box_wv)) # objstr.order, $
            currordr: 0L, $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            res1ymnx: fltarr(2), $
            stdymnx: fltarr(2), $
            size: lonarr(2), $
            help: strarr(50), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            funclist_id: 0L, $
            reject_base_id: 0L, $
            resid1_id: 0L, $
            resid2_id: 0L, $
            stdstr_id: 0L $
          }

  state.resid1[*] = -999.
  state.newflux[*] = 0.

      
  ;;    WIDGET
  base = WIDGET_BASE( title = 'x_echfitstd: Interactive ' + $
                      'Sensitivity Function', /column)
  state.base_id = base
  
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB,LMB -- Delete region',$
             'RMB -- Select order',$
             'W -- Reset screen',$
             's,s -- Delete Region', $
             'x,x -- Undelete Region', $
             'f,f -- Force region to be in fit',  $
             'l -- Set Left Window',$
             'r -- Set Right Window',$
             'b -- Set Bottom Window',$
             't -- Set Top Window',$
             '{} -- Pan',$
             'i,o -- Zoom in/our',$
;             'x -- Remove one Region',$
             'u -- Increase fit order 1',$
             'd -- Decrease fit order 1',$
             'F -- Reset y values on Flux window',$
             'k,K -- Decrease/Increase y-axis of residual window', $
             'R -- Remove the deleted regions of the current order',$
             'q -- Quit and save']
  state.help = strhelp

;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)
;      Error base
;  error_base_id  = widget_base(toolbar, /column)
;  error_btn_id = widget_button(error_base_id, value='Erase Error Msg', $
;                               uvalue='ERRORB')
;          Error Messages
;  state.error_msg_id = widget_text(error_base_id, value='', xsize=40, $
;                                   ysize=2, /scroll)
  
;      Drawing
  d_ysize = round(2.*ysize/3)
  state.size[0] = xsize
  state.size[1] = d_ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)
  
  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(state.draw_base_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
  ;; Resid1
  d_ysize = round(ysize/6)
  resid1b_id = widget_base(base, /column, /base_align_left, $
                                uvalue='RESID1_BASE', frame=2 ) 
  state.resid1_id = widget_draw(resid1b_id, xsize=state.size[0], $
                              ysize=d_ysize, /frame,  uvalue='RESID1')
  stdstrb_id = widget_base(base, /column, /base_align_left, $
                                uvalue='STDSTR_BASE', frame=2 ) 
  state.stdstr_id = widget_draw(stdstrb_id, xsize=state.size[0], $
                              ysize=d_ysize, /frame,  uvalue='STDSTR')
;  resid2b_id = widget_base(base, /column, /base_align_left, $
;                                uvalue='RESID2_BASE', frame=2 ) 
;  state.resid2_id = widget_draw(resid2b_id, xsize=state.size[0], $
;                              ysize=d_ysize, /frame,  uvalue='RESID2')
  
; Realize
  WIDGET_CONTROL, base, /realize
  
  
  ;; Update
  x_echfitstd_Reset, state
  x_echfitstd_UpdateFit, state, /all
  x_echfitstd_UpdatePlot, state, /all
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
  xmanager, 'x_echfitstd', base

  ;; Output
  ordr_fit = objstr.order
  sens_wv = objstr.wave
  nopix = objstr[5<(nord-1)].npix
  save, ordr_fit, tot_fit, nopix, filename=outfil
  ;; Save everything in case you want to fix one of the orders
  IF KEYWORD_SET(SAVEFIL) THEN save, filename = savefil
  return
end
