;+ 
; NAME:
; x_intphotcal   
;  Version 1.2
;
; PURPOSE:
;    Allows the user to interactively perform a photometric solution.
;   At its fullest, the code will calculate a zeropoint, an airmass
;   term and a color term.
;   The routine launches a GUI which allows the deletion 
;   of specific stars.
;
; CALLING SEQUENCE:
;   x_intphotcal, obs, landolt, outfil, XSIZE=, YSIZE=, /NCLR,
;    MIN_NOBS=, MIN_MOBS=, SETAM=
;
; INPUTS:
;   obs - Standard star observations  (stdstruct)
;   landolt - Landolt info (lndltstr)
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL -- ASCII file summary of the photometric solution.
;
; OPTIONAL KEYWORDS:
;   XSIZE = size of gui
;   YSIZE = size of gui
;   MIN_NOBS= Min n value for Landolt star
;   MIN_MOBS= Min m value for Landolt star
;   /NCLR = Solve without color terms
;   SETAM = Value taken for AM term 
;        (Note: You cannot choose a value of 0. exactly!)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_intphotcal, obs, landolt, outfil
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;   07-Jan-2002 Revised by JXP
;   14-Oct-2002 Revised by JXP [setup no AM, no CLR terms]
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Events
;;;;

pro x_iphotcal_event, ev


  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DRAW_Mag' : begin
          case ev.type of
              0 : begin ; Button press
                  state.press = ev.press
                  case ev.press of
                      1 : begin ; Delete/Undelete a star
                          x_iphotcal_DelStr, state, state.tag_Mag 
                          x_iphotcal_Solve, state
;                          x_iphotcal_Reset, state
                      end
                      2 : state.tv_Mag.xymnx = state.tv_Mag.svxymnx    ; Zoom out
                      4 : x_iphotcal_SetZoom, state, 1, state.tag_Mag  ; Zoom 
                  endcase
                  if ev.press NE 4 then x_iphotcal_UpdPlot, state
              end
              1 : begin ; Button release
                  if( state.press EQ 4) then begin
                      x_iphotcal_SetZoom, state, 2, state.tag_Mag
                      x_iphotcal_UpdPlot, state
                  endif
              end
              2 : begin ; Motion event
                  state.tv_Mag.xcurs = ev.x
                  state.tv_Mag.ycurs = ev.y
              end
          endcase
      end
      'DRAW_CLR' : begin
          case ev.type of
              0 : begin ; Button press
                  state.press = ev.press
                  case ev.press of
                      1 : begin ; Delete/Undelete a star
                          x_iphotcal_DelStr, state, state.tag_CLR 
                          x_iphotcal_Solve, state
                          x_iphotcal_Reset, state
                      end
                      2 : state.tv_CLR.xymnx = state.tv_CLR.svxymnx    ; Zoom out
                      4 : x_iphotcal_SetZoom, state, 1, state.tag_CLR  ; Zoom 
                  endcase
                  if ev.press NE 4 then x_iphotcal_UpdPlot, state
              end
              1 : begin ; Button release
                  if( state.press EQ 4) then begin
                      x_iphotcal_SetZoom, state, 2, state.tag_CLR
                      x_iphotcal_UpdPlot, state
                  endif
              end
              2 : begin ; Motion event
                  state.tv_CLR.xcurs = ev.x
                  state.tv_CLR.ycurs = ev.y
              end
          endcase
      end
      'DRAW_AM' : begin
          case ev.type of
              0 : begin ; Button press
                  state.press = ev.press
                  case ev.press of
                      1 : begin ; Delete/Undelete a star
                          x_iphotcal_DelStr, state, state.tag_AM 
                          x_iphotcal_Solve, state
                          x_iphotcal_Reset, state
                      end
                      2 : state.tv_AM.xymnx = state.tv_AM.svxymnx    ; Zoom out
                      4 : x_iphotcal_SetZoom, state, 1, state.tag_AM  ; Zoom 
                  endcase
                  if ev.press NE 4 then x_iphotcal_UpdPlot, state
              end
              1 : begin ; Button release
                  if( state.press EQ 4) then begin
                      x_iphotcal_SetZoom, state, 2, state.tag_AM
                      x_iphotcal_UpdPlot, state
                  endif
              end
              2 : begin ; Motion event
                  state.tv_AM.xcurs = ev.x
                  state.tv_AM.ycurs = ev.y
              end
          endcase
      end
      'NOSV' : begin
          print, 'x_intphotcal: Output will not be saved!'
          widget_control, ev.top, /destroy 
          return
      end
      'DONE' : begin
          x_iphotcal_Output, state
          widget_control, ev.top, /destroy 
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_iphotcal_Reset, state


; Plotting

  ; Mag
  gdmag = where(state.obs.sig_Mag LT 9. AND state.obs.Mag LT 25.)
  state.tv_Mag.svxymnx[0] = min(state.obs[gdmag].Mag,max=mx) - 0.05
  state.tv_Mag.svxymnx[2] = mx + 0.05
  state.tv_Mag.svxymnx[1] = min(state.RES[gdmag], max=mx) - 0.05
  state.tv_Mag.svxymnx[3] = mx + 0.05
  state.tv_Mag.xymnx = state.tv_Mag.svxymnx
  ; Clr
  state.tv_CLR.svxymnx[0] = min(state.obs.CLR,max=mx) - 0.05
  state.tv_CLR.svxymnx[2] = mx + 0.05
  state.tv_CLR.svxymnx[1] = state.tv_Mag.svxymnx[1] 
  state.tv_CLR.svxymnx[3] = state.tv_Mag.svxymnx[3] 
  state.tv_CLR.xymnx = state.tv_CLR.svxymnx
  ; AM
  state.tv_AM.svxymnx[0] = min(state.obs.AM,max=mx) - 0.05
  state.tv_AM.svxymnx[2] = mx + 0.05
  state.tv_AM.svxymnx[1] = state.tv_Mag.svxymnx[1] 
  state.tv_AM.svxymnx[3] = state.tv_Mag.svxymnx[3] 
  state.tv_AM.xymnx = state.tv_AM.svxymnx

end

;;;;;;;;;
; Perform Photo solution


pro x_iphotcal_Solve, state

  widget_control, /hourglass


  tmp_obs = state.obs
  tmpflg = 0

  if state.flg_CLR EQ 0 then tmpflg = tmpflg + 1
  if state.flg_AM EQ 0 then tmpflg = tmpflg + 2


  ; Call photcal
  case tmpflg of
      0: begin  ; CLR + AM
          x_photcal, tmp_obs, state.landolt, fit, sigfit, NCORR=ncorr, CHISQ=chisq,$
            MIN_NOBS=state.min_nobs, MIN_MOBS=state.min_mobs
          state.fit = fit
          state.sigfit = sigfit
      end
      1: begin  ; No CLR
          x_photcal, tmp_obs, state.landolt, fit, sigfit, NCORR=ncorr, $
            MIN_NOBS=state.min_nobs, MIN_MOBS=state.min_mobs, $
            CHISQ=chisq, /NOCLR
          state.fit = fit
          state.sigfit = sigfit
          state.fit[1] = -state.fit[1]  ;; Kludge to make AM positive
      end
      2: begin  ; COLOR but no AM
          x_photcal, tmp_obs, state.landolt, fit, sigfit, NCORR=ncorr, CHISQ=chisq,$
            MIN_NOBS=state.min_nobs, MIN_MOBS=state.min_mobs, SETAM=state.setam
          state.fit[0] = fit[0]
          state.sigfit[0] = sigfit[0]
          state.fit[1] = state.setam
          state.sigfit[1] = 0.
          state.fit[2] = fit[1]
          state.sigfit[2] = sigfit[1]
      end
      3: begin  ; no CLR, no AM
          if state.setam GT 0 then state.setam = -state.setam
          x_photcal, tmp_obs, state.landolt, fit, sigfit, NCORR=ncorr, CHISQ=chisq,$
            MIN_NOBS=state.min_nobs, MIN_MOBS=state.min_mobs, SETAM=state.setam, $
            /NOCLR
          state.fit[0] = fit[0]
          state.sigfit[0] = sigfit[0]
          state.fit[1] = abs(state.setam)
          state.sigfit[1] = 0.
          state.fit[2] = 0.
      end
      else: begin
          print, 'x_intphotcal: Not set up for this'
          stop
      endelse
  endcase

  state.obs.flg_anly = tmp_obs.flg_anly
  state.chisq = chisq
   ; Deal with bad data points
  bad = where(state.obs.flg_anly EQ 0)
  state.nbad = n_elements( temporary(bad) )

; Color and Residual 

  for i=0,n_elements(state.obs)-1 do begin
      istr = where(state.obs[i].Name EQ state.landolt.Name)

      if state.flg_CLR EQ 1 then $
        state.obs[i].CLR = x_lndltclr(state.obs[i].sCLR, state.landolt[istr])

      ; Residual
      state.RES[i] = x_lndltMag(state.obs[i].filter, state.landolt[istr]) - $
        (state.obs[i].Mag + state.fit[0] - state.fit[1]*state.obs[i].AM $
         - state.fit[2]*state.obs[i].CLR)
  endfor
  good = where(state.obs.flg_anly EQ 1, ngd)
  state.rms = sqrt(total( state.Res[good]^2 )/float(ngd-1))


; Update screen

  x_iphotcal_UpdLbl, state

  delvarx, tmp_obs, good

end

;;;;;;;;;;;;;;;;;;;;
;  Update Label for the Fit
;;;;;;;;;;;;;;;;;;;;

pro x_iphotcal_UpdLbl, state

  state.lbl_fit[0] = strjoin([strmid( strtrim(state.fit[0],2), 0, 6), ' +/- ', $
                        strmid( strtrim(state.sigfit[0],2), 0,5) ])
  state.lbl_fit[1] = strjoin([strmid( strtrim(state.fit[1],2), 0, 5), ' +/- ', $
                        strmid( strtrim(state.sigfit[1],2), 0,5) ])
  state.lbl_fit[2] = strjoin([strmid( strtrim(state.fit[2],2), 0, 6), ' +/- ', $
                        strmid( strtrim(state.sigfit[2],2), 0,5) ])

  widget_control, state.lblfit1_id, set_value=state.lbl_fit[0]
  widget_control, state.lblfit2_id, set_value=state.lbl_fit[1]
  widget_control, state.lblfit3_id, set_value=state.lbl_fit[2]
  widget_control, state.chisq_id, set_value=strtrim(state.chisq,2)
  widget_control, state.rms_id, set_value=strtrim(state.rms,2)


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro x_iphotcal_UpdPlot, state
  
;;;;;;;;;;;;;;;;;;
; Mag Window

  widget_control, state.Mag_draw_id, get_value=wind
  	wset, wind

  ; Color
  clr = getcolor(/load)

  good = where(state.obs.flg_anly EQ 1)
  plot, state.obs[good].Mag, state.RES[good], psym=1, $
    position=state.tv_Mag.pos, xrange=[state.tv_Mag.xymnx[0],$
                                       state.tv_Mag.xymnx[2]], $
    yrange=[state.tv_Mag.xymnx[1],state.tv_Mag.xymnx[3]], $
    xstyle=1, ystyle=1, YTITLE='!17RES', TITLE='Mag', $
    background=clr.white, color=clr.black
  ; Deleted
  bad = where(state.obs.flg_anly EQ 0)
  if state.nbad EQ 1 then oplot, extrac(state.obs.Mag,bad[0],1), $
            extrac(state.RES,bad[0],1), psym=2, color=clr.red $
  else oplot, state.obs[bad].Mag, state.RES[bad], psym=2, color=clr.red

;;;;;;;;;;;;;;;;;;
; CLR Window

  if state.flg_CLR then begin
      widget_control, state.CLR_draw_id, get_value=wind
      wset, wind

      good = where(state.obs.flg_anly EQ 1)
      plot, state.obs[good].CLR, state.RES[good], psym=1, $
        position=state.tv_CLR.pos, xrange=[state.tv_CLR.xymnx[0],$
                                           state.tv_CLR.xymnx[2]], $
        yrange=[state.tv_CLR.xymnx[1],state.tv_CLR.xymnx[3]], $
        xstyle=1, ystyle=1, YTITLE='RES', TITLE=state.obs[good[0]].sCLR, $
        background=clr.white, color=clr.black
                                ; Deleted
      bad = where(state.obs.flg_anly EQ 0)
      if state.nbad EQ 1 then oplot, extrac(state.obs.CLR,bad[0],1), $
        extrac(state.RES,bad[0],1), psym=2, color=clr.red $
      else oplot, state.obs[bad].CLR, state.RES[bad], psym=2, color=clr.red
  endif
      
;;;;;;;;;;;;;;;;;;
; AM Window

  if state.flg_AM then begin
      widget_control, state.AM_draw_id, get_value=wind
      wset, wind
      
      good = where(state.obs.flg_anly EQ 1)
      plot, state.obs[good].AM, state.RES[good], psym=1, $
        position=state.tv_AM.pos, xrange=[state.tv_AM.xymnx[0],$
                                      state.tv_AM.xymnx[2]], $
        yrange=[state.tv_AM.xymnx[1],state.tv_AM.xymnx[3]], $
        xstyle=1, ystyle=1, YTITLE='RES', TITLE='AM', $
        background=clr.white, color=clr.black
                                ; Deleted
      bad = where(state.obs.flg_anly EQ 0)
      if state.nbad EQ 1 then oplot, extrac(state.obs.AM,bad[0],1), $
        extrac(state.RES,bad[0],1), psym=2, color=clr.red $
      else oplot, state.obs[bad].AM, state.RES[bad], psym=2, color=clr.red 
  endif
      

end


;;;;;;;;;
;  Set Zoom 
;;;;;;;;;

pro x_iphotcal_SetZoom, state, flg, tvid

  if flg EQ 1 then begin
      state.tmpreg[0] = x_tvx(state.(tvid))
      state.tmpreg[1] = x_tvy(state.(tvid))
  endif else begin
      state.tmpreg[2] = x_tvx(state.(tvid))
      state.tmpreg[3] = x_tvy(state.(tvid))

      if (state.tmpreg[0] NE state.tmpreg[2] AND $
          state.tmpreg[1] NE state.tmpreg[3]) then begin
          state.(tvid).xymnx[0] = state.tmpreg[0] < state.tmpreg[2]
          state.(tvid).xymnx[2] = state.tmpreg[0] > state.tmpreg[2]
          state.(tvid).xymnx[1] = state.tmpreg[1] < state.tmpreg[3]
          state.(tvid).xymnx[3] = state.tmpreg[1] > state.tmpreg[3]
      endif
  endelse

end

;;;;;;;;;
;  Delete/Undelete Stars
;;;;;;;;;

pro x_iphotcal_DelStr, state, tag_tv

  widget_control, /hourglass

  x = x_tvx(state.(tag_tv))
  y = x_tvy(state.(tag_tv))

  case tag_tv of 
      0 : begin  ; Mag tv
          dist = sqrt( (x-state.obs.Mag)^2 + (y-state.RES)^2 )
          mdist = min(dist, jmin)
          state.obs[jmin].flg_anly = abs(state.obs[jmin].flg_anly - 1)
      end
      1 : begin  ; AM tv
          dist = sqrt( (x-state.obs.AM)^2 + (y-state.RES)^2 )
          mdist = min(dist, jmin)
          state.obs[jmin].flg_anly = abs(state.obs[jmin].flg_anly - 1)
      end
      2 : begin  ; CLR tv
          dist = sqrt( (x-state.obs.CLR)^2 + (y-state.RES)^2 )
          mdist = min(dist, jmin)
          state.obs[jmin].flg_anly = abs(state.obs[jmin].flg_anly - 1)
      end
      else :
  endcase
end
  

;;;;
; Output 
;;;;

pro x_iphotcal_Output, state

  widget_control, /hourglass

; Name

; Write
  close, /all
  openw, 1, state.outfil

  printf, 1, state.obs[0].filter
  printf, 1, state.obs[0].sCLR

  printf, 1, state.fit
  printf, 1, state.sigfit

  close, 1
 
return

end

; 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_intphotcal, obs, landolt, outfil, XSIZE=xsize, YSIZE=ysize, NCLR=nclr, $
                  MIN_NOBS=min_nobs, MIN_MOBS=min_mobs, SETAM=setam

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'x_intphotcal, obs, landolt, outfil, XSIZE=, YSIZE=, MIN_NOBS='
      print, '      MIN_MOBS=, SETAM=, /NCLR (v1.1)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1000
  if not keyword_set( YSIZE ) then    ysize = 600
  if not keyword_set( SETAM ) then    setam = 0
  if n_elements( MIN_NOBS ) EQ 0 then min_nobs = 4
  if n_elements( MIN_MOBS ) EQ 0 then min_mobs = 2

;  Some preparsing

  nstrs = n_elements(obs)

;    STATE
  ; tv structure
  tmp = {tvstruct}
  tmp.pos = [0.13, 0.1, 0.95, 0.95]
  

  state = {             $
            tv_Mag: tmp, $              ; Draw Windows
            tv_AM: tmp, $              
            tv_CLR: tmp, $            
            tag_Mag: 0, $
            tag_AM: 1, $
            tag_CLR: 2, $
            nstrs: nstrs, $                  ; Important stuff
            RES: fltarr(nstrs), $
            nbad: 0, $
            obs: obs, $
            landolt: landolt, $              ; Landolt
            min_nobs: min_nobs, $
            min_mobs: min_mobs, $
            tmpreg: fltarr(4), $         
            flg_CLR: 1, $   ; 0 = No color term
            flg_AM: 1, $    ; 0 = No AM 
            setam: 0., $
            fit: fltarr(3), $
            sigfit: fltarr(3), $
            chisq: 0., $
            rms: 0., $
            lbl_fit: strarr(3), $
            outfil: outfil, $
            press: 0, $
            base_id: 0L, $      ; Widgets
            top_base_id: 0L, $
            btm_base_id: 0L, $
            butt_base_id: 0L, $
            lblfit1_id: 0L, $
            lblfit2_id: 0L, $
            lblfit3_id: 0L, $
            chisq_id: 0L, $
            rms_id: 0L, $
            draw_actv_id: 0L, $
            CLR_draw_id: 0L, $
            Mag_draw_id: 0L, $
            AM_draw_id: 0L $
          }

; Color

  if keyword_set( NCLR ) then state.flg_CLR = 0
  if keyword_set( SETAM ) then begin
      state.flg_AM = 0
      state.setam = setam
  endif

;    WIDGET
  base = WIDGET_BASE( title = 'x_intphotcal: ', /column, xsize=xsize,$
                    ysize=ysize, xoffset=300, yoffset=200)
  state.base_id = base

;;;;;;
;    Top Base

  state.top_base_id = widget_base(base, /row, xsize=xsize, $
                                  ysize=ysize/2)
  ; Buttons
  state.butt_base_id = widget_base(state.top_base_id, /row, $
                                   /base_align_center, /align_center, $
                                   xsize=xsize/2)
;        Version + Name
  labelbase = widget_base(state.butt_base_id, /column, /base_align_center, $
                          /frame, /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_intphotcal', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.1)', /align_center)

  titlfit = WIDGET_LABEL(labelbase, value='Fit Values', /align_center)
  fitbase = widget_base(labelbase, /column, /base_align_left, $
                          /frame, /align_left)
  state.lblfit1_id = cw_field(fitbase, value=' ', xsize=17, title='ZP:')
  state.lblfit2_id = cw_field(fitbase, value=' ', xsize=17, title='AM:')
  state.lblfit3_id = cw_field(fitbase, value=' ', xsize=17, title='CR:')
  ; CHISQ AND RMWS
  state.chisq_id = CW_FIELD(fitbase, value=' ', title='CHISQ:', /floating, xsize=6)
  state.rms_id = cw_field(fitbase, value=' ', title='RMS:', /floating, xsize=6)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help & Done
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB -- Delete/Undelete point',$
             'CMB -- Full View',$
             'RMB+drag -- Zoom region', $
             'Red star -- Deleted point' $
            ]
  help_base = widget_base(state.butt_base_id, /column, /align_center)
  help_text_id = widget_text(help_base, value=strhelp, $
                             xsize=30, ysize=10, /scroll)
  done = WIDGET_BUTTON(help_base, value='DONE',uvalue='DONE')
  done = WIDGET_BUTTON(help_base, value='NOSV',uvalue='NOSV')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ; CLR tv
  state.tv_CLR.winsize[0] = xsize/2
  state.tv_CLR.winsize[1] = ysize/2
  state.CLR_draw_id = widget_draw(state.top_base_id, $
                                  xsize=state.tv_CLR.winsize[0], $
                                  ysize=state.tv_CLR.winsize[1], /frame, $
                                 /button_events, /motion_events, $
                                  uvalue='DRAW_CLR',retain=2)
;;;;;;
;    Bottom Base
  state.btm_base_id = widget_base(base, /row, xsize=xsize, $
                                  ysize=ysize/2)
      ; Mag tv
  state.tv_Mag.winsize[0] = xsize/2
  state.tv_Mag.winsize[1] = ysize/2
  state.Mag_draw_id = widget_draw(state.btm_base_id, $
                                  xsize=state.tv_Mag.winsize[0], $
                                  ysize=state.tv_Mag.winsize[1], /frame, $
                                 /button_events, /motion_events, $
                                  uvalue='DRAW_Mag',retain=2)
      ; AM tv
  state.tv_AM.winsize[0] = xsize/2
  state.tv_AM.winsize[1] = ysize/2
  state.AM_draw_id = widget_draw(state.btm_base_id, $
                                  xsize=state.tv_AM.winsize[0], $
                                  ysize=state.tv_AM.winsize[1], /frame, $
                                 /button_events, /motion_events, $
                                  uvalue='DRAW_AM',retain=2)
  
; Realize
  WIDGET_CONTROL, base, /realize

; Color Table 
  loadct, 2, /silent

  x_iphotcal_Solve, state
  x_iphotcal_Reset, state
  x_iphotcal_UpdPlot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_iphotcal', base

  return
end

