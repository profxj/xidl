;+ 
; NAME:
; x_continuum   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to spectroscopic data interactively
;
; CALLING SEQUENCE:
;   
;   x_continuum, [xdat], ydat, ysig, FITSTR= 
;
; INPUTS:
;   xdat       - Values along one dimension [optional]
;   ydat       - Values along the other
;   ysig       - Error in ydat
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   fitstr     - Fit structure
;   INFLG      - 0 = yin, ysin as data arrays (fits allowed)
;                1 = One fits file (flux, sig)
;                2 = One fits file (flux, sig, wave)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_continuum, wav, fx, sigfx
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
;   26-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_continuum_initcommon

;

common x_continuum_fit, fin_fit, fitprm, xfit, svdelpts, fc_val

fitprm = { fitstrct }
fitprm.lsig = 2.5
fitprm.hsig = 3.
fitprm.niter = 1
fitprm.nord = 3
fitprm.flg_rej = 1
fitprm.maxrej = 100
fitprm.func = 'LEGEND'

end

;;;;
; Events
;;;;

function x_continuum_ev, ev
;pro x_continuum_event, ev

common x_continuum_fit

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  flg_plt = 1
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
          x_continuum_UpdateFit, state
      end
      'UP' : begin
          fitprm.nord = fitprm.nord+1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm.nord, format='(i4)')
          ; Update fit
          x_continuum_UpdateFit, state
      end
      'DOWN' : begin
          fitprm.nord = (fitprm.nord-1) > 1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm.nord, format='(i4)')
          ; Update fit
          x_continuum_UpdateFit, state
      end
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return, 1
      end
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Add a Data point with weight 30
                          state.xtot[state.ntot] = xgetx_plt(state,/strct)
                          state.ytot[state.ntot] = xgety_plt(state,/strct)
                          state.wtot[state.ntot] = 1./sqrt(state.med_ivar)/30.
                          ;; Inverse variance
                          state.ivtot[state.ntot] = 100*state.med_ivar 

                          state.gdpix[state.ntot] = 1
                          state.ntot = state.ntot + 1
                          ; Update Fit
                          x_continuum_UpdateFit, state
                      end 
                      4 : x_continuum_Delete, state ; Delete/Undelete a data point
                      else :
                  endcase
              end
              1 : begin ; Button release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return, 1
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return, 1
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          if (state.flg_reg EQ 1 AND eventch NE 's') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another s !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return, 1
          endif
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another z !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return, 1
          endif
          if (state.flg_svcont EQ 1 AND eventch NE 'g') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another g !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return, 1
          endif
          if (state.flg_addcont EQ 1 AND eventch NE 'c') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another c !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return, 1
          endif
          case eventch of
              'u': begin
                  fitprm.nord = fitprm.nord+1
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm.nord, format='(i4)')
                  flg_plt = 0
;                  x_continuum_UpdateFit, state
              end
              'd': begin
                  fitprm.nord = (fitprm.nord-1) > 1
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm.nord, format='(i4)')
;                  x_continuum_UpdateFit, state
                  flg_plt = 0
              end
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'Z': state.xymnx[1] = 0.
              'z': begin  ; Zoom
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif
              end
              'i': x_speczoom, state, 0   ; Zoom out
              'o': x_speczoom, state, 1   ; Zoom out
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              's': begin  ; Region
                  x_continuum_GetReg, state
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif 
              end
              'F': state.svxymnx = state.xymnx ; Save screen area to svxymn
              'f': x_continuum_UpdateFit, state
              'W': begin
                  if state.flg_plot MOD 2 EQ 0 then $
                    state.xymnx = state.svxymnx $ ; Reset the screen
                  else x_continuum_SetResiduals, state
              end
              '!': x_continuum_Reset, state   ; Reset all
              '$': begin ; Reset screen view
                  state.svxymnx = state.svsvxymnx 
                  state.xymnx = state.svxymnx 
              end
              'U': mwrfits, fc_val, state.outfil, /create
              'D': x_continuum_Delete, state ; Delete/Undelete a data point
              'x': x_continuum_DelReg, state        ; Delete one region
              'X': x_continuum_DelReg, state, /all  ; Delete all regions and reset
              'S': x_continuum_SaveCont, state, /all  ; Save entire continuum
              'g': begin
                  x_continuum_SaveCont, state ; Save piece of continuum
                  if state.flg_svcont EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif
              end
              'c': begin ; Add continuum points to the fit
                  x_continuum_AddCont, state 
                  if state.flg_addcont EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif else x_continuum_UpdateFit, state
              end
              'R': x_continuum_SetResiduals, state  ; Plot residuals
              'V': x_continuum_SetValues, state     ; Plot values
              'C': x_continuum_ClearRej, state ; Clear all rejected points
              'q': begin
                  x_continuum_setpnt, state
                  widget_control, ev.top, /destroy
                  return, 0
              end
              else:  ; Nothing
          endcase
      end
;            REJECTION
      'HSIGVAL': begin
          fitprm.hsig = ev.value
          if fitprm.niter EQ 0 then fitprm.niter=3
          x_continuum_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.lsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'LSIGVAL': begin
          fitprm.lsig = ev.value
          if fitprm.niter EQ 0 then fitprm.niter=3
          x_continuum_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.hsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'DONE' : begin
          x_continuum_setpnt, state
          widget_control, ev.top, /destroy
          return, 0
      end
      else:
  endcase

; Update Plot
  if flg_plt EQ 1 then x_continuum_UpdatePlot, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return, 1
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_continuum_UpdatePlot, state
  
common x_continuum_fit

; Plot Data

  widget_control, /hourglass   

  widget_control, state.draw_id, get_value=wind
  wset, wind

  color = getcolor(/load)

  if state.flg_plot MOD 2 EQ 0 then begin  ; NORMAL

;  ORIGINAL DATA NOT IN REGION
      gdorg = where( state.wave[0:state.norg-1] GT state.xymnx[0] $
                    AND state.wave[0:state.norg-1] LT state.xymnx[2], count) 
      if count NE 0 then begin
          plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
            /nodata, xstyle=1, ystyle=1, $
            position=state.pos,  background=color.white, color=color.black 
          oplot, state.wave[gdorg], state.sig[gdorg], psym=10, $
            color=color.orange 
          oplot, state.wave[gdorg], state.fx[gdorg], psym=10, $
            color=color.black 
      endif else $
        plot, [-1.0],  [state.svxymnx[3]+1.], psym=1, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        background=color.white, color=color.black

      ; Overplot regions
      for q=0L, state.nreg-1 do begin
          pts = where(state.wave LT state.reg[q,1] AND $
                      state.wave GT state.reg[q,0], npts) 
          if npts NE 0 then $
            oplot, [state.wave[pts]], [state.fx[pts]], psym=10, $
                    color=color.blue ; DATA IN REGION
      endfor

      ; Overplot rejected/other points
      for q=0,7 do begin
          if q EQ 1 OR q EQ 3 or q EQ 4 OR q EQ 6 OR q EQ 0 OR q EQ 5 then continue
          pts = where(state.gdpix[0:state.norg-1] EQ q AND $
                      state.wave[0:state.norg-1] GT state.xymnx[0] $
                      AND state.wave[0:state.norg-1] LT state.xymnx[2], count)
          if count NE 0 then begin
              case q of
                  0: oplot, [state.wave[pts]], [state.fx[pts]], psym=2, $
                    color=color.green ; DELETED DATA OUT OF REGION
                  2: oplot, [state.wave[pts]], [state.fx[pts]], psym=2, $
                    color=color.purple ; DELETED DATA IN REGION
                  5: oplot, [state.wave[pts]], [state.fx[pts]], psym=5, $
                    color=color.pink ; REJECTED DATA OUT OF REGION
                  7: oplot, [state.wave[pts]], [state.fx[pts]], psym=5, $
                    color=color.red ; REJECTED DATA IN REGION
                  else:
              endcase
          endif
      endfor
      
      
;    Added points
      if(state.ntot GT state.norg) then begin
          oplot, extrac(state.xtot,state.norg,state.ntot-state.norg), $
            extrac(state.ytot,state.norg,state.ntot-state.norg), psym=5, $
            color=color.cyan
      endif

;    Fit
      if state.nreg NE 0 then begin
          pts = where(xfit LT state.xymnx[2] AND xfit GT state.xymnx[0], npts)
          if npts NE 0 then $
            oplot, [xfit[pts]], [fin_fit[pts]], color=color.red ; Best fit Red  
      endif

;    Continuum
      pts = where(state.wave[0:state.norg-1] LT state.xymnx[2] AND $
                  state.wave[0:state.norg-1] GT state.xymnx[0], npts)
      if npts NE 0 then $
        oplot, [state.wave[pts]], [fc_val[pts]], color=color.green ; Best fit Red  

  endif else begin  ; RESIDUALS

;  ORIGINAL DATA NOT IN REGION
      gdorg = where(state.gdpix[0:state.norg-1] EQ 1, count) ; Original data 
      if count NE 0 then begin
          plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
            /nodata, xstyle=1, ystyle=1, $
            position=state.pos,  background=color.white, color=color.black 
          oplot, state.wave[gdorg], state.residuals[gdorg], psym=1, $
            color=color.black 
      endif else $
        plot, [-1.0],  [state.svxymnx[3]+1.], psym=1, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        background=color.white, color=color.black

      for q=0,7 do begin
          if q EQ 1 then continue
          pts = where(state.gdpix[0:state.norg-1] EQ q, count)
          if count NE 0 then begin
              case q of
                  0: oplot, [state.wave[pts]], [state.residuals[pts]], psym=2, $
                    color=color.red ; DELETED DATA OUT OF REGION
                  2: oplot, [state.wave[pts]], [state.residuals[pts]], psym=2, $
                    color=color.purple ; DELETED DATA IN REGION
                  3: oplot, [state.wave[pts]], [state.residuals[pts]], psym=1, $
                    color=color.blue ; GOOD DATA IN REGION
                  4: 
                  5: oplot, [state.wave[pts]], [state.residuals[pts]], psym=5, $
                    color=color.pink ; REJECTED DATA OUT OF REGION
                  6: 
                  7: oplot, [state.wave[pts]], [state.residuals[pts]], psym=5, $
                    color=color.yellow ; REJECTED DATA OUT OF REGION
                  else:
              endcase
          endif
      endfor
      
      
;    Added points
      if(state.ntot GT state.norg) then begin
          oplot, extrac(state.xtot,state.norg,state.ntot-state.norg), $
            extrac(state.residuals,state.norg,state.ntot-state.norg), psym=5, $
            color=color.cyan
      endif
  endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Fit
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_UpdateFit, state

common x_continuum_fit
  widget_control, /hourglass   

; Set Good points

  if state.nreg EQ 0 then return $ ; Require regions
  else begin
      if state.norg LT state.ntot then begin   ; Deal with added points
          gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count) 
          for i=state.norg, state.ntot-1 do gdpt = [gdpt,i]
          count = count + state.ntot - state.norg
      endif else gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count)
  endelse

  mask = lonarr(state.ntot)
  mask[gdpt] = 1
  fitprm.minpt = fitprm.minpt > count/2

; Catch error

  if count EQ 0 then begin
      widget_control, state.error_msg_id, $
        set_value='No good data!  Adjust regions as necessary...'
      widget_control, base_id, set_uvalue=state, /no_copy
      return
  endif

; Error

  if fitprm.nord EQ 0 then begin
      widget_control, state.error_msg_id, set_value='FIT requires nord > 0'
      fitprm.nord = 1
      widget_control, state.lblordr_id, $
        set_value=string(fitprm.nord, format='(i4)')
  endif else begin  ; FIT!
      if fitprm.flg_rej EQ 0 then begin
          fit = x_fit(state.xtot[gdpt], state.ytot[gdpt], $
                      FITSTR=fitprm, $
                      SIG=state.wtot[gdpt], $
                      IVAR=state.ivtot[gdpt] )
          if fit[0] EQ -1 then begin
              fin_fit = fltarr(n_elements(xfit))
              widget_control, state.error_msg_id, $
                set_value='Bad fit! Adjust nord probably'
              widget_control, base_id, set_uvalue=state, /no_copy
              return
          endif
      endif else begin ; FIT with REJECTION!
          rejpt = -1
          fit = x_fitrej(state.xtot[gdpt],$
                         state.ytot[gdpt], $
                         SIG=state.wtot[gdpt], $
                         IVAR=state.ivtot[gdpt], $
                         FITSTR=fitprm, $
                         REJPT=rejpt )
;          fit = x_fitrej(state.xtot[0:state.ntot-1], $
;                         state.ytot[0:state.ntot-1], $
;                         SIG=state.wtot[0:state.ntot-1], $
;                         IVAR=state.ivtot[0:state.ntot-1], $ 
;                         MSK=mask,$
;                         FITSTR=fitprm, $
;                         REJPT=rejpt )
          if rejpt[0] NE -1 then begin
              rejpt = gdpt[rejpt]
              state.gdpix[rejpt] = state.gdpix[rejpt]+4
          endif
      endelse
  endelse
  

  ;; FIT CHECK
  if fit[0] EQ -1 then stop

  fit = x_calcfit(state.xtot[0:state.ntot-1], FITSTR=fitprm)

  ; Calculate values
  fin_fit = fit

  delvarx, mask, fit

  ; RESIDUALS AND RMS
  x_continuum_CalcResid, state

    ; Set Good points

  if state.nreg EQ 0 then $
    gdpt = where(state.gdpix[0:state.ntot-1] EQ 1, count) $
  else begin
      if state.norg LT state.ntot then begin   ; Deal with added points
          gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count) 
          for i=state.norg, state.ntot-1 do gdpt = [gdpt,i]
          count = count + state.ntot - state.norg
      endif else gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count)
  endelse

  fitprm.rms = sqrt(total( state.residuals[gdpt]^2 ) / float(n_elements(gdpt)-1))
  widget_control, state.rms_id, set_value=fitprm.rms

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_Reset, state

common x_continuum_fit

;
  state.xtot[0:n_elements(state.wave)-1] = state.wave
  state.ytot[0:n_elements(state.wave)-1] = state.fx
  if state.flg_wgt EQ 0 then state.wtot[0:n_elements(state.wave)-1] = 1.0 $
  else state.wtot[0:n_elements(state.wave)-1] = state.sig
  state.ivtot[0:n_elements(state.wave)-1] = state.ivar

  state.ntot = n_elements(state.wave)


; gdpix Flag
;   1 = Good
;   2 = In region
;   4 = Rejected by fit
  state.gdpix[0:n_elements(state.wave)-1] = 1

; Regions
  state.nreg = 0

; Plotting
  state.xymnx = state.svsvxymnx
  state.svxymnx = state.svsvxymnx

; Fitting
  fitprm.nord = 3
  fitprm.func = 'LEGEND'

end

;;;;;;;;;;;;;;;;;;;;
;  Delete
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_Delete, state

; Deletes the nearest data point in pixel space

  ; Get into pixel space
  pix_x = xgetxpix_plt(state.xtot[0:state.ntot-1], $
                     state.pos, state.xymnx, state.size)
  if state.flg_plot MOD 2 EQ 0 then $ ; VALUES
    pix_y = xgetypix_plt(state.ytot[0:state.ntot-1], state.pos, $
                         state.xymnx, state.size) $
  else $  ; RESIDUALS
    pix_y = xgetypix_plt(state.residuals[0:state.ntot-1], $
                         state.pos, state.xymnx, state.size)

  ; Minimum distance
  minsep = min((pix_x - state.xcurs)^2 + (pix_y - state.ycurs)^2, jmin)
  delvarx, pix_x, pix_y

  ; Update the points

  if jmin LT state.norg then begin ; ORIGINAL DATA POINT
      if state.gdpix[jmin] MOD 2 EQ 0 then $
        state.gdpix[jmin] = state.gdpix[jmin] + 1 $  ; Add it back
      else state.gdpix[jmin] = state.gdpix[jmin] - 1  ; Delete it
  endif else begin  ; Added data point
      if state.ntot EQ state.norg + 1 then state.ntot = state.norg $
      else begin
          state.gdpix[jmin] = 0 ; Data point deleted!
          gdadd = where(state.gdpix[state.norg:state.ntot-1] NE 0)
          gdadd = gdadd + state.norg
          state.ntot = state.ntot - 1
          state.xtot[state.norg:state.ntot-1] = state.xtot[gdadd]
          state.ytot[state.norg:state.ntot-1] = state.ytot[gdadd]
          state.wtot[state.norg:state.ntot-1] = state.wtot[gdadd]
          state.ivtot[state.norg:state.ntot-1] = state.ivtot[gdadd]
          state.gdpix[state.norg:state.ntot-1] = state.gdpix[gdadd]
      endelse
  endelse

  ; Update Fit
  x_continuum_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  UnDelRej
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_UnDelRej, state

  rej = where(state.gdpix MOD 8 GE 4, count)
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej]-4
  ; Update Fit
  x_continuum_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  x_continuum_GetReg, state
pro x_continuum_GetReg, state
  ; First region?
  if state.flg_reg EQ 0 then begin
      state.flg_reg = 1
      state.reg[state.nreg,0] = xgetx_plt(state, /strct)
  endif else begin
      state.flg_reg = 0
      state.reg[state.nreg,1] = xgetx_plt(state, /strct)
      state.nreg = state.nreg + 1
      x_continuum_SetReg, state
  endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;  Set Regions
;;;;

pro x_continuum_SetReg, state

;  Determines all points in the regions -- Only examines data

;  Require reg[1] > reg[0]

  if state.reg[state.nreg-1,1] LT state.reg[state.nreg-1,0] then begin
      tmp = state.reg[state.nreg-1,1]
      state.reg[state.nreg-1,1] = state.reg[state.nreg-1,0]
      state.reg[state.nreg-1,0] = tmp
  endif

; Require regions dont overlap

  if state.nreg GT 1 then begin
      bla = where(state.reg[0:state.nreg-1,*] LT state.reg[state.nreg-1,1] AND $
                  state.reg[0:state.nreg-1,*] GT state.reg[state.nreg-1,0], count)
      if count NE 0 then begin
          widget_control, state.error_msg_id, set_value='Regions may not overlap!'
          state.nreg = state.nreg - 1
          return
      endif
  endif
                  

; Adjust new pixels

  newpix = where(state.wave LE state.reg[state.nreg-1,1] AND $
                 state.wave GE state.reg[state.nreg-1,0] AND $
                 state.gdpix MOD 4 LT 2, count)
  if count NE 0 then begin
      state.gdpix[newpix] = state.gdpix[newpix] + 2 ; Bitwise
;     Reset Rejected pixels
      rej = where(state.gdpix[newpix] MOD 8 GE 4, count)
      if count NE 0 then $
        state.gdpix[newpix[rej]] = state.gdpix[newpix[rej]] - 4 ; Bitwise
  endif else begin
      widget_control, state.error_msg_id, $
        set_value='Region must contain at least 1 new data point'
      state.nreg = state.nreg - 1
  endelse
  ; Update Fit
  x_continuum_UpdateFit, state

end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
; Delete Region(s)
;;;;

pro x_continuum_DelReg, state, ALL=all

  if state.nreg EQ 0 then return
; Deletes 1 or All regions

  if (keyword_set (ALL) OR state.nreg EQ 1) then begin
      if state.nreg EQ 0 then return $
      else begin
          regpix = where(state.gdpix[0:state.norg-1] MOD 4 GT 1)
          state.gdpix[regpix] = state.gdpix[regpix] - 2
          state.nreg = 0
      endelse
  endif else begin
      xpos = xgetx_plt(state, /strct)
      fndreg = where( xpos GE state.reg[0:state.nreg-1,0] AND $
                      xpos LE state.reg[0:state.nreg-1,1], count, $
                      complement=newreg, ncomplement=nnew)
      if count EQ 0 then begin      ; Not in the region
          widget_control, state.error_msg_id, set_value='No region found!'
          return
      endif
; Reset gdpix
      fndpix = where(state.wave LE state.reg[fndreg,1] AND $
                     state.wave GE state.reg[fndreg,0])
      state.gdpix[fndpix] = state.gdpix[fndpix] - 2
; Reset regions
      for i=0,nnew-1 do begin  
          state.reg[i,0] = state.reg[newreg[i],0]
          state.reg[i,1] = state.reg[newreg[i],1]
      endfor
      state.nreg = state.nreg - 1
  endelse
  ; Update Fit
  x_continuum_UpdateFit, state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_setpnt, state

common x_continuum_fit

; 

  *state.pnt_nreg = state.nreg
  *state.pnt_reg = state.reg
  svdelpts = where(state.gdpix[0:state.norg-1] MOD 2 EQ 0)

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_SetResiduals, state

  ; Set the Flag
  if state.flg_plot MOD 2 NE 1 then state.flg_plot = state.flg_plot + 1

  ; Calculate the residuals
  x_continuum_CalcResid, state

  ; Set svxymnx
  state.res_svxymnx = [min(state.xtot[0:state.ntot-1])-$  ; xmin
                 0.02*abs(max(state.xtot[0:state.ntot-1])-$
                          min(state.xtot[0:state.ntot-1])), $
                 min(state.residuals[0:state.ntot-1])-$
                 0.02*abs(max(state.residuals[0:state.ntot-1])-$
                          min(state.residuals[0:state.ntot-1])), $ ; ymin
                 max(state.xtot[0:state.ntot-1])+$
                 0.02*abs(max(state.xtot[0:state.ntot-1])-$
                          min(state.xtot[0:state.ntot-1])), $
                 max(state.residuals[0:state.ntot-1])+$
                 0.02*abs(max(state.residuals[0:state.ntot-1])-$
                          min(state.residuals[0:state.ntot-1]))]

  state.svxymnx = state.res_svxymnx
  state.xymnx = state.svxymnx

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;

pro x_continuum_CalcResid, state

  common x_continuum_fit

  ; Calculate the residuals
  fit = x_calcfit(state.xtot[0:state.ntot-1], FITSTR=fitprm)
  state.residuals[0:state.ntot-1] = state.ytot[0:state.ntot-1] - fit

  ; Reset screen
  if state.flg_plot MOD 2 EQ 1 then begin
      state.res_svxymnx[1] = min(state.residuals[0:state.ntot-1])-$
        0.02*abs(max(state.residuals[0:state.ntot-1])-$
                 min(state.residuals[0:state.ntot-1])) ; ymin
      state.res_svxymnx[3] = max(state.residuals[0:state.ntot-1])+$
        0.02*abs(max(state.residuals[0:state.ntot-1])-$
                 min(state.residuals[0:state.ntot-1]))
      state.svxymnx = state.res_svxymnx
      state.xymnx = state.svxymnx
  endif

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Reset screen to data
;
pro x_continuum_SetValues, state

  ; Set the Flag
  if state.flg_plot MOD 2 NE 0 then state.flg_plot = state.flg_plot - 1

  ; svxymnx
  state.svxymnx = state.dat_svxymnx
  state.xymnx = state.svxymnx

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Clear all Rejected Points
;
pro x_continuum_ClearRej, state

  ; Rej
  rej = where(state.gdpix[0:state.ntot-1] MOD 8 GE 4, count)

  ; Clear
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej] - 4

  ; Update Fit
  x_continuum_UpdateFit, state

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Save fit to continuum
;
pro x_continuum_SaveCont, state, ALL=all

common x_continuum_fit

  if keyword_set( ALL ) then begin
      fc_val = fin_fit 
      ; Reset the fit
      x_continuum_Reset, state
  endif else begin
      if state.flg_svcont EQ 0 then begin
          state.flg_svcont = 1
          state.tmpxy[0] = xgetx_plt(state, /strct)
      endif else begin
          state.flg_svcont = 0
          x2 = xgetx_plt(state, /strct)
          xmn = x2 < state.tmpxy[0]
          xmx = x2 > state.tmpxy[0]
          ; Find all xdat points
          gdpt = where(xfit LT xmx AND xfit GT xmn, ngd)
          if ngd NE 0 then fc_val[gdpt] = fin_fit[gdpt]
          ; Reset the fit
          x_continuum_Reset, state
      endelse
  endelse
  return
end

;;;;;;;;;;;;;;;;;;;;
; Add continuum to fit
;
pro x_continuum_AddCont, state

common x_continuum_fit

  ; First time?
  if state.flg_addcont EQ 0 then begin
      state.flg_addcont = 1
      state.tmpxy[0] = xgetx_plt(state, /strct)
  endif else begin
      state.flg_addcont = 0
      x2 = xgetx_plt(state, /strct)
      xmn = x2 < state.tmpxy[0]
      xmx = x2 > state.tmpxy[0]
          ; Find all xdat points
      gdpt = where(xfit LT xmx AND xfit GT xmn, ngd)
      if ngd NE 0 then begin 
          state.xtot[state.ntot:state.ntot+ngd-1] = xfit[gdpt]
          state.ytot[state.ntot:state.ntot+ngd-1] = fc_val[gdpt]
          state.wtot[state.ntot:state.ntot+ngd-1] = 1./sqrt(state.med_ivar)/10.
          state.ivtot[state.ntot:state.ntot+ngd-1] = 100*state.med_ivar ; Inverse var
          state.gdpix[state.ntot:state.ntot+ngd-1] = 1
          state.ntot = state.ntot + ngd
      endif
  endelse
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

pro x_continuum, yin, ysin, FITSTR=fitstr, CONTI=conti, LSIG=lsig, $
                 XSIZE=xsize, YSIZE=ysize, INFLG=inflg, OUTFIL=outfil

common x_continuum_fit


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_continuum, ydat, [ysig], FITSTR=, INFLG=, OUTFIL= [v1.0]'
    return
  endif 

  if not keyword_set( LSIG ) then lsig = 2.5
  if not keyword_set( INFLG ) then inflg = 0
  if not keyword_set( OUTFIL ) then begin
      if keyword_set( CONTI ) then outfil = conti else outfil = 'conti.fits'
  endif

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0
  ydat = x_readspec(yin, INFLG=inflg, /dscale, head=head, NPIX=npix, $
                    WAV=xdat, FIL_SIG=ysin, SIG=ysig)

; CONTINUUM
  if keyword_set( CONTI ) then begin
      cdat = x_readspec(conti, /fscale)
      if n_elements(cdat) NE n_elements(xdat) then begin
          print, 'x_continuum: Continuum is the wrong size!'
          return
      endif
      fc_val = cdat
  endif 

; XFIT
  xfit = xdat

; Set inv variance
  gdsig = where( ysig GT 0., ngd)
  if ngd EQ 0 then begin
      ysig[*] = 1.
      gdsig = lindgen(npix)
  endif
  ivar = dblarr(npix)
  ivar[gdsig] = 1./(ysig[gdsig])^2
  med_ivar = median(ivar[gdsig])

; Initialize the common block
  x_continuum_initcommon

  ; Set fit structure as input
  if keyword_set( FITSTR ) then begin
      fitprm = fitstr 
  endif
  fitprm.lsig = lsig

; Screen

  if not keyword_set( XSIZE ) then    xsize = 1300
  if not keyword_set( YSIZE ) then    ysize = 1000


;    Pointer for Final fit
  tmpreg = fltarr(100,2)
  pnt_reg = PTR_NEW( tmpreg )
  nreg = 0
  pnt_nreg = PTR_NEW( nreg )
  
;    STATE
  state = { $
            wave: xdat, $
            fx: ydat, $
            sig: ysig, $
            ivar: ivar, $
            med_ivar: med_ivar, $
            norg: n_elements(xdat), $
            pnt_nreg: pnt_nreg, $
            pnt_reg: pnt_reg, $
            residuals: dblarr(200000), $  ; Residuals
            flg_plot: 0, $
            xtot: dblarr(200000), $ ; Fitting vectors
            ytot: dblarr(200000), $
            wtot: dblarr(200000), $ ; Weighting Factor
            ivtot: dblarr(200000), $ ; Weighting Factor
            flg_wgt: 1, $       ; 0 = No input 1 = Input
            ntot: n_elements(xdat), $
            nreg: 0, $
            reg: tmpreg, $
            outfil: outfil, $  ; Output file
            flg_reg: 0, $
            flg_zoom: 0, $
            flg_svcont: 0, $
            flg_addcont: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            res_svxymnx: fltarr(4), $  ; Residuals
            dat_svxymnx: [min(xdat)-0.01*abs(max(xdat)-min(xdat)), $
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $
                      max(xdat)+0.01*abs(max(xdat)-min(xdat)), $
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $
            svsvxymnx: fltarr(4), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            size: intarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L, $
            reject_base_id: 0L, $
            hsig_id: 0L, $
            lsig_id: 0L, $
            rms_id: 0L, $
            gdpix: intarr(200000) $
          }
; SVXYMNX
  state.svxymnx = state.dat_svxymnx
  state.svsvxymnx = state.dat_svxymnx


;   REG
  if keyword_set( REG ) then begin
      sz = size(reg, /dimensions)
      state.reg[0:sz[0]-1,*] = reg
      state.nreg = sz[0]
  endif
      
;    WIDGET
  base = WIDGET_BASE( title = 'x_continuum: Interactive Mode', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
  strlbl = strarr(10)
  strlbl = ['x_continuum', ' ', 'Ver 1.1']
  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)
;        Function
  funcval = ['POLY', 'LEGEND', 'CHEBY', 'BSPLIN', 'GAUSS']
  state.funclist_id = WIDGET_LIST(toolbar, VALUE=funcval, uvalue='FUNCLIST', $
                                  ysize=4)
  case strtrim(fitprm.func) of 
      'POLY' :  widget_control, state.funclist_id, set_list_select=0
      'LEGEND' :  widget_control, state.funclist_id, set_list_select=1
      'CHEBY' :  widget_control, state.funclist_id, set_list_select=2
      'BSPLIN' :  widget_control, state.funclist_id, set_list_select=3
      'GAUSS' :  widget_control, state.funclist_id, set_list_select=4
  endcase
  
;        Order      
  ordrbase = WIDGET_BASE(toolbar, /column)
  upordr = WIDGET_BUTTON(ordrbase, value='Up',uvalue='UP')
  tmp_string = string(fitprm.nord, format='(i4)')
  state.lblordr_id = WIDGET_LABEL(ordrbase, value=tmp_string, $
                                  /dynamic_resize, /align_center)
  dwnordr = WIDGET_BUTTON(ordrbase, value='Down',uvalue='DOWN')
  
;  RMS
  state.rms_id = cw_field(toolbar, value='0.', /floating, title='RMS')

;        Rejection
  state.reject_base_id = widget_base(toolbar, /column, /align_center,$
                                     /frame)
  hsigbase = widget_base(state.reject_base_id, /row, /align_right)
  state.hsig_id = cw_field(hsigbase, title='High Sig', $
                           value=fitprm.hsig, $
                           xsize=5, ysize=1, /return_events, $
                           uvalue='HSIGVAL')
  lsigbase = widget_base(state.reject_base_id, /row, /align_right)
  state.lsig_id = cw_field(lsigbase, title='Low Sig', $
                           value=fitprm.lsig, $
                           xsize=5, ysize=1, /return_events, $
                           uvalue='LSIGVAL')
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB -- Add data point',$
             'RMB -- (Un)Delete data point',$
             'l,r,b,t -- Set window edge', $
             'w -- Reset screen',$
             '[]{} -- Pan', $
             'F -- Save screen area', $
             '$ -- Reset screen view',$
             '------------------', $
             's,s -- Set Region',$
             'g,g -- Save piece of continuum',$
             'S -- Save entire continuum',$
             'x -- Remove one Region',$
             'D -- Delete data point', $
             'c,c -- Add continuum points to the fit', $
             '--------------------', $
             'f -- Update fit', $
             'u,d -- In(de)crease order 1',$
             'C -- Clear all rejected ponts', $
             '! -- Continuum reset', $
             '--------------------', $
             'X -- Reset fit',$
             'R -- Plot Residuals',$
             'V -- Plot Values',$
             'U -- Write out to file', $
             'q -- Quit and save']
  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
                             /scroll)
;      Error base
  error_base_id  = widget_base(toolbar, /column)
  error_btn_id = widget_button(error_base_id, value='Erase Error Msg', $
                               uvalue='ERRORB')
;          Error Messages
  state.error_msg_id = widget_text(error_base_id, value='', xsize=40, $
                                   ysize=2, /scroll)
  
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)
  
  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
;      Done
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Include all elements at first
  state.gdpix[0:n_elements(xdat)] = 1
  
; Color Table 

; Update
  x_continuum_Reset, state
; Continuum
  if not keyword_set( fc_val ) then fc_val = fltarr(state.norg) + 1.
  x_continuum_UpdateFit, state
  x_continuum_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
; xmanager blows (BLOCK)
  repeat begin
      ev = widget_event(base)
      ans = x_continuum_ev(ev)
  end until ans EQ 0
  
; Passing back
  
;    Fit parameters
  
  hsig = fitprm.hsig
  lsig = fitprm.lsig
  func = fitprm.func
  nord = fitprm.nord
      
;    Regions
  if arg_present( REG ) then begin
      if keyword_set( INTER ) then begin
          nreg = *pnt_nreg
          reg = *pnt_reg
          delvarx, tmpreg
      endif
      if nreg GT 0 then tmpreg = reg[0:nreg-1,*] $
      else begin
          tmpreg = fltarr(1,2)
          tmpreg[0,0] = xdat[0]
          tmpreg[0,1] = xdat[n_elements(xdat)-1]
      endelse
      delvarx, reg
      reg = temporary(tmpreg)
  endif 

; Free the pointers
  if keyword_set( pnt_reg ) then PTR_FREE, pnt_reg
  if keyword_set( pnt_inreg ) then PTR_FREE, pnt_nreg
  if keyword_set( tmpreg ) then delvarx, tmpreg

;  End game

  ; Output
  mwrfits, float(fc_val), outfil, /create

  if arg_present( FITSTR ) then fitstr=fitprm
  if arg_present( FFIT ) then ffit = *fitprm.ffit
  if arg_present( RMS ) then rms = fitprm.rms
  if arg_present( DELPTS ) then delpts = svdelpts
  if arg_present( NRM ) then nrm = fitprm.nrm
  if arg_present( CONTI ) then conti = temporary(fc_val) else delvarx, fc_val


; Memory
  delvarx, fitprm, xfit

  return
end
