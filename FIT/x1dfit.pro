;+ 
; NAME:
; x1dfit   
;   Version 1.3
;
; PURPOSE:
;    GUI for fitting a function to a set of x,y data.  This
;    program does far too much!  Warning:  This was one of my first
;    GUIs.  It works pretty well, but is a bit old.
;
; CALLING SEQUENCE:
;   
;   fit = x1dfit([xdat],ydat,func=,nord=, /inter, xsize=, ysize=)
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
;   fit = x1dfit(x, y, 'POLY', nord=5, /inter)
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

pro x1dfit_initcommon

;

common x1dfit_fit, fin_fit, fitprm, xfit, svdelpts

fitprm = x_setfitstrct()

end

;;;;
; Events
;;;;

function x1dfit_ev, ev
;pro x1dfit_event, ev

common x1dfit_fit 

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
          x1dfit_UpdateFit, state
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
          return, 1
      end
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Add a Data point with weight 10
                          state.xtot[state.ntot] = xgetx_plt(state,/strct)
                          state.ytot[state.ntot] = xgety_plt(state,/strct)
                          state.wtot[state.ntot] = 0.1
                          state.ivtot[state.ntot] = 10
                          state.gdpix[state.ntot] = 1
                          state.ntot = state.ntot + 1
                          x1dfit_UpdateSVXY, state
                          x1dfit_setxfit, state
                          x1dfit_UpdateFit, state
                      end 
                      4 : x1dfit_Delete, state ; Delete/Undelete a data point
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
          case eventch of
              'u': begin
                  fitprm.nord = fitprm.nord+1
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm.nord, format='(i4)')
                  x1dfit_UpdateFit, state
              end
              'd': begin
                  fitprm.nord = (fitprm.nord-1 > 0)
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm.nord, format='(i4)')
                  x1dfit_UpdateFit, state
              end
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'z': begin  ; Zoom
                  ximgd_setzoom, state, /plot
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif
              end
              's': begin  ; Region
                  x1dfit_GetReg, state
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif 
              end
              'W': begin
                  if state.flg_plot MOD 2 EQ 0 then $
                    state.xymnx = state.svxymnx $ ; Reset the screen
                  else x1dfit_SetResiduals, state
              end
              'D': x1dfit_Delete, state ; Delete/Undelete a data point
              'x': x1dfit_DelReg, state        ; Delete one region
              'S': x1dfit_DelReg, state, /all  ; Delete all regions and reset
              'R': x1dfit_SetResiduals, state  ; Plot residuals
              'V': x1dfit_SetValues, state     ; Plot values
              'C': x1dfit_ClearRej, state ; Clear all rejected points
              'q': begin
                  x1dfit_setpnt, state
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
          x1dfit_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.lsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'LSIGVAL': begin
          fitprm.lsig = ev.value
          if fitprm.niter EQ 0 then fitprm.niter=3
          x1dfit_UnDelRej, state
          if ev.value NE 0. then fitprm.flg_rej = 1 else begin
              if fitprm.hsig EQ 0. then fitprm.flg_rej = 0
          endelse
      end
      'DONE' : begin
          x1dfit_setpnt, state
          widget_control, ev.top, /destroy
          return, 0
      end
      else:
  endcase

; Update Plot
  x1dfit_UpdatePlot, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return, 1
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x1dfit_UpdatePlot, state
  
common x1dfit_fit

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  color = getcolor(/load)

  if state.flg_plot MOD 2 EQ 0 then begin  ; NORMAL

;  ORIGINAL DATA NOT IN REGION
      gdorg = where(state.gdpix[0:state.norg-1] EQ 1, count) ; Original data 
      if count NE 0 then begin
          plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
            /nodata, xstyle=1, ystyle=1, $
            position=state.pos,  background=color.black, color=color.white 
          oplot, state.xdat[gdorg], state.ydat[gdorg], psym=1, $
            color=color.green 
      endif else $
        plot, [-1.0],  [state.svxymnx[3]+1.], psym=1, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        background=color.black, color=color.white

      for q=0,7 do begin
          if q EQ 1 then continue
          pts = where(state.gdpix[0:state.norg-1] EQ q, count)
          if count NE 0 then begin
              case q of
                  0: oplot, [state.xdat[pts]], [state.ydat[pts]], psym=2, $
                    color=color.red ; DELETED DATA OUT OF REGION
                  2: oplot, [state.xdat[pts]], [state.ydat[pts]], psym=2, $
                    color=color.purple ; DELETED DATA IN REGION
                  3: oplot, [state.xdat[pts]], [state.ydat[pts]], psym=1, $
                    color=color.blue ; GOOD DATA IN REGION
                  4: 
                  5: oplot, [state.xdat[pts]], [state.ydat[pts]], psym=5, $
                    color=color.pink ; REJECTED DATA OUT OF REGION
                  6: 
                  7: oplot, [state.xdat[pts]], [state.ydat[pts]], psym=5, $
                    color=color.yellow ; REJECTED DATA OUT OF REGION
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
      oplot, xfit, fin_fit, color=color.red ; Best fit Red  

  endif else begin  ; RESIDUALS

;  ORIGINAL DATA NOT IN REGION
      gdorg = where(state.gdpix[0:state.norg-1] EQ 1, count) ; Original data 
      if count NE 0 then begin
          plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
            /nodata, xstyle=1, ystyle=1, $
            position=state.pos,  background=color.black, color=color.white 
          oplot, state.xdat[gdorg], state.residuals[gdorg], psym=1, $
            color=color.green 
      endif else $
        plot, [-1.0],  [state.svxymnx[3]+1.], psym=1, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        background=color.black, color=color.white

      for q=0,7 do begin
          if q EQ 1 then continue
          pts = where(state.gdpix[0:state.norg-1] EQ q, count)
          if count NE 0 then begin
              case q of
                  0: oplot, [state.xdat[pts]], [state.residuals[pts]], psym=2, $
                    color=color.red ; DELETED DATA OUT OF REGION
                  2: oplot, [state.xdat[pts]], [state.residuals[pts]], psym=2, $
                    color=color.purple ; DELETED DATA IN REGION
                  3: oplot, [state.xdat[pts]], [state.residuals[pts]], psym=1, $
                    color=color.blue ; GOOD DATA IN REGION
                  4: 
                  5: oplot, [state.xdat[pts]], [state.residuals[pts]], psym=5, $
                    color=color.pink ; REJECTED DATA OUT OF REGION
                  6: 
                  7: oplot, [state.xdat[pts]], [state.residuals[pts]], psym=5, $
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

pro x1dfit_UpdateFit, state

common x1dfit_fit

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

  mask = lonarr(state.ntot)
  mask[gdpt] = 1
  fitprm.minpt = fitprm.minpt < count/2

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
          fit = x_fit(state.xtot[0:state.ntot-1], state.ytot[0:state.ntot-1], $
                      FITSTR=fitprm, $
                      MSK=mask,$
                      SIG=state.wtot[0:state.ntot-1], $
                      IVAR=state.ivtot[0:state.ntot-1] )
          if fit[0] EQ -1 then begin
              fin_fit = fltarr(n_elements(xfit))
              widget_control, state.error_msg_id, $
                set_value='Bad fit! Adjust nord probably'
              return
          endif
      endif else begin ; FIT with REJECTION!
          ;rejpt = -1
          if state.flg_var EQ 0 then begin
             fit = x_fitrej(state.xtot[0:state.ntot-1], $
                            state.ytot[0:state.ntot-1], $
                            MSK=mask,$
                            FITSTR=fitprm, $
                            REJPT=rejpt )
          endif else begin
             fit = x_fitrej(state.xtot[0:state.ntot-1], $
                            state.ytot[0:state.ntot-1], $
                            SIG=state.wtot[0:state.ntot-1], $
                            IVAR=state.ivtot[0:state.ntot-1], $ 
                            MSK=mask,$
                            FITSTR=fitprm, $
                            REJPT=rejpt )
          endelse
          if keyword_set(REJPT) then begin
             if rejpt[0] NE -1 then state.gdpix[rejpt] = state.gdpix[rejpt]+4
          endif
      endelse
  endelse

  ; Calculate values
  fin_fit = x_calcfit(xfit, FITSTR=fitprm)

  delvarx, mask, fit

  ; RESIDUALS AND RMS
  x1dfit_CalcResid, state

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

  fitprm.rms = sqrt(total( state.residuals[gdpt]^2 ) / $
                    float(n_elements(gdpt)-1))
  widget_control, state.rms_id, set_value=fitprm.rms

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x1dfit_Reset, state

;
  state.xtot[0:n_elements(state.xdat)-1] = state.xdat
  state.ytot[0:n_elements(state.xdat)-1] = state.ydat
  if state.flg_var EQ 0 then begin
      state.wtot[0:n_elements(state.xdat)-1] = 1.0 
      state.ivtot[0:n_elements(state.xdat)-1] = 1.0
  endif else begin
      state.wtot[0:n_elements(state.xdat)-1] = state.sig
      state.ivtot[0:n_elements(state.xdat)-1] = state.ivar
  endelse

  state.ntot = n_elements(state.xdat)


; gdpix Flag
;   1 = Good
;   2 = In region
;   4 = Rejected by fit
  state.gdpix[0:n_elements(state.xdat)-1] = 1

; Mask
  a = where(state.mask[0:n_elements(state.xdat)] EQ 0, na)
  if na NE 0 then state.gdpix[a] = 0

; Regions

  if state.nreg NE 0 then begin
      for i=0,state.nreg-1 do begin
          newpix = where(state.xdat LE state.reg[i,1] AND $
                         state.xdat GE state.reg[i,0], count)
          if count NE 0 then $
            state.gdpix[newpix] = state.gdpix[newpix] + 2 ; Bitwise flag
      endfor
  endif


; Plotting

  state.xymnx = state.svxymnx


end

;;;;;;;;;;;;;;;;;;;;
;  Delete
;;;;;;;;;;;;;;;;;;;;

pro x1dfit_Delete, state

; Deletes the nearest data point in pixel space

  ; Get into pixel space
  pix_x = xgetxpix_plt(state.xtot[0:state.ntot-1], state.pos, state.xymnx, state.size)
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

  ; Update fit
  x1dfit_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;
;  UnDelRej
;;;;;;;;;;;;;;;;;;;;

pro x1dfit_UnDelRej, state

  rej = where(state.gdpix MOD 8 GE 4, count)
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej]-4
  ; Update fit
  x1dfit_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  x1dfit_GetReg, state
pro x1dfit_GetReg, state
  ; First region?
  if state.flg_reg EQ 0 then begin
      state.flg_reg = 1
      state.reg[state.nreg,0] = xgetx_plt(state, /strct)
  endif else begin
      state.flg_reg = 0
      state.reg[state.nreg,1] = xgetx_plt(state, /strct)
      state.nreg = state.nreg + 1
      x1dfit_SetReg, state
  endelse
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;  Set Regions
;;;;

pro x1dfit_SetReg, state

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

  newpix = where(state.xdat LE state.reg[state.nreg-1,1] AND $
                 state.xdat GE state.reg[state.nreg-1,0] AND $
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

; Update fit
  x1dfit_UpdateFit, state

end
  
;;;;
; Delete Region(s)
;;;;

pro x1dfit_DelReg, state, ALL=all

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
      xpos = xgetx_plt(state.xcurs, state.pos, state.xymnx, state.size)
      fndreg = where( xpos GE state.reg[0:state.nreg-1,0] AND $
                      xpos LE state.reg[0:state.nreg-1,1], count, $
                      complement=newreg, ncomplement=nnew)
      if count EQ 0 then begin      ; Not in the region
          widget_control, state.error_msg_id, set_value='No region found!'
          return
      endif
; Reset gdpix
      fndpix = where(state.xdat LE state.reg[fndreg,1] AND $
                     state.xdat GE state.reg[fndreg,0], nfnd)
      if nfnd NE 0 then $
        state.gdpix[fndpix] = state.gdpix[fndpix] - 2
; Reset regions
      for i=0,nnew-1 do begin  
          state.reg[i,0] = state.reg[newreg[i],0]
          state.reg[i,1] = state.reg[newreg[i],1]
      endfor
      state.nreg = state.nreg - 1
  endelse

; Update fit
  x1dfit_UpdateFit, state
end

;;;;;;;;;;;;;;;;;;;;

pro x1dfit_setpnt, state

common x1dfit_fit

; 

  *state.pnt_nreg = state.nreg
  *state.pnt_reg = state.reg
  svdelpts = where((state.gdpix[0:state.norg-1] MOD 2 EQ 0) OR $
                   (state.gdpix[0:state.norg-1] MOD 8 GE 4))
                   

return
end

;;;;;;;;;;;;;;;;;;;;

pro x1dfit_SetResiduals, state

  ; Set the Flag
  if state.flg_plot MOD 2 NE 1 then state.flg_plot = state.flg_plot + 1

  ; Calculate the residuals
  x1dfit_CalcResid, state

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

;;;;;;;;;;;;;;;;;;;;

pro x1dfit_CalcResid, state

  common x1dfit_fit

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
;      state.xymnx = state.svxymnx
  endif

return
end

;;;;;;;;;;;;;;;;;;;;
; Reset screen to data
;
pro x1dfit_SetValues, state

  ; Set the Flag
  if state.flg_plot MOD 2 NE 0 then state.flg_plot = state.flg_plot - 1

  ; svxymnx
  state.svxymnx = state.dat_svxymnx
  state.xymnx = state.svxymnx

return
end

;;;;;;;;;;;;;;;;;;;;
; Clear all Rejected Points
;
pro x1dfit_ClearRej, state

  ; Rej
  rej = where(state.gdpix[0:state.ntot-1] MOD 8 GE 4, count)

  ; Clear
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej] - 4

  ; Update Fit
  x1dfit_UpdateFit, state

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Set xfit

pro x1dfit_setxfit, state
common x1dfit_fit

  xfit = findgen(3*state.norg)*(state.svxymnx[2]-state.svxymnx[0])/ $
    float(3*state.norg)  + state.svxymnx[0] 
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Reset SVXY
pro x1dfit_UpdateSVXY, state
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

function x1dfit, xin, yin, FUNC=func, NORD=nord, INTER=INTER, $
                 XSIZE=xsize, YSIZE=ysize, SIG=sig, REG=reg, FFIT=ffit,$
                 REJ=rej, LSIG=lsig, HSIG=hsig, MINPT=minpt, RMS=rms, $
                 DELPTS=delpts, NRM=nrm, FITSTR=fitstr, IVAR=ivar, MSK=msk

common x1dfit_fit


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'fit = x1dfit([xdat], ydat, FUNC=, NORD=, /INTER, '
    print, '        XSIZE=,YSIZE=, SIG=, REG=, FFIT='
    print, '        REJ=, LSIG=, HSIG=, MINPT=, RMS=, DELPTS=, NRM=, '
    print, '        FITSTR=, MSK=, IVAR= ) [V1.3]'
    return, -1
  endif 

; Set xdat, ydat

  if keyword_set( yin ) then begin
      if n_elements(xin) NE n_elements(yin) then begin
          print, 'x1dfit: Wrong array sizes!'
          return, -1
      endif
      xdat = xin
      ydat = yin
  endif else begin
      xdat = findgen( n_elements(xin) )
      ydat = xin
  endelse

; Inv variance
  if keyword_set(SIG) AND not keyword_set( IVAR ) then begin
      gdsig = where( sig GT 0.)
      ivar = fltarr(n_elements(sig))
      ivar[gdsig] = 1./(sig[gdsig])^2
  endif

; Sigma
  if not keyword_set(SIG) AND keyword_set( IVAR ) then begin
      gd = where( ivar GT 0.)
      sig = fltarr(n_elements(ivar))
      sig[gd] = 1./sqrt(ivar[gd])
  endif
      
; Initialize the common block

  x1dfit_initcommon

  ; Set fit structure as input
  if keyword_set( FITSTR ) then fitprm = fitstr else begin

;  Optional Keywords
      ; FUNC
      if not keyword_set( FUNC ) then fitprm.func = 'POLY' else fitprm.func=func
      ; NORD
      if not keyword_set( NORD ) then fitprm.nord = 3 else fitprm.nord = nord
      if fitprm.nord GE n_elements(xdat)-1 then begin
          print, 'Resetting nord below nelements(xdat)'
          if fitprm.func EQ 'POLY' then fitprm.nord = n_elements(xdat)-2 $
          else fitprm.nord = n_elements(xdat)-2 
      endif
      if keyword_set( REJ ) then begin
          fitprm.lsig = 3.
          fitprm.hsig = 3.
          fitprm.flg_rej = 1
          fitprm.niter = 3
          fitprm.maxrej = 50L
          fitprm.minpt = 10L
      endif
      if not keyword_set( MINPT ) then fitprm.minpt = long(0.5*n_elements(xin))
      if keyword_set( LSIG ) then fitprm.lsig = lsig
      if keyword_set( HSIG ) then fitprm.hsig = hsig
  endelse

;   REJECTION
  if keyword_set(fitprm.lsig) OR keyword_set( fitprm.hsig ) then $ 
    fitprm.flg_rej = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Non-interactive

  if not keyword_set (INTER) then begin

      if fitprm.flg_rej EQ 0 then $
        fit= x_fit(xdat, ydat, FITSTR=fitprm, $
                      SIG=sig, IVAR=ivar, REG=reg) $
      else $
        fit = x_fitrej(xdat, ydat, FITSTR=fitprm, $
                         SIG=sig, REG=reg, IVAR=ivar, RMS=rms) 

  endif else begin
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; INTERACTIVE

  device, get_screen_size=ssz
;  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
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

;    Pointer for Final fit
  tmpreg = fltarr(100,2)
  pnt_reg = PTR_NEW( tmpreg )
  nreg = 0
  pnt_nreg = PTR_NEW( nreg )
  
;    STATE
  state = { $
            xdat: xdat, $
            ydat: ydat, $
            sig: fltarr(n_elements(ydat)), $
            ivar: fltarr(n_elements(ydat)), $
            flg_var: 0, $           ;  0=Nothing, 1=Sig, 2=Var
            norg: n_elements(xdat), $
            pnt_nreg: pnt_nreg, $
            pnt_reg: pnt_reg, $
            residuals: fltarr(200000), $  ; Residuals
            flg_plot: 0, $
            xtot: fltarr(200000), $ ; Fitting vectors
            ytot: fltarr(200000), $
            wtot: fltarr(200000), $ ; Weighting Factor
            ivtot: fltarr(200000), $ ; Weighting Factor
            ntot: n_elements(xdat), $
            nreg: 0, $
            reg: tmpreg, $
            mask: bytarr(200000), $
            flg_reg: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            res_svxymnx: fltarr(4), $  ; Residuals
            dat_svxymnx: [min(xdat)-0.01*abs(max(xdat)-min(xdat)), $
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $
                      max(xdat)+0.01*abs(max(xdat)-min(xdat)), $
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $
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

; Include all elements at first
  state.gdpix[0:n_elements(xdat)] = 1

; SVXYMNX
  state.svxymnx = state.dat_svxymnx

; XFIT
  
  x1dfit_setxfit, state

;   SIG

  if keyword_set( SIG ) then begin
      state.sig = sig
      state.flg_var = state.flg_var + 1
  endif
  if keyword_set( IVAR ) then begin
      state.ivar = ivar
      state.flg_var = state.flg_var + 2
  endif

;   REG

  if keyword_set( REG ) then begin
      sz = size(reg, /dimensions)
      state.reg[0:sz[0]-1,*] = reg
      state.nreg = sz[0]
  endif

; MSK
  if keyword_set( MSK ) then begin
      a = where(msk EQ 1)
      state.mask[a] = 1
  endif else state.mask[*] = 1
      
      
;    WIDGET
  base = WIDGET_BASE( title = 'x1dfit: Interactive Mode', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
  strlbl = strarr(10)
  strlbl = ['x1dfit', ' ', 'Ver 1.1']
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
             'l -- Set Left Window',$
             'r -- Set Right Window',$
             'b -- Set Bottom Window',$
             't -- Set Top Window',$
             'z,z -- Zoom in on area', $
             'W -- Reset screen',$
             's,s -- Set Region',$
             'x -- Remove one Region',$
             'D -- Delete data point', $
             'S -- Remove all Regions',$
             'q -- Quit and save',$
             'u -- Increase order 1',$
             'd -- Decrease order 1',$
             'R -- Plot Residuals',$
             'V -- Plot Values',$
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
  
  
; Color Table 

; Update
  x1dfit_Reset, state
  x1dfit_UpdateFit, state
  x1dfit_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
; xmanager blows (BLOCK)
  repeat begin
      ev = widget_event(base)
      ans = x1dfit_ev(ev)
  end until ans EQ 0
  
endelse

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

  if arg_present( FITSTR ) then fitstr=fitprm
  if arg_present( FFIT ) then ffit = *fitprm.ffit
  if arg_present( RMS ) then rms = fitprm.rms
  if arg_present( DELPTS ) then delpts = svdelpts
  if arg_present( NRM ) then nrm = fitprm.nrm

;  Reset fin_fit to the xin points

  return, x_calcfit(xdat, FITSTR=temporary(fitprm))


  return, -1
end
