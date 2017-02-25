;+ 
; NAME:
; x_spec_zerolvl   
;   Version 1.1
;
; PURPOSE:
;    GUI used to fit the spec_zerolvl of a spectrum (usually quasars).
;    The user can dictate a SPLINE or perform a minimum chi^2 fit to
;    regions with rejection.
;
; CALLING SEQUENCE:
;   x_spec_zerolvl, flux_fil, error_fil, FITSTR=, CONTI=, LSIG=, XSIZE=, YSIZE=,
;   INFLG=, OUTFIL=, /SPLINE, INISPL=
;
; INPUTS:
;   flux_fil   - Values along one dimension [optional]
;   error_fil  - Values along the other
;
; RETURNS:
;
; OUTPUTS:
;  OUTFIL=  -- Name of FITS file to save spec_zerolvl
;
; OPTIONAL KEYWORDS:
;   CONTI=  -- Name of FITS file containing previously saved
;                spec_zerolvl
;   INFLG   -- 0 = yin, ysin as data arrays (fits allowed and
;                expected)
;                1 = One fits file (flux, sig)
;                2 = One fits file (flux, sig, wave)
;                3 = FUSE format
;                5 = SDSS file
;   XSIZE=  -- Window size in pixels
;   YSIZE=  -- Window size in pixels
;   LSIG=   -- Value for sigma rejection on low side [default: 2.5]
;   INISPL= -- IDL file containing a saved SPLINE for the spec_zerolvl
;   FITSTR  -- Fit structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_spec_zerolvl, 'Q0000.fits', OUTFIL='Q0000_c.fits', /SPLINE
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

pro x_spec_zerolvl_initcommon

;

common x_spec_zerolvl_fit, fin_fit, fitprm, xfit, svdelpts, fc_val, $
  spec_2d

tmp = { fitstrct }
fitprm = replicate(tmp, 100)
fitprm.lsig = 2.5
fitprm.hsig = 3.
fitprm.niter = 1
fitprm.nord = 1
fitprm.flg_rej = 1
fitprm.maxrej = 100
fitprm.func = 'LEGEND'

end

;;;;
; Events
;;;;

function x_spec_zerolvl_ev, ev
;pro x_spec_zerolvl_event, ev

common x_spec_zerolvl_fit

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  flg_plt = 1
  case uval of
      'SPLINE': begin
          state.flg_fit = ev.value
          x_spec_zerolvl_UpdateFit, state
      end
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'FUNCLIST' : begin
          nfunc = ev.index
          case nfunc of 
              0 : fitprm[state.curfit].func = 'POLY'
              1 : fitprm[state.curfit].func = 'LEGEND'
              2 : fitprm[state.curfit].func = 'CHEBY'
              3 : fitprm[state.curfit].func = 'BSPLIN'
              4 : fitprm[state.curfit].func = 'GAUSS'
              else:
          endcase
          ; Reset rejected
          rej = where(state.gdpix[0:state.norg-1] MOD 8 GE 4, nrej)
          if nrej NE 0 then state.gdpix[rej] = state.gdpix[rej] - 4
          ; Update fit
          x_spec_zerolvl_UpdateFit, state
      end
      'UP' : begin
          fitprm[state.curfit].nord = fitprm[state.curfit].nord+1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm[state.curfit].nord, format='(i4)')
          ; Update fit
          x_spec_zerolvl_UpdateFit, state
      end
      'DOWN' : begin
          fitprm[state.curfit].nord = (fitprm[state.curfit].nord-1) > 1
          widget_control, state.lblordr_id, $
            set_value=string(fitprm[state.curfit].nord, format='(i4)')
          ; Update fit
          x_spec_zerolvl_UpdateFit, state
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
                          state.ivtot[state.ntot] = 30*state.med_ivar 

                          state.gdpix[state.ntot] = 8
                          state.ntot = state.ntot + 1
                          ; Update Fit
                          x_spec_zerolvl_UpdateFit, state
                      end 
                      4 : x_spec_zerolvl_Delete, state ; Delete/Undelete a data point
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
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
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
                  fitprm[state.curfit].nord = fitprm[state.curfit].nord+1
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm[state.curfit].nord, format='(i4)')
                  flg_plt = 0
              end
              'd': begin
                  fitprm[state.curfit].nord = (fitprm[state.curfit].nord-1) > 1
                  widget_control, state.lblordr_id, $
                    set_value=string(fitprm[state.curfit].nord, format='(i4)')
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
                  x_spec_zerolvl_GetReg, state
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif 
              end
              '0': begin
                 state.zro_lvl[*] = 0. ;; Set zero to 0
                 if keyword_set(FIN_FIT) then fin_fit[*] = 0. $
                 else fin_fit = replicate(0., n_elements(state.fx))
              end
              'A': state.flg_apply = (state.flg_apply + 1) MOD 2
              'F': state.svxymnx = state.xymnx ; Save screen area to svxymn
              'f': x_spec_zerolvl_UpdateFit, state
              'W': begin
                  if state.flg_plot MOD 2 EQ 0 then $
                    state.xymnx = state.svxymnx $ ; Reset the screen
                  else x_spec_zerolvl_SetResiduals, state
              end
              '!': x_spec_zerolvl_Reset, state   ; Reset all
              '$': begin ; Reset screen view
                  state.svxymnx = state.svsvxymnx 
                  state.xymnx = state.svxymnx 
              end
              'U': mwrfits, fc_val, state.outfil, /create
              'D': x_spec_zerolvl_Delete, state ; Delete/Undelete a data point
              'x': x_spec_zerolvl_DelReg, state        ; Delete one region
              'X': x_spec_zerolvl_DelReg, state, /all  ; Delete all regions and reset
              'S': x_spec_zerolvl_SaveCont, state, /all  ; Save entire spec_zerolvl
              'g': begin
                  x_spec_zerolvl_SaveCont, state ; Save piece of spec_zerolvl
                  if state.flg_svcont EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif
              end
              'c': begin ; Add spec_zerolvl points to the fit
                  x_spec_zerolvl_AddCont, state 
                  if state.flg_addcont EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return, 1
                  endif else x_spec_zerolvl_UpdateFit, state
              end
              'R': x_spec_zerolvl_SetResiduals, state  ; Plot residuals
              'V': x_spec_zerolvl_SetValues, state     ; Plot values
              'C': x_spec_zerolvl_ClearRej, state ; Clear all rejected points
              'P': x_spec_zerolvl_print, state ; Print
              'M': x_spec_zerolvl_ClearAdd, state
              '3': begin ;; Add point (Spline)
                  if n_elements(state.zro_lvl) NE 1 then $
                    x_spec_zerolvl_UpdateFit, state, SPFLG=1L
              end
              '4': begin ;; Move a point
                  if n_elements(state.zro_lvl) NE 1 then $
                  x_spec_zerolvl_UpdateFit, state, SPFLG=2L
              end
              'q': begin
                  x_spec_zerolvl_setpnt, state
                  widget_control, ev.top, /destroy
                  return, 0
              end
              else:  ; Nothing
          endcase
      end
;            REJECTION
      'HSIGVAL': begin
          fitprm[state.curfit].hsig = ev.value
          if fitprm[state.curfit].niter EQ 0 then fitprm[state.curfit].niter=3
          x_spec_zerolvl_UnDelRej, state
          if ev.value NE 0. then fitprm[state.curfit].flg_rej = 1 else begin
              if fitprm[state.curfit].lsig EQ 0. then fitprm[state.curfit].flg_rej = 0
          endelse
      end
      'LSIGVAL': begin
          fitprm[state.curfit].lsig = ev.value
          if fitprm[state.curfit].niter EQ 0 then fitprm[state.curfit].niter=3
          x_spec_zerolvl_UnDelRej, state
          if ev.value NE 0. then fitprm[state.curfit].flg_rej = 1 else begin
              if fitprm[state.curfit].hsig EQ 0. then fitprm[state.curfit].flg_rej = 0
          endelse
      end
      ;; FIT
      'CURFIT': begin
          state.flg_fit = 0
          state.curfit = ev.value
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state
          flg_plt = 1
      end
      'UPFIT' : begin
          state.curfit = state.curfit + 1
          state.flg_fit = 0
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state
          flg_plt = 1
      end
      'DOWNFIT' : begin
          state.curfit = (state.curfit - 1) > 0L
          state.flg_fit = 0
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state
          flg_plt = 1
      end
      ;; Spline
      'CURSPL': begin
          state.flg_fit = 1
          state.curspl = ev.value
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state, spflg=3
          flg_plt = 1
      end
      'UPSPL' : begin
          state.curspl = state.curspl + 1
          state.flg_fit = 1
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state, spflg=3
          flg_plt = 1
      end
      'DOWNSPL' : begin
          state.curspl = (state.curspl - 1) > 0L
          state.flg_fit = 1
          x_spec_zerolvl_Reset, state
          x_spec_zerolvl_UpdateFit, state, spflg=3L
          flg_plt = 1
      end
      'SVSPL' : x_spec_zerolvl_splout, state
      'SVALL' : x_spec_zerolvl_idlout, state
      'DONE' : begin
          x_spec_zerolvl_setpnt, state
          widget_control, ev.top, /destroy
          return, 0
      end
      else:
  endcase

; Update Plot
  if flg_plt EQ 1 then x_spec_zerolvl_UpdatePlot, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return, 1
end
  
;;;;;;;;;;;;;;;;;;;;
;  IDL Output
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_splout, state, ERR=err

  zro_lvl = state.zro_lvl
  if tag_exist(state, 'cstr') EQ 1 then cstr = state.cstr else cstr = 0.
  save, cstr, zro_lvl, filename='zro_lvl.idl'

  return
end

;;;;;;;;;;;;;;;;;;;;
;  IDL Output
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_idlout, state, ERR=err

  common x_spec_zerolvl_fit

  widget_control, /hourglass   
  save, xfit, state, fitprm, fc_val, filename=state.fit_all, /compress

  return
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_print, state

  print, 'x_spec_zerolvl:  Printing to idl.ps'
; Device
  device, get_decomposed=svdecomp

  !p.thick = 3
  !p.charthick = 3

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color, /maxs
  state.psfile = 1
  x_spec_zerolvl_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0
  !p.thick = 1
  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_spec_zerolvl_UpdatePlot, state
  
common x_spec_zerolvl_fit

; Plot Data

  widget_control, /hourglass   
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
   endif

  if state.flg_apply then begin 
      case state.flg_fit of 
         0: sfx = state.fx - fin_fit
         1: sfx = state.fx - state.zro_lvl 
         else: stop
      endcase
  endif else sfx = state.fx

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
          oplot, state.wave[gdorg], sfx[gdorg], psym=10, $
            color=color.black 
      endif else $
        plot, [-1.0],  [state.svxymnx[3]+1.], psym=1, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        background=color.white, color=color.black

      oplot, [-1e9, 1e9], [0., 0.], color=color.gray, linesty=1

      ; Overplot regions
      for q=0L, state.nreg[state.curfit]-1 do begin
          pts = where(state.wave LT state.reg[q,1,state.curfit] AND $
                      state.wave GT state.reg[q,0,state.curfit], npts) 
          if npts NE 0 then $
            oplot, [state.wave[pts]], [sfx[pts]], psym=10, $
                    color=color.blue ; DATA IN REGION
      endfor

      ; Overplot rejected/other points
      for q=0,7 do begin
          ;; Nice line!
          if q EQ 1 OR q EQ 3 or q EQ 4 OR q EQ 6 OR q EQ 0 OR q EQ 5 then continue
          pts = where(state.gdpix[0:state.norg-1] EQ q AND $
                      state.wave[0:state.norg-1] GT state.xymnx[0] $
                      AND state.wave[0:state.norg-1] LT state.xymnx[2], count)
          if count NE 0 then begin
              case q of
                  0: oplot, [state.wave[pts]], [sfx[pts]], psym=2, $
                    color=color.green ; DELETED DATA OUT OF REGION
                  2: oplot, [state.wave[pts]], [sfx[pts]], psym=2, $
                    color=color.purple ; DELETED DATA IN REGION
                  5: oplot, [state.wave[pts]], [sfx[pts]], psym=5, $
                    color=color.pink ; REJECTED DATA OUT OF REGION
                  7: oplot, [state.wave[pts]], [sfx[pts]], psym=5, $
                    color=color.red ; REJECTED DATA IN REGION
                  else:
              endcase
          endif
      endfor
      
      
;    Added points
      if(state.ntot GT state.norg) then begin
          gd = where(state.gdpix[0:state.ntot-1] EQ 8, ngd)
          if ngd NE 0 then $
            oplot, [state.xtot[gd]], [state.ytot[gd]], psym=5, color=color.cyan
      endif

      ;;    Continuum
      cpts = where(state.wave[0:state.norg-1] LT state.xymnx[2] AND $
                  state.wave[0:state.norg-1] GT state.xymnx[0], npts)
      if npts NE 0 then $
        oplot, [state.wave[cpts]], [fc_val[cpts]], color=color.green ; Best fit Red  

;    Fit
      case state.flg_fit of 
          0: begin
              if state.nreg[state.curfit] NE 0 then begin
                  pts = where(xfit LT state.xymnx[2] AND $
                              xfit GT state.xymnx[0], npts)
                  if npts NE 0 then $
                    oplot, [xfit[pts]], [fin_fit[pts]], $
                    color=color.red ; Best fit Red  
              endif
          end
          1: begin  ; Spline
              ;; Points
              if state.cstr.npts NE 0 then begin
                  gdc = where(state.cstr.msk EQ 1)
                  oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], $
                    psym=1, color=color.cyan, symsize=5
              endif
              ;; Line
              oplot, state.wave[cpts], state.zro_lvl[cpts], color=color.red
          end
          else: stop
      endcase


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

pro x_spec_zerolvl_UpdateFit, state, SPFLG=spflg

common x_spec_zerolvl_fit
  widget_control, /hourglass   

; Set Good points
  case state.flg_fit of
      0: begin
          if state.nreg[state.curfit] EQ 0 then return ; Require regions
          ;; Deal with added points
          if state.norg LT state.ntot then begin 
              gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count) 
              ;; Extras
              moregd = where(state.gdpix[state.norg:*] EQ 8, nmgd)
              gdpt = [gdpt, state.norg+moregd]
              count = count + nmgd
          endif else gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count)
          
          mask = lonarr(state.ntot)
          mask[gdpt] = 1
          fitprm[state.curfit].minpt = fitprm[state.curfit].minpt > count/2
          
          ;; Catch error
          if count EQ 0 then begin
              widget_control, state.error_msg_id, $
                set_value='No good data!  Adjust regions as necessary...'
              widget_control, base_id, set_uvalue=state, /no_copy
              return
          endif

          ;; Error
          if fitprm[state.curfit].nord EQ 0 then begin
              widget_control, state.error_msg_id, set_value='FIT requires nord > 0'
              fitprm[state.curfit].nord = 1
              widget_control, state.lblordr_id, $
                set_value=string(fitprm[state.curfit].nord, format='(i4)')
          endif else begin      ; FIT!
              tmpfit = fitprm[state.curfit]
;              if tmpfit.func EQ 'BSPLIN' then flg_bsp = 2 else flg_bsp=0
              if fitprm[state.curfit].flg_rej EQ 0 then begin
                  fit = x_fit(state.xtot[gdpt], state.ytot[gdpt], $
                              FITSTR=tmpfit, $
                              SIG=state.wtot[gdpt], $
                              IVAR=state.ivtot[gdpt], $
                              FLG_BSP=flg_bsp)
                  fitprm[state.curfit] = tmpfit
                  if fit[0] EQ -1 then begin
                      fin_fit = fltarr(n_elements(xfit))
                      widget_control, state.error_msg_id, $
                        set_value='Bad fit! Adjust nord probably'
                      widget_control, base_id, set_uvalue=state, /no_copy
                      return
                  endif
              endif else begin  ; FIT with REJECTION!
                  rejpt = -1
                  fit = x_fitrej(state.xtot[gdpt],$
                                 state.ytot[gdpt], $
                                 FITSTR=tmpfit, $
                                 REJPT=rejpt, FLG_BSP=flg_bsp )
;                                 SIG=state.wtot[gdpt], $
;                                 IVAR=state.ivtot[gdpt], $
                  fitprm[state.curfit] = tmpfit
                  if rejpt[0] NE -1 then begin
                      rejpt = gdpt[rejpt]
                      state.gdpix[rejpt] = state.gdpix[rejpt]+4
                  endif
              endelse
          endelse
  

          ;; FIT CHECK
          if fit[0] EQ -1 then stop
          fit = x_calcfit(state.xtot[0:state.ntot-1], $
                          FITSTR=fitprm[state.curfit])

          ;; Calculate values
          fin_fit = fit
          delvarx, mask, fit

          ;; RESIDUALS AND RMS
          x_spec_zerolvl_CalcResid, state

          ;; Set Good points
          if state.nreg[state.curfit] EQ 0 then $
            gdpt = where(state.gdpix[0:state.ntot-1] EQ 1, count) $
          else begin
              if state.norg LT state.ntot then begin ; Deal with added points
                  gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count) 
                  for i=state.norg, state.ntot-1 do gdpt = [gdpt,i]
                  count = count + state.ntot - state.norg
              endif else gdpt = where(state.gdpix[0:state.norg-1] EQ 3, count)
          endelse
          
          fitprm[state.curfit].rms = sqrt(total( state.residuals[gdpt]^2 ) $
                            / float(n_elements(gdpt)-1))
          widget_control, state.rms_id, set_value=fitprm[state.curfit].rms
      end
      1: begin  ; Spline
          x_fitline_spec_zerolvl, state, spflg
      end
      else: stop
  end

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_Reset, state

common x_spec_zerolvl_fit

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
  case state.flg_fit of 
      0: begin
          state.gdpix[0:n_elements(state.wave)-1] = 1
          if state.nreg[state.curfit] NE 0 then begin
              for ii=0L,state.nreg[state.curfit]-1 do begin
                  gd = where(state.wave LE state.reg[ii,1,state.curfit] AND $
                             state.wave GE state.reg[ii,0,state.curfit], ngd)
                  if ngd NE 0 then state.gdpix[gd] = 3
              endfor
          endif 
          widget_control, state.curfit_id, set_value=state.curfit
      end
      1: begin ; spline
          state.cstr = state.svspl[state.curspl]
          if state.cstr.npts EQ 0 then state.zro_lvl[*] = 0.
          widget_control, state.curspl_id, set_value=state.curspl
      end
      else: stop
  endcase

  ;; Set the gui parameters
;  widget_control, state.funclist_id, set_value=fitprm[state.curfit].func
  tmp_string = string(fitprm[state.curfit].nord, format='(i4)')
  widget_control, state.lblordr_id, set_value=tmp_string
  case strtrim(fitprm[state.curfit].func,2) of 
      'POLY' :  widget_control, state.funclist_id, set_list_select=0
      'LEGEND' :  widget_control, state.funclist_id, set_list_select=1
      'CHEBY' :  widget_control, state.funclist_id, set_list_select=2
      'BSPLIN' :  widget_control, state.funclist_id, set_list_select=3
      'GAUSS' :  widget_control, state.funclist_id, set_list_select=4
  endcase

end

;;;;;;;;;;;;;;;;;;;;
;  Clear
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_Clear, state

common x_spec_zerolvl_fit

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
  state.nreg[state.curfit] = 0

; Plotting
  state.xymnx = state.svsvxymnx
  state.svxymnx = state.svsvxymnx

; Fitting
  fitprm[state.curfit].nord = 3
  fitprm[state.curfit].func = 'LEGEND'
  widget_control, state.funclist_id, set_list_select=1
;  widget_control, state.funclist_id, set_value=1
  tmp_string = string(fitprm[state.curfit].nord, format='(i4)')
  widget_control, state.lblordr_id, set_value=tmp_string

; Zero out spline
  state.cstr.npts = 0
  state.zro_lvl[*] = 0.

end

;;;;;;;;;;;;;;;;;;;;
;  Delete
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_Delete, state

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
  x_spec_zerolvl_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  UnDelRej
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_UnDelRej, state

  rej = where(state.gdpix MOD 8 GE 4, count)
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej]-4
  ; Update Fit
  x_spec_zerolvl_UpdateFit, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  x_spec_zerolvl_GetReg, state
pro x_spec_zerolvl_GetReg, state
  ;; 
  if state.flg_fit EQ 1 then begin
      widget_control, state.error_msg_id, set_value='Turn Spline off first!!'
      return
  endif
  ; First region?
  if state.flg_reg EQ 0 then begin
      state.flg_reg = 1
      state.reg[state.nreg[state.curfit],0,state.curfit] = xgetx_plt(state, /strct)
  endif else begin
      state.flg_reg = 0
      state.reg[state.nreg[state.curfit],1,state.curfit] = xgetx_plt(state, /strct)
      state.nreg[state.curfit] = state.nreg[state.curfit] + 1
      x_spec_zerolvl_SetReg, state
  endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;  Set Regions
;;;;

pro x_spec_zerolvl_SetReg, state

;  Determines all points in the regions -- Only examines data

;  Require reg[1] > reg[0]

  if state.reg[state.nreg[state.curfit]-1,1,state.curfit] LT state.reg[state.nreg[state.curfit]-1,0,state.curfit] then begin
      tmp = state.reg[state.nreg[state.curfit]-1,1,state.curfit]
      state.reg[state.nreg[state.curfit]-1,1,state.curfit] = state.reg[state.nreg[state.curfit]-1,0,state.curfit]
      state.reg[state.nreg[state.curfit]-1,0,state.curfit] = tmp
  endif

; Require regions dont overlap

  if state.nreg[state.curfit] GT 1 then begin
      bla = where(state.reg[0:state.nreg[state.curfit]-1,*,state.curfit] LT state.reg[state.nreg[state.curfit]-1,1,state.curfit] AND $
                  state.reg[0:state.nreg[state.curfit]-1,*,state.curfit] GT state.reg[state.nreg[state.curfit]-1,0,state.curfit], count)
      if count NE 0 then begin
          widget_control, state.error_msg_id, set_value='Regions may not overlap!'
          state.nreg[state.curfit] = state.nreg[state.curfit] - 1
          return
      endif
  endif
                  

; Adjust new pixels

  newpix = where(state.wave LE $
                 state.reg[state.nreg[state.curfit]-1,1,state.curfit] AND $
                 state.wave GE $
                 state.reg[state.nreg[state.curfit]-1,0,state.curfit] AND $
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
      state.nreg[state.curfit] = state.nreg[state.curfit] - 1
  endelse
  ; Update Fit
  x_spec_zerolvl_UpdateFit, state

end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
; Delete Region(s)
;;;;

pro x_spec_zerolvl_DelReg, state, ALL=all

  ;; Check spline flag
  if state.flg_fit EQ 1 then begin
      widget_control, state.error_msg_id, set_value='Turn Spline off first!!'
      return
  endif

  ;; 
  if state.nreg[state.curfit] EQ 0 then return
; Deletes 1 or All regions

  if (keyword_set (ALL) OR state.nreg[state.curfit] EQ 1) then begin
      if state.nreg[state.curfit] EQ 0 then return $
      else begin
          regpix = where(state.gdpix[0:state.norg-1] MOD 4 GT 1)
          state.gdpix[regpix] = state.gdpix[regpix] - 2
          state.nreg[state.curfit] = 0
      endelse
  endif else begin
      xpos = xgetx_plt(state, /strct)
      fndreg = where( xpos GE state.reg[0:state.nreg[state.curfit]-1,0, $
                                        state.curfit] AND $
                      xpos LE state.reg[0:state.nreg[state.curfit]-1,1, $
                                        state.curfit], count, $
                      complement=newreg, ncomplement=nnew)
      if count EQ 0 then begin      ; Not in the region
          widget_control, state.error_msg_id, set_value='No region found!'
          return
      endif
; Reset gdpix
      fndpix = where(state.wave LE (state.reg[fndreg,1,state.curfit])[0] AND $
                     state.wave GE (state.reg[fndreg,0,state.curfit])[0])
      state.gdpix[fndpix] = state.gdpix[fndpix] - 2
; Reset regions
      for i=0,nnew-1 do begin  
          state.reg[i,0,state.curfit] = state.reg[newreg[i],0,state.curfit]
          state.reg[i,1,state.curfit] = state.reg[newreg[i],1,state.curfit]
      endfor
      state.nreg[state.curfit] = state.nreg[state.curfit] - 1
  endelse
  ; Update Fit
  x_spec_zerolvl_UpdateFit, state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_setpnt, state

common x_spec_zerolvl_fit

; 

;  *state.pnt_nreg = state.nreg
;  *state.pnt_reg = state.reg
  svdelpts = where(state.gdpix[0:state.norg-1] MOD 2 EQ 0)

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_SetResiduals, state

  ; Set the Flag
  if state.flg_plot MOD 2 NE 1 then state.flg_plot = state.flg_plot + 1

  ; Calculate the residuals
  x_spec_zerolvl_CalcResid, state

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

pro x_spec_zerolvl_CalcResid, state

  common x_spec_zerolvl_fit

  ; Calculate the residuals
  fit = x_calcfit(state.xtot[0:state.ntot-1], FITSTR=fitprm[state.curfit])
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
pro x_spec_zerolvl_SetValues, state

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
pro x_spec_zerolvl_ClearRej, state

  ; Rej
  rej = where(state.gdpix[0:state.ntot-1] MOD 8 GE 4, count)

  ; Clear
  if count NE 0 then state.gdpix[rej] = state.gdpix[rej] - 4

  ; Update Fit
  x_spec_zerolvl_UpdateFit, state

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Clear all added points
;
pro x_spec_zerolvl_ClearAdd, state

  ; Rej
  add = where(state.gdpix[0:state.ntot-1] EQ 8, count)

  ; Clear
  if count NE 0 then state.gdpix[add] = 0

  ; Update Fit
  x_spec_zerolvl_UpdateFit, state

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Save fit to spec_zerolvl
;
pro x_spec_zerolvl_SaveCont, state, ALL=all

common x_spec_zerolvl_fit

  if keyword_set( ALL ) then begin
      if state.flg_fit EQ 0 then fc_val = fin_fit[0:state.norg-1]  $
      else fc_val = state.zro_lvl[0:state.norg-1]
      ;; Current
      case state.flg_fit of
          0: state.curfit = state.curfit + 1
          1: begin
              state.svspl[state.curspl] = state.cstr
              state.curspl = state.curspl + 1
          end
          else: stop
      endcase
      ;; Reset the fit
      x_spec_zerolvl_Reset, state
      x_spec_zerolvl_UpdateFit, state, spflg=3
      return
  endif 

  if state.flg_svcont EQ 0 then begin
      state.flg_svcont = 1
      state.tmpxy[0] = xgetx_plt(state, /strct)
      return
  endif

  ;; Save
  state.flg_svcont = 0
  x2 = xgetx_plt(state, /strct)
  xmn = x2 < state.tmpxy[0]
  xmx = x2 > state.tmpxy[0]
  ;; Find all xdat points
  gdpt = where(xfit LT xmx AND xfit GT xmn, ngd)

  case state.flg_fit of
      0: begin ;; Fit
          if ngd NE 0 then fc_val[gdpt] = fin_fit[gdpt]
      end
      1: begin ;; Spline
          if ngd NE 0 then fc_val[gdpt] = state.zro_lvl[gdpt]
      end
      else: stop
  endcase

  case state.flg_fit of
      0: state.curfit = state.curfit + 1
      1: begin
          state.svspl[state.curspl] = state.cstr
          state.curspl = state.curspl + 1
      end
      else: stop
  endcase

  ;; Reset
  x_spec_zerolvl_Reset, state
  x_spec_zerolvl_UpdateFit, state, spflg=3
  return

end

;;;;;;;;;;;;;;;;;;;;
; Add spec_zerolvl to fit
;
pro x_spec_zerolvl_AddCont, state

common x_spec_zerolvl_fit

  ; First time?
  if state.flg_addcont EQ 0 then begin
      state.flg_addcont = 1
      state.tmpxy[0] = xgetx_plt(state, /strct)
  endif else begin
      state.flg_addcont = 0
      x2 = xgetx_plt(state, /strct)
      xmn = x2 < state.tmpxy[0]
      xmx = x2 > state.tmpxy[0]
      ;; Find all xdat points
      gdpt = where(xfit LT xmx AND xfit GT xmn, ngd)
      ;; Spline?
      if state.flg_fit EQ 1 then begin
          for qq=0L,ngd-1 do begin
              state.xpos = xfit[gdpt[qq]]
              state.ypos = fc_val[gdpt[qq]]
              x_fitline_spec_zerolvl, state, 1
          endfor
      endif else begin
          if ngd NE 0 then begin 
              state.xtot[state.ntot:state.ntot+ngd-1] = xfit[gdpt]
              state.ytot[state.ntot:state.ntot+ngd-1] = fc_val[gdpt]
              state.wtot[state.ntot:state.ntot+ngd-1] = 1./sqrt(state.med_ivar)/10.
              ;; Inverse variance
              state.ivtot[state.ntot:state.ntot+ngd-1] = 30*state.med_ivar 
              state.gdpix[state.ntot:state.ntot+ngd-1] = 8
              state.ntot = state.ntot + ngd
          endif
      endelse
  endelse
  return
end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro x_spec_zerolvl_inispl, state, inispl


  ;; File?
  a = findfile(inispl, count=na)
  if na EQ 0 then begin
      print, 'x_spec_zerolvl: FILE ', inispl, ' does not exist!'
      stop
  endif

  ;; Restore
  restore, inispl

  ;; Continuum
  if keyword_set(zro_lvl) then state.zro_lvl[*] = zro_lvl else state.zro_lvl[*] = 0.

  ;; conti_str
  if keyword_set( CSTR ) then state.cstr = cstr

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

pro x_spec_zerolvl, flux_fil, error_fil, FITSTR=fitstr, ZRO_LVL=zro_lvl, LSIG=lsig, $
                 XSIZE=xsize, YSIZE=ysize, INFLG=inflg, OUTFIL=outfil, $
                 SPLINE=spline, INISPL=inispl, INISTAT=inistat, WAVE=wave, $
                 NOSIG=nosig

common x_spec_zerolvl_fit

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_spec_zerolvl, flux_fil, error_fil, FITSTR=, INFLG=, OUTFIL=, CONTI=, LSIG='
    print, '    XSIZE=, YSIZE=, INFLG=, OUTFIL=, /SPLINE, INISPL=, INISTAT= [v1.1]'
    print, '    /NOSIG '
    return
  endif 

  if not keyword_set( FIT_ALL ) then fit_all = 'zro_lvl_all.idl'
  if not keyword_set( LSIG ) then lsig = 2.5
  if not keyword_set( INFLG ) then inflg = 0
  if not keyword_set( OUTFIL ) then begin
      if keyword_set( CONTI ) then outfil = zro_lvl else outfil = 'zro_lvl.fits'
  endif

  ;; Screen
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

  resolve_routine, 'x_fitline', /no_recompile

; Read in the Data
  x_spec_zerolvl_initcommon
  if keyword_set( INISTAT ) then begin  ; OLD FILE!
      restore, inistat
      ;; Add tags
      ;; flg_2d
;      state.fx = state.fx*1.5
  endif else begin
      if not keyword_set( INFLG ) then inflg = 0
      ydat = x_readspec(flux_fil, INFLG=inflg, head=head, NPIX=npix, $
                        WAV=xdat, FIL_SIG=error_fil, SIG=ysig)
      if keyword_set(WAVE) then xdat = wave
      
      ;; CONTINUUM
      if keyword_set( ZRO_LVL ) then begin
          cdat = x_readspec(zro_lvl, /fscale)
          if n_elements(cdat) NE n_elements(xdat) then begin
              print, 'x_spec_zerolvl: Continuum is the wrong size!'
              print, 'x_spec_zerolvl: Interpolating...'
              fc_val = interpol(cdat, findgen(n_elements(cdat)), findgen(n_elements(xdat)))
;              return
           endif else fc_val = cdat
      endif 
      
      ;; XFIT
      xfit = xdat

      ;; Set inv variance
      gdsig = where( ysig GT 0., ngd)
      if ngd EQ 0 or keyword_set(NOSIG) then begin
          ysig[*] = 1.
          gdsig = lindgen(npix)
      endif
      ivar = dblarr(npix)
      ivar[gdsig] = 1./(ysig[gdsig])^2
      med_ivar = median(ivar[gdsig])
      
      ;; Initialize the common block
;      x_spec_zerolvl_initcommon

  ; Set fit structure as input
      if keyword_set( FITSTR ) then begin
          fitprm = fitstr 
      endif
      fitprm.lsig = lsig
      

      ;;    Pointer for Final fit
      tmpreg = fltarr(100,2,100)
      nreg = 0
      pnt_nreg = PTR_NEW( nreg )
      
      tmp2 = { conti_str, $
               npts: 0L, $
               xval: dblarr(100), $
               yval: dblarr(100), $
               msk: lonarr(100) }
      
      
      ;;    STATE
      state = { $
              wave: xdat, $
              fx: double(ydat), $
              sig: double(ysig), $
              ivar: ivar, $
              med_ivar: med_ivar, $
              norg: n_elements(xdat), $
              residuals: dblarr(200000L), $ ; Residuals
              flg_plot: 0, $
              fit_all: fit_all, $
              flg_apply: 0, $
              xtot: dblarr(200000L), $ ; Fitting vectors
              ytot: dblarr(200000L), $
              wtot: dblarr(200000L), $ ; Weighting Factor
              ivtot: dblarr(200000L), $ ; Weighting Factor
              flg_wgt: 1, $ ; 0 = No input 1 = Input
              zro_lvl: replicate(0.d,npix), $ ; Spline
              psfile: 0, $
              curspl: 0L, $
              cstr: tmp2, $
              svspl: replicate(tmp2,100), $
              nlin: 0, $
              xpos: 0.d, $
              ypos: 0.d, $
              ntot: n_elements(xdat), $
              curfit: 0L, $
              nreg: lonarr(100), $
              reg: tmpreg, $
              outfil: outfil, $ ; Output file
              flg_fit: 0, $ ; 1=Spline
              flg_reg: 0, $
              flg_zoom: 0, $
              flg_svcont: 0, $
              flg_addcont: 0, $
              pos: [0.1,0.1,0.95,0.95], $ ; Plotting
              res_svxymnx: fltarr(4), $ ; Residuals
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
              base_id: 0L, $ ; Widgets
              lblordr_id: 0L, $
              draw_base_id: 0L, $
              draw_id: 0L, $
              text_id: 0L, $
              funclist_id: 0L, $
              error_msg_id: 0L, $
              reject_base_id: 0L, $
              spline_id: 0L, $
              hsig_id: 0L, $
              lsig_id: 0L, $
              curfit_id: 0L, $
              curspl_id: 0L, $
              rms_id: 0L, $
              gdpix: intarr(200000L), $
              flg_2d: 0 $ ;; New tags start here  
              }
      ;; SVXYMNX
      state.xymnx = state.dat_svxymnx
      state.svxymnx = state.dat_svxymnx
      state.svsvxymnx = state.dat_svxymnx
  endelse
      
  if keyword_set( SPLINE ) then state.flg_fit = 1
  
  ;;    WIDGET
  base = WIDGET_BASE( title = 'x_spec_zerolvl: Interactive Mode', /column)
  state.base_id = base
  
  ;;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
  strlbl = strarr(10)
  strlbl = ['x_spec_zerolvl', ' ', 'Ver 1.1']
  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)
;        Function
  funcval = ['POLY', 'LEGEND', 'CHEBY', 'BSPLIN', 'GAUSS']
  state.funclist_id = WIDGET_LIST(toolbar, VALUE=funcval, uvalue='FUNCLIST', $
                                  ysize=4)
  case strtrim(fitprm[state.curfit].func,2) of 
      'POLY' :  widget_control, state.funclist_id, set_list_select=0
      'LEGEND' :  widget_control, state.funclist_id, set_list_select=1
      'CHEBY' :  widget_control, state.funclist_id, set_list_select=2
      'BSPLIN' :  widget_control, state.funclist_id, set_list_select=3
      'GAUSS' :  widget_control, state.funclist_id, set_list_select=4
  endcase
  
;        Order      
  ordrbase = WIDGET_BASE(toolbar, /column)
  upordr = WIDGET_BUTTON(ordrbase, value='Up',uvalue='UP')
  tmp_string = string(fitprm[state.curfit].nord, format='(i4)')
  state.lblordr_id = WIDGET_LABEL(ordrbase, value=tmp_string, $
                                  /dynamic_resize, /align_center)
  dwnordr = WIDGET_BUTTON(ordrbase, value='Down',uvalue='DOWN')
  
;  RMS
  state.rms_id = cw_field(ordrbase, value='0.', /floating, title='RMS', xsize=6)

;        Rejection
  state.reject_base_id = widget_base(toolbar, /column, /align_center,$
                                     /frame)
  hsigbase = widget_base(state.reject_base_id, /row, /align_right)
  state.hsig_id = cw_field(hsigbase, title='High Sig', $
                           value=fitprm[state.curfit].hsig, $
                           xsize=5, ysize=1, /return_events, $
                           uvalue='HSIGVAL')
  lsigbase = widget_base(state.reject_base_id, /row, /align_right)
  state.lsig_id = cw_field(lsigbase, title='Low Sig', $
                           value=fitprm[state.curfit].lsig, $
                           xsize=5, ysize=1, /return_events, $
                           uvalue='LSIGVAL')
  lsigbase = widget_base(state.reject_base_id, /row, /align_right)
  curfitbase = widget_base(toolbar, /column, /align_right)
  state.curfit_id = cw_field(curfitbase, title='Fit#', $
                           value=state.curfit, $
                           xsize=3, ysize=1, /return_events, $
                           uvalue='CURFIT')
  curfitbutbase = widget_base(curfitbase, /row, /align_right)
  upfti = WIDGET_BUTTON(curfitbutbase, value='Up',uvalue='UPFIT')
  dwnfti = WIDGET_BUTTON(curfitbutbase, value='Down',uvalue='DOWNFIT')

  cursplbase = widget_base(toolbar, /column, /align_right)
  state.curspl_id = cw_field(cursplbase, title='Spl#', $
                           value=state.curspl, $
                           xsize=3, ysize=1, /return_events, $
                           uvalue='CURSPL')
  cursplbutbase = widget_base(cursplbase, /row, /align_right)
  upspl = WIDGET_BUTTON(cursplbutbase, value='Up',uvalue='UPSPL')
  dwnspl = WIDGET_BUTTON(cursplbutbase, value='Down',uvalue='DOWNSPL')
  ;; Spline
  crude_err = widget_base(toolbar, /row, /align_center, frame=2)
  state.spline_id = CW_BGROUP(crude_err, ['No', 'Yes'], $
                             label_top='Spline?', $
                             row=2, /exclusive, /no_release, $
                             set_value=state.flg_fit,  uvalue='SPLINE')
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
             'g,g -- Save piece of spec_zerolvl',$
             'S -- Save entire spec_zerolvl',$
             'x -- Remove one Region',$
             'D -- Delete data point', $
             'c,c -- Add spec_zerolvl points to the fit', $
             'Z -- Set zero level to 0', $
             '--------------------', $
             'f -- Update fit', $
             'u,d -- In(de)crease order 1',$
             'C -- Clear all rejected points', $
             'M -- Clear all added points', $
             '! -- Continuum reset', $
             '--------------------', $
             'X -- Reset fit',$
             'R -- Plot Residuals',$
             'V -- Plot Values',$
             'U -- Write out to file', $
             'q -- Quit and save']
  help_text_id = widget_text(toolbar, value=strhelp, xsize=15, ysize=4,$
                             /scroll)
;      Error base
  error_base_id  = widget_base(toolbar, /column)
  error_btn_id = widget_button(error_base_id, value='Erase Error Msg', $
                               uvalue='ERRORB')
;          Error Messages
  state.error_msg_id = widget_text(error_base_id, value='', xsize=20, $
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
  morbutt = widget_base(toolbar, /column, /align_right)
  svspl = WIDGET_BUTTON(morbutt, value='SVSPL',uvalue='SVSPL', /align_right)
  svidl = WIDGET_BUTTON(morbutt, value='SVALL',uvalue='SVALL', /align_right)
  done = WIDGET_BUTTON(morbutt, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Include all elements at first
  state.gdpix[0:n_elements(xdat)] = 1
  

; Update
  x_spec_zerolvl_Reset, state
; Continuum
  if not keyword_set( fc_val ) then fc_val = fltarr(state.norg) 
  ;; Init Fit
  if keyword_set(INISPL) then x_spec_zerolvl_inispl, state, inispl

  x_spec_zerolvl_UpdateFit, state
  x_spec_zerolvl_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
; xmanager blows (BLOCK)
  repeat begin
      ev = widget_event(base)
      ans = x_spec_zerolvl_ev(ev)
  end until ans EQ 0
  
;  End game

  ; Output
  mwrfits, float(fc_val), outfil, /create

;  if arg_present( FITSTR ) then fitstr=fitprm
;  if arg_present( FFIT ) then ffit = *fitprm.ffit
;  if arg_present( RMS ) then rms = fitprm.rms
  if arg_present( DELPTS ) then delpts = svdelpts
;  if arg_present( NRM ) then nrm = fitprm.nrm
  if arg_present( CONTI ) then zro_lvl = temporary(fc_val) else delvarx, fc_val


; Memory
;  delvarx, fitprm, xfit

  return
end
