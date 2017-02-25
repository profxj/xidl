;+ 
; NAME:
; esi_echspecplt   
;   Version 1.0
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   esi_echspecplt, ydat, [head], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           QAL=, ERR=, /GAL
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array
;   ERR        - Error array (fits or image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echspecplt, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Aug-2002 Written by JXP
;   17-Mar-2012 Quick fix to get 'D' option working by MR
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro esi_echspecplt_initcommon
;
common x_specplot_lines, $
   flg_lines, $
   flg_lines2, $
   lines, $
   lines2, $
   zabs, $
   zabs2

  flg_lines = 0
  flg_lines2 = 0
  zabs = 0.
  zabs2 = 0.

common esi_echspecplt_cmmn, esiobj

end

;;;;
; Events
;;;;

pro esi_echspecplt_event, ev

common x_specplot_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          zabs = tmp
          flg_lines = 1
      end
      'ZABS2' : begin
          widget_control, state.zabs2_id, get_value=tmp
          zabs2 = tmp
          flg_lines2 = 1
      end
      'LNLIST' : begin  ; LINE LIST
          state.flg_lines = ev.index + 1
          x_specplot_initLines, state
      end
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
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
                      1 : begin     ; Left button = Zoom
                          x_speczoomreg, state
                          if state.flg_zoom NE 2 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state,  /no_copy
                              return
                          endif else state.flg_zoom = 0
                      end 
                      4 : esi_echspecplt_SetLine, state   ; Set reference line
                      else: 
                  endcase
              end
              1 : begin ; Button Release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state,  /no_copy
                  return
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
                  widget_control, state.xpos_id, set_value=state.xpos
                  widget_control, state.ypos_id, set_value=state.ypos
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
              print, 'Set the other region!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          if (state.flg_EW EQ 1 AND eventch NE 'E') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
              print, 'Set the other side for EW!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'T': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              ; ZOOMING
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'z': begin  ; Region
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
              'w': state.xymnx = state.svxymnx ; Reset the screen
              ; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'H': x_helpwidg, state.help
              ; SMOOTH
              'S': esi_echspecplt_smooth, state
              'U': esi_echspecplt_smooth, state, /reset
              ; Crude Analysis
              'E': esi_echspecplt_EW, state        ; Calc EW
              'N': esi_echspecplt_Colm, state      ; Calc AODM colm
              ' ': print, 'x: '+strtrim(state.xpos,2)+$
                '   y:'+strtrim(state.ypos,2)
              ; LINES
              'L': esi_echspecplt_SetLine, state   ; Set reference line
              'A': begin ; Plot AlIII
                  x_specplot_guess, state, 'A'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'C': begin ; Plot CIV
                  x_specplot_guess, state, 'C'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'I': begin ; Plot SiIV
                  x_specplot_guess, state, 'S'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'M': begin ; Plot MgII
                  x_specplot_guess, state, 'M'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              ;; Damped Lya overplot
              'D': m_specplot_initdla, state
              ;; ORDERS
              '0': esi_echspecplt_Chgordr, state, 0
              '1': esi_echspecplt_Chgordr, state, 1
              '2': esi_echspecplt_Chgordr, state, 2
              '3': esi_echspecplt_Chgordr, state, 3
              '4': esi_echspecplt_Chgordr, state, 4
              '5': esi_echspecplt_Chgordr, state, 5
              '6': esi_echspecplt_Chgordr, state, 6
              '7': esi_echspecplt_Chgordr, state, 7
              '8': esi_echspecplt_Chgordr, state, 8
              '9': esi_echspecplt_Chgordr, state, 9
              '=': esi_echspecplt_Chgordr, state, /add
              '-': esi_echspecplt_Chgordr, state, /sub
              ; Postscript
              'P': esi_echspecplt_psfile, state  
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'esi_echspecplt: Not a valid key!' ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  esi_echspecplt_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro esi_echspecplt_UpdatePlot, state
  
common x_specplot_lines

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

  if state.flg_smooth EQ 0 then $
    plot, state.wave[0:state.npix-1], state.fx[0:state.npix-1], $
    psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.7, $
    ycharsize=1.7 $
  else $
    plot, state.wave[0:state.npix-1], state.smooth[0:state.npix-1], $
    psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.7, $
    ycharsize=1.7 

  ; Plot Error array
  if state.flg_sig EQ 1 then $
    oplot, state.wave[0:state.npix-1], state.sig[0:state.npix-1], $
    psym=state.psym, color=clr.red

  ; Plot Lines as required
  if flg_lines EQ 1 then esi_echspecplt_PltLines, state
  
  ; EW
  if state.flg_EW NE 0 then esi_echspecplt_PltEW, state

  ; Colm
  if state.flg_Colm NE 0 then esi_echspecplt_PltClm, state

  ;; DLA
  if state.flg_DLA NE 0 then m_specplot_PltDLA, state

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro esi_echspecplt_Reset, state, FLUX=flux

common esi_echspecplt_cmmn

  if state.flg_fspec EQ 0 then begin
      ;; Set wave, fx
      a = where(esiobj[state.ordr].fx NE 0., na, COMPLEMENT=b, NCOMPLEMENT=nb)
      amin = (a[0] > 0)
      if amin GT 100 then amin = 0L
      bmax = b[nb-1]
      if bmax LT 1000 then bmax = 4999
      npix = bmax-amin+1
      state.npix = npix
      
      state.wave[0:npix-1] = esiobj[state.ordr].wave[amin:bmax]
      if keyword_set( FLUX ) then begin
          state.fx[0:npix-1] = esiobj[state.ordr].flux[amin:bmax] 
          state.sig[0:npix-1] = esiobj[state.ordr].sig[amin:bmax]
      endif else begin
          state.fx[0:npix-1] = esiobj[state.ordr].fx[amin:bmax] 
          state.sig[0:npix-1] = sqrt(esiobj[state.ordr].var[amin:bmax] > 0.)
       endelse
      ;; Sky?
      if state.flag_sky EQ 1 then state.fx[0:npix-1] = esiobj[state.ordr].sky[amin:bmax] 
  endif else begin
      ;; Set wave, fx
      a = where(esiobj.fx[*,state.ordr] NE 0., na, $
                COMPLEMENT=b, NCOMPLEMENT=nb)
      amin = (a[0] > 0)
      if amin GT 100 then amin = 0L
      bmax = b[nb-1]
      if bmax LT 1000 then bmax = 4999
      npix = bmax-amin+1
      state.npix = npix

      state.wave[0:npix-1] = esiobj.wave[amin:bmax,state.ordr]
      if state.flag_sky EQ 1 then state.fx[0:npix-1] = esiobj.sky[amin:bmax,state.ordr] $
         else state.fx[0:npix-1] = esiobj.fx[amin:bmax,state.ordr] 
      state.sig[0:npix-1] = $
        sqrt(esiobj.var[amin:bmax,state.ordr] > 0.)
   endelse

  ;; Plotting
  srt = sort(state.fx[0:npix-1])
  ymx = state.fx[srt[0.9*npix]] * 1.2
  state.svxymnx = [state.wave[0], 0., state.wave[npix-1], ymx]
  state.xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sets the line and redshift with a gui
pro esi_echspecplt_SetLine, state

common x_specplot_lines

  flg_lines = 1
  ; Set Line with gui
  case state.flg_lines of 
      1: setwave = x_slctline(lines, /ISM)
      2: setwave = x_slctline(lines, /GAL)
      3: setwave = x_slctline(lines)
      4: setwave = x_slctline(lines,/ISM)
  endcase

  ; Set redshift
  diff = abs(lines.wave - setwave)
  mndiff = min(diff, imin)
  mnwav = lines[imin].wave

  mrkwav = xgetx_plt(state, /strct)
  zabs = mrkwav / mnwav - 1.d

  widget_control, state.zabs_id, set_value=strmid(strtrim(zabs,2),0,10)

  print, 'zabs = '+strtrim(zabs,2)

  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echspecplt_PltLines, state, IMG=img

common x_specplot_lines

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.blue
      rclr = clr.red
  endif else pclr = 3

  ; Find Min and Max in region
  wvmin = state.xymnx[0]/(zabs+1.)
  wvmax = state.xymnx[2]/(zabs+1.)

  ; Parse list
  allwv = where(lines.wave GE wvmin AND lines.wave LE wvmax, cntwv)

  ; Plot
  ymax = state.xymnx[1] + 0.02*(state.xymnx[3]-state.xymnx[1])
  ymax2 = state.xymnx[1] + 0.8*(state.xymnx[3]-state.xymnx[1])

  for q=0L,cntwv-1 do begin
      xplt = lines[allwv[q]].wave*(1.+zabs)
      ; Name
      case state.flg_lines of
          1: xyouts, xplt, ymax, $ ; QAL
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          2: xyouts, xplt, ymax2, $ ; GAL
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          3: xyouts, xplt, ymax, $
            strtrim(lines[allwv[q]].name,2)+$
            string(lines[allwv[q]].wave,format='(f7.1)'),$
            color=pclr, orientation=90., $
            charsize=1.5
          4: xyouts, xplt, ymax, $ ; LLS
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
	else: stop
      endcase
      ; Dotted line
      oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
        color=pclr, linestyle=1
   endfor

  if flg_lines2 then begin
     ;; Find Min and Max in region
     wvmin = state.xymnx[0]/(zabs2+1.)
     wvmax = state.xymnx[2]/(zabs2+1.)
     
     ;; Parse list
     allwv = where(lines.wave GE wvmin AND lines.wave LE wvmax, cntwv)
     
     ;; Plot
     ymax = state.xymnx[1] + 0.02*(state.xymnx[3]-state.xymnx[1])
     ymax2 = state.xymnx[1] + 0.8*(state.xymnx[3]-state.xymnx[1])
     
     for q=0L,cntwv-1 do begin
        xplt = lines[allwv[q]].wave*(1.+zabs2)
                                ; Name
        case state.flg_lines of
           1: xyouts, xplt, ymax, $ ; QAL
                      strtrim(lines[allwv[q]].name,2), color=rclr, $
                      orientation=90., charsize=1.5
           2: xyouts, xplt, ymax2, $ ; GAL
                      strtrim(lines[allwv[q]].name,2), color=rclr, $
                      orientation=90., charsize=1.5 
           3: xyouts, xplt, ymax, $
                      strtrim(lines[allwv[q]].name,2)+$
                      string(lines[allwv[q]].wave,format='(f7.1)'),$
                      color=rclr, orientation=90., $
                      charsize=1.5
           4: xyouts, xplt, ymax, $ ; LLS
                      strtrim(lines[allwv[q]].name,2), color=rclr, $
                      orientation=90., charsize=1.5
           else: stop
        endcase
                                ; Dotted line
        oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
               color=rclr, linestyle=1
     endfor
  end

  ; Mark

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc EW
;  EW has [pts, x/y]

pro esi_echspecplt_EW, state

common x_specplot_lines

  ; Set the flag
  if state.flg_EW MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_EW = 1 
      state.EW_lmt[0,0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.EW_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.EW_lmt[0] then begin
          state.EW_lmt[1,0] = tmp 
          state.EW_lmt[1,1] = xgety_plt(state, /strct)
      endif else begin
          state.EW_lmt[1,0] = state.EW_lmt[0,0]
          state.EW_lmt[1,1] = state.EW_lmt[0,1]
          state.EW_lmt[0,0] = tmp
          state.EW_lmt[0,1] = xgety_plt(state, /strct)
      endelse

      ; Calc EW
      sumpix = where(state.wave GE state.EW_lmt[0,0] AND $
                     state.wave LE state.EW_lmt[1,0], npix)
      cntm = (state.EW_lmt[0,1]+state.EW_lmt[1,1])/2.
      if npix NE 0 then begin
          state.EW = int_tabulated(state.wave[sumpix], 1.-state.fx[sumpix]/cntm,$
                                   /double)/(1.+zabs)
;          state.EW = total(cntm - state.fx[sumpix], /double)/(1.+zabs)
          ; Multiply by delta lambda
;          state.EW = state.EW * state.xd
          print, 'Rest EW = ', state.EW*1e3, ' mA'
      endif
      
      ; Reset the flag
      state.flg_EW = 2
  endelse

  ; EW xlimit

return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot EW stuff to the screen

pro esi_echspecplt_PltEW, state

  clr = getcolor(/load)
  ; flg
  case state.flg_EW of 
      1: oplot, [state.EW_lmt[0,0]], [state.EW_lmt[0,1]], $
        psym=2, color=clr.red
      2: begin
          oplot, [state.EW_lmt[0,0],state.EW_lmt[1,0]], $
            [state.EW_lmt[0,1],state.EW_lmt[1,1]], $
            psym=-2, color=clr.red, linestyle=2
;          xyouts, 0.5, 0.97, 'EW = '+string(state.EW)+' Ang', $
;            /normal, charsize=1.5, alignment=0.5
          xyouts, 0.5, 0.97, 'Rest EW = '+string(state.EW*1.e3)+' mA', $
            /normal, charsize=1.5, alignment=0.5, color=clr.red
      end
      else :
  endcase
  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc AODM colm
;  EW has [pts, x/y]

pro esi_echspecplt_Colm, state

common x_specplot_lines

  ; Check for error array
  if state.flg_sig NE 1 then begin
      print, 'esi_echspecplt_Colm: Need to specify the error array!'
      return
  endif

  ; Set the flag
  if state.flg_Colm EQ 0 then begin
      ; Set one limit
      state.flg_Colm = 1 
      state.Colm_lmt[0]= xgetx_plt(state, /strct)
      state.EW_lmt[0,1]= xgety_plt(state, /strct)
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.Colm_lmt[0] then begin
          state.Colm_lmt[1] = tmp 
          state.EW_lmt[1,1] = xgety_plt(state, /strct)
      endif else begin
          state.Colm_lmt[1] = state.Colm_lmt[0]
          state.Colm_lmt[0] = tmp
          state.EW_lmt[1,1] = state.EW_lmt[0,1]
          state.EW_lmt[0,1] = xgety_plt(state, /strct)
      endelse

      ;; continuum
      cntm = (state.EW_lmt[0,1]+state.EW_lmt[1,1])/2.

      ; Set pixels
      mn = min(abs(state.Colm_lmt[0]-state.wave),pxmin)
      mn = min(abs(state.Colm_lmt[1]-state.wave),pxmax)

      ; Choose atomic line
      if zabs NE 0. then begin
          wvmin = state.Colm_lmt[0]/(1.+zabs)
          wvmax = state.Colm_lmt[1]/(1.+zabs)
          gdlin = where(lines.wave LT wvmax AND lines.wave GT wvmin, count)
          case count of
              0: begin
                  print, 'No line in region so choose your own!'
                  case state.flg_lines of
                      1: gdwave = x_slctline(lines, /ISM) 
                      2: gdwave = x_slctline(lines, /GAL) 
                      3: gdwave = x_slctline(lines, /GAL) 
                      4: gdwave = x_slctline(lines, /ISM) 
                      else: gdwave = x_slctline(lines)
                  endcase
              end
              1: begin
                  print, 'Using '+strtrim(lines[gdlin].name,2)
                  gdwave = lines[gdlin].wave
              end
              else: begin
                  print, 'More than one line!'
                  print, 'Taking: '+strtrim(lines[gdlin[0]].wave,2)
                  gdwave = lines[gdlin[0]].wave
              end
          endcase
      endif else begin
          print, 'Choose a line'
          case state.flg_lines of
              1: gdwave = x_slctline(lines, /ISM) 
              2: gdwave = x_slctline(lines, /GAL) 
              3: gdwave = x_slctline(lines, /GAL) 
              4: gdwave = x_slctline(lines, /ISM) 
              else: gdwave = x_slctline(lines)
          endcase
      endelse
                  
      ; Calc AODM Colm
      x_aodm, state.wave[pxmin:pxmax], $
        (state.fx[pxmin:pxmax]/cntm), $
        (state.sig[pxmin:pxmax]/cntm), $
        gdwave, clm, sig_clm, /LOG
      if clm GT 0. then begin
          state.colm = clm
          state.sig_colm = sig_clm
      endif else begin
          state.colm = alog10(sig_clm*3.)
          state.sig_colm = -9.99
      endelse
      state.flg_Colm = 2
  endelse

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot Colm stuff to the screen

pro esi_echspecplt_PltClm, state

  clr = getcolor(/load)
  ; flg
  case state.flg_Colm of 
      1: oplot, [state.Colm_lmt[0],state.Colm_lmt[0]], $
        [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
      2: begin
          oplot, [state.Colm_lmt[0],state.Colm_lmt[0]], $
            [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
          oplot, [state.Colm_lmt[1],state.Colm_lmt[1]], $
            [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
          xyouts, 0.2, 0.97, 'Colm = '+strtrim(state.Colm,2)+$
            ' Err = '+strtrim(state.sig_colm,2), $
            /normal, charsize=1.5, alignment=0.5, color=clr.green
          state.flg_colm = 0
      end
      else :
  endcase
  return
end
          
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Plot Guess of lines

pro esi_echspecplt_guess, state, val, IMG=img

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.orange
  endif else pclr = 4

  ; Color
  clr = getcolor(/load)

  ; Plot symbol
  plotsym, 2, color=pclr

  ; Set ypt
  scrn = state.xymnx[3]-state.xymnx[1]
  ypt = state.xymnx[1] + 0.1*scrn
  xpt = xgetx_plt(state, /strct)

  ; Value
  case val of
      'A': begin  ; AlIII
          oplot, [ xpt*1854.7164/1862.7895, $
                   xpt, $
                   xpt*1862.7895/1854.7164 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'Al III', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1854.7164 - 1),2)+', '+$
            strtrim((xpt/1862.7164 -1),2)+']'
      end
      'C': begin  ; CIV
          oplot, [ xpt*1548.195/1550.770, $
                   xpt, $
                   xpt*1548.195/1550.770 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'C IV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1548.195 - 1),2)+', '+$
            strtrim((xpt/1550.770 -1),2)+']'
      end
      'CaHK': begin  ; CaHK
          OII = xpt*3729./3934.79
          CaH = xpt*3969.61/3934.79
          oplot, [ OII, xpt, CaH ], [ypt, ypt, ypt], psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'CaK', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, OII, ypt - 0.05*scrn, 'O[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, CaH, ypt - 0.05*scrn, 'CaH', $
            charsize=1.5, color=pclr, alignment=0.5
          print, 'zgal = [ '+strtrim((xpt/3934.79 - 1),2)+']'
          xyouts, 0.5, 0.9, 'zgal = [ '+strtrim((xpt/3934.79 - 1),2)+']', /normal
      end
      'Ha': begin  ; Halpha
          NIIa = xpt*6549.91/6564.63
          NIIb = xpt*6585.42/6564.63
          oplot, [ NIIa, xpt, NIIb ], [ypt, ypt, ypt], psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'Ha', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, NIIa, ypt - 0.05*scrn, 'N[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, NIIb, ypt - 0.05*scrn, 'N[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          print, 'zgal = [ '+strtrim((xpt/6564.63 - 1),2)+']'
          xyouts, 0.5, 0.9, 'zgal = [ '+strtrim((xpt/6564.63 - 1),2)+']', /normal
      end
      else:
  endcase
  return
end


;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro esi_echspecplt_psfile, state

; Device
  device, get_decomposed=svdecomp

  !p.thick = 3
  !p.charthick = 3

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  esi_echspecplt_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0
  !p.thick = 1
  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;
;  SMOOTH
;;;;;;;;;;;;;;;;;;;;

pro esi_echspecplt_smooth, state, RESET=reset

  if not keyword_set(RESET) then begin
      ; Flag
      state.flg_smooth = state.flg_smooth + 1
      ; Smooth
      state.smooth = smooth(state.fx, 2*state.flg_smooth+1)
  endif else begin
      ; Flag
      state.flg_smooth = 0
      ; Smooth
      state.smooth = state.fx
  endelse

end

;;;;;;;;;;;;;;;;;;;;
;  Init Lines
;;;;;;;;;;;;;;;;;;;;

pro m_specplot_initLines, state

  common x_specplot_lines

  ;; Grab the lines
  case state.flg_lines of
      1: begin ; DLA
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/dla.lst'
          x_specplot_initQAL, llist
      end
      2: begin ; GAL
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/gal.lst'
          x_specplot_initGAL, llist
      end
      3: begin ; QSO
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qso.lst'
          X_specplot_initGAL, llist
      end
      4: begin ; LLS
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
          X_specplot_initQAL, llist
      end
      5: begin ; H2
          x_specplot_initH2, llist
      end
      6: begin ; GRB
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/grb.lst'
          x_specplot_initQAL, llist
       end
      7: begin ; LBG
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lbg.lst'
          x_specplot_initGAL, llist
       end         
      8: begin ;MOLECULES in MM (MF)
         llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/mol.lst'
         x_specplot_initMOL, llist
      end         
      else:
  endcase

  ;; Reset z
  zabs = 0.
  flg_lines = 0
end

;;;;;;;;;;;;;;;;;;;;
;  Init DLA
;;;;;;;;;;;;;;;;;;;;

pro m_specplot_initdla, state, BETA=beta

  common x_specplot_lines

  ;; Grab x,y pos
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ;; Setup HI
  if not keyword_set(BETA) then wrest = 1215.6701d else wrest=1025.7223d
  tmp = x_setline(wrest)

  ;; Get z
  state.dla_line = tmp
  state.dla_line.zabs = (xpt / wrest) - 1.
  state.dla_line.N = 20.3
  state.dla_line.b = 30.0

  ;; Calculate
  xmin = state.wave[0] > (xpt - 40.*(1+state.dla_line.zabs))
  xmax = state.wave[state.npix-1] < (xpt + 40.*(1+state.dla_line.zabs))

  mn = min(abs(state.wave-xmin), imn)
  mx = min(abs(state.wave-xmax), imx)

  state.dla_fx = 1.
  state.dla_fx[imn:imx] = x_voigt(state.wave[imn:imx], state.dla_line, $
                                  FWHM=state.FWHM) 
  state.dla_fx = state.dla_fx * ypt
  state.flg_DLA = 1

  ;; Overplot
  state.flg_lines = 4
  m_specplot_initLines, state
  flg_lines = 1
  zabs = state.dla_line.zabs
  widget_control, state.zabs_id, set_value=zabs
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

  return
end


;;;;;;;;;;;;;;;;;;;;
;  Change ORDR
;;;;;;;;;;;;;;;;;;;;

pro esi_echspecplt_Chgordr, state, newordr, ADD=add, SUB=sub

  ;; Add+Subtract
  if keyword_set( ADD ) or keyword_set( SUB ) then begin
      if keyword_set( ADD ) then begin
          newordr = state.ordr + 1
          if newordr EQ 10 then newordr = 0
          state.ordr = newordr
      endif
      if keyword_set( SUB ) then begin
          newordr = state.ordr - 1
          if newordr LT 0 then newordr = 9
          state.ordr = newordr
      endif
  endif else state.ordr = newordr
  
  ;; RESET
  widget_control, state.ordr_id, set_value=newordr
  esi_echspecplt_Reset, state, FLUX=state.flux

end

;;;;;;;;;;;;;;;;;;;;
;  Plot DLA
;;;;;;;;;;;;;;;;;;;;

pro m_specplot_PltDLA, state

  clr = getcolor(/load)

  oplot, state.wave, state.dla_fx, color=clr.green

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

pro esi_echspecplt, spec_fil, obj_id, XSIZE=xsize, YSIZE=ysize, $
                WAVE=wave, LLIST=llist, QAL=QAL, ERR=err, GAL=gal, INFLG=inflg,$
                BLOCK=block, ZIN=zin, NRM=nrm, FSPEC=fspec, FLUX=flux, SKY=sky

common x_specplot_lines
common esi_echspecplt_cmmn

;
;  print,'Syntax - ' + $
;    'esi_echspecplt, spec_fil, obj_id, XSIZE=,YSIZE=, /NOCTB, WAVE=, '
;  print, '            ERR=, /QAL, /GAL, ZIN=, /BLOCK, /NRM) [v1.0]'

;  Optional Keywords

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
  if not keyword_set( OBJ_ID ) then    obj_id = 'a'

; Input the Obj structure
  if keyword_set( SKY ) then flag_sky = 1 else flag_sky = 0
  if not keyword_set( FSPEC ) then begin
      if keyword_set( spec_fil ) then begin
         objstr = xmrdfits(spec_fil, 1, structyp='dblsobjstrct', /silent)
         ydat = objstr.fx
      endif else begin
          fils = findfile('Extract/Obj*fits*', count=nfil)
          if nfil EQ 0 then return 
          obj_fil = x_guilist(fils)
          objstr = xmrdfits(obj_fil, 1, structyp='dblsobjstrct', /silent)
         ydat = objstr.fx
      endelse 
      gdobj = where(objstr.obj_id EQ obj_id) 
  endif else begin
      if not keyword_set( spec_fil ) then begin
          fils = findfile('FSpec/*_ech.fits*', count=nfil)
          if nfil EQ 0 then return 
          spec_fil = x_guilist(fils)
      endif
      x_wrechfspec, objstr, spec_fil, /read
      gdobj = [0L]
      ydat = objstr[*,0]
   endelse

  tmp1 = { newabslinstrct }

; Init common
  esi_echspecplt_initcommon
  esiobj = objstr[gdobj]
  delvarx, objstr

; ZABS
  if keyword_set( ZIN ) then begin
      flg_lines = 1
      zabs = zin
  endif

;    STATE
  state = { fx: fltarr(5000), $
            wave: dblarr(5000), $
            sig: fltarr(5000), $
            npix: 0L, $
            ordr: 4L, $
            flg_smooth: 0, $   ; Smoothing
            smooth: fltarr(5000), $
            flg_sig: 1, $
            flg_zoom: 0, $
            flg_flux: 0, $
            flg_fspec: 0, $
            flux: 0L, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_EW: 0, $ ; EW flag
            EW_lmt: dblarr(2,2), $
            EW: 0.d, $    
            FWHM: 4., $   ; FWHM of instrument (pix)
            flg_DLA: 0, $ ; DLA flag
            dla_line: tmp1, $ ; DLA flag
            dla_fx: fltarr(n_elements(ydat)) + 1., $
            flg_Colm: 0, $ ; Colm flag
            Colm_lmt: dblarr(2), $
            Colm: 0., $
            sig_Colm: 0., $
            xpos: 0.d, $
            ypos: 0.d, $
            psfile: 0, $ ; Postscript
            flag_sky: flag_sky, $  
            flg_lines: 0, $  ; QAL=1, GAL=2
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            zabs2_id: 0L, $
            ordr_id: 0L, $
            lines_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }

; FSPEC
  if keyword_set( FSPEC ) then state.flg_fspec = 1
  if keyword_set( FLUX ) then state.flux = 1

; LINELIST

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  if keyword_set( QAL ) then state.flg_lines = 1
  if keyword_set( GAL ) then state.flg_lines = 2
  if keyword_set( QSO ) then state.flg_lines = 3
  if keyword_set( LLS ) then state.flg_lines = 4
  if state.flg_lines NE 0 then m_specplot_initLines, state
  
;    WIDGET
  base = WIDGET_BASE( title = 'esi_echspecplt', /column, $
                      xoffset=100L, yoffset=100L)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = '-/= -- Change order number'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'T -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = 'A -- Al III doublet'
  state.help[19] = 'C -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'N -- AODM'
  state.help[22] = 'q -- Quit '


;;;;;;;;;
;  Toolbar

; zabs

  state.zabs_id = cw_field(toolbar, title='zabs', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='ZABS')
  state.zabs2_id = cw_field(toolbar, title='zabs2', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='ZABS2')



; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Objname

  objord_id= widget_base(toolbar, /column, /align_center, frame=2)
  objnm_id = cw_field(objord_id, title='Object: ', value=obj_id, /row, $
                      xsize=2)
  state.ordr_id=cw_field(objord_id, title='Order: ', value=state.ordr, /row, $
                      xsize=2, /long, /return_events, uvalue='ORDR')
  
  
;;;;;;;;;;;;
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
;      Lines
  state.lines_id = WIDGET_LIST(toolbar, $
                             VALUE=['QAL','GAL','QSO','LLS'], $
                             uvalue='LNLIST', ysize = 3)
  widget_control, state.lines_id, set_list_select=state.flg_lines-1
;      Done
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Update
  esi_echspecplt_Reset, state, FLUX=state.flux
  esi_echspecplt_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'esi_echspecplt', base, /no_block $
  else xmanager, 'esi_echspecplt', base

return
end
