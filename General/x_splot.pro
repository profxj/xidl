;+ 
; NAME:
; x_splot   
;     Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   x_splot, [xdat], ydat, XSIZE=, YSIZE=, TITLE=, XTWO=, YTWO=, PSYM_Y2=,
;      /BLOCK
;
; INPUTS:
;   [xdat]     - x values (optional)
;   ydat       - Values 
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   TITLE
;   YTWO
;   PSYM_y2
;   BLOCK
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_splot, y
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP (modified from x1dfit)
;   24-Nov-2001 Added ytwo
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_splot_initcomm

common x_splot_common, xdat, ydat, y2, x2, y3, x3

end

;;;;
; Events
;;;;

pro x_splot_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
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
                      1 : state.xymnx[0] = xgetx_plt(state, /strct) ; left
                      2 : state.xymnx[3] = xgety_plt(state, /strct) ; top
                      4 : state.xymnx[2] = xgetx_plt(state, /strct) ; right
                      else :
                  endcase
              end
              1 : ; Button release
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
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another s !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'N': state.xymnx[3] = 1.1
              'h': state.psym = 10
              'H': x_helpwidg, state.help
              'z': begin  ; Zoom
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, $
                        set_uvalue = state, /no_copy
                      return
                  endif
              end
              'M': begin  ; Median/mean stats
                  xsplot_stats, state
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              '}': xsplot_pan, state, /noy
              '{': xsplot_pan, state, /left, /noy
              ']': xsplot_pan, state
              '[': xsplot_pan, state, /left
              ' ': state.showval = 1 ; Show value
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'W': state.xymnx = state.svxymnx ; Reset the screen
              'P': xsplot_psfile, state  ; Send Current screen to ps file
              'q': begin
                  if keyword_set( X2 ) then delvarx, x2
                  if keyword_set( Y2 ) then delvarx, y2
                  if keyword_set( XDAT ) then delvarx, xdat
                  if keyword_set( YDAT ) then delvarx, ydat
                  widget_control, ev.top, /destroy
                  return
              end
              else:  ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  xsplot_UpdatePlot, state
;
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xsplot_UpdatePlot, state
  
common x_splot_common

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  color = getcolor(/load)

  ;; Where
  a = where(xdat LE state.xymnx[2] AND xdat GE state.xymnx[0], na)
  plot, [0], [0], psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    title=state.title, background=color.white, $
    xcharsize=1.5, $
    ycharsize=1.5, $
    color=color.black, /nodata

  if na NE 0 then $
    oplot, xdat[a], ydat[a], psym=state.psym, color=color.black

; y2

  if state.flg_y2 then begin
      oplot, x2, y2, color=color.red, psym=state.psym2
  endif

; y3

  if state.flg_y3 then begin
      oplot, x3, y3, color=color.blue, psym=state.psym3
  endif

; Show val

  if state.showval EQ 1 then begin
      xyouts, 0.7, 0.97, $
        'Value = '+string(xgetx_plt(state, /strct), $
                          xgety_plt(state, /strct), $
                          format='(2g14.7)'), /NORMAL, charsize=1.5, $
        color=color.black
      state.showval = 0
  endif
      
end

;;;;;;;;;;;;;;;;;;;;
;  Pan
;;;;;;;;;;;;;;;;;;;;

pro xsplot_pan, state, LEFT=left, NOY=noy

common x_splot_common

  ; Size of screen
  sz = state.xymnx[2]-state.xymnx[0]

  ; Right or left
  if not keyword_set( LEFT ) then begin ; RIGHT
      state.xymnx[0] = state.xymnx[2] - 0.1*sz
      state.xymnx[2] = state.xymnx[2] + 0.9*sz
  endif else begin
      state.xymnx[2] = state.xymnx[0] + 0.1*sz
      state.xymnx[0] = state.xymnx[0] - 0.9*sz
  endelse

  ; Reset top and bottom
  if not keyword_set( NOY ) then begin
      gdwv = where(xdat GT state.xymnx[0] AND xdat LT state.xymnx[2], ngd)
      if ngd GT 0 then begin
          mn = min(ydat[gdwv], max=mx)
          state.xymnx[1] = mn - (mx-mn)*0.05
          state.xymnx[3] = mx + (mx-mn)*0.05
      endif
  endif

  return
end

      

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;
pro xsplot_Reset, state
common x_splot_common
  state.ntot = n_elements(ydat)
; Plotting
  state.xymnx = state.svxymnx
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro xsplot_psfile, state

; Device
  device, get_decomposed=svdecomp

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  xsplot_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  Stats
;;;;;;;;;;;;;;;;;;;;
pro xsplot_stats, state
common x_splot_common
  if state.flg_stats EQ 0 then begin
      state.sv_statx = round(xgetx_plt(state,/strct))
      state.flg_stats = 1
  endif else begin
      x2 = round(xgetx_plt(state,/strct))
      x1 = x2 < state.sv_statx
      x2 = x2 > state.sv_statx
      print, 'xsplot: Stats -- ;'
      print, ' median = ', median(ydat[x1:x2])
      print, ' mean = ', mean(ydat[x1:x2])
      state.flg_stats = 0
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

pro x_splot, xin, yin, XSIZE=xsize, YSIZE=ysize, TITLE=title, YTWO=ytwo, $
             PSYM2=psym2, BLOCK=block, XTWO=xtwo, PSYM1=psym1, XYOFF=xyoff, $
             YMNX=ymnx, XTHR=xthr, YTHR=ythr, PSYM3=psym3


common x_splot_common
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_splot, [xdat], ydat, YTWO=ytwo, XSIZE=, YSIZE=, TITLE= '
    print, '            XTWO=, /BLOCK, PSYM_Y2=, PSYM1=, YMNX= [v1.1]'
    return
  endif 

; Set xdat, ydat

  x_splot_initcomm

  if keyword_set( yin ) then begin
      if n_elements(xin) NE n_elements(yin) then begin
          print, 'x_splot: Wrong array sizes!'
          return
      endif
      xdat = xin
      ydat = yin
  endif else begin
      xdat = findgen( n_elements(xin) )
      ydat = xin
  endelse

  if keyword_set( YTWO ) then begin
      y2 = ytwo 
      if keyword_set( XTWO ) then x2 = xtwo else x2 = xdat
  endif

  if keyword_set( YTHR ) then begin
      y3 = ythr 
      if keyword_set( XTHR ) then x3 = xthr else x3 = xdat
  endif
  
;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1200
  if not keyword_set( YSIZE ) then    ysize = 800
  if not keyword_set( PSYM2 ) then    psym2 = 0
  if not keyword_set( PSYM3 ) then    psym3 = 1
  if not keyword_set( PSYM1 ) then    psym1 = 0
  if not keyword_set( XYOFF ) then    xyoff=[300L,400L]

;    STATE
  state = { ntot: n_elements(ydat), $
            reg: fltarr(100,2), $
            flg_y2: 0, $
            flg_y3: 0, $
            psym2: psym2, $
            psym3: psym3, $
            flg_zoom: 0, $
            flg_flip: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            svxymnx: [ min(xdat)-0.01*abs(max(xdat)-min(xdat)), $ ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      max(xdat)+0.01*abs(max(xdat)-min(xdat)), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $  ; ymax
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            flg_stats: 0, $
            sv_statx: 0L, $
            help: strarr(20), $
            showval: 0, $
            psym: psym1, $
            psfile: 0, $
            title: '', $
            size: intarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }
  
; YMNX
  if keyword_set( YMNX ) then begin
      state.svxymnx[1] = ymnx[0]
      state.svxymnx[3] = ymnx[1]
  endif
      
; Flag

  if keyword_set( YTWO ) then state.flg_y2 = 1
  if keyword_set( YTHR ) then state.flg_y3 = 1

; Set xvxymnx[0,2]

  state.svxymnx[0] = min(xdat)
  state.svxymnx[2] = max(xdat)

;    Title
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_splot', /column, xoffset=xyoff[0], $
                    yoffset=xyoff[1])
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)

;        Help
  state.help[0] = '     Help Menu   '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'N -- Set ymax to 1.1'
  state.help[9] = 'h -- Switch to histogram mode'
  state.help[10] = 'H -- Show this screen'
  state.help[11] = 'W -- Reset screen'
  state.help[12] = 'q -- Quit '

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
;      Done
;  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
  
; Update
  xsplot_Reset, state
  xsplot_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set( BLOCK ) then xmanager, 'x_splot', base, /no_block $
    else xmanager, 'x_splot', base


return
end
