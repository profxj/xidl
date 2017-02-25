;+ 
; NAME:
; x_prspeaks   
;     Version 1.1
;
; PURPOSE:
;    Launches a GUI to enable the user to fiddle with peaks by hand
;
; CALLING SEQUENCE:
;   
;   x_prspeaks, [xdat], ydat, XSIZE=, YSIZE=, TITLE=, XTWO=, YTWO=, PSYM_Y2=,
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
;   x_prspeaks, y
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

pro x_prspeaks_initcomm

common x_prspeaks_common, xdat, ydat, y2, x2, msk

end

;;;;
; Events
;;;;

pro x_prspeaks_event, ev

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
                set_value='Expecting another z !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
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
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              '}': xprspeaks_pan, state, /noy
              '{': xprspeaks_pan, state, /left, /noy
              ']': xprspeaks_pan, state
              '[': xprspeaks_pan, state, /left
              ' ': state.showval = 1 ; Show value
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'W': state.xymnx = state.svxymnx ; Reset the screen
              'P': xprspeaks_psfile, state  ; Send Current screen to ps file
              'd': xprspeaks_delpeak, state  ; Delete peak
              'u': xprspeaks_undelpeak, state
              'n': xprspeaks_newpeak, state  ; Add new peak
              'N': xprspeaks_newpeak, state  ; Add new peak
              'q': begin
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
  xprspeaks_UpdatePlot, state
;
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xprspeaks_UpdatePlot, state
  
common x_prspeaks_common

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  color = getcolor(/load)

  ;; Where
  a = where(xdat LE state.xymnx[2] AND xdat GT state.xymnx[0], na)

  plot, [0], [0], psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    title=state.title, background=color.white, $
    xcharsize=1.5, $
    ycharsize=1.5, $
    color=color.black, /nodata

  if na NE 0 then $
    oplot, xdat[a], ydat[a], psym=state.psym, color=color.black

; Peaks

  oplot, x2[where(msk EQ 1)], y2[where(msk EQ 1)], $
    color=color.blue, psym=state.psym_y2

  del = where(msk EQ 0, ndel)
  if ndel NE 0 then $
    oplot, x2[del], y2[del], color=color.red, psym=1
  

; Show val

  if state.showval EQ 1 then begin
      xyouts, 0.7, 0.97, $
        'Value = '+string(xgetx_plt(state, /strct), $
                          format='(g14.7)'), /NORMAL, charsize=1.5, $
        color=color.black
      state.showval = 0
  endif
      
end

;;;;;;;;;;;;;;;;;;;;
;  Pan
;;;;;;;;;;;;;;;;;;;;

pro xprspeaks_pan, state, LEFT=left, NOY=noy

common x_prspeaks_common

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
;  Delete Peak
;;;;;;;;;;;;;;;;;;;;

pro xprspeaks_newpeak, state

common x_prspeaks_common

  widget_control, /hourglass   
  ;; FIT
  autofit = x_fitrej(xdat, ydat, 'BSPLIN', 31L, $
                     hsigma=2., lsigma=4.)
  ;; xval
  xval = xgetx_plt(state, /strct)

  npix = n_elements(xdat)
  ;; Center up
  peak = round(xval)
  pmn = (peak-4L) > 0
  pmx = (peak+4L) < (npix-1)
  ;; Center
  center = x_centspln(xdat[pmn:pmx], $
                      ydat[pmn:pmx]-autofit[pmn:pmx], 0.3, /SILENT, /FORCE)

  x2 = [x2,center]
  msk = [msk,1]

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Delete Peak
;;;;;;;;;;;;;;;;;;;;

pro xprspeaks_delpeak, state

common x_prspeaks_common

  ;; xval
  xval = xgetx_plt(state, /strct)
  ;; Find nearest peak
  gdpk = where(msk EQ 1)
  mn = min(abs(x2[gdpk]-xval),imn)
  msk[gdpk[imn]] = 0

  return
end

;;;;;;;;;;;;;;;;;;;;
;  UnDelete Peak
;;;;;;;;;;;;;;;;;;;;

pro xprspeaks_undelpeak, state

common x_prspeaks_common

  ;; xval
  xval = xgetx_plt(state, /strct)

  ;; Find nearest peak
  delpk = where(msk EQ 0, ndel)
  if ndel EQ 0 then return
  mn = min(abs(x2[delpk]-xval),imn)
  msk[delpk[imn]] = 1

  return
end

      

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;
pro xprspeaks_Reset, state
common x_prspeaks_common
  state.ntot = n_elements(ydat)
; Plotting
  state.xymnx = state.svxymnx
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro xprspeaks_psfile, state

; Device
  device, get_decomposed=svdecomp

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  xprspeaks_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0

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

pro x_prspeaks, yin, xtwo, outmsk, XIN=xin, XSIZE=xsize, YSIZE=ysize, TITLE=title, $
                YTWO=ytwo, PSYM_Y2=psym_y2, BLOCK=block, PSYM1=psym1, XYOFF=xyoff, $
                YMNX=ymnx


common x_prspeaks_common
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_prspeaks, ydat, peaks, [msk], YTWO=ytwo, XSIZE=, YSIZE=, TITLE= '
    print, '            XTWO=, /BLOCK, PSYM_Y2=, PSYM1=, YMNX= [v1.1]'
    return
  endif 

; Set xdat, ydat

  x_prspeaks_initcomm

  ydat = yin
  if keyword_set( xin ) then begin
      if n_elements(xin) NE n_elements(yin) then begin
          print, 'x_prspeaks: Wrong array sizes!'
          return
      endif
      xdat = xin
  endif else xdat = findgen( n_elements(yin) )

  x2 = xtwo
  npk = n_elements(xtwo)
  if keyword_set( YTWO ) then y2 = ytwo else y2 = fltarr(npk)

  msk = bytarr(npk) + 1B
  
;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1000
  if not keyword_set( YSIZE ) then    ysize = 600
  if not keyword_set( PSYM_Y2 ) then    psym_y2 = 2
  if not keyword_set( PSYM1 ) then    psym1 = 10
  if not keyword_set( XYOFF ) then    xyoff=[300L,400L]

;    STATE
  state = { ntot: n_elements(ydat), $
            reg: fltarr(100,2), $
            flg_y2: 0, $
            psym_y2: psym_y2, $
            npk: npk, $
            flg_zoom: 0, $
            flg_flip: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            svxymnx: [ min(xdat)-0.01*abs(max(xdat)-min(xdat)), $ ; xmin
                      -10., $ ; ymin
                      max(xdat)+0.01*abs(max(xdat)-min(xdat)), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $  ; ymax
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
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

; Set xvxymnx[0,2]

  state.svxymnx[0] = min(xdat)
  state.svxymnx[2] = max(xdat)

;    Title
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_prspeaks', /column, xoffset=xyoff[0], $
                    yoffset=xyoff[1])
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)

;        Help
  state.help[0] = '     Help Menu   '
  state.help[1] = 'LMB -- Set left'
  state.help[2] = 'RMB -- Set right'
  state.help[3] = 'CMB -- Set top'
  state.help[4] = 'l -- Set Left '
  state.help[5] = 'r -- Set Right '
  state.help[6] = 'b -- Set Bottom '
  state.help[7] = 't -- Set Top '
  state.help[8] = 'z/z -- Zoom'
  state.help[9] = 'N -- Set ymax to 1.1'
  state.help[10] = 'i/o -- Zoom in/out'
  state.help[11] = '{[]} -- Pan left, right'
  state.help[12] = 'd -- Delete peak'
  state.help[13] = 'u -- UnDelete peak'
  state.help[14] = 'N -- Add new peak'
  state.help[15] = 'H -- Show this screen'
  state.help[16] = 'W or w -- Reset screen'
  state.help[17] = 'q -- Quit '
  state.help[18] = 'h -- Switch to histogram mode'

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
;  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
  
; Update
  xprspeaks_Reset, state
  xprspeaks_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set( BLOCK ) then xmanager, 'x_prspeaks', base, /no_block $
    else xmanager, 'x_prspeaks', base

  delvarx, y2, xdat, ydat
  xtwo = temporary(x2)
  outmsk = temporary(msk)

return
end
