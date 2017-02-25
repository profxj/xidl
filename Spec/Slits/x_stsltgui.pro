;+ 
; NAME:
; x_stsltgui   
;    Version 1.0
;
; PURPOSE:
;    Allows the user to interactively (GUI) identify edges of 
;     inidividual slits in the undistorted image
;
; CALLING SEQUENCE:
;   
;   x_stsltgui, img, slitstr, XSIZE=, YSIZE=
;
; INPUTS:
;   img        - Undistorted image
;   slitstr    - Slit structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 700)
;   YSIZE      - Size of gui in screen y-pixels (default = 700)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_stsltgui, img, slitstr
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

pro x_stsltgui_initcommon

 ; colors
 common xcommon_color, r_vector, g_vector, b_vector
 xicmm_colors
 ; images
 common x_stsltgui_images, $
   main_image, $
   img_size, $ 
   display_image, $
   tv_image
 common x_stsltgui_output, flg_out, yedg_flt
 flg_out = 0

end

;;;;
; Events
;;;;

pro x_stsltgui_event, ev

common x_stsltgui_images
common xcommon_color
common x_stsltgui_output

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  if not keyword_set( uval ) then begin
      uval = WIDGET_INFO(ev.id, /uname )
      if uval NE 'BASE' then stop
  endif

  case uval of
      'BASE': begin     
          c = where(tag_names(ev) EQ 'ENTER', count)
          if (count EQ 0) then begin ; resize event
              ximgd_resize, state
              x_stsltgui_ReDisplay, state
          endif
      end
      'INVERT': begin
          r_vector = reverse(r_vector)
          g_vector = reverse(g_vector)
          b_vector = reverse(b_vector)
          ximgd_stretchct, state
          x_stsltgui_PlotImg, state
      end
      'PMIN': begin
          state.pltmin = ev.value
          x_stsltgui_ReDisplay, state, /update
      end
      'PMAX': begin
          state.pltmax = ev.value
          x_stsltgui_ReDisplay, state, /update
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
              0 : begin         ; Button press
                  state.press = ev.press
                  case state.press of
                      1 : x_stsltgui_Contrast, state ; CONTRAST/BRIGHTNESS
                      2 : 
                      4 : begin ; LMB = Center
                          ximgd_zoom, state, 'none', /recenter
                          x_stsltgui_ReDisplay, state
                      end 
                  endcase
              end
              1 : state.zoom.flg = 0 ; Button release
              2 : begin         ; Motion event
                  state.tv.xcurs = ev.x
                  state.tv.ycurs = ev.y
                  if state.zoom.flg EQ 2 then x_stsltgui_Contrast, state $
                  else begin
                      state.xpos = x_tvx(state.tv, /intg)
                      state.ypos = x_tvy(state.tv, /intg)
                      widget_control, state.xpos_id, set_value=state.xpos
                      widget_control, state.ypos_id, set_value=state.ypos
                  endelse
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          case eventch of
              ; ZOOM
              'z': begin
                  ximgd_zoom, state, 'in', /recenter
                  x_stsltgui_ReDisplay, state
              end
              'Z': begin
                  ximgd_zoom, state, 'out', /recenter
                  x_stsltgui_ReDisplay, state
              end
              ; CENTER
              'c': begin
                  ximgd_zoom, state, 'none', /recenter
                  x_stsltgui_ReDisplay, state
              end
              ; BIG SHIFT
              's': begin  ; Key on Top
                  x_stsltgui_Shift, state
                  x_stsltgui_ReDisplay, state
              end
              'S': begin  ; Key on Bottom
                  x_stsltgui_Shift, state, /bottom
                  x_stsltgui_ReDisplay, state
              end
              ; Single shift
              'b': begin
                  x_stsltgui_SngShift, state, /bottom
                  x_stsltgui_ReDisplay, state
              end
              't': begin
                  x_stsltgui_SngShift, state
                  x_stsltgui_ReDisplay, state
              end
              else: 
          endcase
      end
      'CENTER' : begin
          state.zoom.centerpix = round(state.sz_img / 2.)
          ximgd_getoffset, state
          x_stsltgui_ReDisplay, state
      end
      'DONE' : begin
          flg_out = 1
          yedg_flt = state.slitstr.yedg_flt
          widget_control, ev.top, /destroy
          return
      end
      'NOSV' : begin
          flg_out = 0
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Update display and plot
;;;;;;;;;;;;;;;;;;;;


pro x_stsltgui_ReDisplay, state, UPDATE=update
 
  ; 
  ximgd_stretchct, state
  ; Dont always update display image
  if keyword_set( UPDATE ) then x_stsltgui_UpdDisplay, state
  x_stsltgui_UpdTVImg, state
  x_stsltgui_PlotImg, state
  return
end


;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_stsltgui_PlotImg, state
  
common x_stsltgui_images

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  tv, tv_image


; Slits

  xrange=[state.tv.xymnx[0], state.tv.xymnx[2]]
  yrange=[state.tv.xymnx[1], state.tv.xymnx[3]]
  
  clr = getcolor(/load)

  ; TOP
  plot, [0., 10000.], [state.slitstr[0].yedg_flt[1],$
                       state.slitstr[0].yedg_flt[1]], $
    position=[0.,0.,1., 1.], $
    xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /noerase, $
    color=clr.red
  for i=1L,state.nslit-1 do $
    oplot, [0., 10000.], [state.slitstr[i].yedg_flt[1],$
                          state.slitstr[i].yedg_flt[1]], $
    color=clr.red

  ; BOTTOM
  for i=0L,state.nslit-1 do $
    oplot, [0., 10000.], [state.slitstr[i].yedg_flt[0],$
                          state.slitstr[i].yedg_flt[0]], $
    color=clr.blue

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_stsltgui_Reset, state

; Plotting

  state.brightness = 0.5
  state.contrast = 0.5
  ximgd_stretchct, state

  ; ZOOM
  state.zoom.centerpix = round(state.sz_img / 2.)

  factor = (float(state.tv.winsize[0])/state.sz_img[0]) < $
    (float(state.tv.winsize[1])/state.sz_img[1])

  if factor LT 1. then state.zoom.level = -round(2.*(1./factor - 1.)) $
    else stop
  ximgd_zoom, state, 'none'
  ximgd_getoffset, state

end

;;;;;;;;;;;;;;;;;;;;
;  InitImage
;;;;;;;;;;;;;;;;;;;;

pro x_stsltgui_InitImg, state

common x_stsltgui_images

  widget_control, /hourglass

  ; Size
  img_size = size(main_image, /dimensions)
  state.sz_img = img_size

  ; Set min max
  state.imgmin = min(main_image, max=fmax)
  state.imgmax = fmax

  med = median(main_image, /even)
  sig = stddev(main_image)

  state.pltmax = (med + (2 * sig)) < state.imgmax
  state.pltmin = (med - (2 * sig))  > state.imgmin

  widget_control, state.pmin_id, set_value=state.pltmin
  widget_control, state.pmax_id, set_value=state.pltmax

  ; ZOOMREG
  state.zoom.reg = [0L, 0, state.sz_img[0]-1, state.sz_img[1]-1]
;  state.tv.xymnx = float(state.zoom.reg)
;  x_setgridxy, state, state.zoom.reg, /FILL

; Display Image
;  display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax)
;  tv_image = congrid(display_image, state.tv.winsize[0], state.tv.winsize[1])

end


;;;;;;;;;
;  Update tv image
;;;;;;;;;

pro x_stsltgui_UpdTVimg, state

common x_stsltgui_images

  widget_control, /hourglass

; Reset Display

  tv_image = bytarr(state.tv.winsize[0], state.tv.winsize[1])

  view_min = round(state.zoom.centerpix - $
                   (0.5 * state.tv.winsize / state.zoom.factor))
  view_max = round(view_min + state.tv.winsize / state.zoom.factor)
  
  view_min = (0 > view_min < (state.sz_img - 1)) 
  view_max = (0 > view_max < (state.sz_img - 1)) 
  
  newsize = round( (view_max - view_min + 1) * state.zoom.factor) > 1
  startpos = abs( round(state.zoom.offset * state.zoom.factor) < 0)

  tmp_image = congrid(display_image[view_min[0]:view_max[0], $
                                    view_min[1]:view_max[1]], $
                      newsize[0], newsize[1])
  ; xymnx
  dum = [0., 0.]
  state.tv.xymnx[0:1] = (dum / state.zoom.factor) + state.zoom.offset
  dum = [state.tv.winsize[0], state.tv.winsize[1]]
  state.tv.xymnx[2:3] = (dum / state.zoom.factor) + state.zoom.offset

  xmax = newsize[0] < (state.tv.winsize[0] - startpos[0])
  ymax = newsize[1] < (state.tv.winsize[1] - startpos[1])

  tv_image[startpos[0],startpos[1]] = tmp_image[0:xmax-1, 0:ymax-1]
  delvarx, tmp_image

end

;;;;;;;;;
;  Update Display image
;;;;;;;;;

pro x_stsltgui_UpdDisplay, state

common x_stsltgui_images

  widget_control, /hourglass

  display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax)

end

;;;;
; Contrast/Brightness
;;;;

pro x_stsltgui_Contrast, state

  ; Set value
  state.contrast = state.tv.xcurs/$
    float(state.tv.winsize[0])
  state.brightness = state.tv.ycurs/$
    float(state.tv.winsize[1])
  ximgd_stretchct, state
  x_stsltgui_PlotImg, state
  state.zoom.flg = 2

  return

end
  
  
;;;;
; Shift 
;;;;

pro x_stsltgui_Shift, state, BOTTOM=bottom

  ; Move all lower edges down by same amount to align on top/bottom

  if not keyword_set( BOTTOM ) then begin
      ; Find closest top
      mn = min(abs(state.slitstr.yedg_flt[1]-state.ypos), imn)
      ; Shift
      shft = state.slitstr[imn].yedg_flt[1]-state.ypos
      ; Move em all
      lwr = where(state.slitstr.yedg_flt[1] LE state.slitstr[imn].yedg_flt[1], nlwr)
      if nlwr NE -1 then state.slitstr[lwr].yedg_flt = $
        state.slitstr[lwr].yedg_flt-shft
  endif else begin
      ; Find closest top
      mn = min(abs(state.slitstr.yedg_flt[0]-state.ypos), imn)
      ; Shift
      shft = state.slitstr[imn].yedg_flt[0]-state.ypos
      ; Move em all
      upr = where(state.slitstr.yedg_flt[0] LE state.slitstr[imn].yedg_flt[0], nupr)
      if nupr NE -1 then $
        state.slitstr[upr].yedg_flt = state.slitstr[upr].yedg_flt-shft
  endelse

  return

end

;;;;
; Single Shift 
;;;;

pro x_stsltgui_SngShift, state, BOTTOM=bottom

  ; Move all lower edges down by same amount to align on top/bottom

  if not keyword_set( BOTTOM ) then begin
      ; Find closest top
      mn = min(abs(state.slitstr.yedg_flt[1]-state.ypos), imn)
      ; Shift
      shft = state.slitstr[imn].yedg_flt[1]-state.ypos
      state.slitstr[imn].yedg_flt = state.slitstr[imn].yedg_flt-shft
  endif else begin
      ; Find closest top
      mn = min(abs(state.slitstr.yedg_flt[0]-state.ypos), imn)
      ; Shift
      shft = state.slitstr[imn].yedg_flt[0]-state.ypos
      state.slitstr[imn].yedg_flt = state.slitstr[imn].yedg_flt-shft
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

pro x_stsltgui, img, slitstr, XSIZE=xsize, YSIZE=ysize

common xcommon_color
common x_stsltgui_images
common x_stsltgui_output

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'x_stsltgui, img, slitstr, XSIZE=, YSIZE=,'
    print, '      (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 700
  if not keyword_set( YSIZE ) then    ysize = 700

;    STATE

  ; tv structure
  tmp = {tvstruct}
  tmp.pos = [0., 0., 1., 1.]
  tmp2 = {zoomstrct}

  state = {             $
            img: img, $                  ; Important stuff
            slitstr: slitstr, $          ; Slit structure
            nslit: n_elements(slitstr), $
            title: '', $
            imgmin: 0., $                ; Image stuff
            imgmax: 0., $
            sz_img: lonarr(2), $
            zoom: tmp2, $
            tmpreg: lonarr(4), $         
            tv: tmp, $              ; TV structure
            xpos: 0., $
            ypos: 0., $
            ncolors: 0, $
            brightness: 0.0, $
            contrast: 0.0, $
            pltmax: 0.0, $
            pltmin: 0.0, $
            press: 0, $
            base_id: 0L, $      ; Widgets
            lbltrc_id: 0L, $
            drag_contr_id: 0L, $
            drag_brght_id: 0L, $
            drawbase_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            pmin_id: 0L, $
            pmax_id: 0L, $
            name_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            error_msg_id: 0L, $
            help_text_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'x_stsltgui: Check traces', /column, $
                    UNAME='BASE', /tlb_size_events)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_stsltgui', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Img max/min

  pltmnxbase = widget_base(toolbar, /column, /base_align_center, /align_center, $
                        /frame)
  state.pmin_id = cw_field(pltmnxbase, value=0., /floating, title='pmin', $
                          /return_events, xsize=8, UVALUE='PMIN') 
  state.pmax_id = cw_field(pltmnxbase, value=0., /floating, title='pmax', $
                          /return_events, xsize=8, UVALUE='PMAX') 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  state.drawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='DRAW_BASE', frame=2)

  state.tv.winsize[0] = xsize
  state.tv.winsize[1] = ysize
  state.draw_id = widget_draw(state.drawbase_id, xsize=state.tv.winsize[0], $
                              ysize=state.tv.winsize[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; TEXT
  state.text_id = widget_text(state.drawbase_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0L, xsize=5, /long)
  state.ypos_id = cw_field(xy_id, title='y:', value=0L, xsize=5, /long)

;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
  center = WIDGET_BUTTON(butbase, value='CENTER',uvalue='CENTER')
  nosv = WIDGET_BUTTON(butbase, value='NOSV',uvalue='NOSV')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB - Contrast/Brightness', $
             'RMB - center', $ 
             'z/Z - Zoom/Unzoom', $ 
             'c - center', $
             's - Shift all (key on top of slit)', $
             'S - Shift all (key on bottom of slit)', $
             'b - Shift one slit (key on bottom of slit)', $
             't - Shift one slit (key on top of slit)' $
             ]
  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
                             /scroll)

; Realize
  WIDGET_CONTROL, base, /realize

  ; COLORS 
  device, get_decomposed=val_decomp
  device, decompose=0
;  device, true=8
  loadct, 0, /silent
  if (!d.table_size LT 12) then begin
      message, 'x_stsltgui: Too few colors available for color table'
      stop
  endif
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  ximgd_getct, state, 0



  ; IMAGE
  tv_image = bytarr(xsize, ysize)
  main_image = img

  ; PLOT

  x_stsltgui_InitImg, state
  x_stsltgui_Reset, state
  x_stsltgui_ReDisplay, state, /update
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_stsltgui', base

  delvarx, main_image, tv_image, display_image, trc_img
  device, decomposed=val_decomp


  ; OUTPUT
  if flg_out EQ 1 then slitstr.yedg_flt = temporary(yedg_flt)

  return
end


