;+ 
; NAME:
; x_setobjgui
;    Version 1.0
;
; PURPOSE:
;    Allows the user to interactively examine the trace to slit
;     edges used to trace the y-distortion
;
; CALLING SEQUENCE:
;   
;   x_setobjgui, img, tracestrct, [newmnx], XSIZE=, YSIZE=
;
; INPUTS:
;   img        - Traced Image 
;   trace      - x,y positions of traces fltarr(N,M,2)
;   flgs       - Flags describing the traces
;
; RETURNS:
;
; OUTPUTS:
;  newmnx      - The new min and max values
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_setobjgui, img, tracestrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

pro x_setobjgui_initcommon

;

common x_setobjgui_color, r_vector, g_vector, b_vector
common x_setobjgui_images, $
  main_image, $
  img_size, $ 
  slit_img, $
  slit_size, $ 
  display_image, $
  tv_image
common x_setobjgui_cmmn, sobj_slit, sobj_obj

end

;;;;
; Events
;;;;

pro x_setobjgui_event, ev

common x_setobjgui_images
common x_setobjgui_color

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
              x_setobjgui_UpdDisplay, state
              x_setobjgui_stretchct, state
              x_setobjgui_PlotImg, state
          endif
      end
      'INVERT': begin
          r_vector = reverse(r_vector)
          g_vector = reverse(g_vector)
          b_vector = reverse(b_vector)
          x_setobjgui_stretchct, state
          x_setobjgui_PlotImg, state
      end
      'PMIN': begin
          state.pltmin = ev.value
          x_setobjgui_UpdDisplay, state
          x_setobjgui_PlotImg, state
      end
      'PMAX': begin
          state.pltmax = ev.value
          x_setobjgui_UpdDisplay, state
          x_setobjgui_PlotImg, state
      end
      'DRAW' : begin
          case ev.type of
              0 : begin         ; Button press
                  state.press = ev.press
                  if state.zoom EQ 1 then begin ; Zoom overrides everything
                      x_setobjgui_SetZoom, state, 2
                      x_setobjgui_UpdDisplay, state
                      x_setobjgui_PlotImg, state
                      state.zoom = 0
                  endif else begin
                      case state.press of
                          1 : begin ; LMB = Cut
                              x_setobjgui_Cutext, state
                              x_setobjgui_stretchct, state
                              x_setobjgui_PlotImg, state
                          end 
                          2 : x_setobjgui_SetZoom, state, 1 ; Set Zoom region
                          4 : x_setobjgui_Contrast, state ; CONTRAST/BRIGHTNESS
                      endcase
                  endelse
              end
              1 : if state.zoom EQ 2 then state.zoom = 0 ; Button release
              2 : begin         ; Motion event
                  state.tv.xcurs = ev.x
                  state.tv.ycurs = ev.y
                  if state.zoom EQ 2 then x_setobjgui_Contrast, state 
              end
          endcase
      end
      'UNZOOM' : begin
          state.zoomreg = [0,0,slit_size[0]-1,slit_size[1]-1]
          state.tv.xymnx = float(state.zoomreg)
          x_setgridxy, state, state.zoomreg, /FILL
          x_setobjgui_UpdDisplay, state
          x_setobjgui_PlotImg, state
      end
      'NEXT' : begin
          state.curslit = state.curslit + 1
          if state.curslit GT state.nslit-1 then state.curslit = 0
          widget_control, state.lblslit_id, set_value=state.curslit
          x_setobjgui_NewSlitImg, state
          x_setobjgui_stretchct, state
          x_setobjgui_PlotImg, state
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_setobjgui_PlotImg, state
  
common x_setobjgui_images
common x_setobjgui_cmmn

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  tv, tv_image


; TRACE

  xrange=[state.tv.xymnx[0], state.tv.xymnx[2]]
  yrange=[state.tv.xymnx[1], state.tv.xymnx[3]]
  
  clr = getcolor(/load)

; SCIENCE OBJECT
  sci = state.obj[where(sobj_obj[state.obj[0:state.nobj-1]].obj_id EQ 'a')]
  ; Point
  plot, [sobj_obj[sci].xcen], [sobj_obj[sci].ycen]-state.ydelt, $
    position=[0.,0.,1., 1.], $
    xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /noerase, $
    psym=1, color=clr.green
  ; Trace
  oplot, findgen(state.sz[0]), sobj_obj[sci].trace[0:state.sz[0]-1]$
    -state.ydelt, psym=3, color=clr.green

  ; EXCLUDED
; if state.trcstr.minx[state.curtrc] GT 0 then $
;   oplot, $
;   state.trcstr.xcen[0:state.trcstr.minx[state.curtrc],$
;                     state.curtrc], $
;   state.trcstr.ycen[0:state.trcstr.minx[state.curtrc],$
;                     state.curtrc]-state.ydelt, $
;   psym=3, color=clr.green
;    psym=3, color=2
; if state.trcstr.maxx[state.curtrc] LT state.sz[0] then $
;   oplot, state.trcstr.xcen[state.trcstr.maxx[state.curtrc]:state.sz[0]-1,$
;                            state.curtrc], $
;   state.trcstr.ycen[state.trcstr.maxx[state.curtrc]:state.sz[0]-1,$
;                     state.curtrc]-state.ydelt, $
;   psym=3, color=clr.green
;    psym=3, color=2
  

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_setobjgui_Reset, state

; Plotting

  state.brightness = 0.5
  state.contrast = 0.5
  x_setobjgui_stretchct, state

end

;;;;;;;;;;;;;;;;;;;;
;  InitImage
;;;;;;;;;;;;;;;;;;;;

pro x_setobjgui_InitImg, state

common x_setobjgui_images

  widget_control, /hourglass

  img_size = size(main_image, /dimensions)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Update trace image
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_setobjgui_NewSlitImg, state

common x_setobjgui_images
common x_setobjgui_cmmn

; Set Region

  mn = min(sobj_slit[state.curslit].yedg_orig[0:state.sz[0]-1,0])
  mx = min(sobj_slit[state.curslit].yedg_orig[0:state.sz[0]-1,1])
  
  state.slitreg = [0, mn, state.sz[0]-1, mx]

  slit_img = main_image[state.slitreg[0]:state.slitreg[2],$
                        state.slitreg[1]:state.slitreg[3]] 

  ; Set min max
  slit_size = size(slit_img, /dimensions)

  state.imgmin = min(slit_img, max=fmax)
  state.imgmax = fmax
  med = median(slit_img, /even)
  state.pltmax = med*1.1 + 300
  state.pltmin = med*0.8 - 200

  widget_control, state.pmin_id, set_value=state.pltmin
  widget_control, state.pmax_id, set_value=state.pltmax

;  state.pltmax = (med + 200) < state.imgmax
;  state.pltmin = (med - 100) > state.imgmin

  ; ZOOMREG
  state.zoomreg = [0L, 0, slit_size[0]-1, slit_size[1]-1]
  state.tv.xymnx = float(state.zoomreg)
  x_setgridxy, state, state.zoomreg, /FILL

; Display Image
  display_image = bytscl(slit_img, min=state.pltmin, max=state.pltmax)
  tv_image[*] = congrid(display_image, state.tv.winsize[0], state.tv.winsize[1])

; Set ydelt

  state.ydelt = state.slitreg[1]

; Name
  widget_control, state.lblslit_id, set_value=state.curslit

; Objects
  obj = where(sobj_obj.slit_id EQ sobj_slit[state.curslit].id, nobj)
  state.nobj = nobj
  state.obj[0:nobj-1] = obj

end


;;;;;;;;;
;  Update Display
;;;;;;;;;

pro x_setobjgui_UpdDisplay, state

common x_setobjgui_images

  widget_control, /hourglass

; Reset Display

  delvarx, display_image
  display_image = bytscl(slit_img[state.zoomreg[0]:state.zoomreg[2], $
                                 state.zoomreg[1]:state.zoomreg[3]], $
                                    min=state.pltmin, max=state.pltmax)
  tv_image = congrid(display_image, state.tv.winsize[0], state.tv.winsize[1])
end

;;;;;;;;;
;  Set Zoom 
;;;;;;;;;

pro x_setobjgui_SetZoom, state, flg

common x_setobjgui_images

  if flg EQ 1 then begin
      state.tmpreg[0] = x_tvx(state.tv, /intg)
      state.tmpreg[1] = x_tvy(state.tv, /intg)
      state.zoom = 1
  endif else begin
      state.tmpreg[2] = x_tvx(state.tv, /intg)
      state.tmpreg[3] = x_tvy(state.tv, /intg)

      if (state.tmpreg[0] NE state.tmpreg[2] AND $
          state.tmpreg[1] NE state.tmpreg[3]) then begin
          state.zoomreg[0] = state.tmpreg[0] < state.tmpreg[2]
          state.zoomreg[2] = state.tmpreg[0] > state.tmpreg[2]
          state.zoomreg[1] = state.tmpreg[1] < state.tmpreg[3]
          state.zoomreg[3] = state.tmpreg[1] > state.tmpreg[3]

          state.zoomreg[0] = state.zoomreg[0] > 0
          state.zoomreg[2] = state.zoomreg[2] > 0
          state.zoomreg[1] = state.zoomreg[1] < slit_size[0]-1
          state.zoomreg[3] = state.zoomreg[3] < slit_size[1]-1
      endif
      x_setgridxy, state, state.zoomreg, /FILL
  endelse
  state.tv.xymnx = float(state.zoomreg)

end

;;;;
; Contrast/Brightness
;;;;

pro x_setobjgui_Contrast, state

  ; Set value
  state.contrast = state.tv.xcurs/$
    float(state.tv.winsize[0])
  state.brightness = state.tv.ycurs/$
    float(state.tv.winsize[1])
  x_setobjgui_stretchct, state
  x_setobjgui_PlotImg, state
  state.zoom = 2

  return

end
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_setobjgui_initcolors, state

common x_setobjgui_color

; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.

;rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
;gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
;btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
;tvlct, 255*rtiny, 255*gtiny, 255*btiny

tvlct, [255],[255],[255], !d.table_size-1

end

;--------------------------------------------------------------------

pro x_setobjgui_getct, state, tablenum

common x_setobjgui_color

; Read in a pre-defined color table, and invert if necessary.

loadct, tablenum, /silent, bottom=8
tvlct, r, g, b, /get
;tvlct, r, g, b, 8, /get

x_setobjgui_initcolors, state

r = r[0:state.ncolors-2]
g = g[0:state.ncolors-2]
b = b[0:state.ncolors-2]

;r = reverse(r)
;g = reverse(g)
;b = reverse(b)

r_vector = r
g_vector = g
b_vector = b

x_setobjgui_stretchct, state
;if (state.bitdepth EQ 24 AND (n_elements(pan_image) GT 10) ) then $
;  atv_refresh

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_setobjgui_stretchct, state

; routine to change color stretch for given values of 
; brightness and contrast.
; Complete rewrite 2000-Sep-21 - Doug Finkbeiner
; This routine is now shorter and easier to understand.  

; if GETMOUSE then assume mouse positoin passed; otherwise ignore
; inputs

common x_setobjgui_color

x = state.brightness*(state.ncolors-1)
y = state.contrast*(state.ncolors-1) > 2   ; Minor change by AJB 
high = x+y & low = x-y
diff = (high-low) > 1

slope = float(state.ncolors-1)/diff ;Scale to range of 0 : nc-1
intercept = -slope*low
p = long(findgen(state.ncolors)*slope+intercept) ;subscripts to select
tvlct, r_vector[p], g_vector[p], b_vector[p], 8

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

pro x_setobjgui, img, slitstr, objstr, XSIZE=xsize, YSIZE=ysize

common x_setobjgui_color
common x_setobjgui_images
common x_setobjgui_cmmn

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'x_setobjgui, img, slitstr, objstr, XSIZE=, YSIZE=,'
    print, '      (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1200
  if not keyword_set( YSIZE ) then    ysize = 600

;    STATE
  
  x_setobjgui_initcommon
  sobj_slit = slitstr
  sobj_obj = objstr

  ; tv structure
  tmp = {tvstruct}
  tmp.pos = [0., 0., 1., 1.]

  state = {             $
            img: img, $                  ; Important stuff
            sz: size(img, /dimensions), $
            nslit: n_elements(sobj_slit), $
            curslit: 0L, $
            nobj: 0, $
            obj: lonarr(100), $
            title: '', $
            imgmin: 0., $                ; Image stuff
            imgmax: 0., $
            yoff: 15L, $
            posneg: 0, $
            ydelt: 0., $  ; Offset between display and image
            zoom: 0, $
            zoomreg: lonarr(4), $         
            slitreg: lonarr(4), $         
            tmpreg: lonarr(4), $         
            tv: tmp, $              ; TV structure
            ncolors: 0, $
            brightness: 0.0, $
            contrast: 0.0, $
            pltmax: 0.0, $
            pltmin: 0.0, $
            press: 0, $
            base_id: 0L, $      ; Widgets
            lblslit_id: 0L, $
            drag_contr_id: 0L, $
            drag_brght_id: 0L, $
            draw_id: 0L, $
            drawbase_id: 0L, $
            xmin_id: 0L, $
            xmax_id: 0L, $
            pmin_id: 0L, $
            pmax_id: 0L, $
            name_id: 0L, $
            error_msg_id: 0L, $
            help_text_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'x_setobjgui: Check traces', /column, $
                    UNAME='BASE', /tlb_size_events)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_setobjgui', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)
  nslitlbl = WIDGET_LABEL(labelbase, value=strtrim(state.nslit-1,2), /align_center)
  state.lblslit_id = cw_field(labelbase, value=state.curslit, $
                             /long, /return_events, xsize=7,$
                             title='Trace')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Trace ends

  endsbase = widget_base(toolbar, /column, /base_align_center, /align_center, $
                        /frame)
  state.xmin_id = cw_field(endsbase, value=0L, /long, title='xmin', $
                          /return_events, xsize=8, UVALUE='XMIN') 
  state.xmax_id = cw_field(endsbase, value=state.sz[0]-1, /long, title='xmax', $
                          /return_events, xsize=8, UVALUE='XMAX') 

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
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.tv.winsize[0] = xsize
  state.tv.winsize[1] = ysize
  state.draw_id = widget_draw(state.drawbase_id, xsize=state.tv.winsize[0], $
                              ysize=state.tv.winsize[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')

;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
  zoombt = WIDGET_BUTTON(butbase, value='UNZOOM',uvalue='UNZOOM')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB - Truncate/Extend trace', $ 
             'RMB - Contrast/Brightness', $
             'CMB/CMB - Zoom' $ 
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
      message, 'x_setobjgui: Too few colors available for color table'
      stop
  endif
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  x_setobjgui_getct, state, 0

  x_setobjgui_Reset, state


  ; IMAGE
  tv_image = bytarr(xsize, ysize)
  main_image = img

  ; PLOT
  x_setobjgui_InitImg, state
  x_setobjgui_NewSlitImg, state
  x_setobjgui_PlotImg, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_setobjgui', base

  delvarx, main_image, tv_image, display_image, slit_img
  device, decomposed=val_decomp

  if arg_present(newmnx) then newmnx = mnx

  ; Update 

  objstr = temporary(sobj_obj)

  return
end

