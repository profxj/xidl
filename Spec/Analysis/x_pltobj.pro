;+ 
; NAME:
; x_pltobj
;    Version 1.1
;
; PURPOSE:
;   Sophisticated GUI for plotting 1D and 2D spectra together.
;
; CALLING SEQUENCE:
;  x_pltobj, wave, fx, sig, infx, inwv, [imgwv], YSIZE=
;             XSIZE=, PMNX=, YSIZE=, OBJNM=, ZIN=, XMAX=, LLIST=
;
; INPUTS:
; wave -- 1D Wavelength array
; fx   -- 1D flux array
; sig  -- 1D sigma array
; infx -- 2D flux image
; inwv -- 2D wavelength image
; [imgwv] -- 1D wavelength array giving approximate wave of the 2D
;            image
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XSIZE= -- Size of gui in screen x-pixels [default = ssz-300]
;  YSIZE= -- Size of gui in screen y-pixels [default = ssz-400]
;  PMNX=  -- Plot range of the 1D spectrum
;  ZIN=   -- Input redshift for the object
;  LLIST= -- Line list to use upon initialization (1: QAL, 2: GAL, 3:
;            QSO)
;  OBJNM= -- Title for the plot
;  XMAX=  -- Maximum x value for the 1D plot
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_pltobj, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

pro x_pltobj_icmmn


  ; COLOR
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  ; Images
  common x_pltobj_images, $
    subfx, $
    subwv, $
    img_size, $ 
    main_image, $
    display_image, $
    tv_image, $
    pan_image

  common x_specplot_lines, $
    flg_lines, $
    lines, $
    zabs
  flg_lines = 0

end

;;;;
; Events
;;;;

pro x_pltobj_event, ev

  common xcommon_color
  common x_pltobj_images
  common x_specplot_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'PMIN': begin
          state.pltmin = ev.value
          x_pltobj_ReDisplay, state
      end
      'PMAX': begin
          state.pltmax = ev.value
          x_pltobj_ReDisplay, state
      end
      'INVERT': begin
          r_vector = reverse(r_vector)
          g_vector = reverse(g_vector)
          b_vector = reverse(b_vector)
          ximgd_stretchct, state
          x_pltobj_Display, state
      end
      'IDRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.itext_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'IDRAW' : begin
          widget_control, state.itext_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin         ; Button press
                  state.ipress = ev.press
                  case state.ipress of
                      1 : begin
                          ximgd_contrast, state ; CONTRAST/BRIGHTNESS
                          x_pltobj_Display, state
                      end
                      2 : 
                      4 : begin ; RMB = Center
                          x_pltobj_Plot, state
                          x_pltobj_pltreflin, state
                      end 
                  endcase
              end
              1 : state.zoom.flg = 0 ; Button release
              2 : begin         ; Motion event
                  state.tv.xcurs = ev.x
                  state.tv.ycurs = ev.y
                  if state.zoom.flg EQ 2 then begin
                      ximgd_contrast, state
                      x_pltobj_Display, state
                  endif else begin
                      state.xpos = x_tvx(state.tv, /intg)
                      state.ypos = x_tvy(state.tv, /intg)
                      x_pltobj_track, state
                  endelse
              end
          endcase
      end
      'ITEXT' : begin
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
              ; Single shift
              'b': begin
                  yval = long(xgety_plt(state.tv, /strct))
                  state.zoomreg[1] = yval
                  x_pltobj_ReDisplay, state
              end
              't': begin
                  yval = long(xgety_plt(state.tv, /strct))
                  state.zoomreg[3] = yval
                  x_pltobj_ReDisplay, state
              end
              'R': begin
                  state.zoomreg[1] = 0L
                  state.zoomreg[3] = img_size[1]-1
                  x_pltobj_ReDisplay, state
              end
              else: 
          endcase
      end
;;;;;;;;; SPEC ;;;;;;;;;;
      'XMAX' : begin
          widget_control, state.xmax_id, get_value=tmp
          state.svxymnx[3] = tmp
          state.xymnx[3] = tmp
          x_pltobj_Plot, state
      end
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          zabs = tmp
          flg_lines = 1
          x_pltobj_Plot, state
      end
      'SDRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.stext_id, sensitive = 0
          endif
      end
      'SDRAW' : begin
          widget_control, state.stext_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = Zoom
                          x_speczoomreg, state
                          if state.flg_zoom EQ 2 then begin
                              x_pltobj_update, state
                              state.flg_zoom = 0
                          endif
                      end 
                      4 : begin
                          x_specplot_SetLine, state ; Set reference line
                          x_pltobj_Plot, state
                      end
                  endcase
              end
              1 : ; Button Release
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
                  widget_control, state.swvval_id, set_value=state.xpos
;                  widget_control, state.ypos_id, set_value=state.ypos
              end
          endcase
      end
      'STEXT' : begin
          eventch = string(ev.ch)
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
          endif
          if (state.flg_EW EQ 1 AND eventch NE 'E') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
              print, 'Set the other side for EW!'
;              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
;              return
          endif
          tmpxy = state.xymnx
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'N': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              'z': begin  ; Region
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'L': x_specplot_SetLine, state   ; Set reference line
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'H': x_helpwidg, state.help
              ; Crude Analysis
              'E': x_specplot_EW, state        ; Calc EW
              'C': x_specplot_Colm, state      ; Calc AODM colm
              ' ': print, 'x: '+strtrim(state.xpos,2)+$
                '   y:'+strtrim(state.ypos,2)
              ; SMOOTH
              'S': begin
                  x_specplot_smooth, state
                  x_pltobj_Plot, state
              end
              'U': begin
                  x_specplot_smooth, state, /reset
                  x_pltobj_Plot, state
              end
              ; LINES
              '1': begin ; Plot CaH+K and OII (key on CaK)
                  x_pltobj_Plot, state
                  x_specplot_guess, state, 'CaHK', /IMG
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              '2': begin ; Plot Ha+N[II]
                  x_pltobj_Plot, state
                  x_specplot_guess, state, 'Ha', /IMG
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              ; Postscript
              'P': x_pltobj_specpsfile, state
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  ; Nothing
          endcase
          b = where(tmpxy NE state.xymnx,nb) 
          if nb NE 0 then x_pltobj_update, state
      end
      'LNLIST' : begin  ; LINE LIST
          state.flg_lines = ev.index + 1
          x_specplot_initLines, state
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
;  SPEC Plot
;;;;;;;;;;;;;;;;;;;;


pro x_pltobj_Plot, state
  
  common xcommon_color
  common x_specplot_lines

  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.sdraw_id, get_value=wind
      wset, wind
  endif

  ; Plot
  if state.flg_smooth EQ 0 then $
    plot, state.wave, state.fx, $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
    charsize=1.2, psym=10, background=7, color=0, $
    xtitle='!17Wavelength', xmargin=[0,0], ymargin=[0,0], xstyle=1,$
    ystyle=1 $
  else $
    plot, state.wave, state.smooth, $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
    charsize=1.2, psym=10, background=7, color=0, $
    xtitle='!17Wavelength', xmargin=[0,0], ymargin=[0,0], xstyle=1,$
    ystyle=1 

  ; Plot Lines as required
  if flg_lines EQ 1 then x_specplot_PltLines, state, /IMG

  ; SIGMA
  oplot, state.wave, state.sig, psym=10, color=1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Overplot wave value
;;;;;;;;;;;;;;;;;;;;


pro x_pltobj_pltreflin, state

  common x_pltobj_images
  
  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.sdraw_id, get_value=wind
      wset, wind
  endif

  xval = state.zoomreg[2] < $
    long(xgetx_plt(state.tv.xcurs, state.pos, state.tv.xymnx, state.size)) $
    > state.zoomreg[0]
  yval = long(xgety_plt(state.tv, /strct))
  wav = subwv[xval,yval]

  ; Plot
  oplot, [wav, wav], [-1e15,1e15], color=2

end

;;;;;;
; UPDATE PLOT+IMG

pro x_pltobj_update, state
  x_pltobj_Plot, state
  x_pltobj_ReDisplay, state
  return
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro x_pltobj_specpsfile, state

; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
  !p.charthick = 3

  device, decompose=0
  ps_open, file='spec.ps', font=1, /color
  state.psfile = 1
  x_pltobj_plot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0
  !p.thick = 1
  !p.charthick = 1

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; IMAGES ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;
;  InitImage
;;;;;;;;;;;;;;;;;;;;

pro x_pltobj_InitImg, state

common x_pltobj_images

  widget_control, /hourglass

  ; Size
  img_size = size(subfx, /dimensions)

  ; Set min max
  state.imgmin = min(subfx, max=fmax)
  state.imgmax = fmax

  ; pmin, pmax
  if state.pltmax EQ 0. then begin
;      med = median(subfx, /even)
;      sig = stddev(subfx)
;      state.pltmax = (med + (2 * sig)) < state.imgmax
;      state.pltmin = (med - (2 * sig))  > state.imgmin
      state.pltmax = 300.
      state.pltmin = -50.
      widget_control, state.pmin_id, set_value=state.pltmin
      widget_control, state.pmax_id, set_value=state.pltmax
  endif

;;;;;;;;;;;;;;;;;;
; Set Region

  img_left = min( abs(state.xymnx[0]-state.imgwv), iml)
  img_right = min( abs(state.xymnx[2]-state.imgwv), imr)

  state.sz_img[0] = imr-iml+1L
  state.sz_img[1] = img_size[1]

  ; ZOOMREG
  state.zoomreg = [iml, 0L, imr, state.sz_img[1]-1]
  state.tv.xymnx = float(state.zoomreg)
  x_setgridxy, state, state.zoomreg, /FILL

; Display Image
  display_image = bytscl(subfx[iml:imr,0:img_size[1]-1], $
                         min=state.pltmin, max=state.pltmax, /nan, $
                         top=state.ncolors-1) + 8B
  tv_image[*] = congrid(display_image, state.tv.winsize[0], state.tv.winsize[1])

end

;;;;;;;;;;;;;;;;;;;;
;  Display
;;;;;;;;;;;;;;;;;;;;

pro x_pltobj_Display, state
  
common x_pltobj_images

  widget_control, state.idraw_id, get_value=wind
  wset, wind

  ; TV
  tv, tv_image, round(state.size[0]*state.pos[0]), 0L

;  yrange=[state.tv.xymnx[1], state.tv.xymnx[3]]
  
end

;;;;;;;;;
;  ReDisplay
;;;;;;;;;

pro x_pltobj_ReDisplay, state

common x_pltobj_images

  widget_control, /hourglass

  ; Reset Zoom
  img_left = min( abs(state.xymnx[0]-state.imgwv), iml)
  img_right = min( abs(state.xymnx[2]-state.imgwv), imr)
  state.zoomreg[0] = iml
  state.zoomreg[2] = imr
  state.sz_img[0] = imr-iml+1L
  state.sz_img[1] = state.zoomreg[3]-state.zoomreg[1]+1
  state.tv.xymnx = float(state.zoomreg)

  ; Reset Display
  delvarx, display_image
  display_image = bytscl(subfx[state.zoomreg[0]:state.zoomreg[2], $
                               state.zoomreg[1]:state.zoomreg[3]], $
                         min=state.pltmin, max=state.pltmax, /nan, $
                        top=state.ncolors-1) + 8B
  tv_image = congrid(display_image, state.tv.winsize[0], state.tv.winsize[1])
  x_pltobj_Display, state
end

;;;;;;;;;
;  Track
;;;;;;;;;

pro x_pltobj_track, state

common x_pltobj_images

  ; Coordinates
  xval = state.zoomreg[2] < $
    long(xgetx_plt(state.tv.xcurs, state.pos, state.tv.xymnx, state.size)) $
    > state.zoomreg[0]
  yval = long(xgety_plt(state.tv, /strct))

  zcenter = [xval, yval] 

  track = bytarr(13,13)
  boxsize=6
  xmin = 0 > (zcenter[0] - boxsize - state.zoomreg[0])
  xmax = (zcenter[0] + boxsize - state.zoomreg[0]) < (state.sz_img[0] - 1) 
  ymin = 0 > (zcenter[1] - boxsize - state.zoomreg[1]) 
  ymax = (zcenter[1] + boxsize - state.zoomreg[1]) < (state.sz_img[1] - 1)

  startx = abs( (zcenter[0] - boxsize) < 0 )
  starty = abs( (zcenter[1] - boxsize) < 0 ) 

  track[startx,starty] = display_image[xmin:xmax,ymin:ymax]
  track_image = rebin(track, $
                      state.track_window_size, state.track_window_size, $
                      /sample)

  widget_control, state.tdraw_id, get_value=wind
  wset, wind
  tv, track_image

; Overplot an X on the central pixel in the track window, to show the
; current mouse position
  plots, [0.46, 0.54], [0.46, 0.54], /normal, color = 2, psym=0
  plots, [0.46, 0.54], [0.54, 0.46], /normal, color = 2, psym=0

; LABELS
  fx_string = string(subfx[xval,yval], format = '(g13.5)' )
  widget_control, state.fxval_id, set_value = fx_string
  wav_string = string(subwv[xval,yval], format = '(g13.6)' )
  widget_control, state.iwvval_id, set_value = wav_string

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; MAIN PROGRAM ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_pltobj, wave, fx, sig, infx, inwv, imgwv, $
              XSIZE=xsize, PMNX=pmnx, $
              YSIZE=ysize, OBJNM=objnm, ZIN=zin, XMAX=xmax, LLIST=llist

common x_pltobj_images
common xcommon_color
common x_specplot_lines

;
  if  N_params() LT 5  then begin 
    print,'Syntax - ' + $
      'x_pltobj, wave, fx, sig, subfx, subwv, [imgwv], XSIZE=, YSIZE= '
    print, '   PMNX=, OBJNM=, ZIN=, XMAX=, LLIST=  [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-300
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-400
  if not keyword_set( I_YSIZE ) then i_ysize = ysize/4
  if not keyword_set( S_YSIZE ) then s_ysize = 3*ysize/4

; Initialize the common blcok
  x_pltobj_icmmn
  subfx = infx
  subwv = inwv

; imgwv

  sz = size(subwv, /dimensions)
  if not keyword_set( IMGWV ) then imgwv = subwv[*, sz[1]/2]

; STATE

  ; tv structure
  tmp = {tvstruct}
  tmp.pos = [0.1, 0., 0.95, 1.]
  ; zoom structure
  tmp2 = {zoomstrct}

  state = {             $
            npix: n_elements(wave), $
            fx: fx, $
            wave: wave, $
            var: fltarr(n_elements(wave)), $
            sig: sig, $
            flg_smooth: 0, $   ; Smoothing
            smooth: fltarr(n_elements(wave)), $
            imgmin: 0., $                ; Image stuff
            imgmax: 0., $
            imgwv: imgwv, $
            sz_img: lonarr(2), $
            zoom: tmp2, $
            zoomreg: lonarr(4), $         
            tmpreg: lonarr(4), $         
            pltmax: 0.0, $
            pltmin: 0.0, $
            xpos: 0.0, $
            ypos: 0.0, $
            ipress: 0L, $
            track_window_size: 169L, $ ; size of tracking window
            tv: tmp, $              ; TV structure
            ncolors: 0, $
            brightness: 0.5, $
            contrast: 0.5, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_zoom: 0, $
            flg_EW: 0, $
            psfile: 0, $
            help: strarr(50), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            xcurs: 0., $
            ycurs: 0., $
            size: lonarr(2), $
            flg_lines: 0, $  ; QAL=1, GAL=2, QSO=3
            base_id: 0L, $      ; Widgets
            tdraw_id: 0L, $    ; Pan
            tdrawbase_id: 0L, $
            idraw_id: 0L, $    ; Image draw window
            idrawbase_id: 0L, $
            itext_id: 0L, $
            fxval_id: 0L, $
            iwvval_id: 0L, $
            sdraw_id: 0L, $       ; Spec Window
            sdrawbase_id: 0L, $
            stext_id: 0L, $
            swvval_id: 0L, $
            zabs_id: 0L, $
            xmax_id: 0L, $
            name_id: 0L, $
            nspec_id: 0L, $
            pmin_id: 0L, $
            pmax_id: 0L, $
            lines_id: 0L, $
            help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES
;  x_specplot_initcommon
  resolve_routine, 'x_specplot', /NO_RECOMPILE
  if not keyword_set( LLIST ) then state.flg_lines = 2 $
    else state.flg_lines = llist
  x_specplot_initLines, state

; ZABS
  if keyword_set( ZIN ) then begin
      zabs = zin
      flg_lines = 1
  endif else zin = 0.


; PMIN, PMAX
  if keyword_set(PMNX) then begin
      state.pltmin = pmnx[0]
      state.pltmax = pmnx[1]
  endif

;;;;
; XYMNX

  if not keyword_set( XYMNX ) then begin
      gdvar = where(state.sig GT 0., COMPLEMENT=badpix)
      dumwv = state.wave[gdvar]
      mdwv = where(dumwv[gdvar] GT 5000. AND $
                   dumwv[gdvar] LT 7000, nmdwv)
      if nmdwv EQ 0 then begin
          gdvar = lindgen(state.npix)
          dumwv = state.wave
          mdwv = where(dumwv GT 5000. AND $
                       dumwv LT 7000, nmdwv)
      endif

      mdfx = 2*median(state.fx[gdvar[mdwv]])
      
      state.xymnx = [min(dumwv), 0.0, max(dumwv), mdfx]
      state.svxymnx = state.xymnx
      ; SET bad pix
      if badpix[0] NE -1 then state.sig[badpix] = 1.e5
  endif

  ; XMAX keyword
  if keyword_set( XMAX ) then begin
      state.xymnx[3] = xmax
      state.svxymnx[3] = xmax
  endif

;    WIDGET
  base = WIDGET_BASE( title = 'x_pltobj: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_pltobj', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
  if not keyword_set(objnm) then objnm = ' '
  state.name_id = WIDGET_LABEL(labelbase, value='Obj: '+objnm, /align_center)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      IMAGE DRAW
  state.idrawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='IDRAW_BASE', frame=2)

  state.tv.winsize[0] = round(xsize*(state.pos[2]-state.pos[0]))
  state.tv.winsize[1] = i_ysize

  state.idraw_id = widget_draw(state.idrawbase_id, xsize=xsize, $
                              ysize=state.tv.winsize[1], /frame, $
                              /button_events, /motion_events, uvalue='IDRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; IMAGE TEXT
  state.itext_id = widget_text(state.idrawbase_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'ITEXT', $
                              value = '')
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Img max/min

  pltmnxbase = widget_base(toolbar, /column, /base_align_center, /align_center, $
                        /frame)
  state.pmin_id = cw_field(pltmnxbase, value=state.pltmin, /floating, $
                           title='imgmin', $
                          /return_events, xsize=8, UVALUE='PMIN') 
  state.pmax_id = cw_field(pltmnxbase, value=state.pltmax, /floating, $
                           title='imgmax', $
                          /return_events, xsize=8, UVALUE='PMAX') 

; XY position
  fw_id = widget_base(toolbar, /column, /align_center, frame=2)
  imglbl = WIDGET_LABEL(fw_id, value='Image', /align_center)
  state.fxval_id = cw_field(fw_id, title='fx:', value=0., xsize=13)
  state.iwvval_id = cw_field(fw_id, title='wv:', value=0., xsize=13)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      SPEC DRAW
  state.sdrawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='SDRAW_BASE', frame=2)

  state.size[0] = xsize
  state.size[1] = s_ysize

  state.sdraw_id = widget_draw(state.sdrawbase_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='SDRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; SPEC TEXT
  state.stext_id = widget_text(state.sdrawbase_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'STEXT', $
                              value = '')
; XY position
  spec_wv_id = widget_base(toolbar, /column, /align_center, frame=2)
;  state.fxval_id = cw_field(fw_id, title='fx:', value=0., xsize=13)
  imglbl = WIDGET_LABEL(spec_wv_id, value='Spec', /align_center)
  state.swvval_id = cw_field(spec_wv_id, title='wv:', value=0., xsize=13)
  state.zabs_id = cw_field(spec_wv_id, title='z:', value=zin, /floating, $
                           xsize=10, /return_events, uvalue='ZABS')
  state.xmax_id = cw_field(spec_wv_id, title='xmax:', value=state.xymnx[3], $
                           /floating, $
                           xsize=12, /return_events, uvalue='XMAX')


;      Lines
  state.lines_id = WIDGET_LIST(toolbar, $
                             VALUE=['QAL','GAL','QSO'], $
                             uvalue='LNLIST', ysize = 3)
  widget_control, state.lines_id, set_list_select=state.flg_lines-1
;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
;  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
;  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      PAN DRAW
  state.tdrawbase_id = $
    WIDGET_BASE( toolbar, /row, /base_align_center,/align_center, $
               uvalue='PDRAW_BASE', frame=2)

  state.tdraw_id = widget_draw(state.tdrawbase_id, xsize=169L, $
                              ysize=169L, frame=2, $
                              /button_events, /motion_events, uvalue='TDRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB - Truncate/Extend trace', $ 
             'RMB - Contrast/Brightness', $
             'CMB/CMB - Zoom' $ 
             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; COLORS
  device, get_decomposed=val_decomp
  device, decompose=0
  loadct, 0, /silent
  if (!d.table_size LT 12) then begin
      message, 'x_stsltgui: Too few colors available for color table'
      stop
  endif
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  ximgd_getct, state, 0, /CLR

; INVERT
;  r_vector = reverse(r_vector)
;  g_vector = reverse(g_vector)
;  b_vector = reverse(b_vector)
;  ximgd_stretchct, state

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; IMAGE
  tv_image = bytarr(state.tv.winsize[0], state.tv.winsize[1])
;  main_image = subfx
  x_pltobj_InitImg, state
;  x_pltobj_Reset, state
  x_pltobj_Display, state

  ; PLOT
  x_pltobj_Plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_pltobj', base

  delvarx, main_image, tv_image, display_image, subfx, subwv
  device, decomposed=val_decomp

  return
end
	
