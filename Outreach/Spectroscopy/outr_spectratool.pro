;+ 
; NAME:
; outr_spectratool 
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_spectratool_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           /QAL, /GAL, ZIN=, /BLOCK, /NRM, /LLS
;
; INPUTS:
;   flux  - Flux array (or FITS file)
;   [ysin]  - Sigma array (or FITS file)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave=      - wavelength array
;   INFLG=     - Specifies the type of input files (or arrays)
;   /LLS       - Use Lyman limit line list
;   /QSO       - Use Quasar line list
;   /GAL       - Use galaxy line list
;   /QAL       - Use quasar absorption line list
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   outr_spectratool_ 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;   19-Dec-2001 Added Error array, Colm
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;
pro outr_spectratool_common
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  return
end


;;;;
; Events
;;;;

pro outr_spectratool_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'OVERLAY_GROUP': 
      'CLRBAR' : 
      'LINE_LIST': begin
          tmp = state.lines
          outr_spectra_lines_modify, ev.value, tmp
          state.lines = tmp
      end
      'SLIDER1': begin
          widget_control, state.slider1.slider_id, get_value=tmpv
          case state.slider1.curslide of
              0: begin ;; Redshift
                  state.zabs = tmpv
                  ;; Save
                  state.slider1.value[state.slider1.curslide] = tmpv
              end
              1: begin ;; Temperature
                  ;; Overlayed Star
                  mt = where(strmatch(state.overlays.name, 'Star'))
                  mt = mt[0]
                  state.overlays[mt].param[0] = tmpv
                  tmp = state.overlays[mt]
                  outr_overlays, tmp, /fiddle, wvmnx=[state.xymnx[0]/10., $
                                                      state.xymnx[2]*10.]
                  state.overlays[mt] = tmp
              end
              else: stop
          endcase
      end
      'SLIDER1_MIN': begin
          widget_control, state.slider1.min_id, get_value=tmpv
          state.slider1.min[state.slider1.curslide] = tmpv
          x_widget_slider, state.slider1, /reinit
      end
      'SLIDER1_MAX': begin
          widget_control, state.slider1.max_id, get_value=tmpv
          state.slider1.max[state.slider1.curslide] = tmpv
          x_widget_slider, state.slider1, /reinit
      end
      'SLIDER1_DROP': begin
          state.slider1.curslide = ev.index
          x_widget_slider, state.slider1, /reinit
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
                      4 : if state.flg_lines NE 0 then $
                            outr_spectratool_SetLine, state ; Set reference line
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
              ; Crude Analysis
              'E': outr_spectratool_EW, state        ; Calc EW
              'N': outr_spectratool_Colm, state      ; Calc AODM colm
              'G': outr_spectratool_Gauss, state      ; Fit a Gaussian
              ' ': print, 'x: '+strtrim(state.xpos,2)+$
                '   y:'+strtrim(state.ypos,2)
              ;; Postscript
              'P': outr_spectratool_psfile, state  
              ;; Overplot a Quasar spectrum
              'Q': outr_spectratool_qsotempl, state  
              ;; Overplot a galaxy spectrum
              'g': outr_spectratool_galtempl, state  
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'outr_spectratool_ Not a valid key!' ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; Bound the wavelengths (Angstroms)
  state.xymnx[0] = state.xymnx[0] > 1e-7
  state.xymnx[2] = state.xymnx[2] < 1e20

; Update Plot
  outr_spectratool_UpdatePlot, state
  outr_spectratool_colorbar, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;
pro outr_spectratool_UpdatePlot, state
  
;common outr_spectratool_lines

; Plot Data
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  lthick = 2.

  ;; Change to Log?
  if (state.xymnx[2]-state.xymnx[0]) GT 1e4 or  $
    (state.xymnx[2]/state.xymnx[0]) GT 1e3 then begin
      if state.flg_logx EQ 0 then state.xymnx[0] = state.xymnx[2]/1e3
      state.flg_logx = 1 
  endif else state.flg_logx = 0

  clr = getcolor(/load)
  csize = 2.1

  plot, state.spectrum.wave*(1+state.zabs), state.spectrum.fx, psym=state.psym, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        xtitle='!17Wavelength'+state.wave_lbl, ytitle='Flux', thick=lthick, $
        title=state.title, background=clr.white, $
        color=clr.black, charsize=csize, $
        xlog=state.flg_logx
  
  ;; Input spectrum
  xyouts, 0.7, 0.9, 'Input: '+state.spectrum.name, $
          /normal, alignment=0., color=clr.black, charsize=csize
          

  ;;;;;;;;;;;;;;;;;;;
  ;; Overlays
  widget_control, state.overlay_id, get_value=over_indx
  gd = where(over_indx, nover)
  if nover GT 0 then begin
      clrs = (x_setclrs(nover+1))[1:*]
      for jj=0L,nover-1 do begin
          oplot, state.overlays[gd[jj]].wave, $
                 state.overlays[gd[jj]].fx*state.xymnx[3]*0.95, $
                 color=clrs[jj], thick=lthick
          xyouts, 0.7, 0.9-0.06*(jj+1), state.overlays[gd[jj]].name, $
                  /normal, alignment=0., color=clrs[jj], charsize=csize
      endfor
  endif

  ;;;;;;;;;;;;;;;;;;;
  ;; Spectral lines
  gd = where(state.lines.j GT 0, ngd)
  for ss=0L,ngd-1 do begin
      ;; Emission
      mn = min(abs(state.lines[gd[ss]].wrest-$
                   state.spectrum.wave*(1+state.zabs)), imn)

      if state.lines[gd[ss]].A LT 0 then begin
          plotsym, 1, 1.5, thick=3 
          off = (state.xymnx[3]-state.xymnx[1])*0.1
          ali = 0.
          lclr = clr.blue
      endif else begin
          plotsym, 2, 1.5, thick=3 
          off = -1*(state.xymnx[3]-state.xymnx[1])*0.1
          ali = 1.
          lclr = clr.red
      endelse 

      ;; Plot
      oplot, [state.lines[gd[ss]].wrest], [state.spectrum.fx[imn]+off], $
             psym=8, color=lclr
      xyouts, state.lines[gd[ss]].wrest, state.spectrum.fx[imn]+off, $
              state.lines[gd[ss]].ion, charsize=csize/2., color=lclr, $
              orientation=90., alignment=ali
  endfor
      
      
  return

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro outr_spectratool_Reset, state


; Plotting
  state.xymnx = state.svxymnx

; Sort wave
  srt = sort(state.spectrum.wave)
  state.spectrum.wave = state.spectrum.wave[srt]
  state.spectrum.fx = state.spectrum.fx[srt]
  
  
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;
pro outr_spectratool_psfile, state

; Device
  device, get_decomposed=svdecomp

  x_psopen, 'idl.ps', /maxs
  state.psfile = 1
  outr_spectratool_UpdatePlot, state
  x_psclose
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  Color Bar
;;;;;;;;;;;;;;;;;;;;
pro outr_spectratool_colorbar, state

  widget_control, state.clrbar_id, get_value=wind
  wset, wind

  state.contrast=0.5
  state.brightness=0.5
  ximgd_getct, state, 13, /CLR
  xsize = (widget_info(state.clrbar_id, /geometry)).xsize
  ysize = (widget_info(state.clrbar_id, /geometry)).ysize

  ;; Map for 4000A to 7000A
  nwav = 3000L
  cmap = congrid( findgen(state.ncolors), nwav) + 8
  wvm = 4000. + findgen(nwav)
  

  ;; Get endpoints (full screen)
  dum = state.xcurs
  state.xcurs = 0.
  x0 = xgetx_plt(state, /strct) ;; Ang
  state.xcurs = xsize
  x1 = xgetx_plt(state, /strct) ;; Ang
  state.xcurs = dum
  if state.flg_logx EQ 0 then $
    xv = x0 + findgen(xsize)*(x1-x0)/float(xsize-1) $
  else xv = 10^(alog10(x0) + $
                findgen(xsize)*(alog10(x1)-alog10(x0))/float(xsize-1))

  ;; Interpolate colors
  b = interpol(cmap, wvm, xv) > cmap[0]
  cuttop = where(xv GT 8000., ncut)
  if ncut NE 0 then b[cuttop] = cmap[0]
  c = replicate(1, ysize)
  a = b # c

  tv, a

  ;; Label
  plot, [0.], [0.], position=state.pos, $
        xrange=[state.xymnx[0], state.xymnx[2]], $
        yrange=[0., 1.], xstyle=5, ystyle=5, /nodata, $
        /noerase, xlog=state.flg_logx
  clr = getcolor(/load)
  wavebands_x = [0.01, 10., 1000., 5600., 50000., 1e10]
  wavebands_y = [0.4, 0.6, 0.1, 0.6, 0.1, 0.4]
  wavebands_nm = ['Gamma-rays', 'X-rays', 'Ultraviolet', 'Optical', 'Infrared', 'Radio']
  nwave = n_elements(wavebands_x)

  for ss=0L,nwave-1 do $
         xyouts, wavebands_x[ss], wavebands_y[ss], $
                 wavebands_nm[ss], color=clr.white, charsize=2., $
                 alignment=0.5
  

  widget_control, state.draw_id, get_value=wind
  wset, wind


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

pro outr_spectratool, spectrum, XSIZE=xsize, YSIZE=ysize, TITLE=title $
                , BLOCK=block, IMGFIL=imgfil, DEFAULT=default, SPIRAL=spiral

  common xcommon_color

;
  ;; Default
  if keyword_set(DEFAULT) then begin
;      wv= 4000 + findgen(3000L)
      wv= 100 * 10.^(3*findgen(7000L)/7000)
      spectrum = {$
                 name: 'Sun (BB)', $
                 wave: wv, $
                 fx: blackbody(wv, 5500.) $
                 }
      imgfil = getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Images/sun.jpg'
  endif

  if keyword_set(SPIRAL) then begin
      dat = xmrdfits(getenv('XIDL_DIR')+ $
                     'Outreach/Spectroscopy/Spectra/sdss_spiral.fit',0,head)
      npix = n_elements(dat[*,0])
      wv= 10.d^(sxpar(head,'COEFF0') + dindgen(npix)*sxpar(head,'COEFF1'))
      spectrum = {$
                 name: 'Spiral galaxy', $
                 wave: wv/(1+sxpar(head,'Z')), $
                 fx: smooth(dat[*,0],3) $
                 }
      imgfil = getenv('XIDL_DIR')+'Outreach/Spectroscopy/Images/sdss_spiral.jpg'
  endif

  if not keyword_set(spectrum) then begin
      print,'Syntax - ' + $
            'outr_spectratool, spectrum, WAVE=, DUNIT='
      print, '            /BLOCK, /DEFAULT) [v1.0]'
      return
  endif 

                 
  outr_spectratool_common

  ;;  Optional Keywords
  device, get_screen_size=ssz

  ;;  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-100
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100

  fonts = x_widget_setfont(ssz[0])

  tmp1 = { newabslinstrct }

;  outr_spectratool_initcommon

  ;; Plot range (initial)
  xmin = min(spectrum.wave, max=xmax)
  ymin = min(spectrum.fx, max=ymax)

  ;; Spectral lines
  resolve_routine, 'outr_spectra_lines', /compile_full_file
  outr_spectra_lines, pd, lines

;    STATE
  state = { spectrum: spectrum, $
            zabs: 0., $
            flg_zoom: 0, $
            flg_logx: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_GS: 0, $ ; Gaussian stuff
            GS_lmt: dblarr(2,2), $
            GS_fit: fltarr(1000L), $
            GS_xpix: lonarr(2), $
            Gauss: fltarr(4), $    
            xpos: 0.d, $
            ypos: 0.d, $
            psfile: 0, $ ; Postscript
            svxymnx: double([xmin, ymin, xmax, ymax+0.10*(ymax-ymin)]), $
            xymnx: dblarr(4), $
            wave_lbl: ' (Angstroms)', $
            lines: lines, $  ;; Lines
            pd_lines: pd, $
            fonts: fonts, $
            tmpxy: dblarr(4), $
            xlog: 0, $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            ncolors: 0L, $      ; Image stuff
            brightness: 0.5, $
            contrast: 0.5, $
            overlays: replicate({overlaystrct},10), $
            slider1: {wsliderstrct}, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            clrbar_base_id: 0L, $
            clrbar_id: 0L, $
            draw_base_id: 0L, $
            overlay_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            lines_id: 0L, $
            img_draw_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }

  state.size[0] = xsize
  state.size[1] = ysize

  
  ;;    Title
;  if size(flux, /type) EQ 7 then state.title = flux
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'outr_spectratool -- Version 1.0', /column, $
                    xoffset=50, yoffset=50)
  state.base_id = base
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color bar
  state.clrbar_base_id = widget_base(base, /column, /base_align_left, $
                                     uvalue='CLRBAR_BASE', frame=1 )
  state.clrbar_id = widget_draw(state.clrbar_base_id, xsize=state.size[0], $
                              ysize=round(state.size[1]/18), retain=2, $
                              uvalue='CLRBAR')
  
  ;;;;;;;;;;;;
  ;;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=1, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=round(2.*state.size[1]/3), /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
  state.size[0] = xsize
  state.size[1] = 2*ysize/3.

  ;;;;;;;;;
  ;;  Stuff
          
  imgsz = round(xsize/6.) < round(ysize/6.)
  stuff = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center, ysize=imgsz)

; zabs

  explorebar = WIDGET_BASE( stuff, /row, /frame, /base_align_center,$
                           /align_left, xsize=xsize-imgsz)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; OVERLAYS
  overlays = WIDGET_BASE( explorebar, /row, /frame, /base_align_left,$
                           /align_center)

  state.overlays[0:3].name = ['Star', 'Human', 'Na Lamp', 'QSO']
  gd = where(strlen(state.overlays.name) GT 0, ngd)
  ;; Initialize
  tmp = state.overlays
  outr_overlays, tmp, /INIT
  state.overlays = tmp
  ;; Create button group
  state.overlay_id = cw_bgroup(overlays, state.overlays[gd].name, $
                               /nonexclusive, row=6, $
                               UVALUE='OVERLAY_GROUP',frame=1, $
                               LABEL_TOP='Overlays (Up to 3)',$
                               FONT=state.fonts.big_font)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; TOOLS
  tools = WIDGET_BASE( explorebar, /row, /frame, /base_align_left,$
                       /align_center)
  ;;;;;;;;;;
  ;; Slider
  state.slider1.base_id = WIDGET_BASE( tools, /column, /frame, /base_align_left,$
                                       /align_center)
  ;; Set the values
  state.slider1.max = 1.
  state.slider1.xsize = 50.
  state.slider1.uname = 'SLIDER1'
  values = ['Redshift', 'Star Temperature']
  state.slider1.nslide = 2
  state.slider1.value[0:state.slider1.nslide-1] = [0., 5500.]
  state.slider1.droplist[0:state.slider1.nslide-1]=values
  state.slider1.min[0:state.slider1.nslide-1] = [0., 1000.]
  state.slider1.max[0:state.slider1.nslide-1] = [10., 100000.]
  state.slider1.fonts = state.fonts

  ;; Make the widget
  tmp = state.slider1
  x_widget_slider, tmp, /init, /DRAG
  state.slider1 = tmp
  

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; xy position

  xy_id = widget_base(explorebar, /column, /align_center, frame=2, $
                      ysize=round(ysize/6.))
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13, $
                           font=state.fonts.small_font)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13, $
                           font=state.fonts.small_font)

  
  ;; Objname
  if keyword_set( head ) then begin
      objnm_id = cw_field(explorebar, title='Object: ', value=spectrum.name, $
                          /column, $
                          xsize=strlen(spectrum.name))
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Spectral lines
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  list_id = cw_pdmenu(explorebar, state.pd_lines, $
                      font=state.fonts.big_font, $
                      /help, $
                      /return_full_name, $
                      uvalue = 'LINE_LIST')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Done
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  done = WIDGET_BUTTON(explorebar, value='Done',uvalue='DONE', /align_right, $
                      font=state.fonts.big_font)

  ;; IMAGE 
  imagewin = WIDGET_BASE( stuff, /row, /frame, /base_align_center,$
                          /align_left, xsize=imgsz)
  state.img_draw_id = widget_draw(imagewin, xsize=imgsz, $
                              ysize=imgsz, retain=2, $
                              uvalue='IMGDRAW')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Realize
  WIDGET_CONTROL, base, /realize

  !p.font = 1
  device, set_FONT='Times', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Update spectrum
  outr_spectratool_Reset, state
  outr_spectratool_UpdatePlot, state

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color Bar
  loadct, 0, /silent
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  ximgd_getct, state, 13, /CLR
  outr_spectratool_colorbar, state

  if keyword_set(IMGFIL) then begin
      ;; Check extension
      ipos = strpos(imgfil, '.', /reverse_sear)
      case strmid(imgfil,ipos+1) of
          'jpg': read_jpeg, imgfil, img;, ctable, COLORS=state.ncolors
          else: stop
      endcase
      widget_control, state.img_draw_id, get_value=wind
      wset, wind 
      xsz = round((widget_info(state.img_draw_id, /geometry)).xsize)

      nw_img = congrid(img, 3, xsz, xsz)
      tv, nw_img, /true
;      tv, img, /true
;      tvlct, ctable
  endif

  

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'outr_spectratool', base, /no_block $
  else xmanager, 'outr_spectratool', base

return
end
