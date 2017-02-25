;+ 
; NAME:
; outr_speclinetool
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_speclinetool_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
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
;   dunit=     - Extension in the FITS file
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
;   outr_speclinetool_ 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   30-Mar-2008 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;
pro outr_speclinetool_common
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  return
end


;;;;
; Events
;;;;

pro outr_speclinetool_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SLIDER1': begin
          widget_control, state.slider1.slider_id, get_value=tmpv
          case state.slider1.curslide of
              0: begin ;; Velocity
                  state.zabs = tmpv/3e5
                  ;; Save
                  state.slider1.value[state.slider1.curslide] = tmpv
              end
              1: begin ;; Redshift
                  state.zabs = tmpv
                  ;; Save
                  state.slider1.value[state.slider1.curslide] = tmpv
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
      'ABS_EMISS' : begin
          widget_control, state.abse_id, get_value=abse_indx
          state.flg_emiss = (where(abse_indx))[0]
      end
      'DRAW_BASE' : begin
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'DRAW' : begin
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = upper level
                          state.xcurs = ev.x
                          state.ycurs = ev.y
                          state.xclick = xgetx_plt(state, /strct)
                          state.yclick = xgety_plt(state, /strct)
                          r = sqrt(state.xclick^2 + state.yclick^2)
                          state.nupper = ((state.nlower+1) > $
                                         round(sqrt(r/state.rad))) < 6
                      end 
                      4 : begin     ; Right button = lower level
                          state.xcurs = ev.x
                          state.ycurs = ev.y
                          state.xclick = xgetx_plt(state, /strct)
                          state.yclick = xgety_plt(state, /strct)
                          r = sqrt(state.xclick^2 + state.yclick^2)
                          state.nlower = ((state.nupper-1) < $
                                         round(sqrt(r/state.rad))) > 1
                      end 
                      else: 
                  endcase
                  outr_speclinetool_newwave, state
              end
              1 : begin ; Button Release 
                  state.flg_click = 0
                  WIDGET_CONTROL, state.base_id, set_uvalue = state,  /no_copy
                  return
              end
              2 : begin ; Motion event
                  if state.flg_click EQ 0 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; Update Plot
  outr_speclinetool_UpdatePlot, state
  outr_speclinetool_colorbar, state

  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;
pro outr_speclinetool_UpdatePlot, state
  
;common outr_speclinetool_lines

; Plot Data
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  lthick = 3.
  csize = 2.2
  lsz = 2.0
  lsz_big = 2.5

  clr = getcolor(/load)

  ;; 
  plot, [0],  [0], $
        position=state.pos, $
        xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1], state.xymnx[3]], xstyle=4, ystyle=4, $
        thick=lthick, background=clr.white, $
        color=clr.black, charsize=csize, /nodata

  ;; Plot all orbits
  for nn=1L,6 do begin
      if nn EQ state.nlower or nn EQ state.nupper then clrl = clr.black $
      else clrl = clr.gray
      x_oplotcirc, state.rad*nn^2, color=clrl, linestyle=2, thick=3
      xyouts, 0., state.rad*(nn^2+0.1), 'n='+strtrim(nn,2), $
              color=clrl, charsize=csize*0.8, alignment=0.5
  endfor

  ;; Plot current
  x_oplotcirc, state.rad*state.nlower^2, color=clr.black, thick=3
  x_oplotcirc, state.rad*state.nupper^2, color=clr.black, thick=3

  ;; Proton
  plotsym, 0, 0.5, /fill
  oplot, [0], [0], psym=8, color=clr.blue

  ;; Label
  xyouts, 0., 0.98, 'Hydrogen Atom', color=clr.black, charsize=csize, $
          alignment=0.5

  xyouts, 0., -0.99, 'Click to change where the electron starts/ends', $
          color=clr.black, charsize=csize, alignment=0.5

  ;; Arrow
  if state.flg_emiss then begin
      n1 = state.nupper
      n2 = state.nlower
  endif else begin
      n2 = state.nupper
      n1 = state.nlower
  endelse
  clra = clr.brown
  arrow, state.rad*n1^2, 0., state.rad*n2^2, 0., $
         color=clra, thick=3, /solid, /data
  if state.flg_emiss then lbl = 'Emission' else lbl='Absorption'
  xyouts, mean(state.rad*[n1^2,n2^2]), 0.03, lbl, color=clra, charsize=csize,$
          alignment=0.5

  return

end

;;;;;;;;;;;;;;;;;;;;
;  New wave range
;;;;;;;;;;;;;;;;;;;;

pro outr_speclinetool_newwave, state

  ;; Calculate the wavelength (Rest)
  state.wavelength = state.const / abs(1./state.nlower^2 - 1./state.nupper^2) 

  ;; Reset wavelength region
  state.wvmin = (state.wavelength-300) < 3300.
  state.wvmax = (state.wavelength+500) > 7500.

  return
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;
pro outr_speclinetool_psfile, state

; Device
  device, get_decomposed=svdecomp

  x_psopen, 'idl.ps', /maxs
  state.psfile = 1
  outr_speclinetool_UpdatePlot, state
  x_psclose
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  Reset 
;;;;;;;;;;;;;;;;;;;;

pro outr_speclinetool_Reset, state, INWAV=inwav

  if keyword_set(INWAV) then state.wavelength = inwav 

  state.xymnx = state.svxymnx
  outr_speclinetool_newwave, state

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Color Bar
;;;;;;;;;;;;;;;;;;;;
pro outr_speclinetool_colorbar, state

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
  x0 = state.wvmin
  x1 = state.wvmax
  xv = x0 + findgen(xsize)*(x1-x0)/float(xsize-1) 

  ;; Interpolate colors
  b = interpol(cmap, wvm, xv) > cmap[0]
  cuttop = where(xv GT 7000., ncut)
  if ncut NE 0 then b[cuttop] = b[cuttop[0]]

  cutlow = where(xv LT 4200., ncut)
  if ncut NE 0 then b[cutlow] = b[cutlow[ncut-1]+20]

  ;; Emission or Absorption?
  mn = min(abs(xv-state.wavelength*(1+state.zabs)), imn)
  i0 = (imn-3) > 0
  i1 = (imn+3) < (xsize-2)
  if state.flg_emiss then begin
      b[0:i0] = cmap[0]
      b[i1:*] = cmap[0]
  endif else b[i0:i1] = cmap[0]
      

  c = replicate(1, ysize)
  a = b # c

  tv, a

  ;; Label
  clr = getcolor(/load)
  plot, [0.], [0.], position=[0., 0.15, 1., 0.95], $
        xrange=[state.wvmin, state.wvmax], $
        yrange=[0., 1.], xstyle=1, ystyle=5, /nodata, $
        /noerase, background=clr.black, color=clr.white, $
        xtitle='Wavelength (Angstroms)', charsize=1.5
  wavebands_x = [0.01, 10., 1500., 5600., 12000., 1e10]
  wavebands_y = [0.4, 0.1, 0.1, 0.1, 0.1, 0.1]
  wavebands_nm = ['Gamma-rays', 'X-rays', 'Ultraviolet', 'Optical', 'Infrared', 'Radio']
  nwave = n_elements(wavebands_x)

  csize = 2.
  for ss=0L,nwave-1 do $
         xyouts, wavebands_x[ss], wavebands_y[ss], $
                 wavebands_nm[ss], color=clr.white, charsize=csize, $
                 alignment=0.5

  xyouts, mean([state.wvmin, state.wvmax]), 0.8, 'Spectral Line at '+$
          strtrim(round(state.wavelength*(1+state.zabs)),2)+' Angstroms', $
          color=clr.white, charsize=csize, alignment=0.5
  

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

pro outr_speclinetool, XSIZE=xsize, YSIZE=ysize, TITLE=title, BLOCK=block

  common xcommon_color

  outr_speclinetool_common

  ;;  Optional Keywords
  device, get_screen_size=ssz

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

  ;; Get fonts
  fonts = x_widget_setfont(ssz[0])
;	print, fonts.small_font

  npix = 10000L
  c = x_constants()
  const = c.h*c.c/c.Ryd * 1e8 ;; Ang

;    STATE
  state = { wavelength: 5500., $ ;; Angstroms
            zabs: 0., $
            pos: [0.05,0.05,0.95,0.95], $ ; Plotting
            nlower: 2L, $
            nupper: 3L, $
            rad: 0.025, $
            const: const, $
            psfile: 0, $ ; Postscript
            svxymnx: [-1, -1., 1, 1.], $
            xymnx: dblarr(4), $
            tmpxy: dblarr(4), $
            wvmin: 3300., $
            wvmax: 7500., $
            flg_emiss: 1L, $
            psym: 10, $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            xclick: 0., $
            yclick: 0., $
            slider1: {wsliderstrct}, $
            wv_click: 0., $
            flg_click: 0, $
            ncolors: 0L, $      ; Image stuff
            brightness: 0.5, $
            contrast: 0.5, $
            fonts: fonts, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            clrbar_base_id: 0L, $
            clrbar_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            img_draw_id: 0L, $
            abse_id: 0L, $
            error_msg_id: 0L $
          }

  ;;    Title
;  if size(flux, /type) EQ 7 then state.title = flux
  if keyword_set( TITLE ) then state.title = title

;  device, set_FONT='Times', /TT_FONT

;    WIDGET
  base = WIDGET_BASE( title = 'outr_speclinetool -- Version 1.0', /row, $
                    xoffset=30, yoffset=30)
  state.base_id = base

  ;;;;;;;;;;;;
  ;;   Orbital
  state.size[0] = 2*ysize/3
  state.size[1] = 2*ysize/3
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=1, $
                                   /tracking_events)
  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[1], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')

  ;;;;;;; Right
  right_base = widget_base(base, /column, /base_align_left, $
                      uvalue='DRAW_BASE', frame=1)
  
  ;;;;;;;;;;;
  ;; Text
  fil = getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Text/spectral_line.txt'
  nlines = MIN([FILE_LINES(fil), 10000])
  OPENR, unit, FIL, /GET_LUN
  a = strarr(nlines)
  readf, unit, a
  CATCH, /CANCEL
  FREE_LUN, unit
  text_id = WIDGET_TEXT(right_base, scr_xsize=xsize-state.size[1],$
		UVALUE='TEXT', /scroll, VALUE = a, ysize=12, $
		FONT = state.fonts.small_font, /WRAP)

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color bar
  state.clrbar_base_id = widget_base(right_base, /column, /base_align_left, $
                                     uvalue='CLRBAR_BASE', frame=1 )
  state.clrbar_id = widget_draw(state.clrbar_base_id, $
                                xsize=xsize-state.size[1],$
                                ysize=round(state.size[1]/2), retain=2, $
                                uvalue='CLRBAR')


  ;; Twiddle base
  fiddle_base = widget_base(right_base, /row, /base_align_left, $
                      uvalue='ABSE_BASE')
  ;; Emission/Absorption
  state.abse_id = cw_bgroup(fiddle_base, ['Emission', 'Absorption'], $
                               /exclusive, row=2, $
                               UVALUE='ABS_EMISS',frame=1, $
                               LABEL_TOP='Process',$
                               FONT=state.fonts.big_font)
  widget_control, state.abse_id, set_value=0


  ;;;;;;;;;;
  ;; Slider
  state.slider1.base_id = WIDGET_BASE( fiddle_base, /column, $
                                       /frame, /base_align_left,$
                                       /align_left)
  ;; Set the values
  state.slider1.max = 1.
  state.slider1.xsize = 100.
  state.slider1.uname = 'SLIDER1'
  values = ['Speed of Atom (km/s)', 'Redshift']
  state.slider1.nslide = 2
  state.slider1.value[0:state.slider1.nslide-1] = [0., 0.]
  state.slider1.droplist[0:state.slider1.nslide-1]=values
  state.slider1.min[0:state.slider1.nslide-1] = [-30000., 0.]
  state.slider1.max[0:state.slider1.nslide-1] = [30000., 1.]
  state.slider1.fonts = state.fonts

  ;; Make the widget
  tmp = state.slider1
  x_widget_slider, tmp, /init, /DRAG
  state.slider1 = tmp
  widget_control, state.slider1.slider_id, set_value=0.
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Done
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  done = WIDGET_BUTTON(right_base, value='Done',uvalue='DONE', /align_right, $
                      FONT=state.fonts.big_font)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Realize
  WIDGET_CONTROL, base, /realize

  !p.font = 1
  device, set_FONT='Times', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Update spectrum
  outr_speclinetool_Reset, state
  outr_speclinetool_UpdatePlot, state

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color Bar
  loadct, 0, /silent
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)
  ximgd_getct, state, 13, /CLR
  
  outr_speclinetool_colorbar, state

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'outr_speclinetool', base, /no_block $
  else xmanager, 'outr_speclinetool', base

return
end
