;+ 
; NAME:
; outr_basicspectra 
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_basicspectra_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
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
;   outr_basicspectra_ 'spec.fits'
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
pro outr_basicspectra_common
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  return
end


;;;;
; Events
;;;;

pro outr_basicspectra_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'MIRROR_POS': begin
          mirror = widget_info( state.mirrorpos_id, /droplist_select)
          prism = widget_info( state.prismpos_id, /droplist_select)
          if prism EQ mirror then widget_control, state.prismpos_id, $
            set_droplist_select=0
          state.mirror_angle = 0
          widget_control, state.slider_id, set_value=0.
      end
      'PRISM_POS': begin
          mirror = widget_info( state.mirrorpos_id, /droplist_select)
          prism = widget_info( state.prismpos_id, /droplist_select)
          if prism EQ mirror then widget_control, state.mirrorpos_id, $
            set_droplist_select=0
      end
      'EXERCISES': 
      'OVERLAYS': 
      'MIRROR': begin
          state.fiddle = 1
          state.mirror_angle = ev.value
          if abs(ev.value-45) LT 3 then $
            widget_control, state.mode_id, sensitive=1
      end
      'SOURCE': begin
          case state.sources[ev.index] of
              'Sun': begin 
                  tmp = state.spectrum
                  tmp.name = 'Sun'
                  outr_overlays, tmp, /init, wvmnx=state.wvmnx
                  state.spectrum = tmp
                  state.imgfil = $
                    getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/sun.jpg'
              end
              'Spiral Galaxy': begin 
                  dat = xmrdfits(getenv('OUTR_DIR')+ $
                                 'Outreach/Spectroscopy/Spectra/sdss_spiral.fit',0,head,/silent)
                  npix = n_elements(dat[*,0])
                  state.spectrum.npix = npix
                  wv= 10.d^(sxpar(head,'COEFF0') + $
                            dindgen(npix)*sxpar(head,'COEFF1'))
                  state.spectrum.name= 'Spiral Galaxy'
                  state.spectrum.wave[0:npix-1]= wv/(1+sxpar(head,'Z'))
                  state.spectrum.wave[npix:*]= 0.
                  state.spectrum.fx[0:npix-1]= smooth(dat[*,0],3) 
                  state.imgfil = $
                    getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/sdss_spiral.jpg'
              end
              'Sodium Lamp': begin 
                  tmp = state.spectrum
                  tmp.name = 'Na Lamp'
                  outr_overlays, tmp, /init, WVMNX=state.wvmnx
                  state.spectrum = tmp
                  state.spectrum.name = 'Sodium Lamp'
                  state.imgfil = $
                    getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/sodium_lamp.jpg'
              end
              'Quasar': begin 
                  tmp = state.spectrum
                  tmp.name = 'QSO'
                  outr_overlays, tmp, /init, WVMNX=state.wvmnx
                  state.spectrum = tmp
                  state.spectrum.name = 'Quasar (distant)'
                  state.spectrum.wave = state.spectrum.wave * 2 
                  state.imgfil = $
                    getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/quasar.jpg'
              end
              else: stop
          endcase
          state.spectrum.fx = state.spectrum.fx / max(state.spectrum.fx)
          outr_basicspectra_PltIMG, state
      end
      'MODE': begin
          state.mode = ev.index
          outr_basicspectra_setmode, state
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; Bound the wavelengths (Angstroms)
  state.xymnx[0] = state.xymnx[0] > 1e-7
  state.xymnx[2] = state.xymnx[2] < 1e20

  ;; Update Plot
  outr_basicspectra_UpdatePlot, state

  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end

;;;;;;;;;;;;;;;;;;;;;
;; Set mode
pro outr_basicspectra_setmode, state
  case state.mode of 
      0: begin  ;; Experiment
          widget_control, state.overlay_id, set_value=[0,0]
          widget_control, state.overlay_id, sensitive=0
          widget_control, state.spec_id, sensitive=0
          widget_control, state.prismpos_id, sensitive=1
          widget_control, state.mirrorpos_id, sensitive=1
      end
      1: begin  ;; Spectrum
          widget_control, state.overlay_id, sensitive=1
          widget_control, state.overlay_id, set_value=[1,0]
          widget_control, state.spec_id, sensitive=1
          widget_control, state.slider_id, sensitive=0
          widget_control, state.prismpos_id, sensitive=0
          widget_control, state.mirrorpos_id, sensitive=0
      end
      else:
  endcase
  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;
pro outr_basicspectra_UpdatePlot, state
  
;common outr_basicspectra_lines

  ;; Experiment
  if state.mode EQ 0 then begin
      outr_basicspectra_PltExperiment, state
      return
  endif

  !p.font = 1
  device, set_FONT='Times Bold', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  widget_control, state.draw_id, get_value=wind
  wset, wind

  widget_control, state.overlay_id, get_value=overlays

  lthick = 2.
  clr = getcolor(/load)
  csize = 2.1


  plot, [0], [0], psym=state.psym, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[0., 1.], xstyle=9, ystyle=9, $
        xtitle='Wavelength (Angstroms)', $
        ytitle='Brightness', thick=lthick, $
        background=clr.white, $
        color=clr.black, charsize=csize, $
        xlog=state.flg_logx

  ;; Histogram?
  if overlays[0] then begin
      rainbow = x_getrainbow()
      nclr = n_elements(rainbow)
      ;; Get y-values
      dwv = abs(state.spectrum.wave - shift(state.spectrum.wave,1))
      dwv[0]=dwv[1]
      yval = fltarr(nclr)
;      if strmatch(state.spectrum.name, 'Spiral Galaxy') then stop
      for ss=0L,nclr-1 do begin
          pix = where(state.spectrum.wave GT rainbow[ss].wvmnx[0] and $
                      state.spectrum.wave LT rainbow[ss].wvmnx[1], gdp)
          if gdp NE 0 then yval[ss] = total(state.spectrum.fx[pix]*dwv[pix])
      endfor
      mx = max(yval)
      yval = yval / mx
      ;; Plot
      for ss=0L,nclr-1 do begin
          if yval[ss] LT 1e-5 then continue
          x_curvefill, rainbow[ss].wvmnx, [0, 0], [yval[ss], yval[ss]], $
                       color=rainbow[ss].color, outthick=lthick
      endfor
  endif

  if overlays[1] then begin
      pix = lindgen(state.spectrum.npix)
      oplot, state.spectrum.wave[pix], state.spectrum.fx[pix], color=clr.darkgray, thick=2
  endif
  
  ;; Input spectrum
  xyouts, 3150., 1.0, 'Source: '+state.spectrum.name, $
          alignment=0., color=clr.black, charsize=csize
          
  xyouts, 0.05, 0.9, 'More Light', $
          /normal, alignment=0., color=clr.black, charsize=csize*0.8
  xyouts, 0.05, 0.2, 'Less Light', $
          /normal, alignment=0., color=clr.black, charsize=csize*0.8
      
      
  return

end

;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;

pro outr_basicspectra_PltExperiment, state

  common xcommon_color

  atoler = 2.

  ;; Grab the window
  widget_control, state.draw_id, get_value=wind
  wset, wind

  ;; Clear
  clr = getcolor(/load)
  plot, [0], [0], psym=state.psym, $
        position=[0., 0., 1., 1.], $
        xrange=[0., 1.], $
        yrange=[0., 1.], xstyle=4, ystyle=4, $
        thick=lthick, background=clr.black, $
        charsize=csize

  totx = round((widget_info(state.draw_id, /geometry)).xsize)
  xsz = round((widget_info(state.draw_id, /geometry)).xsize)/10.
  ysz = round((widget_info(state.draw_id, /geometry)).ysize)

  ;; Mirror
  mirror = widget_info( state.mirrorpos_id, /droplist_select)
  if mirror GT 0 then begin
      if mirror EQ 1 then xpos = [0.05, 0.25] else xpos = [0.4,0.6]
      clr = getcolor(/load)
      frac = totx/ysz
      yplt = 0.2*totx/ysz
      plot, [0], [0], psym=state.psym, $
            position=[xpos[0], 0.5-(0.1*totx/ysz), xpos[1], $
                      0.5+(0.1*totx/ysz)], $ ; Scale y,x
        xrange=[0., 1.], $
        yrange=[0., 1.], xstyle=4, ystyle=4, $
        thick=lthick, title=state.title, background=clr.black, $
        charsize=csize, /noerase
      ;; Black out
      x_curvefill, [0., 1], [0., 0], [1,1], color=clr.black
      half_len = (totx/15.) / (0.2*totx)
      x1 = 0.5 - half_len*cos(state.mirror_angle*!pi/180)
      y1 = 0.5 + half_len*sin(state.mirror_angle*!pi/180)
      x2 = 0.5 + half_len*cos(state.mirror_angle*!pi/180)
      y2 = 0.5 - half_len*sin(state.mirror_angle*!pi/180)
      oplot, [x1,x2], [y1,y2], color=clr.gray, thick=6
      xyouts, mean(xpos), 0.1, 'MIRROR', $
              color=clr.gray, charsiz=2.5, alignm=0.5, /normal
  endif

  ;; Prism  0.3 to 0.6
  prism = widget_info( state.prismpos_id, /droplist_select)
  if prism GT 0 then begin
      if prism EQ 2 then begin
          xpos = [0.3, 0.6] 
          ypos = [0., 1.]
      endif else begin
          xpos = [0.05, 0.25]
          ypos = [0., 0.8]
      endelse
      width = totx * (xpos[1]-xpos[0])
      nrm_side = 0.8
      side = nrm_side * width
      nrm_height = nrm_side / 2. * tan(60*!pi/180.) * (width/ysz)
      height = nrm_height * totx* (xpos[1]-xpos[0])
      plot, [0], [0], psym=state.psym, $
            position=[xpos[0], ypos[0], xpos[1], ypos[1]], xrange=[0., 1.], $
            yrange=[0., 1.], xstyle=5, ystyle=5, $
            thick=lthick, title=state.title, background=clr.black, $
            charsize=csize, /noerase
      x_curvefill, [0., 1], [0., 0], [1,1], color=clr.black
      ybase = (1.-height/ysz)/2
      x_curvefill, [0.1, 0.5], replicate(ybase,2), $
                   [ybase, ybase + nrm_height], color=clr.gray
      x_curvefill, [0.5, 0.9], replicate(ybase,2), $
             [ybase + nrm_height, ybase], color=clr.gray
      xyouts, mean(xpos), 0.1, 'PRISM', $
              color=clr.gray, charsiz=2.5, alignm=0.5, /normal
  endif


  ;; Sun [x = 0 to 0.25]
  nw_img = congrid(state.sunimg, 3, xsz, xsz)
  tv, nw_img, 0.10*totx, ysz-xsz, /true
      
  clr = getcolor(/load)
  xyouts, 0.05, 0.9, 'SUN', color=clr.gray, $
          charsiz=2.5, alignm=0.5, /normal


  ;; Rainbow arrows
  if MIRROR EQ 1 and PRISM EQ 2 then begin
      x0 = 0.7+0.01
      y0 = ybase + nrm_height/2.
      x1 = 1.3
      if abs(state.mirror_angle-45) LT atoler then begin
          arrow, x0, y0, x1, 0.7, color=clr.red, /data, /solid, thick=5, hthick=2.
          arrow, x0, y0, x1, y0, color=clr.green, /data, /solid, thick=5, hthick=2.
          arrow, x0, y0, x1, 0.4, color=clr.blue, /data, /solid, thick=5, hthick=2.
      endif 
  endif else begin
      arrow, x0, y0, x1, 0.7, color=clr.black, /data, /solid, thick=5, hthick=2.
      arrow, x0, y0, x1, y0, color=clr.black, /data, /solid, thick=5, hthick=2.
      arrow, x0, y0, x1, 0.4, color=clr.black, /data, /solid, thick=5, hthick=2.
  endelse

  ;; Light [0.25 to 0.3]
  plot, [0], [0], psym=state.psym, $
        position=[0.22, 0., 0.33, 1.], xrange=[0., 1.], $
        yrange=[0., 1.], xstyle=4, ystyle=4, $
        thick=lthick, title=state.title, background=clr.black, $
        charsize=csize, /noerase
  if abs(state.mirror_angle-45) LT atoler and MIRROR EQ 1 then $
    aclr = clr.white else aclr=clr.black
  arrow, 0.1, 0.6, 1., 0.6, /data, color=aclr, /solid, thick=5, hthick=2.
  arrow, 0.1, 0.5, 1., 0.5, /data, color=aclr, /solid, thick=5, hthick=2.
  arrow, 0.1, 0.4, 1., 0.4, /data, color=aclr, /solid, thick=5, hthick=2.

  
  ;; Spectrum [0.6 to 1.]
  
;;  off = 0.02
  
  device, decomposed=0
  loadct, 0, /silent
  state.ncolors = !d.table_size - 9
  state.contrast=0.5
  state.brightness=0.5
  ximgd_getct, state, 13, /CLR

  xwid = 0.2 * totx 
  yclr = round(0.8 * ysz)

  ;; Map for 4000A to 7000A
  nwav = 3000L
  cmap = congrid( findgen(state.ncolors), nwav) + 8
  wvm = 4000. + findgen(nwav)

  ;; Get endpoints (full screen)
  dum = state.xcurs
  state.xcurs = 0.
  x0 = 3500.
  x1 = 7500.
  state.xcurs = dum
  xv = x0 + findgen(yclr)*(x1-x0)/float(yclr-1) 

  ;; Interpolate colors
  b = interpol(cmap, wvm, xv) > cmap[0]
  cuttop = where(xv GT 7500., ncut)
  if ncut NE 0 then b[cuttop] = cmap[0]
  c = replicate(1, xwid)
  a = c # b
  if abs(state.mirror_angle-45) GT atoler OR $
    (MIRROR NE 1) OR (PRISM NE 2) then a[*] = cmap[0]

  tv, a, 0.75*totx, 0.1*ysz
  device, decomposed=1

  clr = getcolor(/load)
  off = 0.02
  plot, [0], [0], psym=state.psym, $
        position=[0.75-off,0.2-off,0.95+off,0.9+off], xrange=[0., 1.], $
        yrange=[0., 1.], xstyle=5, ystyle=5, $
        thick=lthick, title=state.title, background=clr.white, $
            charsize=csize, /noerase
  oplot, [0., 0., 1., 1.,0.], [0., 1., 1., 0.,0.], color=clr.gray, thick=4
  xyouts, 0.85, 0.1, 'SPECTRUM', color=clr.gray, charsiz=2.5, alignm=0.5,$
          /normal

  state.fiddle = 0
  if abs(state.mirror_angle-45) GT atoler AND $
    MIRROR EQ 1 AND PRISM EQ 2 then widget_control, state.overlay_id, sensitive=1$
  else widget_control, state.overlay_id, sensitive=0

  return
end

pro outr_basicspectra_PltIMG, state

  loadct, 0, /silent
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  ximgd_getct, state, 13, /CLR
  ;; Check extension
  ipos = strpos(state.imgfil, '.', /reverse_sear)
  case strmid(state.imgfil,ipos+1) of
      'jpg': read_jpeg, state.imgfil, img ;, ctable, COLORS=state.ncolors
      else: stop
  endcase
  widget_control, state.img_draw_id, get_value=wind
  wset, wind 
  xsz = round((widget_info(state.img_draw_id, /geometry)).xsize)
  
  nw_img = congrid(img, 3, xsz, xsz)
  tv, nw_img, /true

 return
end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro outr_basicspectra_Reset, state

  ;; Set the mode
  outr_basicspectra_setmode, state

; Plotting
  state.xymnx = state.svxymnx

; Sort wave
  srt = sort(state.spectrum.wave)
  state.spectrum.wave = state.spectrum.wave[srt]
  state.spectrum.fx = state.spectrum.fx[srt]
  
  
end

;;;;;;;;;;;;;;;;;;;;
;  Text
;;;;;;;;;;;;;;;;;;;;

pro outr_basicspectra_text, state, txt

  case state.mode of
      0: begin ;; Experiment
          txt = 'The experiment above shows a simple "spectrometer", a'+$
                ' device which splits light apart into smaller pieces (color).'+$
                ' In this example, we use a prism which is simple piece of glass.'+$
                ' Once the sunlight is directed through the prism, one may add'+$
                ' up how much light is in each color.  This is called ' + $
                'a spectrum, a product'+$
                ' that serves many purposes in astronomy.'
      end
      else:
  endcase
;  widget_control, state.txt, set_value=txt

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

pro outr_basicspectra, XSIZE=xsize, YSIZE=ysize, TITLE=title $
  , NOBLOCK=noblock, IMGFIL=imgfil, DEFAULT=default, SPIRAL=spiral, $
  SPECTRUM=spectrum

  common xcommon_color

;
  ;; Default
  solar_spec = {overlaystrct}
  solar_spec.name = 'Sun'
  outr_overlays, solar_spec, /init, wvmnx=[3000., 9000]
  imgfil = getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/sun.jpg'
  solar_spec.fx = solar_spec.fx / max(solar_spec.fx)

  outr_basicspectra_common

  ;; Sun
  imgfil = getenv('OUTR_DIR')+'Outreach/Spectroscopy/Images/sun.jpg'
  ipos = strpos(imgfil, '.', /reverse_sear)
  case strmid(imgfil,ipos+1) of
      'jpg': read_jpeg, imgfil, img ;, ctable, COLORS=state.ncolors
      else: stop
  endcase

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

;  outr_basicspectra_initcommon

  ;; Plot range (initial)
  xmin = 3000.
  xmax = 7500.
  ymin = 0.
  ymax = 1.

;    STATE
  state = { spectrum: solar_spec, $
            mode: keyword_set(SPECTRUM), $
            flg_zoom: 0, $
            flg_logx: 0, $
            imgfil: imgfil, $
            mirror_angle: 0., $
            sunimg: img, $
            wvmnx: [3000., 9000], $
            fiddle: 0, $
            sources: ['Sun', 'Spiral Galaxy', $
                      'Sodium Lamp', 'Quasar'], $
            pos: [0.2,0.2,0.95,0.90], $ ; Plotting
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
            mode_id: 0L, $
            spec_id: 0L, $
            mirrorpos_id: 0L, $
            prismpos_id: 0L, $
            slider_id: 0L, $
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
  base = WIDGET_BASE( title = 'outr_basicspectra -- Version 1.0', /column, $
                    xoffset=50, yoffset=50)
  state.base_id = base
  
  ;;;;;;;;;;;;
  ;;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=1)
  state.size[0] = xsize
  state.size[1] = (3./5)*ysize

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=round(state.size[1]), /frame, $
                              retain=2, $
                              uvalue='DRAW')

  ;;;;;;;;;
  ;;  Stuff
          
  imgsz = round(xsize/3.) < round(ysize/3.)
  stuff = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center, ysize=imgsz)

; zabs

  explorebar = WIDGET_BASE( stuff, /row, /frame, /base_align_center,$
                           /align_left, xsize=xsize-imgsz)
  main_base = WIDGET_BASE( explorebar, /column, /frame, /base_align_center,$
                           /align_center)
  mode_base = WIDGET_BASE( main_base, /column, /frame, /base_align_center,$
                           /align_center)

  mode_tit = widget_label(mode_base, value='Mode:',$
                      uvalue='TMODE', FONT=state.fonts.big_font)
  state.mode_id = widget_droplist(mode_base, $
                                  value=['Experiment', 'Spectrum'], $
                                  uvalue='MODE', FONT=state.fonts.big_font)
  widget_control, state.mode_id, sensitive=0
  prism_tit = widget_label(mode_base, value='Prism Position:',$
                      uvalue='TMIRR', FONT=state.fonts.big_font)
  state.prismpos_id = widget_droplist(mode_base, $
                                  value=['None','Position 1', 'Position 2'], $
                                  uvalue='PRISM_POS', FONT=state.fonts.big_font)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Slider
  tools_base = WIDGET_BASE( explorebar, /column, /frame, /base_align_center,$
                           /align_center)

  mirror_tit = widget_label(tools_base, value='Mirror Position:',$
                      uvalue='TMIRR', FONT=state.fonts.big_font)
  state.mirrorpos_id = widget_droplist(tools_base, $
                                  value=['None','Position 1', 'Position 2'], $
                                  uvalue='MIRROR_POS', FONT=state.fonts.big_font)
  state.slider_id = cw_fslider(tools_base, MAX=90., MIN=0., $
;                               xsize=strct.xsize[indx], $
;                               /DOUBLE, $
                               /DRAG, $
                               title='Mirror Angle',$
                               VALUE=0., $
                               UVALUE='MIRROR')
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Overlays
  source_base = WIDGET_BASE( explorebar, /column, /frame, /base_align_left,$
                           /align_center)
  src_tit = widget_label(source_base, value='Source:',$
                      uvalue='TSRC', FONT=state.fonts.big_font)
  state.spec_id = widget_droplist(source_base, sensitive=0, $
                                  value=state.sources, $
                                  uvalue='SOURCE', FONT=state.fonts.big_font)

  overlay = ['Histogram', 'Exact Spectrum']
  gd = where(strlen(state.overlays.name) GT 0, ngd)
  ;; Create button group
  state.overlay_id = cw_bgroup(source_base, overlay, $
                               /nonexclusive, row=3, $
                               UVALUE='OVERLAYS',frame=1, $
                               LABEL_TOP='Spectral Plot',$
                               FONT=state.fonts.big_font)
  widget_control, state.overlay_id, sensitive=0


  ;; Objname
  if keyword_set( head ) then begin
      objnm_id = cw_field(explorebar, title='Object: ', value=spectrum.name, $
                          /column, $
                          xsize=strlen(spectrum.name))
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Done
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  buttons = WIDGET_BASE( explorebar, /column, /frame, /base_align_left,$
                           /align_center)
  exercises = WIDGET_BUTTON(buttons, value='Exercises',uvalue='EXERCISES', $
                            /align_right, $
                            font=state.fonts.big_font)
  done = WIDGET_BUTTON(buttons, value='Done',uvalue='DONE', /align_right, $
                      font=state.fonts.big_font)

  outr_basicspectra_text, state, txt
  state.text_id = widget_text(explorebar, value=txt, xsize=50,$
                                    ysize=10, /scroll, /wrap, $
                                    FONT=state.fonts.big_font)

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
  device, set_FONT='Times Bold', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Update spectrum
  outr_basicspectra_Reset, state
  outr_basicspectra_UpdatePlot, state

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color Bar

  if keyword_set(IMGFIL) then outr_basicspectra_PltIMG, state

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Send to the xmanager
  if not keyword_set(BLOCK) then $
    xmanager, 'outr_basicspectra', base, no_block=noblock $
  else xmanager, 'outr_basicspectra', base

return
end
