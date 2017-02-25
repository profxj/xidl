;+ 
; NAME:
; outr_wavelngthtool 
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_wavelngthtool_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
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
;   outr_wavelngthtool_ 'spec.fits'
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
pro outr_wavelngthtool_common
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  return
end


;;;;
; Events
;;;;

pro outr_wavelngthtool_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'PHOTON_COUNTER_TEXT' : x_widget_displayfile, getenv('XIDL_DIR')+ $
        '/Outreach/Spectroscopy/Text/photons.txt', /WRAP
      'WAVEB_DROP' : begin
          state.curwaveb = ev.index
          state.wavelength = alog10(mean(state.wavebstrct[state.curwaveb].wv_mnx))
          outr_wavelngthtool_newwave, state
      end
      'DRAW_BASE' : begin
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'DRAW' : begin
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = Start stretching
                          state.xcurs = ev.x
                          state.ycurs = ev.y
;                          state.xclick = xgetx_plt(state, /strct)
                          state.xclick = ev.x
;                          state.yclick = xgety_plt(state, /strct)
                          state.wv_click = state.wavelength
                          state.flg_click = 1
                      end 
                      else: 
                  endcase
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

                  ;; Modify the wavelength
                  dx = (ev.x-state.xclick)/state.size[0]
                  state.wavelength = state.wv_click + dx
;                  state.wavelength = state.wv_click * 10.^(dx)
                  outr_wavelngthtool_newwave, state
                  state.ycurs = ev.y
              end
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; Bound the wavelengths (Angstroms)
;  state.xymnx[0] = state.xymnx[0] > 1e-10
  state.xymnx[2] = state.xymnx[2] < 1e20

  ;; Update Plot
  outr_wavelngthtool_UpdatePlot, state

  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;
pro outr_wavelngthtool_UpdatePlot, state
  
;common outr_wavelngthtool_lines

; Plot Data
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  lthick = 3.
  csize = 2.1
  lsz = 2.0
  lsz_big = 2.5

  clr = getcolor(/load)

  state.pow_wave = 10.^state.wavelength

  ;; 
  plot, [0],  [0], $
        position=state.pos, $
        xrange=[state.xymnx[0],state.xymnx[2]]/ $
        state.sizestrct[state.cursize].wv_scl, $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=5, ystyle=9, $
        xtitle='Length (Size)', ytitle='Power (Energy per second)', $
        thick=lthick, background=clr.white, $
        color=clr.black, charsize=csize

  oplot, state.xplt/state.sizestrct[state.cursize].wv_scl, $
         state.yplt, color=clr.black, thick=lthick

  ;;;;;;;;;;;;;;
  ;; Label wavelength 'sign'
  sign_clr = clr.royalblue
  yval = (max(state.yplt)+(state.xymnx[3]*0.15)) < state.xymnx[3]*0.95
  xval = [0.25,1.25]*state.pow_wave/state.sizestrct[state.cursize].wv_scl
  xmn = mean(xval)*0.75
  plotsym, 0, 2., /fill
  oplot, xval, replicate(yval, 2), color=sign_clr,  thick=5, psym=-8
  xyouts, xmn, yval+0.03*state.xymnx[3], 'Wavelength', $
          color=sign_clr,  charsiz=lsz, alignment=0.5
  xyouts, xmn, yval-0.07*state.xymnx[3], $
          string(state.pow_wave/state.sizestrct[state.cursize].wv_scl,$
          format='(f7.1)')+$
          ' '+state.wave_lbl, $
          color=sign_clr,  charsiz=lsz, alignment=0.5
  

  ;;;;;;;;;;;;;;;;;
  ;; x-axis
  xrng=[state.xymnx[0],state.xymnx[2]]/state.sizestrct[state.cursize].wv_scl
  plot, [0],  [0],  psym=psym, xrange=xrng, $
        yrange=[state.xymnx[1],state.xymnx[3]], color=clr.black, $
        position=state.pos, charsiz=csize, $
        thick=lthick, $
        background=clr.white, ystyle=5, xstyle=5, $
        /nodata, /noerase
  axis, xaxis=0, xrange=xrng, charsize=csize, xstyle=1, color=clr.black, $
        xtitle='!17Length'
  xyouts, state.pos[2], state.pos[1]-0.10, state.sizestrct[state.cursize].wv_lbl,$
          charsize=lsz, color=clr.black, /normal, alignment=0.5
  
  ;;;;;;;;;;;;;;;;;
  ;; Waveband
  xyouts, state.pos[2]*0.8, state.pos[3]*1.10, state.wavebstrct[state.curwaveb].name,$
          charsize=lsz_big, color=clr.black, /normal, alignment=0.5
  
  !p.font = 1
  device, set_FONT='Times', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  return

end

;;;;;;;;;;;;;;;;;;;;
;  Make Sine wave
;;;;;;;;;;;;;;;;;;;;

pro outr_wavelngthtool_mkwaves, state

  ;; Simple Sine wave
  state.yplt = (sin(2*!pi*state.xplt/(10.^state.wavelength)) + 1.)/2.

end

;;;;;;;;;;;;;;;;;;;;
;  Reset x-values
;;;;;;;;;;;;;;;;;;;;

pro outr_wavelngthtool_reset_xplt, state

  state.pow_wave = 10.^state.wavelength

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find the Size Struct
  nsize = n_elements(state.sizestrct)
  state.cursize = -1
  if state.pow_wave GT state.sizestrct[nsize-1].wv_min then $
    state.cursize = nsize-1
  if state.pow_wave LT state.sizestrct[0].wv_min then $
    state.cursize = 0
  if state.cursize LT 0 then begin
      gd = where(state.pow_wave GT state.sizestrct.wv_min AND $
                 state.pow_wave LT shift(state.sizestrct.wv_min,-1), ngd)
      state.cursize = gd[ngd-1]
  endif
  state.xymnx[2] = state.sizestrct[state.cursize].wv_axis
  state.wave_lbl = state.sizestrct[state.cursize].wv_lbl

  ;;  Reset position for fiddling
 ; state.xclick = xgetx_plt(state,/strct)
 ; state.wv_click = state.pow_wave

  ;; Make the x values
  state.xplt = findgen(state.npix)*state.xymnx[2]/float(state.npix - 1)

  ;; Upload the Image
  imgfil = state.sizestrct[state.cursize].size_imgfil
  ipos = strpos(imgfil, '.', /reverse_sear)
  case strmid(imgfil,ipos+1) of
      'jpg': read_jpeg, imgfil, img ;, ctable, COLORS=state.ncolors
      else: stop
  endcase
  widget_control, state.size_draw_id, get_value=wind
  wset, wind 
  xsz = round((widget_info(state.size_draw_id, /geometry)).xsize)
  
  nw_img = congrid(img, 3, xsz, xsz)
  tv, nw_img, /true
  close, /all

  ;; Label
  widget_control, state.size_lbl_id, $
                  set_value=state.sizestrct[state.cursize].size_imglbl

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find the Telescope Struct
  ntele = n_elements(state.telestrct)
  state.curtele = -1
  if state.pow_wave GT state.telestrct[ntele-1].wv_mnx[0] then $
    state.curtele = ntele-1
  if state.pow_wave LT state.telestrct[0].wv_mnx[1] then $
    state.curtele = 0
  if state.curtele LT 0 then begin
      gd = where(state.pow_wave GT state.telestrct.wv_mnx[0] AND $
                 state.pow_wave LT shift(state.telestrct.wv_mnx[1],-1), ngd)
      state.curtele = gd[ngd-1]
  endif

  ;; Upload the Image
  imgfil = state.telestrct[state.curtele].tele_imgfil
  ipos = strpos(imgfil, '.', /reverse_sear)
  case strmid(imgfil,ipos+1) of
      'jpg': read_jpeg, imgfil, img ;, ctable, COLORS=state.ncolors
      else: img = bytarr(3,100,100)
  endcase
  widget_control, state.telescope_draw_id, get_value=wind
  wset, wind 
  xsz = round((widget_info(state.telescope_draw_id, /geometry)).xsize)
  ysz = round((widget_info(state.telescope_draw_id, /geometry)).ysize)
  
  nw_img = congrid(img, 3, xsz, ysz)
  tv, nw_img, /true

  ;; Label
  widget_control, state.telescope_lbl_id, $
                  set_value=state.telestrct[state.curtele].telescope
  ;; Text
  text = x_readtxtfil(state.telestrct[state.curtele].tele_textfil)

  widget_control, state.telescope_notes_id, set_value=text
                  


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Waveband
  nwaveb = n_elements(state.wavebstrct)
  state.curwaveb = -1
  if state.pow_wave GT state.wavebstrct[nwaveb-1].wv_mnx[0] then $
    state.curwaveb = nwaveb-1
  if state.pow_wave LT state.wavebstrct[0].wv_mnx[1] then $
    state.curwaveb = 0
  if state.curwaveb LT 0 then begin
      gd = where(state.pow_wave GT state.wavebstrct.wv_mnx[0] AND $
                 state.pow_wave LT shift(state.wavebstrct.wv_mnx[1],-1), ngd)
      state.curwaveb = gd[ngd-1]
  endif
  widget_control, state.waveb_drop_id, set_combobox_select=state.curwaveb

end

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;
  pro outr_wavelngthtool_emspec, state

  ;; Upload the Image
  imgfil = getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Images/em_spectrum.jpg'
  read_jpeg, imgfil, img ;, ctable, COLORS=state.ncolors
  widget_control, state.EM_draw_id, get_value=wind
  wset, wind 
  xsz = round((widget_info(state.EM_draw_id, /geometry)).xsize)
  ysz = round((widget_info(state.EM_draw_id, /geometry)).ysize)
  
  nw_img = congrid(img, 3, xsz, ysz)
  tv, nw_img, /true
  close, /all

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;
;  Wavelength has changed so check about changing other things
;;;;;;;;;;;;;;;;;;;;;;;;;

pro outr_wavelngthtool_newwave, state

  state.pow_wave = 10.^state.wavelength

  ;; Photon energy
  state.E_photon = 2d-15 / (state.pow_wave)  ; E in Joules, Wave in Angstroms
  Nphoton = state.luminosity / state.E_photon  ;; Photons per second
  widget_control, state.photon_id, set_value=Nphoton

  ;; Minimum
  state.wavelength = state.wavelength > state.min_wave
  ;; Check for new axis size
  if state.pow_wave LT state.sizestrct[state.cursize].wv_min $
    OR state.pow_wave GT $
    state.sizestrct[(state.cursize+1)< $
                    (n_elements(state.sizestrct)-1)].wv_min then begin
      outr_wavelngthtool_reset_xplt, state
;      state.flg_click = 0
  endif

  ;; Check for new telescope
  if state.pow_wave LT state.telestrct[state.curtele].wv_mnx[0] $
    or state.pow_wave GT state.telestrct[state.curtele].wv_mnx[1]  then begin
      outr_wavelngthtool_reset_xplt, state
;      state.flg_click = 0
  endif

  ;; Check for new waveband
  if state.pow_wave LT state.wavebstrct[state.curwaveb].wv_mnx[0] $
    or state.pow_wave GT state.wavebstrct[state.curwaveb].wv_mnx[1]  then begin
      outr_wavelngthtool_reset_xplt, state
  endif

  ;; Remake the wave
  outr_wavelngthtool_mkwaves, state

  return
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;
pro outr_wavelngthtool_psfile, state

; Device
  device, get_decomposed=svdecomp

  x_psopen, 'idl.ps', /maxs
  state.psfile = 1
  outr_wavelngthtool_UpdatePlot, state
  x_psclose
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  Size Struct
;;;;;;;;;;;;;;;;;;;;

function outr_wavelngthtool_sizestrct

  tmp = {sizestrct, $
         wv_axis: 0., $
         wv_min: 0., $
         wv_scl: 0., $
         wv_lbl: '', $
         human_lbl: '', $
         size_imgfil: '', $
         size_imglbl: '' $
        }

  ;; Start it up
  readcol, getenv('XIDL_DIR')+'/Outreach/Spectroscopy/outr_sizestrct.dat', $
           all_wvaxis, all_wvmin, all_wvlbl, all_hlbl, all_wvscl, $
           all_sizeimgfil, all_sizeimglbl, $
           FORMAT='F,F,A,A,F,A,A', DELIMITER=','
;  all_wvaxis = [0.4, 2.,  8,   40., 200, 800., 4000., 20000,     1e5, 5e5, $
;                2e6, 1e7, 5e7, 2e8, 1e9, 5e9, 2e10]
;  all_wvmin  = [ 0.04, 0.2, 0.8, 4., 20., 80.,  400.,  2000., 10000., 50000., $
;               2e5, 1e6, 5e6, 2e7, 1e8, 5e8, 2e9]
;  all_wvlbl  = ['Angstroms', 'Angstroms', 'Angstroms', 'Angstroms', $
;                'Angstroms', 'Angstroms', $
;                'Angstroms', 'Angstroms', 'Microns', 'Microns', $
;                'Microns', 'Microns', 'Micros', 'Microns', 'cm', 'cm', 'cm' $
;               ]
;  all_hlbl  =  ['Atoms', 'Atoms', 'Atoms', 'Atoms', 'Atoms', 'Atoms', $
;                'Angstroms', $
;                'Angstroms', 'Hair Strands','TBD', 'TBD', 'TBD', 'TBD', $
;                'TBD', 'TBD', 'TBD', 'TBD']
;  all_wvscl  = [ 1., 1., 1., 1., 1., 1., 1., 1., 1e4, 1e4, 1e4, 1e4, 1e4, $
;                 1e4, 1e8, 1e8, 1e8]
  nsize = n_elements(all_wvaxis)
;  all_sizeimgfil  =  replicate('fluvirus.jpg',nsize)
;  all_sizeimglbl  =  replicate('Flu Virus:  Diameter = 1000 Angstroms',nsize)

  ;; Make the real structure
  all_size = replicate(tmp, nsize)
  all_size.wv_axis = all_wvaxis
  all_size.wv_min = all_wvmin
  all_size.wv_lbl = all_wvlbl
  all_size.wv_scl = all_wvscl
  all_size.human_lbl = all_hlbl
  all_size.size_imgfil = $
    getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Images/'+strtrim(all_sizeimgfil,2)
  all_size.size_imgLbl = all_sizeimglbl

  return, all_size

end

;;;;;;;;;;;;;;;;;;;;
;  Telescope Struct
;;;;;;;;;;;;;;;;;;;;

function outr_wavelngthtool_telestrct

  tmp = {telestrct, $
         telescope: '', $
         wv_mnx: fltarr(2), $
         name: '', $
         tele_imgfil: '', $
         tele_textfil: '' $
        }

  ;; Start it up
  name = ['Glast', 'Swift', 'Chandra', 'SOHO', 'Galex', 'Keck', $
          'Spitzer', 'JCMT', 'ALMA', 'Arecibo']
  all_wvmnx = [ $
              [4e-13,  0.1],$    ; Glast
              [0.1, 1], $       ; Swift
              [1, 50], $        ; Chandra
              [50, 1000], $     ; SOHO
              [1000, 3200], $   ; Galex
              [3200, 20000], $  ; Keck
              [20000, 3000000.], $ ; Spitzer
              [3e6, 3e7], $     ; JCMT
              [3e7, 3e8], $     ; ALMA
              [3e8, 7e9] $      ; Arecibo
                  ]
  
  nsize = n_elements(name)
  all_teleimgfil  =  ['glast.jpg', $
                      'swift.jpg', $
                      'chandra2.jpg', $
                      'soho.jpg', $
                      'galex.jpg', $
                      'keck.jpg', $
                      'spitzer.jpg', $ 
                      'jcmt.jpg', $
                      'alma.jpg', $
                      'arecibo.jpg']
  all_teletext  =  ['glast.txt', $
                    'swift.txt', $
                    'chandra.txt', $
                    'soho.txt', $
                    'galex.txt', $
                    'keck.txt', $
                    'spitzer.txt', $ 
                    'jcmt.txt', $
                    'alma.txt', $
                    'arecibo.txt']

  ;; Make the real structure
  all_size = replicate(tmp, nsize)
  all_size.telescope = name
  all_size.wv_mnx = all_wvmnx
  all_size.tele_imgfil = $
    getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Telescopes/'+all_teleimgfil
  all_size.tele_textfil = $
    getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Telescopes/'+all_teletext
;  all_size.tele_imglbl = all_sizeimglbl

  return, all_size

end

;;;;;;;;;;;;;;;;;;;;
;  Waveband Struct
;;;;;;;;;;;;;;;;;;;;

function outr_wavelngthtool_wavebstrct

  tmp = {wavebstrct, $
         name: '', $
         wv_mnx: fltarr(2), $
         waveb_imgfil: '', $
         waveb_textfil: '' $
        }

  ;; Start it up
  name = ['Gamma-rays', 'X-rays', 'Ultraviolet', 'Optical', $
          'Infrared', 'Submm', 'Millimeter', 'Radio']
  all_wvmnx = [ $
              [4e-13,  0.1],$   ; Gamma
              [0.1, 70], $     ; X-ray
              [70, 3500], $    ; Ultraviolet
              [3500, 8000], $   ; Optical
              [8000, 3e6], $    ; IR
              [3e6, 1e7], $     ; submm
              [1e7, 1e8], $     ; mm
              [3e8, 7e9] $      ; Radio
                  ]
  
  nsize = n_elements(name)
  ;; Make the real structure
  all_size = replicate(tmp, nsize)
  all_size.name = name
  all_size.wv_mnx = all_wvmnx
;  all_size.tele_imgfil = $
;    getenv('XIDL_DIR')+'/Outreach/Spectroscopy/Images/'+all_sizeimgfil
;  all_size.tele_imglbl = all_sizeimglbl

  return, all_size

end

;;;;;;;;;;;;;;;;;;;;
;  Reset 
;;;;;;;;;;;;;;;;;;;;

pro outr_wavelngthtool_Reset, state, INWAV=inwav

  if keyword_set(INWAV) then state.wavelength = inwav 

  ;; 
  state.xymnx[3] = 1.2
      
  outr_wavelngthtool_newwave, state
  outr_wavelngthtool_reset_xplt, state
;  outr_wavelength_reset_yplt, state
  outr_wavelngthtool_mkwaves, state

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

pro outr_wavelngthtool, XSIZE=xsize, YSIZE=ysize, TITLE=title $
                , DUNIT = dunit, BLOCK=block, IMGFIL=imgfil, DEFAULT=default

  common xcommon_color

  outr_wavelngthtool_common

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

;    STATE
  state = { zabs: 0., $
            wavelength: alog10(5500.), $ ;; Angstroms
            pow_wave: 5500., $
            luminosity: 100.d , $  ;; Watts (J/s)
            E_photon: 0.d, $       ;; Joules
            min_wave: -4., $ ;; Angstroms (Log10)
            npix: npix, $
            xplt: findgen(npix), $
            yplt: findgen(npix), $
            flg_zoom: 0, $
            pos: [0.1,0.15,0.95,0.85], $ ; Plotting
            psfile: 0, $ ; Postscript
            svxymnx: [0., 0., 20000., 2.], $
            xymnx: dblarr(4), $
            sizestrct: outr_wavelngthtool_sizestrct(), $  ; Structures
            telestrct: outr_wavelngthtool_telestrct(), $
            wavebstrct: outr_wavelngthtool_wavebstrct(), $
            cursize: 0, $
            curtele: 0, $
            curwaveb: 0, $
            wave_lbl: '', $
            flux_lbl: '', $
            tmpxy: dblarr(4), $
            flg_logx: 0, $
            psym: 10, $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            xclick: 0., $
            yclick: 0., $
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
            EM_draw_id: 0L, $
            size_draw_id: 0L, $
            size_lbl_id: 0L, $
            overlay_id: 0L, $
            draw_id: 0L, $
            telescope_draw_id: 0L, $
            telescope_notes_id: 0L, $
            telescope_lbl_id: 0L, $
            photon_id: 0L, $
            waveb_drop_id: 0L, $
            img_draw_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }

  state.size[0] = xsize
  state.size[1] = ysize

  ;;    Title
;  if size(flux, /type) EQ 7 then state.title = flux
  if keyword_set( TITLE ) then state.title = title

;  device, set_FONT='Times', /TT_FONT

;    WIDGET
  base = WIDGET_BASE( title = 'outr_wavelngthtool -- Version 1.0', /row, $
                    xoffset=10, yoffset=10)
  state.base_id = base
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color bar
  lhs_base = widget_base(base, /column, /base_align_left, $
                         frame=1, xsize=state.size[0]*0.8)
  rhs_base = widget_base(base, /column, /base_align_left, $
                         frame=1, xsize=state.size[0]*0.2)
  ;;;;;;;;;;;;
  ;;      Drawing
  state.draw_base_id = widget_base(lhs_base, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=1)

  state.draw_id = widget_draw(state.draw_base_id, /frame, retain=2, $
                              ysize=state.size[1]*0.7, $
                              xsize=state.size[0]*0.8, $
                              /button_events, /motion_events, uvalue='DRAW')

  ;; Size Image
  size_base = widget_base(rhs_base, /column, /base_align_center, $
                          /align_center, $
                          frame=1, xsize=state.size[0]*0.2)
  state.size_draw_id = widget_draw(size_base, retain=2, $
                                   xsize=state.size[0]*0.2, $
                                   ysize=state.size[0]*0.2, $
                                   uvalue='SIZE_DRAW')
  state.size_lbl_id = widget_text(size_base, value='', $
                                  xsize=state.size[0]*0.18, $
                                  FONT=state.fonts.small_font)


  ;;;;;;;;;;;;;;;;;;;;;
  ;; Waveband
  state.waveb_drop_id = widget_combobox(rhs_base, $
                                  VALUE=state.wavebstrct.name, $
                                  UVALUE='WAVEB_DROP', FONT=state.fonts.big_font)
  ;;;;;;;;;
  ;; EM spectrum
  EM_base = widget_base(rhs_base, /column, /base_align_center, $
                          /align_center, $
                          frame=1, xsize=state.size[0]*0.2)
  state.EM_draw_id = widget_draw(EM_base, retain=2, $
                                   xsize=state.size[0]*0.2, $
                                   ysize=state.size[1]*0.5, $
                                   uvalue='EM_DRAW')

  ;;;;;;;;;
  ;;  Lower bar
  explorebar = WIDGET_BASE( lhs_base, /row, /frame, /base_align_center,$
                           /align_left, ysize=ysize*0.3, xsize=xsize*0.8)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Telescope
  telescope_base = WIDGET_BASE( explorebar, /row, /frame, /base_align_left,$
                                /align_center)

  tele2_base = WIDGET_BASE( telescope_base, /column, /frame, /base_align_center,$
                                /align_center)
  state.telescope_lbl_id = widget_text(tele2_base, value='', xsize=20, $
                                       FONT=state.fonts.big_font)
  state.telescope_draw_id = widget_draw(tele2_base, xsize=xsize*0.3, $
                                        ysize=ysize*0.2, retain=2, $
                                        uvalue='TELESCOPE_DRAW')
  state.telescope_notes_id = widget_text(telescope_base, xsize=50, $
                                         value='',$
                                         ysize=10, /scroll, /wrap, $
                                        FONT=state.fonts.small_font)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Photon Counter
  photon_base = WIDGET_BASE( explorebar, /column, /frame, /base_align_center,$
                                /align_center, xsize=xsize*0.2)

  txt = 'PHOTON COUNTER'
  photon_title = widget_text(photon_base, value=txt, xsize=strlen(txt), $
                                       FONT=state.fonts.big_font)
  state.photon_id = cw_field(photon_base, $
                             title='Photons per Second =', value=0.d, $
                             xsize=13, FONT=state.fonts.small_font)
  photon_text = WIDGET_BUTTON(photon_base, value='How does this counter work?',$
                              uvalue='PHOTON_COUNTER_TEXT', /align_right)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; TOOLS

  ;;;;;;;;;;
  ;; Fiddle
  fiddle_base = WIDGET_BASE( explorebar, /row, /frame, /base_align_left,$
                             /align_center, xsize=xsize*0.2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Done
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  done = WIDGET_BUTTON(fiddle_base, value='Done',uvalue='DONE', /align_right)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Realize
  WIDGET_CONTROL, base, /realize

  !p.font = 1
  device, set_FONT='Times', /TT_FONT, SET_CHARACTER_SIZE=[10,10]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Update spectrum
  outr_wavelngthtool_Reset, state
  outr_wavelngthtool_UpdatePlot, state
  outr_wavelngthtool_emspec, state

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Color Bar
  loadct, 0, /silent
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)
  ximgd_getct, state, 13, /CLR
  

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'outr_wavelngthtool', base, /no_block $
  else xmanager, 'outr_wavelngthtool', base

return
end
