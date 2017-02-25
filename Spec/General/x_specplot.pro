;+ 
; NAME:
; x_specplot   
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   x_specplot, flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           /QAL, /GAL, ZIN=, /BLOCK, /NRM, /LLS, /QSO, /LBG
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
;   /LBG       - use the Lyman Break Galaxy line list (Shapley
;                et al '03)
;   /STAR      - use the star list (derived from gal.lst)
;   XRANGE=    - Opening plot x-axis (and default)
;   YRANGE=    - Opening plot y-axis (and default)
;   /GUI       - Choose from a list of FITS files (current directory)
;   /ASYM_SIG  - Plot asymmetric error bars (e.g. Poisson data from HST/COS)
;  /DISP       - Print dispersion of the spectrum per pixel (km/s)
;   linid_fil= - CSV file of zabs,wrest,label to label in GUI
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_specplot, 'spec.fits'
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
;   21-May-2008 Added XRANGE option, KLC
;   16-Oct-2008 Enable flux, ysin to be arrays, KLC
;;  29-Apr-2010 Added Shapley et al (2003) LBG composite, fixed 'g',
;;              updated Help table, KLC
;;   4-Feb-2014 Added /stars, KLC
;   21-Feb-2014 Handle 2+ portrait monitors, KLC
;   25-Feb-2014 Enable linid_fil to always mark identified lines, KLC
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_specplot_initcommon
;
common x_specplot_lines, $
   flg_lines, $
   lines, $
   linid, $
   zabs

  flg_lines = 0
  zabs = 0.
end

;;;;
; Events
;;;;

pro x_specplot_event, ev

common x_specplot_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          if state.flg_lines EQ 0 then begin
              state.flg_lines = 4 ;; LLS is the default
              x_specplot_initLines, state
          endif
          zabs = tmp
          flg_lines = 1
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
                      4 : if state.flg_lines NE 0 then $
                            x_specplot_SetLine, state ; Set reference line
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
              'S': x_specplot_smooth, state
              'U': x_specplot_smooth, state, /reset
              's': x_specplot_smooth, state, /YTWO
              'u': x_specplot_smooth, state, /reset, /YTWO
              ; Crude Analysis
              'E': x_specplot_EW, state        ; Calc EW
              'N': x_specplot_Colm, state      ; Calc AODM colm
              'n': x_specplot_SN, state        ; S/N
              'G': x_specplot_Gauss, state      ; Fit a Gaussian
              ' ': print, 'x: '+strtrim(state.xpos,2)+$
                '   y:'+strtrim(state.ypos,2)
              ; LINES
              '*': begin ; Set as MgII or [OIII] (for z~1 QSOs)
                 mrkwav = xgetx_plt(state, /strct)
                 zabs = mrkwav / 2799.4 - 1.d
                 widget_control, state.zabs_id, set_value=strmid(strtrim(zabs,2),0,10)
                 flg_lines = 1
              end
              '&': begin ; Set as MgII or [OIII] (for z~1 QSOs)
                 mrkwav = xgetx_plt(state, /strct)
                 zabs = mrkwav / 5008.24- 1.d
                 widget_control, state.zabs_id, set_value=strmid(strtrim(zabs,2),0,10)
                 flg_lines = 1
              end
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
              'M': begin ; Plot MgII
                  x_specplot_guess, state, 'M'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'V': begin ; Plot SiIV
                 x_specplot_guess, state, 'V'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'O': begin ; Plot OVI 
                  x_specplot_guess, state, 'O'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'L': begin ; Plot LLS
                  x_specplot_guess, state, 'L'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'K': begin ; Plot CaHK
                  x_specplot_guess, state, 'CaHK'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              ;; Damped Lya overplot
              'D': x_specplot_initdla, state
              ;; Beta Overplot
              'B': x_specplot_initdla, state, /beta
              ;; Super-LLS overplot
              'R': x_specplot_initslls, state, /esi
              ;; Postscript
              'P': x_specplot_psfile, state  
              ;; Overplot a Quasar spectrum
              'Q': x_specplot_qsotempl, state  
              ;; Overplot a Atmospheric transmission
              'I': x_specplot_atm, state  
              ;; Overplot a galaxy spectrum
              'g': x_specplot_galtempl, state  
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'x_specplot: Not a valid key!' ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  xspecplot_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xspecplot_UpdatePlot, state
  
common x_specplot_lines

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

  if state.flg_smooth EQ 0 then $
     plot, state.wave, state.fx, psym=state.psym, $
           position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
           yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
           xtitle='!17Wavelength', ytitle='Flux', $
           title=state.title, $
           background=clr.white, $
           color=clr.black, $
           xcharsize=1.9, $
           ycharsize=1.9 $
  else $
     plot, state.wave, state.smooth, psym=state.psym, $
           position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
           yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
           xtitle='!17Wavelength', ytitle='Flux', $
           title=state.title, $
           background=clr.white, $
           color=clr.black, $
           xcharsize=1.9, $
           ycharsize=1.9 

  ;; YWO
  if state.flg_ytwo EQ 1 AND state.flg_smooth2 EQ 0 then $
    oplot, state.wave_two, state.ytwo, psym = state.psym2, color = clr.purple $
  ELSE if state.flg_ytwo EQ 1 AND state.flg_smooth2 GT 0 THEN $
    oplot, state.wave_two, state.smooth2, psym = state.psym2, color = clr.purple 

  ;; YTHREE
  if state.flg_ythree EQ 1 THEN $
     oplot, state.wave_three, state.ythree, psym = state.psym3, color = clr.blue

  ;; YFOUR
  if state.flg_yfour EQ 1 THEN $
     oplot, state.wave_four, state.yfour, psym = state.psym4, color = clr.cyan

  ;; Plot Error array
  case state.flg_sig of
     1: oplot, state.wave, state.sig, psym=state.psym, color=clr.red
     2: begin
        ;; Error array
        oploterror, state.wave, state.fx, state.asig1, /hibar, errcolor=clr.red, psym=1
        oploterror, state.wave, state.fx, state.asig2, /lobar, errcolor=clr.red, psym=1
     end
     else: 
  endcase


  ;; Line ID
  if state.flg_linid eq 1 then x_specplot_pltlinid, state

  ;; Plot Lines as required
  if flg_lines EQ 1 then x_specplot_PltLines, state

  ;; EW
  if state.flg_EW NE 0 then x_specplot_PltEW, state

  ;; 
  if state.flg_GS NE 0 then x_specplot_PltGS, state

;  if state.flg_EW EQ 1 then oplot, [state.EW_lmt[0,0]], [state.EW_lmt[0,1]], $
;    psym=2, color=getcolor('red')

  ;; Colm
  if state.flg_Colm NE 0 then x_specplot_PltClm, state

  ;; DLA
  if state.flg_DLA NE 0 then x_specplot_PltDLA, state

  ;; QSO
  if state.flg_qsotempl NE 0 then x_specplot_pltqsot, state

  ;; IR Atmospheric transmission
  if state.flg_atm NE 0 then x_specplot_pltatm, state
  ;; GAL
  if state.flg_galtempl NE 0 then x_specplot_pltgalt, state
end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro xspecplot_Reset, state


; Plotting
  state.xymnx = state.svxymnx

; Sort wave
  srt = sort(state.wave)
  state.wave = state.wave[srt]
  state.fx = state.fx[srt]
  state.sig = state.sig[srt]
  IF KEYWORD_SET(state.flg_ytwo) THEN BEGIN
     srt = sort(state.wave_two)
     state.wave_two = state.wave_two[srt]
     state.ytwo = state.ytwo[srt]
  ENDIF
  IF KEYWORD_SET(state.flg_ythree) THEN BEGIN
     srt = sort(state.wave_three)
     state.wave_three = state.wave_three[srt]
     state.ythree = state.ythree[srt]
  ENDIF
  
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the QAL line list

pro x_specplot_initQAL, llist

common x_specplot_lines
  lines = x_setllst(llist, 0)

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the H2 line list

pro x_specplot_initH2, llist

  common x_specplot_lines
  h2list = fuse_h2lin()
  tmp = { lliststrct }
  lines = replicate(tmp, n_elements(h2list))
  lines.wave = h2list.wrest
  lines.name = h2list.label

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the GAL line list
pro x_specplot_initGAL, llist

common x_specplot_lines
  lines = x_setllst(llist, 1)

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the MOL line list
pro x_specplot_initMOL, llist

common x_specplot_lines
  lines = x_setllst(llist, 2)

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the TEL line list
pro x_specplot_initTEL, llist

common x_specplot_lines
  lines = x_setllst(llist, 2)

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the user's line identification file

pro x_specplot_initlinid, state
  
  common x_specplot_lines
  tmp = file_info(state.linid_fil)    ;  since (re-)reading
  if tmp.mtime eq state.linid_mtime then return ; nothing new

  state.linid_mtime = tmp.mtime ; reset file modification time
  
  readcol, state.linid_fil, zlin, wrest, label, /silent, delimiter=',', $
           format='(f, f, a)'
  nid = (size(zlin,/dim))[0] > 1 
  tmp = {wave:0., zabs:0., wrest:0., label:''}
  linid = replicate(tmp,nid)
  
  linid.wave = wrest*(1+zlin)
  linid.zabs = zlin
  linid.wrest = wrest
  linid.label = strtrim(label,2) + ' z='+string(zlin,format='(f6.4)')
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sets the line and redshift with a gui
pro x_specplot_SetLine, state

common x_specplot_lines

  flg_lines = 1
  ; Set Line with gui
  case state.flg_lines of 
      1: setwave = x_slctline(lines, /ISM)
      2: setwave = x_slctline(lines, /GAL)
      3: setwave = x_slctline(lines, /GAL)
      4: setwave = x_slctline(lines, /ISM)
      5: setwave = x_slctline(lines, /ISM)
      6: setwave = x_slctline(lines, /ISM)
      7: setwave = x_slctline(lines, /GAL)
     10: setwave = x_slctline(lines, /GAL)
      else: setwave = x_slctline(lines)
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

pro x_specplot_PltLines, state, IMG=img

common x_specplot_lines

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.blue
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
          3: xyouts, xplt, ymax2, $ ; QSO
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          4: xyouts, xplt, ymax, $ ; LLS
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          5: xyouts, xplt, ymax, $ ; LLS
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          6: xyouts, xplt, ymax, $ ; GRB
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          7: xyouts, xplt, ymax2, $ ; LBG
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          10: xyouts, xplt, ymax2, $ ; QSO2
                      strtrim(lines[allwv[q]].name, 2), color = pclr, $
                      orientation = 90., charsize = 1.5
          else: xyouts, xplt, ymax, $
            strtrim(lines[allwv[q]].name,2)+$
            string(lines[allwv[q]].wave,format='(f7.1)'),$
            color=pclr, orientation=90., $
            charsize=1.5
      endcase
      ; Dotted line
      oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
        color=pclr, linestyle=1
  endfor

  ; Mark

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_specplot_pltlinid, state, IMG=img

common x_specplot_lines

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.limegreen
  endif else pclr = 3

 ;; Test whether need to reload linid_fil in the init call
  x_specplot_initlinid, state

  ; Parse list
  allwv = where(linid.wave GE state.xymnx[0] AND $
                linid.wave LE state.xymnx[2], cntwv)

  ; Plot
  ymax = state.xymnx[1] + 0.02*(state.xymnx[3]-state.xymnx[1])
  ymax2 = state.xymnx[1] + 0.8*(state.xymnx[3]-state.xymnx[1])

  for q=0L,cntwv-1 do begin
      xplt = linid[allwv[q]].wave
      ; Name
      xyouts, xplt, ymax2, linid[allwv[q]].label, color=pclr, $
              orientation=90., charsize=1.5
      ;; Dashed line
      oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
        color=pclr, linestyle=2 
  endfor

  ; Mark

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_specplot_qsotempl, state

common x_specplot_lines

  ;; Find Min and Max in region
  state.qsot_nrm = [state.xpos, state.ypos]
  state.flg_qsotempl = 1L

  return
end

pro x_specplot_pltqsot, state

common x_specplot_lines

  fx = x_readspec(getenv('XIDL_DIR')+'/Spec/SDSS/vanden.fits', wav=wav, $
                  inflg=2)
  owv = wav*(1+zabs)
  mn = min(abs(owv-state.qsot_nrm[0]),imn)
  nrm = state.qsot_nrm[1]/fx[imn]

  clr = getcolor(/load)
  oplot, wav*(1+zabs), fx*nrm, color=clr.blue, psym = state.psym ; 10

  return
end

pro x_specplot_atm, state

common x_specplot_lines

  ;; Find Min and Max in region
  state.atm_nrm = [state.xpos, state.ypos]
  state.flg_atm = 1L

  return
end

pro x_specplot_pltatm, state

common x_specplot_lines

  atm_file =  getenv('LONGSLIT_DIR') + '/calib/extinction/atm_trans_am1.0.dat'
  rdfloat, atm_file, wav,fx, skip = 2
  wav=wav*1d4
  mn = min(abs(wav-state.qsot_nrm[0]),imn)
  nrm = state.atm_nrm[1]/fx[imn]

  clr = getcolor(/load)
  oplot, wav, fx*nrm, color=clr.magenta
  return
end


pro x_specplot_galtempl, state

common x_specplot_lines

  ;; Find Min and Max in region
  state.galt_nrm = [state.xpos, state.ypos]
  state.flg_galtempl = 1L

  return
end

pro x_specplot_pltgalt, state

  common x_specplot_lines

  if strtrim(state.galtempl_fil,2) eq '' then $
     state.galtempl_fil = getenv('XIDL_DIR')+'/Spec/LBG/lbg_abs.dat'
  rdfloat, state.galtempl_fil, wav, fx, skip = 2
  widget_control, state.zabs_id, get_value=zabs ; retrieve current zabs
  owv = wav*(1+zabs)
  fx_rebin = 0*state.wave
  igd = where(state.wave GE min(owv) AND state.wave LE max(owv), ngd)

  IF ngd GT 0 THEN BEGIN
     fx_rebin[igd] = interpol(fx, owv, state.wave[igd])
     mn = min(abs(state.wave[igd]-state.galt_nrm[0]), imn)
     nrm = state.galt_nrm[1]/fx_rebin[igd[imn]]
  ENDIF ELSE nrm = 1.0d

  clr = getcolor(/load)

  oplot, state.wave, fx_rebin*nrm, color = clr.green, psym = state.psym ;10

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc EW
;  EW has [pts, x/y]

pro x_specplot_EW, state

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
      cntm = 2 / (state.EW_lmt[0,1]+state.EW_lmt[1,1]) ; local conti
      if npix NE 0 then begin
;          state.EW = int_tabulated(state.wave[sumpix], $
;                                   1-state.fx[sumpix]/cntm,$
;                                   /double)/(1.+zabs)
          zinv = 1./(1.+zabs)    ; just save and use
          dwv = state.wave[sumpix[npix-1]+1] - state.wave[sumpix[npix-1]]
          dwv = abs(dwv)
          state.EW = total(1.-(state.fx[sumpix]>0.)*cntm)*dwv * zinv
          ;; ERROR
          sumvar = total(state.sig[sumpix]^2)
          state.sigEW = sqrt(sumvar)*dwv * zinv * cntm 
          wvcent = 0.5*(state.EW_lmt[0,0] + state.EW_lmt[1,0]) * zinv

          print, 'Rest EW('+string(wvcent,format='(f7.2)')+') = ',$
                 strtrim(state.EW*1.e3,2), $
            ' +/- ', state.sigEW*1.d3, ' mA'
      endif
      
      ; Reset the flag
      state.flg_EW = 2
  endelse

  ; EW xlimit

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit a Gaussian

pro x_specplot_Gauss, state

common x_specplot_lines

  ; Set the flag
  if state.flg_GS MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_GS = 1 
      state.GS_lmt[0,0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.GS_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.GS_lmt[0] then begin
          state.GS_lmt[1,0] = tmp 
          state.GS_lmt[1,1] = xgety_plt(state, /strct)
      endif else begin
          state.GS_lmt[1,0] = state.GS_lmt[0,0]
          state.GS_lmt[1,1] = state.GS_lmt[0,1]
          state.GS_lmt[0,0] = tmp
          state.GS_lmt[0,1] = xgety_plt(state, /strct)
      endelse

      ;; Calc continuum
      m_con = (state.GS_lmt[1,1]-state.GS_lmt[0,1])/$
        (state.GS_lmt[1,0]-state.GS_lmt[0,0])
      b_con = state.GS_lmt[1,1] - m_con*state.GS_lmt[1,0]
      
      ;; Subtract continuum
      mn = min(abs(state.wave-state.GS_lmt[0,0]),ipx)
      mn = min(abs(state.wave-state.GS_lmt[1,0]),fpx)
      px = ipx + lindgen(fpx-ipx+1)
      state.GS_xpix = [ipx,fpx]

      gprof = state.fx[px] / (state.wave[px]*m_con + b_con)

      ;; Tip over as necessary
      if total(1.-gprof) GT 0. then begin
          gprof = 1 - gprof
          state.flg_GS = -2
      endif else state.flg_GS = 2

      ;; Wave centroid
      wcen = total(gprof*state.wave[px])/total(gprof)
      dwv = abs(state.wave[px[1]]-state.wave[px[0]])
;      sig_g = n_elements(px)/4. * (state.wave[px[1]]-state.wave[px[0]])
      gsssig = abs(min(state.wave[px], max=mxwv) - mxwv) / 6.
      
      ;; FIT
      yfit = gaussfit(state.wave[px], gprof, acoeff, $
                      estimates=[max(gprof), wcen, gsssig], $
                      sigma=sigma, nterms=3)
;                      estimates=[max(gprof), wcen, gsssig, 0.], $ for
;                      nterms=4
      if state.flg_GS EQ (-2) then yfit = 1. - yfit
      state.GS_fit[0:n_elements(yfit)-1] = yfit * $
        (state.wave[px]*m_con + b_con)
      print, 'x_specplot: Gaussian = ', acoeff
      print, 'x_specplot: sigGaussian = ', sigma
      print, 'x_specplot: FWHM (Ang, km/s, pix) = ', $
             2* sqrt(2.*alog(2.))* $
             [acoeff[2], acoeff[2]/acoeff[1]*3e5, acoeff[2]/dwv]
      ;; The following is wrong for EW
      print, 'x_specplot: Rest EW (Ang) = ', $
             acoeff[0] * acoeff[2] * sqrt(!pi*2.) / (1+zabs) ; Rest EW (Ang)
      state.gauss = acoeff

      ;; Report the EW
      
  endelse

  ; EW xlimit

return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot EW stuff to the screen

pro x_specplot_PltEW, state

  clr = getcolor(/load)
  ; flg
  case state.flg_EW of 
      1: oplot, [state.EW_lmt[0,0]], [state.EW_lmt[0,1]], $
        psym=2, color=clr.red
      2: begin
          oplot, [state.EW_lmt[0,0],state.EW_lmt[1,0]], $
            [state.EW_lmt[0,1],state.EW_lmt[1,1]], $
            psym=-2, color=clr.red, linestyle=2
          xyouts, 0.5, 0.97, 'Rest EW = '+string(state.EW*1.e3)+ $
            ' +/- '+strtrim(state.sigEW*1.e3,2)+' mA', $
            /normal, charsize=1.5, alignment=0.5, color=clr.red
      end
      else :
  endcase
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot Gauss stuff to the screen

pro x_specplot_PltGS, state

  clr = getcolor(/load)
  ; flg
  case abs(state.flg_GS) of 
      1: oplot, [state.GS_lmt[0,0]], [state.GS_lmt[0,1]], $
        psym=2, color=clr.red
      2: begin
          oplot, [state.GS_lmt[0,0],state.GS_lmt[1,0]], $
            [state.GS_lmt[0,1],state.GS_lmt[1,1]], $
            psym=-2, color=clr.red, linestyle=2
          xyouts, 0.1, 0.97, 'Gauss = '+$
            string(state.gauss, FORMAT='(4f12.4)'), $
            /normal, charsize=1.5, alignment=0.0, color=clr.red
          ;; Gaussian
          oplot, state.wave[state.GS_xpix[0]:state.GS_xpix[1]], $
            state.GS_fit[0:state.GS_xpix[1]-state.GS_xpix[0]+1], $
            color=clr.blue
      end
      else :
  endcase
  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc AODM colm
;  EW has [pts, x/y]

pro x_specplot_Colm, state

common x_specplot_lines

  ; Check for error array
  if state.flg_sig NE 1 then begin
      print, 'x_specplot_Colm: Need to specify the error array!'
      return
  endif

  ; Set the flag
  if state.flg_Colm EQ 0 then begin
      ; Set one limit
      state.flg_Colm = 1 
      state.Colm_lmt[0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.EW_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
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
                      7: gdwave = x_slctline(lines, /GAL) 
                      10: gdwave = x_slctline(lines, /GAL) 
                      else: gdwave = x_slctline(lines)
                  endcase
              end
              1: begin
                  print, 'Using '+strtrim(lines[gdlin].name,2)
                  gdwave = lines[gdlin].wave
              end
              else: begin
                  ;; 
                  gwv = x_guilist(strtrim(lines[gdlin].wave), indx=imn, MAXY=40)
                  print, 'Taking: '+strtrim(lines[gdlin[imn]].wave,2)
                  gdwave = lines[gdlin[imn]].wave
              end
          endcase
      endif else begin
          print, 'Choose a line'
          case state.flg_lines of
              1: gdwave = x_slctline(lines, /ISM) 
              2: gdwave = x_slctline(lines, /GAL) 
              3: gdwave = x_slctline(lines, /GAL) 
              4: gdwave = x_slctline(lines, /ISM) 
              7: gdwave = x_slctline(lines, /GAL) 
             10: gdwave = x_slctline(lines, /GAL) 
              else: gdwave = x_slctline(lines)
          endcase
      endelse
                  
      ; Calc AODM Colm
      x_aodm, state.wave[pxmin:pxmax], (state.fx[pxmin:pxmax]/cntm), $
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

pro x_specplot_PltClm, state

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
          print, 'Colm = '+strtrim(state.Colm,2)+$
            ' Err = '+strtrim(state.sig_colm,2)
      end
      else :
  endcase
  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate S/N

pro x_specplot_SN, state

common x_specplot_lines

  ;; Pixels to smooth over
  wave= xgetx_plt(state.xcurs,state.pos,state.xymnx, state.size) 
  mn = min(abs(state.wave-wave),imn)
  gdp = 0 > (imn+lindgen(51)-25) < (n_elements(state.wave)-1)

  ;; Noise
  noise = median(state.sig[gdp])

  ;; Signal
  signal = xgety_plt(state.ycurs,state.pos,state.xymnx,state.size) 

  print, 'x_specplot: S/N = ', signal/noise

return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Plot Guess of lines

pro x_specplot_guess, state, val, IMG=img

  common x_specplot_lines

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
;  ypt = state.xymnx[1] + 0.1*scrn
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ; Value
  case val of
      'L': begin  ; LLS
          zlls = xpt/911.8 - 1
          oplot, [xpt,xpt],  [-1e20, 1e20],  color=pclr
          oplot, [xpt*1215.6701/914.039,xpt*1215.6701/914.039], $
            [-1e20, 1e20],  color=pclr
          oplot, [xpt*1548.195/914.039,xpt*1548.195/914.039], $
            [-1e20, 1e20],  color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'LLS', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim(zlls,2)+' ]'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim(zlls,2)+' ]', $
            /normal, color=pclr, charsize=2.
          ;; Set zabs etc.
          state.flg_lines = 4
          x_specplot_initLines, state
          zabs = zlls
          flg_lines = 1
          widget_control, state.zabs_id, set_value=zlls
          widget_control, state.lines_id, set_list_select=state.flg_lines-1
      end
      'A': begin  ; AlIII
          oplot, [ xpt*1854.7164/1862.7895, $
                   xpt, $
                   xpt*1862.7895/1854.7164 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'Al III', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1854.7164 - 1),2)+', '+$
            strtrim((xpt/1862.7164 -1),2)+']'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim((xpt/1854.7164 - 1),2)+', '+$
            strtrim((xpt/1862.7164 -1),2)+']', $
            /normal, color=pclr, charsize=2.
      end
      'C': begin  ; CIV
          oplot, [ xpt*1548.195/1550.770, $
                   xpt, $
                   xpt*1550.770/1548.195 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'C IV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1548.195 - 1),2)+', '+$
            strtrim((xpt/1550.770 -1),2)+']'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim((xpt/1548.195 - 1),2)+', '+$
            strtrim((xpt/1550.770 -1),2)+']',  /normal, color=pclr, charsize=2.
      end
      'V': begin  ; SiIV
          oplot, [ xpt*1393.755/1402.770, $
                   xpt, $
                   xpt*1402.770/1393.755 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'SiIV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1393.755 - 1),2)+', '+$
            strtrim((xpt/1402.770 -1),2)+']'
      end
      'O': begin  ; OVI
          oplot, [ xpt*1031.9261/1037.6167,$
                   xpt, $
                   xpt*1037.6167/1031.9261], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'OIV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1031.9261 - 1),2)+', '+$
            strtrim((xpt/1037.6167 -1),2)+']'
      end
      'M': begin  ; MgII
          oplot, [ xpt*2796.352/2803.531, $
                   xpt, $
                   xpt*2803.531/2796.352 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'MgII', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/2796.352 - 1),2)+', '+$
            strtrim((xpt/2803.531 -1),2)+']'
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
          xyouts, 0.5, 0.9, 'zgal = [ '+strtrim((xpt/3934.79 - 1),2)+']', /normal, color=clr.red, charsize=2.5
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

pro x_specplot_psfile, state

; Device
  device, get_decomposed=svdecomp

  x_psopen, 'idl.ps', /maxs
  state.psfile = 1
  xspecplot_UpdatePlot, state
  x_psclose
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  SMOOTH
;;;;;;;;;;;;;;;;;;;;

pro x_specplot_smooth, state, RESET = reset, YTWO = YTWO

  if not keyword_set(RESET) then begin
      IF KEYWORD_SET(YTWO) THEN BEGIN
          state.flg_smooth2 = state.flg_smooth2 + 1
                                ; Smooth
          state.smooth2 = smooth(state.YTWO, 2*state.flg_smooth2+1, /NAN)
      ENDIF ELSE BEGIN
                                ; Flag
          state.flg_smooth = state.flg_smooth + 1
                                ; Smooth
          state.smooth = smooth(state.fx, 2*state.flg_smooth+1, /NAN)
      ENDELSE 
  endif else begin
      IF KEYWORD_SET(YTWO) THEN BEGIN
                                ; Flag
          state.flg_smooth2 = 0
      ; Smooth
          state.smooth2 = state.ytwo
ENDIF ELSE BEGIN
      ; Flag
          state.flg_smooth = 0
      ; Smooth
          state.smooth = state.fx
      ENDELSE
  endelse
end


;;;;;;;;;;;;;;;;;;;;
;  Init Lines
;;;;;;;;;;;;;;;;;;;;

pro x_specplot_initLines, state

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
      9: begin ;telluric lines HIRES (MF)
         llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/tell.lst'
         x_specplot_initTEL, llist
      end    
      10: begin                 ;telluric lines HIRES (MF)
         llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qso2.lst'
         x_specplot_initGAL, llist
      end
      11: begin ; stars, taken from gal.lst (so vacuum wavelengths)
         llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/star.lst'
         x_specplot_initGAL, llist
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

pro x_specplot_initdla, state, BETA=beta

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
  x_specplot_initLines, state
  flg_lines = 1
  zabs = state.dla_line.zabs
  widget_control, state.zabs_id, set_value=zabs
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Init LLS
;;;;;;;;;;;;;;;;;;;;

pro x_specplot_initslls, state, ESI=esi

  common x_specplot_lines

  ;; Grab x,y pos
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ;; Setup HI
  tmp = x_setline(1215.670d)

  ;; Get z
  state.dla_line = tmp
  state.dla_line.zabs = (xpt / 1215.6701) - 1.
  if keyword_set(ESI) then state.dla_line.N = 19.3 $
  else state.dla_line.N = 19.0
  state.dla_line.b = 30.0

  ;; Calculate
  xmin = state.wave[0] > (xpt - 20.*(1+state.dla_line.zabs))
  xmax = state.wave[state.npix-1] < (xpt + 20.*(1+state.dla_line.zabs))

  mn = min(abs(state.wave-xmin), imn)
  mx = min(abs(state.wave-xmax), imx)

  state.dla_fx = 1.
  state.dla_fx[imn:imx] = x_voigt(state.wave[imn:imx], state.dla_line, $
                                  FWHM=state.FWHM) 
  state.dla_fx = state.dla_fx * ypt
  state.flg_DLA = 1

  ;; Overplot
  state.flg_lines = 4
  x_specplot_initLines, state
  flg_lines = 1
  zabs = state.dla_line.zabs
  widget_control, state.zabs_id, set_value=zabs
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Plot DLA
;;;;;;;;;;;;;;;;;;;;

pro x_specplot_PltDLA, state

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

pro x_specplot, flux, ysin, XSIZE=xsize, YSIZE=ysize, TITLE=title $
                , WAVE = wave, QAL = QAL, GAL = gal, INFLG = inflg $
                , DUNIT = dunit, BLOCK=block, ZIN=zin, NIR = NIR $
                , NRM=nrm, QSO=qso, GRB=grb, LBG=lbg $
                , YTWO = ytwo, TWO_WAVE = WAVE_TWO, PSYM2 = PSYM2, ZOUT = zout $
                , YTHREE = YTHREE, THREE_WAVE = WAVE_THREE $
                , YFOUR = YFOUR, FOUR_WAVE = WAVE_FOUR $
                , PSYM3 = PSYM3, PSYM4 = PSYM4 $
                , TWO_IFLG = two_iflg, FIL_GALTEMP = galtempl_fil $
                , LLS = lls, AIR = air, YRANGE = YRANGE, XRANGE = XRANGE $
                , AUTO = auto, GUI = gui, ASYM_SIG = asym_sig, DISP = disp $
                , T2QSO = T2QSO, STAR = STAR, LINID_FIL = LINID_FIL $
                , XOFFSET = xoffset, YOFFSET = yoffset, PSN=PSN $
                , D2UNIT=d2unit, _EXTRA = EXTRA


common x_specplot_lines

;
  if  N_params() LT 1 and not keyword_set(GUI) then begin 
    print,'Syntax - ' + $
             'x_specplot, fx, ys_in XSIZE=,YSIZE=, TITLE=, WAVE=, DUNIT='
    print, '            /QAL, /GAL, /QSO, ZIN=, /BLOCK, /NRM, /LLS, /QSO, /AIR, /AUTO, /NIR, YTWO=, /ASYM_SIG ) [v1.2]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( TWO_IFLG ) and keyword_set(INFLG) then two_iflg = inflg
  if not keyword_set( DUNIT) then dunit=0L
  
  device, get_screen_size=ssz
  ;if szz GT 
;  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
       endif
      if abs(ssz[1]/float(ssz[0]) - 1.6) lt 1e-3 then begin ; in case of portrait monitors
         ;; assume at least two
         ;; Force aspect ratio 1.6
         ssz[0] = 2*ssz[0] - 200 ; and another 200 taken off
         ssz[1] = ssz[0]/1.6
      endif 
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0

  if keyword_set(GUI) then begin
      files =findfile('*.fits*', count=nfil)
      if nfil EQ 0 then begin
          print, 'x_specplot: No fits files in the directory. Returning..'
          return
      endif
      flux = x_guilist(files,'FITS files')
  endif

  if size(flux,/type) eq 7 then begin
     ydat = x_readspec(flux, dunit, INFLG=inflg, head=head, NPIX=npix, $
                       WAV=xdat, FIL_SIG=ysin, SIG=ysig, FLG_FIL=flg_fil, $
                       AUTO=auto, ASIG1=asig1, ASIG2=asig2, _EXTRA=EXTRA)
     if flg_fil EQ 0L then return
  endif else begin
     ;; Allow flux and ysin to be arrays (KLC)
     if keyword_set(ASYM_SIG) then stop ;; Read it in somehow
     ydat = flux
     if keyword_set(ysin) then ysig = ysin
     npix = n_elements(ydat)
  endelse

  if not keyword_set(YSIG) then ysig = fltarr(npix)
  if n_elements(ydat) EQ 1 AND ydat[0] EQ -1 then begin
      print, 'x_specplot: Returning...'
      help, flux, ydat
      return
  endif

  ;;IF PSN plot the signal to noise instead
  if keyword_set(PSN) then begin
     ;;find non zeros
     nzr=where(ysig gt 0)
     ydat[nzr] = ydat[nzr]/ysig[nzr]
     ysig = ysig*0.0
  endif


  ;; WAVE
  if keyword_set(WAVE) then xdat = wave
  If keyword_set(ytwo) THEN BEGIN
      IF NOT keyword_set(d2unit) then d2unit=0L
      IF size(ytwo, /type) EQ 7 THEN $
        ydat2 = x_readspec(ytwo, d2unit, INFLG = two_iflg, head = head2, NPIX = npix2 $
                           , WAV = wave_two, FIL_SIG = ysin, SIG = ysig2 $
                           , FLG_FIL = flg_fil, _EXTRA=EXTRA) $
      ELSE ydat2 = ytwo
      IF keyword_set(wave_two) THEN xdat2 = wave_two ELSE xdat2 = xdat
      flg_ytwo = 1
  ENDIF ELSE BEGIN
      ydat2 = fltarr(npix)
      xdat2 = fltarr(npix)
      flg_ytwo = 0
   ENDELSE

  If keyword_set(ythree) THEN BEGIN
     IF size(ythree, /type) EQ 7 THEN $
        ydat3 = x_readspec(ythree, INFLG = inflg, head = head2, NPIX = npix3 $
                           , WAV = wave_three, FIL_SIG = ysin, SIG = ysig3 $
                           , FLG_FIL = flg_fil, _EXTRA=EXTRA) $
     ELSE ydat3 = ythree
     IF keyword_set(wave_three) THEN xdat3 = wave_three ELSE xdat3 = xdat
     flg_ythree = 1
  ENDIF ELSE BEGIN
     ydat3 = fltarr(npix)
     xdat3 = fltarr(npix)
     flg_ythree = 0
  ENDELSE

  If keyword_set(yfour) THEN BEGIN
     IF size(yfour, /type) EQ 7 THEN $
        ydat4 = x_readspec(yfour, INFLG = inflg, head = head2, NPIX = npix3 $
                           , WAV = wave_four, FIL_SIG = ysin, SIG = ysig3 $
                           , FLG_FIL = flg_fil, _EXTRA=EXTRA) $
     ELSE ydat4 = yfour
     IF keyword_set(wave_four) THEN xdat4 = wave_four ELSE xdat4 = xdat
     flg_yfour = 1
  ENDIF ELSE BEGIN
     ydat4 = fltarr(npix)
     xdat4 = fltarr(npix)
     flg_yfour = 0
  ENDELSE

  
  ;; Symbols
  if not keyword_set(psym2) then psym2 = 10 $
  else if psym2 eq -3 then psym2 = 0 ; instead of dots with lines, just smooth line
  if not keyword_set(psym3) then psym3 = 10 $
  else if psym3 eq -3 then psym3 = 0
  if not keyword_set(psym4) then psym4 = 10 $
  else if psym4 eq -3 then psym4 = 0


  if not keyword_set(ASIG1) then asig1 = fltarr(npix)
  if not keyword_set(ASIG2) then asig2 = fltarr(npix)
  
;  tmp1 = { abslinstrct }
  tmp1 = { newabslinstrct }

; Init common

  x_specplot_initcommon
  
  IF KEYWORD_SET(XRANGE) THEN BEGIN
      XMIN = XRANGE[0]
      XMAX = XRANGE[1]
  ENDIF ELSE BEGIN
      xmin = min(xdat,max=xmax)
  ENDELSE 
  IF KEYWORD_SET(YRANGE) THEN BEGIN
      YMIN = YRANGE[0]
      YMAX = YRANGE[1]
  ENDIF ELSE BEGIN
      gd = where(xdat GE XMIN AND xdat LE XMAX, ngd)
      IF ngd EQ 0 THEN gd = lindgen(n_elements(ydat))
      ymin = min(djs_median(ydat[gd], width = 5, boundary = 'reflect')) $
        - 0.01*abs(max(ydat[gd])-min(ydat[gd]))
      ymax = max(djs_median(ydat[gd], width = 5, boundary = 'reflect')) $
        +0.01*abs(max(ydat[gd])-min(ydat[gd]))
   ENDELSE

  ;; Dispersion
  if keyword_set(DISP) then begin
     dwv = xdat - shift(xdat,1)
     med_disp = abs(median(dwv/xdat)) * 3e5
     print, 'x_specplot: Median dispersion ', med_disp, 'km /s'
  endif
     
;    STATE
;state = { fx: ((ydat >
   ;(-1.0d7)) <  1.0d7), $  -- JXP :: Don't do this again!!
  state = { fx: ydat, $
            wave: xdat, $
            wave_two: xdat2, $
            wave_three: xdat3, $
            wave_four: xdat4, $
            sig: ysig, $
            asig1: asig1, $
            asig2: asig2, $
            ytwo: ((ydat2 > (-1.0d7)) <  1.0d7), $
            flg_ytwo: flg_ytwo, $
            ythree: ((ydat3 > (-1.0d7)) <  1.0d7), $
            flg_ythree: flg_ythree, $
            yfour: ((ydat4 > (-1.0d7)) <  1.0d7), $
            flg_yfour: flg_yfour, $
            npix: npix, $
            flg_smooth: 0, $    ; Smoothing
            flg_smooth2:0, $
            smooth: fltarr(n_elements(ydat)), $
            smooth2: fltarr(n_elements(ydat2)), $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_EW: 0, $ ; EW flag
            EW_lmt: dblarr(2,2), $
            EW: 0.d, $    
            flg_GS: 0, $ ; Gaussian stuff
            GS_lmt: dblarr(2,2), $
            GS_fit: fltarr(1000L), $
            GS_xpix: lonarr(2), $
            Gauss: fltarr(4), $    
            sigEW: 0.d, $    
            FWHM: 4., $   ; FWHM of instrument (pix)
            flg_DLA: 0, $ ; DLA flag
            flg_qsotempl: 0, $ ; QSO template flag
            flg_atm: 0, $ ; IR atmospheric transmission template flag
            flg_galtempl: 0, $  ; Gal template flag
            galtempl_fil: '', $ ; Gal template file
            qsot_nrm: fltarr(2), $ ; QSO template flag
            atm_nrm: fltarr(2), $ ; atmospheric transmission flag
            galt_nrm: fltarr(2), $ ; Gal template flag
            dla_line: tmp1, $ ; DLA flag
            dla_fx: fltarr(n_elements(ydat)) + 1., $
            flg_Colm: 0, $ ; Colm flag
            Colm_lmt: dblarr(2), $
            Colm: 0., $
            sig_Colm: 0., $
            xpos: 0.d, $
            ypos: 0.d, $
            psfile: 0, $ ; Postscript
            flg_lines: 0, $  ; QAL=1, GAL=2; LLS=4, LBG=7, star=11
            flg_linid: keyword_set(linid_fil), $ ; CSV of zabs,wrest,label
            linid_fil:'', $
            linid_mtime:long64(0), $             ; (file_info(linid_fil)).mtime
            svxymnx: [xmin, ymin, xmax, ymax], $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            psym: 10, $
            psym2: psym2, $
            psym3: psym3, $
            psym4: psym4, $
            title: '', $
            help: strarr(50), $
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
            lines_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }


  
; NRM
  if keyword_set(NRM) then begin
      state.svxymnx[1] = -0.09
      state.svxymnx[3] = 1.3
  endif
  if keyword_set(linid_fil) then begin
     state.linid_fil = linid_fil
     x_specplot_initlinid, state ; load linid structure in common block
  endif 

; WAVE
  if keyword_set(WAVE) then state.wave = wave

; YTWO
;  if keyword_set(YTWO) then begin
;      IF KEYWORD_SET(WAVE_TWO) THEN state.wave_two = wave_two $
;      ELSE IF KEYWORD_SET(WAVE) THEN state.WAVE_TWO = wave 
;      state.ytwo = ytwo
;      state.flg_ytwo = 1
;  endif
              
  if keyword_set( YSIG ) then state.flg_sig = 1
  if keyword_set(ASYM_SIG) then begin
     state.flg_sig = 2
     state.psym = 2
  endif

; LINELIST

  if keyword_set( QAL ) then state.flg_lines = 1
  if keyword_set( GAL ) then state.flg_lines = 2
  if keyword_set( QSO ) then state.flg_lines = 3
  if keyword_set( LLS ) then state.flg_lines = 4
  if keyword_set( GRB ) then state.flg_lines = 6
  if keyword_set( LBG ) then state.flg_lines = 7
  if keyword_set( MOL ) then state.flg_lines = 8
  if keyword_set( TEL ) then state.flg_lines = 9
  if keyword_set( T2QSO ) then state.flg_lines = 10
  if keyword_set( STAR ) then state.flg_lines = 11
  if keyword_set( ZIN ) and state.flg_lines EQ 0 then state.flg_lines = 1
  
;; 
  if keyword_set( NIR ) THEN BEGIN
     state.atm_nrm = [state.xpos, 0.95*ymax]
     state.flg_atm = 1
  ENDIF

  

  if state.flg_lines NE 0 then x_specplot_initLines, state

  if keyword_set( ZIN ) then begin
      zabs = zin
      flg_lines = 1
   endif

  if keyword_set( GALTEMPL_FIL ) then $
     state.galtempl_fil = galtempl_fil

  ;; Air wavelengths (why!?)
  if keyword_set(AIR) then begin
      tmp = state.wave
      airtovac, tmp
      state.wave = tmp
  endif

; Set svxymnx[0,2]

  IF NOT KEYWORD_SET(XRANGE) THEN BEGIN
      state.svxymnx[0] = min(state.wave)
      state.svxymnx[2] = max(state.wave)
  ENDIF 

;    Title
  if size(flux, /type) EQ 7 then state.title = flux
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_specplot', /column, $
                    xoffset=xoffset, yoffset=yoffset)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['x_specplot', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  ii = 0L
  state.help[ii] = '  :::Help Menu::: '
  ii=ii+1
  state.help[ii] = 'LMB/LMB -- Set region'
;  ii=ii+1
;  state.help[ii] = 's/s -- Set region' ; KLC: this contradicts
;                           's -- ytwo smooth'
  ii=ii+1
  state.help[ii] = '  -- Print x, y'
  ii=ii+1
  state.help[ii] = 'l -- Set Left '
  ii=ii+1
  state.help[ii] = 'r -- Set Right '
  ii=ii+1
  state.help[ii] = 'b -- Set Bottom '
  ii=ii+1
  state.help[ii] = 't -- Set Top '
  ii=ii+1
  state.help[ii] = 'Z -- Set ymin to 0.'
;  state.help[ii] = 'Z -- Set redshift by hand'
  ii=ii+1
  state.help[ii] = 'T -- Set ymax to 1.1'
  ii=ii+1
  state.help[ii] = 'w -- Reset the screen'
  ii=ii+1
  state.help[ii] = 'L -- Set redshift with a line'
  ii=ii+1
  state.help[ii] = 'i -- Zoom in'
  ii=ii+1
  state.help[ii] = 'o -- Zoom out'
  ii=ii+1
  state.help[ii] = 'z -- Zoom region'
;  state.help[ii] = 'z -- Set ymin to 0.'
  ii=ii+1
  state.help[ii] = '[ -- Pan left'
  ii=ii+1
  state.help[ii] = '] -- Pan right'
  ii=ii+1
  state.help[ii] = '{ -- Pan left, no y-rescale'
  ii=ii+1
  state.help[ii] = '} -- Pan right, no y-rescale'
  ii=ii+1
  state.help[ii] = 'H -- Show this screen'
  ii=ii+1
  state.help[ii] = 'A -- Al III doublet'
  ii=ii+1
  state.help[ii] = 'C -- C IV doublet'
  ii=ii+1
  state.help[ii] = 'M -- Mg II doublet'
  ii=ii+1
  state.help[ii] = 'V -- Si IV doublet'
  ii=ii+1
  state.help[ii] = 'O -- O VI doublet'
  ii=ii+1
  state.help[ii] = 'L -- LLS'
  ii=ii+1
  state.help[ii] = 'K -- Ca H+K doublet'
  ii=ii+1
  state.help[ii] = 'D -- DLA overplot'
  ii=ii+1
  state.help[ii] = 'B -- DLB overplot'
  ii=ii+1
  state.help[ii] = 'R -- sLLS overplot'
  ii=ii+1
  state.help[ii] = 'Q -- QSO overplot'
  ii=ii+1
  state.help[ii] = 'I -- Atomsphere'
  ii=ii+1
  state.help[ii] = 'g -- Galaxy overplot'
  ii=ii+1
  state.help[ii] = 'E/E -- EW measurement'
  ii=ii+1
  state.help[ii] = 'G/G -- fit Gaussian'
  ii=ii+1
  state.help[ii] = 'N/N -- AODM'
  ii=ii+1
  state.help[ii] = 'S -- Smooth'
  ii=ii+1
  state.help[ii] = 'U -- UnSmooth'
  ii=ii+1
  state.help[ii] = 's -- Smooth ytwo'
  ii=ii+1
  state.help[ii] = 'u -- Unsmooth ytwo'
  ii=ii+1
  state.help[ii] = 'P -- Print PS'
  ii=ii+1
  state.help[ii] = 'q -- Quit '
;  print,'Number of help strings',ii ; KLC: 40 currently, v1.70 29 Apr 2010

;;;;;;;;;
;  Toolbar

; zabs

  state.zabs_id = cw_field(toolbar, title="zabs", value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='ZABS')


; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Objname

  if keyword_set( head ) then begin
      objnm = strtrim(sxpar(head, 'TARGNAME'),2)
      objnm_id = cw_field(toolbar, title='Object: ', value=objnm, /column, $
                         xsize=strlen(objnm))
  endif
  
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
                               VALUE = ['DLA', 'GAL', 'QSO', 'LLS', 'H2', 'GRB', 'LBG', 'MOL', 'TEL', 'T2QSO','STAR'], $
                             uvalue='LNLIST', ysize = 4)
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

;      Done
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Update
  xspecplot_Reset, state
  xspecplot_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'x_specplot', base, /no_block $
  else xmanager, 'x_specplot', base
  zout = zabs

return
end
