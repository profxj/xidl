;+ 
; NAME:
; sdss_llsfit
;   Version 1.11
;
; PURPOSE:
;    GUI used to fit DLA profiles interactively.  The user determines
;    the continuum at the same time.
;
; CALLING SEQUENCE:
;  sdss_llsfit, flux_fil, [err_fil], XSIZE=, YSIZE=, TITLE=, 
;               WAVE=, /BLOCK, FWHM=, INIFIT=, INFLG=
;
; INPUTS:
;   flux_fil   - FITS file (or array) containing flux
;   [err_fil]  - FITS error array
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array 
;   FWHM=      - Resoltuion (FWHM in pixels)  [default: 2]
;   INFLG=     - Usual flag for reading the FITS file (see x_readspec)
;   INIFIT=    - IDL file containing a saved version of the fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_llsfit, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_llsfit_icmmn

  common sdss_llsfit_cmm, $
    llssearch, $
    template, $
    allnh

  nnhi = 20
  nh1 = xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',0)
  npnh = n_elements(nh1)
  allnh = fltarr(npnh,nnhi)
  for i=0L,nnhi-1 do allnh[*,i] = xmrdfits(getenv('XIDL_DIR')+ $
                                           '/SDSS/LLS/nhi16_19b30.fits',i)

  return
end


;;;;
; Events
;;;;

pro sdss_llsfit_event, ev

  common sdss_llsfit_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval
  flg_plt = 0

  case uval of
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          llssearch.zlls = tmp
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
                      1 : begin     ; Center on LLS
                          llssearch.zlls = (state.xpos / 911.7633) - 1.
                          widget_control, state.zabs_id, $
                            set_value=llssearch.zlls
                      end 
                      4 : begin  ;; Center on Lya
                          llssearch.zlls = (state.xpos / 1215.67) - 1.
                          widget_control, state.zabs_id, $
                            set_value=llssearch.zlls
                      end
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
          if eventch EQ  'q' then begin
              widget_control, ev.top, /destroy
              return
          end
          case eventch of
              'b': state.xymnx[1] = state.ypos
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = state.xpos
              'r': state.xymnx[2] = state.xpos
              't': state.xymnx[3] = state.ypos
              'T': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              ;; ZOOMING
              'i': x_speczoom, state, 0 ; Zoom in
              'o': x_speczoom, state, 1 ; Zoom out
              'z': begin
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      flg_plt = 0
                      return
                  endif
              end
              ;; Smoothing
              'S': state.smooth = state.smooth + 1
              'R': state.smooth = 1
              ;; Reset
              'w': state.xymnx = state.svxymnx ; Reset the screen
              ;; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              ;; HELP
              'H': begin 
                  x_helpwidg, state.help
                  flg_plt = 0
              end
              ;; Recenter
              'L': begin
                  llssearch.zlls = state.xpos/911.7633 - 1
                  widget_control, state.zabs_id, set_value=llssearch.zlls
              end
              'A': begin
                  llssearch.zlls = state.xpos/1215.6701 - 1
                  widget_control, state.zabs_id, set_value=llssearch.zlls
              end
              ;; Continuum
              'C': begin ;; Normalize to the cursor
                  mn = min(abs(state.wave-state.xpos),imn)
                  state.scale = state.ypos / state.conti[imn]
                  llssearch.conti_scale = state.scale
;                  state.conti = state.conti * scl
              end
              ;; Colm
              'n': begin
                  llssearch.taulls = (llssearch.taulls - 1) > 0
                  widget_control, state.Ncolm_id, $
                                  set_value=16+llssearch.taulls*0.2
              end
              'N': begin
                  llssearch.taulls = (llssearch.taulls + 1) < 19
                  widget_control, state.Ncolm_id, $
                                  set_value=16+llssearch.taulls*0.2
              end
              ;; Colm
              'm': begin
                  llssearch.sigtau = (llssearch.sigtau - 0.2) > 0.2
                  widget_control, state.sigN_id, set_value=llssearch.sigtau
              end
              'M': begin
                  llssearch.sigtau = llssearch.sigtau + 0.2
                  widget_control, state.sigN_id, set_value=llssearch.sigtau
              end
              ;; Output
              'P': x_fitline_psfile, state  
                                ; QUIT
              else:  print, 'x_fitline: Not a valid key!' ; Nothing
          endcase
      end
      ;; BUTTONS
      'IDLOUT' : x_fitline_idlout, state, /ERR
      'SPECPLT' : begin
          zin = llssearch.zlls
          x_specplot, state.fx, state.sig, wave=state.wave, inflg=4, $
            zin=zin, /lls, /block, ZOUT=zout
;          llssearch.zlls = zout
;          widget_control, state.zabs_id, set_value=llssearch.zlls
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; Update Plot
  sdss_llsfit_updfit, state
  sdss_llsfit_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_llsfit_UpdatePlot, state
  
  common sdss_llsfit_cmm
; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

;  if state.flg_dum EQ 1 then stop
;  state.flg_dum = 0
;  stop
  plot, state.wave, $
        smooth(state.fx,2*state.smooth+1), $
        psym=state.psym, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        xtitle='!17Wavelength', ytitle='Flux', $
        title=state.title, $
        background=clr.white, $
        color=clr.black, $
        xcharsize=1.7, $
        ycharsize=1.7, /nodat


  ;; Error interval
  oplot, [-1e9,1e9], [0., 0], color=clr.green, linestyl=2, thick=2

  ;; Error interval
  x_curvefill, state.wave, state.scale*state.conti*state.low_fit, $
               state.scale*state.conti*state.hi_fit, color=clr.cyan

  ; Plot Error array
  oplot, state.wave, state.sig, psym=state.psym, color=clr.orange

  ;; Data
  oplot, state.wave, smooth(state.fx,2*state.smooth+1), psym=state.psym,$
        color=clr.black

  ;; Plot smoothed points
  npts = round(state.npix/float(state.nbin)) - 1
  for ii=0L,npts-1 do begin
      p1 = ii*state.nbin
      p2 = (ii+1)*state.nbin - 1
      med_sig = median(state.sig[p1:p2])
      mean_fx = mean(state.fx[p1:p2])
      ;; Plot
      mnwv = mean(state.wave[p1:p2])
      oploterror, mnwv, mean_fx, mnwv-state.wave[p1], $
                  med_sig/sqrt(state.nbin), $
                  color=clr.red, errcolor=clr.red, thick=3
  endfor

  ;; FIT
  oplot, state.wave, state.scale*state.conti*state.fit, color=clr.blue, thick=3

; CONTINUUM
  ;; Line
  oplot, state.wave, state.scale*state.conti, color=clr.purple, linestyle=1

  ;; Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.orange, symsize=5
  endif

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro sdss_llsfit_Reset, state


  ;; Plotting
  state.xymnx = state.svxymnx
  state.xymnx[2] = 1220. * (1+state.speczem)
  state.xymnx[1] = -2*median(state.sig)
  state.old_xymnx = state.xymnx

end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro sdss_llsfit_updfit, state

  common sdss_llsfit_cmm

  ;; NHI absorption (rest-frame)
  model = allnh[*,llssearch.taulls]

  npnh = n_elements(model)
  wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
  mn = min(abs(wv_mod-911.7633),mmn)

  mn = min(abs((1+llssearch.zlls)*911.7633-state.wave),imn)  

  i1 = mmn-imn
  i2 = (npnh-1) < (i1+n_elements(state.fit)-1) 
  state.fit = model[i1:i2]

  ;; Sigma low/hi
  model = allnh[*,(llssearch.taulls-round(llssearch.sigtau/0.2))>0]
  state.low_fit = model[i1:i2]

  model = allnh[*,(llssearch.taulls+round(llssearch.sigtau/0.2))<19]
  state.hi_fit = model[i1:i2]

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

pro sdss_llsfit, plt_fib, i_llssearch, i_template, XSIZE=xsize, YSIZE=ysize, $
                 BLOCK=block, INIFIT=inifit, DR5=dr5, OUTSTR=outstr, $
                 SEARCH_FIL=search_fil, TMPL_FIL=tmpl_fil, INDIV=indiv, $
                 MOCK_DATA=mock_data

  common sdss_llsfit_cmm

  if  N_params() LT 3  then begin 
      if ((not keyword_set(SEARCH_FIL)) or (not keyword_set(TMPL_FIL))) $
        AND not keyword_set(INDIV) then begin
          print,'Syntax - ' + $
                'sdss_llsfit, llssearch, lls_fstrct, template, XSIZE=,YSIZE=, '
          print, '            INIFIT=, INFLG=, FWHM=, /BLOCK, ZIN=, /BETA) [v1.1]'
          return
      endif else begin
          if keyword_set(INDIV) then begin
              ;; LLS search
              i_llssearch = {sdssllsstrct}
              ;; Grab QSO redshift
              sdss_objinf, plt_fib, zem=zem
              if not keyword_set(ZEM) then return
              zem = fix(zem*10)
              ;; Grab the template file
              temp_files = findfile(getenv('SDSSPATH')+ $
                                    '/DR7_QSO/LLS/hiztempl*fits', count=nfil)
              if nfil EQ 0 then return
              flg_tmpl = 0
              for qq=0L,nfil-1 do begin
                  ;; Parse
                  prs = strsplit(temp_files[qq],'_',/extrac)
                  np = n_elements(prs)
                  z2 = long(strmid(prs[np-1],0,2))
                  z1 = long(prs[np-2])
                  if zem LT z2 and zem GE z1 then tfil = temp_files[qq]
              endfor
              if not keyword_set(TFIL) then return
              i_template = xmrdfits(tfil,/silen)
          endif else begin
              i_template = xmrdfits(tmpl_fil,/silen)
              search = xmrdfits(search_fil, 1, /silen)
              mtch = where(search.plate EQ plt_fib[0] and $
                           search.fiber EQ plt_fib[1],nmtch)
              if nmtch NE 1 then begin
                  print, 'Object not in the search file.  Please try again'
                  return
              endif
              i_llssearch = search[mtch]
          endelse
      endelse
  endif 

;  Optional Keywords

  ;; common
  sdss_llsfit_icmmn
  llssearch = i_llssearch
  template = i_template

  ;; Error
  if llssearch.sigtau LT 1e-5 then llssearch.sigtau = 0.2
  if llssearch.conti_scale LT 1e-5 then llssearch.conti_scale = 1.

  ;; Screen
  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( FWHM ) then    fwhm = 2.

  ;; Read in the Data
  if not keyword_set(MOCK_DATA) then begin
      sdss_objinf, plt_fib, DR5=dr5, filnm=filnm
      parse_sdss, filnm, fx, wave, sig=sig, npix=npix, HEAD=hdr

      ;; Generate the continuum
      specwav0 = sxpar(hdr,"COEFF0")
      speczem = sxpar(hdr,"Z")
  endif else begin
      fx = mock_data.spec
      sig = mock_data.sig
      npix = n_elements(fx)
      bad = where(sig GT 999.,nbad)
      if nbad GT 0 then sig[bad] = 0.
      specwav0 = mock_data.w0
      speczem = mock_data.zem
  endelse
  logwav = findgen(npix)*0.0001 + specwav0
  wave = 10.^logwav
  sub_conti = sdss_qsotempl(template, 2.8, logwav, $
                            fx, sig, speczem)

  conti = fltarr(npix)
  conti[0:n_elements(sub_conti)-1] = sub_conti
  conti[n_elements(sub_conti):*] = conti[n_elements(sub_conti)-1]
  

  tmp2 = { conti_str, $
           npts: 0L, $
           xval: fltarr(100), $
           yval: fltarr(100), $
           msk: lonarr(100) }

;    STATE
  state = { fx: fx, $
            speczem: speczem, $
            wave: wave, $
            sig: sig, $
            npix: npix, $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            fit: fltarr(npix) + 1., $
            low_fit: fltarr(npix) + 1., $
            hi_fit: fltarr(npix) + 1., $
            conti: conti, $
            scale: llssearch.conti_scale, $
            cstr: tmp2, $
            curlin: 0L, $
            smooth: 0, $
            nbin: 100L, $
            nset: -1, $
            xpos: 0.d, $
            ypos: 0.d, $
            flg_dum: 1, $
            psfile: 0, $ ; Postscript
            svxymnx: [0., $     ; xmin
                      min(fx)-0.01*abs(max(fx)-min(fx)), $ ; ymin
                      float(n_elements(fx)-1), $ ; xmax
                      max(fx)+0.01*abs(max(fx)-min(fx))], $ ; ymax
            xymnx: fltarr(4), $
            old_xymnx:fltarr(4), $
            tmpxy: fltarr(4), $
            xpmnx: lonarr(2), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
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
            bval_id: 0L, $
            Ncolm_id: 0L, $
            sigN_id: 0L, $
            lines_id: 0L, $
            crude_id: 0L, $
            exact_id: 0L, $
            crudehi_id: 0L, $
            crudelo_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }

; WAVE
  if keyword_set(BETA) then begin
      state.wrest = 1025.7223d
      state.beta = 1
  endif
  if keyword_set(WAVE) then state.wave = wave
  if keyword_set( YSIG ) then state.flg_sig = 1

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE
; Set svxymnx[0,2]

  state.svxymnx[0] = min(state.wave)
  state.svxymnx[2] = max(state.wave)

;    Title
  if size(yin, /type) EQ 7 then state.title = yin
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_llsfit', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['sdss_llsfit', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'T -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = 'A -- Al III doublet'
  state.help[19] = 'C -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'N -- AODM'
  state.help[22] = 'q -- Quit '


;;;;;;;;;
;  Toolbar

; zabs

  linbar = WIDGET_BASE( toolbar, /row, /frame, /base_align_center,$
                           /align_center)
  state.zabs_id = cw_field(linbar, title='zabs', value=llssearch.zlls, $
                           /floating, $
                           /column, /return_events, xsize=10, uvalue='ZABS')
  state.Ncolm_id = cw_field(linbar, title='Ncolm', $
                            value=(16+llssearch.taulls*0.2), /floating, $
                            /column, xsize=10, /return_events, uvalue='NCOLM')
  state.sigN_id = cw_field(toolbar, title='sigN', value=llssearch.sigtau, /floating, $
                           /column, xsize=12, /return_events, uvalue='SIGN')


; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Crude error
;  crude_err = widget_base(toolbar, /row, /align_center, frame=2)
;  state.crude_id = CW_BGROUP(crude_err, ['No', 'Yes'], $
;                              row=2, /exclusive, /no_release, $
;                              set_value=0,  uvalue='CRUDE')
;  crude_val = widget_base(crude_err, /column, /align_center, frame=2)
;  state.crudehi_id = CW_FIELD(crude_val, value=state.crude_val[1], xsize=5, $
;                              title='Crude High', /floating, $
;                              /return_events, uvalue='CRUDEHI')
;  state.crudelo_id = CW_FIELD(crude_val, value=state.crude_val[0], xsize=5, $
;                              title='Crude Low', /floating, $
;                              /return_events, uvalue='CRUDELO')
  
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

;      Done
  buttbase = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  specplt = WIDGET_BUTTON(buttbase, value='Specplt',uvalue='SPECPLT', $
                          /align_right)
;  idlout = WIDGET_BUTTON(buttbase, value='IDL Out',uvalue='IDLOUT', /align_right)
  buttbase2 = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  done = WIDGET_BUTTON(buttbase2, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
;  if not keyword_set( NOCTB ) then loadct, 2, /silent

; HELP
  i = 0
  state.help[0] = '  :::Help Menu for sdss_llsfit::: '
  state.help[i++] = 'RMB -- Recenter line'
  state.help[i++] = 's/s -- Set region'
  state.help[i++] = 'lrbt -- Set Left Right Bottom Top'
  state.help[i++] = 'T -- Set ymax to 1.1'
  state.help[i++] = 'zz -- Zoom corners'
  state.help[i++] = 'io -- Zoom in/out'
  state.help[i++] = '{}[] -- Pan'
  state.help[i++] = 'w -- Reset the screen'
  state.help[i++] = 'H -- This widget'
  state.help[i++] = 'c -- Add new Lya line'
  state.help[i++] = 'L -- Add new LLS'
  state.help[i++] = 'B -- Add new Beta'
  state.help[i++] = 'd -- Delete current line'
  state.help[i++] = 'C -- Set continnuum' 
  state.help[i++] = 'nN -- Adjust colm'
  state.help[i++] = 'vV -- Adjust bvalue'
  state.help[i++] = '=- -- Loop through lines'
  state.help[i++] = 'P -- Print to postrscipt'
  state.help[i++] = 'I -- IDL output'
  state.help[i++] = '3,4 -- Add/move continuum point'

; Update
  sdss_llsfit_Reset, state

  ;; Zin
  if keyword_set( ZIN ) then begin
      state.xpos = (1.+zin)*state.wrest
      x_fitline_newlin, state, BETA=state.beta
  endif

  ;; Init Fit
;  if keyword_set(INIFIT) then sdss_llsfit_inifit, state, inifit
  ;; Set pmnx
  state.xpmnx = x_getxpmnx(state)
  ;; Plot
  sdss_llsfit_updfit, state
  sdss_llsfit_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
;  if not keyword_set(BLOCK) then xmanager, 'sdss_llsfit', base, /no_block $
;  else xmanager, 'sdss_llsfit', base
  xmanager, 'sdss_llsfit', base

  i_llssearch = llssearch
  outstr=llssearch

return
end
