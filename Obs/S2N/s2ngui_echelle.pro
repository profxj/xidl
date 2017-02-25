;+ 
; NAME:
; hires_s2ngui
;    Version 1.0
;
; PURPOSE:
;   This launches a GUI which allows a user to fiddle around and
;   estimate S2N for HIRES observations
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
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
;   sdss_chkciv, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
pro hires_s2ngui_initcommon
;
common hires_s2ngui_cmmn, $
  s2n_fstrct
end

;;;;
; Events
;;;;

pro hires_s2ngui_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'HELP': x_helpwidg, state.help, YSIZE=40
      'REPLOT': begin
          hires_s2ngui_gets2n, state
          hires_s2ngui_updplot, state
      end
      'READNO': begin
          widget_control, state.rn_id, get_value=tmp
          state.str_instr.readno = tmp
          widget_control, state.rnd_id, set_value=tmp
      end
      'DARK': begin
          widget_control, state.dark_id, get_value=tmp
          state.str_instr.dark = tmp
          widget_control, state.darkd_id, set_value=tmp
      end
      'BINS': begin
          widget_control, state.bins_id, get_value=tmp
          state.str_instr.bins = tmp
          widget_control, state.binsd_id, set_value=tmp
      end
      'BIND': begin
          widget_control, state.bind_id, get_value=tmp
          state.str_instr.bind = tmp
          widget_control, state.bindd_id, set_value=tmp
      end
      'EXPTIME': begin
          widget_control, state.expt_id, get_value=tmp
          state.str_obs.exptime = tmp
          widget_control, state.exptd_id, set_value=tmp
      end
      'MAGNITUDE': begin
          widget_control, state.mag_id, get_value=tmp
          state.str_obs.mstar = tmp
          widget_control, state.magd_id, set_value=tmp
      end
      'SEEING': begin
          widget_control, state.seeing_id, get_value=tmp
          state.str_obs.seeing = tmp
          widget_control, state.seeingd_id, set_value=tmp
      end
      'MOON': begin
          widget_control, state.moon_id, get_value=tmp
          state.str_obs.mphase = tmp
          widget_control, state.moond_id, set_value=tmp
      end
      'AIRMASS': begin
          widget_control, state.airm_id, get_value=tmp
          state.str_obs.airmass = tmp
          widget_control, state.airmd_id, set_value=tmp
      end
      'SWIDTH': begin
          widget_control, state.swidth_id, get_value=tmp
          state.str_instr.swidth = tmp
          widget_control, state.swidthd_id, set_value=tmp
          ;; Reset the resolution
          resolu = state.str_instr.R / ( state.str_instr.swidth / $
                                         state.str_tel.plate_scale / $
                                         (state.str_instr.pixel_size/100) )
          widget_control, state.resd_id, set_value=resolu
      end
      'PLOT': begin
          state.flg_plot = ev.index
          hires_s2ngui_gets2n, state
          hires_s2ngui_updplot, state
      end
      'INSTR': begin
          state.flg_instr = ev.index + 1
          x_inithires, str_instr, state.flg_instr, STR_TEL=str_tel, $
            DECKER=state.deckers[state.deckidx]
          state.str_instr = str_instr
          state.str_tel = str_tel
          hires_s2ngui_reinit, state
      end
      'DWV': begin
          widget_control, state.dwv_id, get_value=tmp
          state.dwv = tmp
          hires_s2ngui_setwv, state
          widget_control, state.dwvd_id, set_value=tmp
      end
      'WVMN': begin
          widget_control, state.wvmn_id, get_value=tmp
          state.wvmn = tmp
          hires_s2ngui_setwv, state
          widget_control, state.wvmnd_id, set_value=tmp
      end
      'WVMX': begin
          widget_control, state.wvmx_id, get_value=tmp
          state.wvmx = tmp
          hires_s2ngui_setwv, state
          widget_control, state.wvmxd_id, set_value=tmp
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
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;

pro hires_s2ngui_setwv, state

 ;; Wave array
 state.nwv = long((state.wvmx-state.wvmn)/state.dwv) + 1
 state.wave[0:state.nwv-1] = state.wvmn + findgen(state.nwv)*state.dwv

 return
end

pro hires_s2ngui_gets2n, state

common hires_s2ngui_cmmn

  ;; Calc
  hires_calcs2n, state.wave[0:state.nwv-1], state.flg_instr, $
    STATE=state, /nopr, S2N=s2n,  IORDER=iorder, FSTRCT=s2n_fstrct
  state.s2n[0:state.nwv-1] = s2n
  state.iorder[0:state.nwv-1] = iorder

  ;; Update slit throughput
  widget_control, state.slit_id, set_value=(1.-s2n_fstrct.slit0)*100

  return
end
 
  
pro hires_s2ngui_updplot, state

common hires_s2ngui_cmmn
  ;; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

  wave= state.wave[0:state.nwv-1]
  mxx = max(wave, min=mnx)


  case state.flg_plot of
      0: begin ; S2N

          sn= state.s2n[0:state.nwv-1]
          mxy = max(sn, min=mny)
          state.xymnx = [mnx, mny, mxx, mxy]
          
          binc = state.str_instr.bind
          pixel    = binc*(3e5/state.str_instr.R)
          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.7,$
            xmargin=[10,2], ymargin=[5,2], xtitle='!17Wavelength (Ang)', $
            ytitle='S/N (per '+string(pixel,'(f4.2)')+'km/s pix)', $
            yrange=[mny*0.9,mxy*1.1], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata
          
          ;; Blue
          iorder = state.iorder[0:state.nwv-1]
          blu = where(iorder EQ 2, nblu) 
          if nblu NE 0 then oplot, [wave[blu]], [sn[blu]], $
            psym=-1, color=clr.blue
          
          ;; Red
          red = where(iorder EQ 1, nred) 
          if nred NE 0 then oplot, [wave[red]], [sn[red]], $
            psym=-2, color=clr.red
          
          ;; Label
          xpos = mnx + (mxx-mnx)*0.95
          dy = mxy-mny
      end
      1: begin ; Efficiency
          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.7,$
            xmargin=[10,2], ymargin=[5,3], xtitle='!17Wavelength (Ang)', $
            ytitle='Efficiency', $
            yrange=[0., 0.5], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata
          oplot, wave, s2n_fstrct.thru, psym=-1, color=clr.black, thick=3
      end
      2: begin ; Noise
          mxy = max([s2n_fstrct.noise^2, s2n_fstrct.star, s2n_fstrct.sky])
                    
          state.xymnx = [mnx, 0., mxx, mxy]

          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.4,$
            xmargin=[14,4], ymargin=[5,3], xtitle='!17Wavelength (Ang)', $
            ytitle='Variance', $
            yrange=[1.,mxy*1.05], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata

          ;; Sky + Obj
          oplot, wave, s2n_fstrct.star, psym=-1, color=clr.black, thick=3
          oplot, wave, s2n_fstrct.sky, psym=-1, color=clr.blue, thick=3
          ;; Readno + dark
          oplot, [-9e9,9e9], replicate(s2n_fstrct.noise^2,2), linestyle=2, $
            color=clr.green, thick=2
          oplot, [-9e9,9e9], replicate(s2n_fstrct.ndark,2), linestyle=2, $
            color=clr.brown, thick=2

          ;; Label
          xyouts, mnx*0.95, mxy*0.9, 'Object', color=clr.black, chars=2.1
          xyouts, mnx*0.95, mxy*0.8, 'Sky', color=clr.blue, chars=2.1
          xyouts, mnx*0.95, mxy*0.7, 'Dark', color=clr.brown, chars=2.1
          xyouts, mnx*0.95, mxy*0.6, 'ReadNo', color=clr.green, chars=2.1

      end
      else: stop
  endcase

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
pro hires_s2ngui_reinit, state
  ;; READ, DARK
  widget_control, state.rnd_id, set_value=state.str_instr.readno 
  widget_control, state.rn_id, set_value=state.str_instr.readno 
  widget_control, state.dark_id, set_value=state.str_instr.dark
  widget_control, state.darkd_id, set_value=state.str_instr.dark

  ;; Binning
  widget_control, state.bins_id, set_value=state.str_instr.bins 
  widget_control, state.binsd_id, set_value=state.str_instr.bins 
  widget_control, state.bind_id, set_value=state.str_instr.bind
  widget_control, state.bindd_id, set_value=state.str_instr.bind

  ;; Pixel size
  widget_control, state.pixsz_id, set_value=state.str_instr.PIXEL_SIZE

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; MAIN PROGRAM ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
pro hires_s2ngui, flg_instr, infil=infil, XSIZE=xsize,  YSIZE=ysize

common hires_s2ngui_cmmn

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'hires_s2ngui, flg (1=OHIRES,2=NHIRES,3=MTHR), INFIL= [v1.0]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100

  hires_s2ngui_initcommon

  ;; Initialize
  x_inithires, str_instr, flg_instr, INFIL=infil, STR_TEL=str_tel, $
    DECKER=decker
  str_obs = x_obsinit(infil)

; STATE

  state = {             $
            nwv: 0L, $
            dwv: 100., $
            wave: fltarr(1000), $
            s2n: fltarr(1000), $
            iorder: intarr(1000), $
            wvmn: 4000., $
            wvmx: 8000., $
            infil: '', $
            deckers: ['B5','C1','C5','D1'], $
            deckidx: 0L, $
            instr: ['OHIRES','NHIRES', 'MTHR'], $
            pixel: 0., $
            flg_plot: 0, $
            str_instr: str_instr, $       ; PLOTTING LINES
            str_tel: str_tel, $
            str_obs: str_obs, $
            flg_instr: flg_instr, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            psfile: 0, $
            help: strarr(50), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            size: lonarr(2), $
            base_id: 0L, $      ; Widgets
            ldraw_id: 0L, $    ; Lya
            ltext_id: 0L, $    ; 
            ldrawbase_id: 0L, $
            fxval_id: 0L, $
            iwvval_id: 0L, $
            mdraw_id: 0L, $       ; Spec Window
            mdrawbase_id: 0L, $
            swvval_id: 0L, $
            swidth_id: 0L, $
            swidthd_id: 0L, $
            xmax_id: 0L, $
            instr_id: 0L, $
            expt_id: 0L, $
            exptd_id: 0L, $
            moon_id: 0L, $
            moond_id: 0L, $
            seeing_id: 0L, $
            seeingd_id: 0L, $
            airm_id: 0L, $
            airmd_id: 0L, $
            top_id: 0L, $
            draw_id: 0L, $
            bottom_id: 0L, $
            rhs_id: 0L, $
            info_id: 0L, $
            dwv_id: 0L, $
            dwvd_id: 0L, $
            scr1_id: 0L, $
            scr2_id: 0L, $
            hits_id: 0L, $
            NHI_id: 0L, $
            NHIb_id: 0L, $
            mtl_id: 0L, $
            stat_id: 0L, $
            wvmn_id: 0L, $
            wvmnd_id: 0L, $
            wvmx_id: 0L, $
            wvmxd_id: 0L, $
            bind_id: 0L, $
            bindd_id: 0L, $
            bins_id: 0L, $
            binsd_id: 0L, $
            slit_id: 0L, $
            rn_id: 0L, $
            rnd_id: 0L, $
            dark_id: 0L, $
            darkd_id: 0L, $
            plot_id: 0L, $
            mag_id: 0L, $
            magd_id: 0L, $
            magtype_id: 0L, $
            pixsz_id: 0L, $
            resd_id: 0L, $
            help_text_id: 0L, $
	    ew: 0L$
          }

;;;;;;;;;;;;;;
; SETUP LINES

; Other setup

;    WIDGET
  base = WIDGET_BASE( title = 'hires_s2ngui: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=20L,/align_center)
  state.base_id = base
  

  state.top_id = WIDGET_BASE( state.base_id, /row, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(ysize/3.), $
                              uvalue='TOP_BASE', frame=2)
  state.bottom_id = WIDGET_BASE( state.base_id, /row, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(2*ysize/3.), $
                              uvalue='BOT_BASE', frame=2)
  state.draw_id = widget_draw(state.bottom_id, xsize=xsize, $
                              ysize=round(3*ysize/4), /frame, retain=2, $
                              uvalue='DRAW')

;;;;;; Info window ;;;;;;;;;;;

  ;; Slit
  instrinf = widget_base(state.top_id, /row, /align_center, frame=2)
  state.instr_id = WIDGET_LIST(instrinf, VALUE=state.instr, $
                                uvalue='INSTR', ysize=5)
  widget_control, state.instr_id, set_list_select=(state.flg_instr-1)
  slitb = widget_base(instrinf, /column, /align_left, frame=2)
  swidthb = widget_base(slitb, /row, /align_left, frame=2)
  state.swidth_id = CW_FIELD(swidthb, value=state.str_instr.swidth, xsize=5, $
                           ysize=1, title='SWIDTH', /return_events, $
                           uvalue='SWIDTH')
  state.swidthd_id = CW_FIELD(swidthb, value=state.str_instr.swidth, $
                              xsize=5, ysize=1, title=' ', uvalue='SWIDTHD')
  resolu = state.str_instr.R / ( state.str_instr.swidth / $
                                 state.str_tel.plate_scale / $
                                 (state.str_instr.pixel_size/100) )
  state.resd_id = CW_FIELD(slitb, value=resolu, $
                           xsize=6, ysize=1, title='RESOLUTION', $
                           uvalue='RESO')
;  mtch = where(strtrim(decker,2) EQ state.deckers, nmt)
;  mtch = where(strtrim(decker,2) EQ state.deckers, nmt)
;  if nmt EQ 0 then stop else begin
;      widget_control, state.decker_id, set_list_select=mtch[0] 
;      state.deckidx = mtch[0]
;  endelse

  ;; CCD
  binning = widget_base(instrinf, /column, /align_left, frame=2)
  ;;
  bindbase = widget_base(binning, /row, /align_left)
  state.bind_id = CW_FIELD(bindbase, value=state.str_instr.bind, xsize=2, $
                           ysize=1, title='BINWAVE', /return_events, $
                           uvalue='BIND')
  state.bindd_id = CW_FIELD(bindbase, value=state.str_instr.bind, xsize=2, $
                           ysize=1, title=' ', uvalue='BINDD')
  ;;
  binsbase = widget_base(binning, /row, /align_left)
  state.bins_id = CW_FIELD(binsbase, value=state.str_instr.bins, xsize=2, $
                           ysize=1, title='BINSPAT', /return_events, $
                           uvalue='BINS')
  state.binsd_id = CW_FIELD(binsbase, value=state.str_instr.bins, xsize=2, $
                           ysize=1, title=' ', uvalue='BINSD')
  ;;
  rnbase = widget_base(binning, /row, /align_left)
  state.rn_id = CW_FIELD(rnbase, value=state.str_instr.readno, xsize=4, $
                           ysize=1, title='READNO', /return_events, $
                           uvalue='READNO')
  state.rnd_id = CW_FIELD(rnbase, value=state.str_instr.readno, xsize=4, $
                           ysize=1, title=' ', uvalue='RND')
  ;;
  darkbase = widget_base(binning, /row, /align_left)
  state.dark_id = CW_FIELD(darkbase, value=state.str_instr.dark, xsize=4, $
                           ysize=1, title='DARK', /return_events, $
                           uvalue='DARK')
  state.darkd_id = CW_FIELD(darkbase, value=state.str_instr.dark, xsize=4, $
                           ysize=1, title=' ', uvalue='DARKD')
  state.pixsz_id = CW_FIELD(binning, value=state.str_instr.PIXEL_SIZE, xsize=4, $
                    ysize=1, title='PIXEL SIZE (mic)', $
                    uvalue='PIXSZ')

  ;; Wavelength
  waveinf = widget_base(state.top_id, /column, /align_center, frame=2)
  dwvbase = widget_base(waveinf, /row, /align_left)
  state.dwv_id = CW_FIELD(dwvbase, value=state.dwv, xsize=5, ysize=1, $
                            title='DELTA WAVE', /return_events, uvalue='DWV')
  state.dwvd_id = CW_FIELD(dwvbase, value=state.dwv, xsize=5, $
                           ysize=1, title=' ', uvalue='DWVD')
  ;;
  wvmnbase = widget_base(waveinf, /row, /align_left)
  state.wvmn_id = CW_FIELD(wvmnbase, value=state.wvmn, xsize=6, ysize=1, $
                            title='WAVE MIN', /return_events, uvalue='WVMN')
  state.wvmnd_id = CW_FIELD(wvmnbase, value=state.wvmn, xsize=6, ysize=1, $
                            title=' ', uvalue='WVMND')
  ;;
  wvmxbase = widget_base(waveinf, /row, /align_left)
  state.wvmx_id = CW_FIELD(wvmxbase, value=state.wvmx, xsize=6, ysize=1, $
                            title='WAVE MAX', /return_events, uvalue='WVMX')
  state.wvmxd_id = CW_FIELD(wvmxbase, value=state.wvmx, xsize=6, ysize=1, $
                            title=' ', uvalue='WVMXD')

  ;; Obs
  obsinf = widget_base(state.top_id, /column, /align_center, frame=2)
  expbase = widget_base(obsinf, /row, /align_left)
  state.expt_id = CW_FIELD(expbase, value=state.str_obs.exptime, $
                           xsize=7, ysize=1, $
                           title='EXPTIME (s)', /return_events, $
                           uvalue='EXPTIME')
  state.exptd_id = CW_FIELD(expbase, value=state.str_obs.exptime, $
                           xsize=7, ysize=1, $
                           title=' ', uvalue='EXPTIMED')
  ;;
  seeingbase = widget_base(obsinf, /row, /align_left)
  state.seeing_id = CW_FIELD(seeingbase, value=state.str_obs.seeing, $
                             xsize=4, ysize=1, $
                             title='SEEING (FWHM)', /return_events, $
                             uvalue='SEEING')
  state.seeingd_id = CW_FIELD(seeingbase, value=state.str_obs.seeing, $
                             xsize=4, ysize=1, $
                             title=' ', uvalue='SEEINGD')
  state.slit_id = CW_FIELD(obsinf, value=0., xsize=4, ysize=1, $
                             title='SLIT LOSS (%)', $
                             uvalue='SLIT')

  moonbase = widget_base(obsinf, /row, /align_left)
  state.moon_id = CW_FIELD(moonbase, value=state.str_obs.mphase, $
                           xsize=2, ysize=1, $
                           title='MOON PHASE (days)', /return_events, $
                           uvalue='MOON')
  state.moond_id = CW_FIELD(moonbase, value=state.str_obs.mphase, $
                           xsize=2, ysize=1, $
                           title=' ',  uvalue='MOOND')

  ;; OBS 2
  obsinf2 = widget_base(state.top_id, /column, /align_center, frame=2)
  airbase = widget_base(obsinf2, /row, /align_left)
  state.airm_id = CW_FIELD(airbase, value=state.str_obs.airmass, $
                           xsize=4, ysize=1, $
                            title='AIRMASS', /return_events, uvalue='AIRMASS')
  state.airmd_id = CW_FIELD(airbase, value=state.str_obs.airmass, $
                           xsize=4, ysize=1, $
                            title=' ', uvalue='AIRMASSD')
  ;;
  magbase = widget_base(obsinf2, /row, /align_left)
  state.mag_id = CW_FIELD(magbase, value=state.str_obs.mstar, $
                           xsize=7, ysize=1, $
                            title='MAG', /return_events, uvalue='MAGNITUDE')
  state.magd_id = CW_FIELD(magbase, value=state.str_obs.mstar, $
                           xsize=7, ysize=1, $
                            title=' ', uvalue='MAGNITUDED')
  case state.str_obs.mtype of
      1: stype = 'Johnson'
      else: stype = 'AB'
  endcase
  state.magtype_id = CW_FIELD(obsinf2, value=stype, $
                             xsize=7, ysize=1, $
                            title='MAGTYPE', /return_events, uvalue='MAGTYPE')


  ;; PLOT
  state.plot_id = WIDGET_LIST(instrinf, VALUE=['S2N', 'Eff', 'Noise'], $
                                uvalue='PLOT', ysize=3)
  widget_control, state.plot_id, set_list_select=0L
  state.flg_plot = 0

  ;;BUTTONS
  butbase = widget_base(state.top_id, /column, /align_center, frame=2)
  help = WIDGET_BUTTON(butbase, value='HELP',uvalue='HELP')
  plot = WIDGET_BUTTON(butbase, value='REPLOT',uvalue='REPLOT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
;  good = WIDGET_BUTTON(butbase, value='GOOD',uvalue='GOOD')
;  maybe = WIDGET_BUTTON(butbase, value='MAYBE',uvalue='MAYBE')
;  bad  = WIDGET_BUTTON(butbase, value='BAD', uvalue='BAD')
;  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
;  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
;  splt = WIDGET_BUTTON(butbase2, value='SPLT',uvalue='SPLT')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  ii = 0L
  state.help[ii] = '  :::Help Menu::: '
  ii=ii+1
  state.help[ii] = '---------------------'
  ii=ii+1
  state.help[ii] = 'Basics:  Change inputs in left window'
  ii=ii+1
  state.help[ii] = 'Hit return to udpate the value.'
  ii=ii+1
  state.help[ii] = 'Click on REPLOT to recalculate and plot'
  ii=ii+1
  state.help[ii] = '---------------------'
  ii=ii+1
  state.help[ii] = 'Keywords:'
  ii=ii+1
  state.help[ii] = 'SWIDTH = slit width (")'
  ii=ii+1
  state.help[ii] = 'BINWAVE = Binning in Spectral dimension'
  ii=ii+1
  state.help[ii] = 'BINSPAT = Binning in Spatial dimension'

;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ;; Set qso and qal
;  if qalstr[0].zabs[0] EQ 0. then sdss_chkciv_next, state

  ;; Load data
  hires_s2ngui_setwv, state
  hires_s2ngui_gets2n, state

  ; PLOT
  hires_s2ngui_updplot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'hires_s2ngui', base

  return
end
	
