;+ 
; NAME:
; esi_s2ngui
;    Version 1.0
;
; PURPOSE:
;   This launches a GUI which allows a user to fiddle around and
;   estimate S2N for ESI observations
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
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   March-2011  JXP     v1.0    Original version
;   May-2011    GDW     v1.2    Change formula for "resolu" to match esi_ions2n;
;                               Change minimum wavelength to 4000
;                               (else get bad results);
;                               use !x.crange for plotting
;   2011-Oct-12 GDW     v1.3    Undo changes to resolution; now uses
;                               same method as HIRES
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
pro esi_s2ngui_initcommon
;
common esi_s2ngui_cmmn, $
  s2n_fstrct
end

;;;;
; Events
;;;;

pro esi_s2ngui_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'HELP': x_helpwidg, state.help, YSIZE=40
      'BINNING': begin
          widget_control, state.binning_id, get_value=tmp
          state.str_instr.bins = long(strmid(state.binning[ev.index],0,1))
          state.str_instr.bind = long(strmid(state.binning[ev.index],2,1))
      end
      'EXPTIME': begin
          widget_control, state.expt_id, get_value=tmp
          state.str_obs.exptime = float(tmp)>1
          widget_control, state.exptd_id, set_value=state.str_obs.exptime
      end
      'MAGNITUDE': begin
          widget_control, state.mag_id, get_value=tmp
          state.str_obs.mstar = float(tmp)
          widget_control, state.magd_id, set_value=tmp
      end
      'SEEING': begin
          widget_control, state.seeing_id, get_value=tmp
          state.str_obs.seeing = float(tmp) > 0.4
          widget_control, state.seeingd_id, set_value=state.str_obs.seeing
      end
      'MOON': begin
          widget_control, state.moon_id, get_value=tmp
          state.str_obs.mphase = 0 > long(tmp) < 14
          widget_control, state.moond_id, set_value=state.str_obs.mphase
      end
      'AIRMASS': begin
          widget_control, state.airm_id, get_value=tmp
          state.str_obs.airmass = (float(tmp)>1.)
          widget_control, state.airmd_id, set_value=state.str_obs.airmass
      end
      'SWIDTH': begin
          widget_control, state.swidth_id, get_value=tmp
          state.str_instr.swidth = float(tmp[ev.index])
          resolu = state.str_instr.R / ( state.str_instr.swidth / $
                                         state.str_instr.scale_para)
          widget_control, state.resd_id, set_value=resolu
      end
      'PLOT': begin
          state.flg_plot = ev.index
          esi_s2ngui_gets2n, state
          esi_s2ngui_updplot, state
      end
      'DWV': begin
          widget_control, state.dwv_id, get_value=tmp
          state.dwv = tmp
          esi_s2ngui_setwv, state
          widget_control, state.dwvd_id, set_value=tmp
      end
      'WVMN': begin
          widget_control, state.wvmn_id, get_value=tmp
          state.wvmn = tmp > 2000.
          esi_s2ngui_setwv, state
          widget_control, state.wvmnd_id, set_value=state.wvmn
      end
      'WVMX': begin
          widget_control, state.wvmx_id, get_value=tmp
          state.wvmx = 11000. < tmp > (state.wvmn+100)
          esi_s2ngui_setwv, state
          widget_control, state.wvmxd_id, set_value=state.wvmx
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  esi_s2ngui_gets2n, state
  esi_s2ngui_updplot, state
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;

pro esi_s2ngui_setwv, state

 ;; Wave array
 state.nwv = long((state.wvmx-state.wvmn)/state.dwv) + 1
 state.wave[0:state.nwv-1] = state.wvmn + findgen(state.nwv)*state.dwv

 return
end

pro esi_s2ngui_gets2n, state

common esi_s2ngui_cmmn

  ;; Calc
  keck_calcs2n, state.wave[0:state.nwv-1], 3, $
    STATE=state, /nopr, S2N=s2n,  FSTRCT=s2n_fstrct
  state.s2n[0:state.nwv-1] = s2n

  ;; Update slit throughput
  widget_control, state.slit_id, set_value=(1.-s2n_fstrct.slit0)*100

  return
end
 
  
pro esi_s2ngui_updplot, state

common esi_s2ngui_cmmn
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
            ytitle='S/N (per '+string(round(pixel),'(i2)')+'km/s pix)', $
            yrange=[mny*0.9,mxy*1.1], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata
          
          oplot, [wave], [sn],  psym=-1, color=clr.blue
          
          ;; Label
          xpos = mnx + (mxx-mnx)*0.95
          dy = mxy-mny
      end
      1: begin ; Efficiency
          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.7,$
            xmargin=[10,2], ymargin=[5,3], xtitle='!17Wavelength (Ang)', $
            ytitle='Efficiency Curve', $
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
          oplot, !x.crange, replicate(s2n_fstrct.noise^2,2), linestyle=2, $
            color=clr.green, thick=2
          oplot, !x.crange, replicate(s2n_fstrct.ndark,2), linestyle=2, $
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
pro esi_s2ngui_reinit, state
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
pro esi_s2ngui, infil=infil, XSIZE=xsize,  YSIZE=ysize

common esi_s2ngui_cmmn

;
;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-100
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100

  esi_s2ngui_initcommon

  ;; Initialize
  x_initesi, str_instr, INFIL=infil, STR_TEL=str_tel, SLIT=slit
  str_obs = x_obsinit(infil)

; STATE

  state = {             $
            nwv: 0L, $
            dwv: 10., $
            wave: fltarr(1000), $
            s2n: fltarr(1000), $
            wvmn: 4000., $
            wvmx: 10000., $
            infil: '', $
            slits: ['0.3', '0.5', '0.75', '1.0'], $
            binning: ['1x1', '2x1', '3x1', '2x2'], $
          sky: ['NewMoon','FullMoon'], $
            deckidx: 0L, $
            instr: ['ESI'], $
            pixel: 0., $
            flg_plot: 0, $
            str_instr: str_instr, $       ; PLOTTING LINES
            str_tel: str_tel, $
            str_obs: str_obs, $
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
          sky_id: 0L, $
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
            binning_id: 0L, $
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
  base = WIDGET_BASE( title = 'esi_s2ngui: Check spectra', /column, $
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

  dumi = widget_text(state.top_id, value='ESI S/N GUI v1.0')

  ;;BUTTONS
  butbase = widget_base(state.top_id, /column, /align_center, frame=2)
  help = WIDGET_BUTTON(butbase, value='HELP',uvalue='HELP')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  state.plot_id = WIDGET_DROPLIST(butbase, VALUE=['S2N', 'Eff', 'Noise'], $
                                title='PLOT', uvalue='PLOT')
  widget_control, state.plot_id, set_list_select=0L
  state.flg_plot = 0


  ;; Slit
  instrinf = widget_base(state.top_id, /row, /align_center, frame=2)
  slitb = widget_base(instrinf, /column, /align_left, frame=2)
  swidthb = widget_base(slitb, /row, /align_left, frame=2)
  state.swidth_id = widget_droplist(swidthb, value=state.slits, $
                                    title='SlitWidth(")', uvalue='SWIDTH')
  widget_control, state.swidth_id, set_droplist_select=2L
  state.str_instr.swidth = float(state.slits[2])
  resolu = state.str_instr.R / ( state.str_instr.swidth / $
                                 state.str_instr.scale_para)
  state.resd_id = CW_FIELD(slitb, value=resolu, $
                           xsize=6, ysize=1, title='RESOLUTION', $
                           uvalue='RESO', /NOEDIT)
  state.slit_id = CW_FIELD(slitb, value=0., xsize=4, ysize=1, $
                             title='SLIT LOSS (%)', /NOEDIT, $
                             uvalue='SLIT')

  binb = widget_base(instrinf, /column, /align_center, frame=2)
  state.binning_id = widget_droplist(binb, value=state.binning, $
                                    title='Binning', uvalue='BINNING')
  state.str_instr.bins = long(strmid(state.binning[0],0,1))
  state.str_instr.bind = long(strmid(state.binning[0],2,1))
  state.pixsz_id = CW_FIELD(binb, value=state.str_instr.PIXEL_SIZE, xsize=4, $
                    ysize=1, title='NATIVE PIXEL (microns)', $
                    uvalue='PIXSZ', /noedit)

  ;; Wavelength
  options = widget_base(state.top_id, /row, /align_center, frame=2)
  dwvbase = widget_base(options, /column, /align_left)
  state.dwv_id = CW_FIELD(dwvbase, value=state.dwv, xsize=5, ysize=1, $
                            title='DELTA WAVE', /return_events, uvalue='DWV', /column)
  state.dwvd_id = CW_FIELD(dwvbase, value=state.dwv, xsize=5, $
                           ysize=1, title=' ', uvalue='DWVD', /noedit)
  ;;
  wvmnbase = widget_base(options, /column, /align_center)
  state.wvmn_id = CW_FIELD(wvmnbase, value=state.wvmn, xsize=6, ysize=1, $
                            title='WAVE MIN', /return_events, uvalue='WVMN', /column)
  state.wvmnd_id = CW_FIELD(wvmnbase, value=state.wvmn, xsize=6, ysize=1, $
                            title=' ', uvalue='WVMND', /noedit)
  ;;
  wvmxbase = widget_base(options, /column, /align_center)
  state.wvmx_id = CW_FIELD(wvmxbase, value=state.wvmx, xsize=6, ysize=1, $
                            title='WAVE MAX', /return_events, uvalue='WVMX', /column)
  state.wvmxd_id = CW_FIELD(wvmxbase, value=state.wvmx, xsize=6, ysize=1, $
                            title=' ', uvalue='WVMXD', /noedit)

  ;; Obs
  expbase = widget_base(options, /column, /align_center)
  state.expt_id = CW_FIELD(expbase, value=state.str_obs.exptime, $
                           xsize=8, ysize=1, $
                           title='EXPTIME (s)', /return_events, $
                           uvalue='EXPTIME', /column)
  state.exptd_id = CW_FIELD(expbase, value=state.str_obs.exptime, $
                           xsize=7, ysize=1, $
                           title=' ', uvalue='EXPTIMED', /noedit)
  ;;
  seeingbase = widget_base(options, /column, /align_center)
  state.seeing_id = CW_FIELD(seeingbase, value=state.str_obs.seeing, $
                             xsize=4, ysize=1, $
                             title='SEEING (FWHM)', /return_events, $
                             uvalue='SEEING', /column)
  state.seeingd_id = CW_FIELD(seeingbase, value=state.str_obs.seeing, $
                             xsize=4, ysize=1, $
                             title=' ', uvalue='SEEINGD', /noedit)

  moonbase = widget_base(options, /column, /align_center)
  state.sky_id = widget_droplist(moonbase, value=state.sky, $
                                   title='Sky:', uvalue='SKY')
;  state.moon_id = CW_FIELD(moonbase, value=state.str_obs.mphase, $
;                           xsize=2, ysize=1, $
;                           title='MOON (days)', /return_events, $
;                           uvalue='MOON', /column)
;  state.moond_id = CW_FIELD(moonbase, value=state.str_obs.mphase, $
;                           xsize=2, ysize=1, $
;                           title=' ',  uvalue='MOOND')

  ;; OBS 2
  airbase = widget_base(options, /column, /align_center)
  state.airm_id = CW_FIELD(airbase, value=state.str_obs.airmass, $
                           xsize=4, ysize=1, $
                            title='AIRMASS', /return_events, uvalue='AIRMASS', /column)
  state.airmd_id = CW_FIELD(airbase, value=state.str_obs.airmass, $
                           xsize=4, ysize=1, $
                            title=' ', uvalue='AIRMASSD',/noedit)
  ;;
  magbase = widget_base(options, /column, /align_center)
  state.mag_id = CW_FIELD(magbase, value=state.str_obs.mstar, $
                           xsize=7, ysize=1, /column, $
                            title='AB MAG', /return_events, uvalue='MAGNITUDE')
  state.magd_id = CW_FIELD(magbase, value=state.str_obs.mstar, $
                           xsize=7, ysize=1, $
                            title=' ', uvalue='MAGNITUDED',/noedit)
  state.str_obs.mtype = 2 ;; AB



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  ii = 0L
  state.help[ii] = '  :::Help Menu::: '
  ii=ii+1
  state.help[ii] = '---------------------'
  ii=ii+1
  state.help[ii] = 'Basics:  Change inputs as you wish'
  ii=ii+1
  state.help[ii] = '---------------------'
  ii=ii+1
  state.help[ii] = 'Keywords:'
  ii=ii+1
  state.help[ii] = 'In S2N, blue curve is for HIRESb and red is HIRESr'
  ii=ii+1

;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ;; Set qso and qal
;  if qalstr[0].zabs[0] EQ 0. then sdss_chkciv_next, state

  ;; Load data
  esi_s2ngui_setwv, state
  esi_s2ngui_gets2n, state

  ; PLOT
  esi_s2ngui_updplot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'esi_s2ngui', base

  return
end
	
