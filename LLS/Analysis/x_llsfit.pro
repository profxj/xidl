;+ 
; NAME:
; x_llsfit
;    Version 1.1
;
; PURPOSE:
;   Implements a GUI to perform a N(HI) fit to a LLS.  I wonder
;  if this really works right now.
;
; CALLING SEQUENCE:
;   x_llsfit, yin, zin, INFLG=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_llsfit, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   10-Jun-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_llsfit_icmmn

  common x_specplot_lines, $
    flg_lines, $
    lines, $
    zabs
  flg_lines = 0

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro x_llsfit_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  flg_allplt = 0
  flg_sngplt = 0

  case uval of
      'SNGLIST' : begin
          state.cursngl = ev.index
          x_llsfit_resetsng, state
          flg_sngplt = 1
      end
      'BVAL' : begin
          widget_control, state.bval_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = tmp
          x_fitline_updfit, state, FLG_PLT=flg_sngplt
          flg_allplt = flg_sngplt
      end
      'NCOLM' : begin
          widget_control, state.Ncolm_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = tmp
          x_fitline_updfit, state, FLG_PLT=flg_sngplt
          flg_allplt = flg_sngplt
      end
      'VMIN' : begin
          widget_control, state.vmin_id, get_value=tmp
          state.vmnx[0] = tmp
          x_llsfit_setvelo, state
          flg_allplt = 1
      end
      'VMAX' : begin
          widget_control, state.vmax_id, get_value=tmp
          state.vmnx[1] = tmp
          x_llsfit_setvelo, state
          flg_allplt = 1
      end
      'SNGDRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
      end
      'SNGLDRAW' : begin
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
                      4 : begin
                          gdset = where(state.lines[0:state.nlin-1].set $
                                        EQ state.curlin,ngd)
                          state.lines[gdset].zabs = $
                            (state.xpos / 1215.6701) - 1.
                          x_fitline_updfit, state, FLG_PLT=flg_sngplt
                          flg_allplt = flg_sngplt
                          widget_control, state.zabs_id, $
                            set_value=state.lines[gdset[0]].zabs
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
;                  widget_control, state.xpos_id, set_value=state.xpos
;                  widget_control, state.ypos_id, set_value=state.ypos
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
          ;; FIT LINE
          x_fitline, state, eventch, FLG_PLT=flg_sngplt
          flg_allplt = flg_sngplt
      end
      'PRINT' : x_llsfit_print, state
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  ;; PXMNX
  if state.xymnx[0] NE state.old_xymnx[0] OR $
    state.xymnx[2] NE state.old_xymnx[2] then begin
      state.old_xymnx = state.xymnx
      state.xpmnx = x_getxpmnx(state)
  endif

  ;; PLOT
  if flg_allplt EQ 1 then x_llsfit_allplt, state
  if flg_sngplt EQ 1 then x_llsfit_sngplt, state

  ;; RETURN
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SET VELO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_llsfit_setvelo, state

  ;; Call x_allvelo
  state.all_velo[*,0:state.ntrans-1] = x_allvelo(state.wave, state.all_zabs, $
                                                 state.velplt[0:state.ntrans-1].wrest,$
                                                 state.vmnx, $
                                                 all_pmnx=all_pmnx, NPIX=5000L)
  state.all_pmnx[*,0:state.ntrans-1] = all_pmnx

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro x_llsfit_allplt, state
  
common x_specplot_lines

  ;; GD LIN
;  gdlin = where(lines.flg NE 0)

  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.alldraw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  !P.MULTI= [0,state.nx,state.hplt,0,1]

  pmax = state.nx*state.hplt - 1
  for i=0L, pmax do begin
      ;; Checks
      if i EQ state.ntrans then break
      if state.all_pmnx[2,0] EQ 0 then continue

      ;; Get pxmin, pxmax + velo array
;      x_pixminmax, state.wave, lines[gdlin[i]].wave, state.zabs, $
;        state.svxymnx[0], state.svxymnx[2], PIXMIN=pixmin, PIXMAX=pixmax, $
;        VELO=velo
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]
      
      ;; Plot
      if (i NE state.ntrans-1) AND ((i+1) MOD state.hplt NE 0) then begin
          spaces = replicate('!17 ',30)
          plot, state.all_velo[0:state.all_pmnx[2,i]], $
            state.fx[pixmin:pixmax], xrange=state.vmnx, $
            yrange=state.velplt[i].ymnx, xtickn=spaces, xmargin=[9,3], $
            ymargin=[0,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.all_velo[0:state.all_pmnx[2,i]], $
            state.fx[pixmin:pixmax], xrange=state.vmnx, $
            yrange=state.velplt[i].ymnx, xmargin=[9,3], ymargin=[3,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse

      ;; VOIGT
      if state.nlin NE 0 then begin
          oplot, state.all_velo[0:state.all_pmnx[2,i]], $
            state.fit[pixmin:pixmax], color=clr.green
      endif
      
      ;; Labels
;      xyouts, 0.07*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], $
;        mn+(mx-mn)*0.05, strtrim(lines[gdlin[i]].name,2), $
;        color=clr.red, charsize=1.5
      
      ;; Lines
      oplot, [0., 0.], state.velplt[i].ymnx, color=clr.blue, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.green, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  endfor
  
  !P.MULTI= [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro x_llsfit_sngplt, state
  
  ;; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.sngldraw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  plot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.fx[state.xpmnx[0]:state.xpmnx[1]], psym=10, $
    position=state.pos, $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.7, $
    ycharsize=1.7

  ;; VOIGT
  if state.nlin NE 0 then begin
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.fit[state.xpmnx[0]:state.xpmnx[1]], color=clr.green
      ;; Mark all lines
      for i=0L,state.nlin-1 do begin
          oplot, replicate( (state.lines[i].zabs+1.)*$
                            state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.blue, linestyle=1
      endfor
      ;; Mark current set
      gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin, ngd)
      for i=0L,ngd-1 do begin
          oplot, replicate( (state.lines[gdlin[i]].zabs+1.)*$
                            state.lines[gdlin[i]].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.red, linestyle=1
      endfor
  endif
      
  ;; Lines
  oplot, [-10000., 10000.], [0.,0.], color=clr.green, linestyle=3
  oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  
end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro x_llsfit_print, state

  widget_control, /hourglass   
; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
;  !p.charthick = 1

  device, decompose=0
  ; Get file name
  ipos = strpos(state.obj_fil, '.fits')
  psfil = strmid(state.obj_fil, 0, ipos)+'.ps'
  ps_open, file=psfil, font=1, /color
  state.psfile = 1
  for qq=0L,state.npg-1 do begin
      state.curpg = qq
      x_llsfit_plot, state
  endfor
  ps_close, /noprint, /noid
  spawn, 'gzip -f '+psfil
  device, decomposed=svdecomp
  state.psfile = 0
;  !p.thick = 1
;  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;
;  Reset SNGL
;;;;;;;;;;;;;;;;;;;;

pro x_llsfit_resetsng, state

  ;; Current single
  cen = state.velplt[state.cursngl].wrest*(state.all_zabs+1.)

  ;; Ranges
  state.xymnx[0] = cen - 50.
  state.xymnx[2] = cen + 50.

  state.xymnx[1] = state.velplt[state.cursngl].ymnx[0]
  state.xymnx[3] = state.velplt[state.cursngl].ymnx[1]
  
  ;; XPMNX
  state.xpmnx = x_getxpmnx(state)

  ;; SVXYMNX, OLD
  state.svxymnx = state.xymnx
  state.old_xymnx = state.xymnx

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_llsfit, yin, zin, ysin=ysin, INFLG=inflg, NPLT=nplt, NRM=nrm

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'x_llsfit, spec, zin, INFLG=, NPLT=, /NRM  [v1.1]'
    return
  endif 

;  Optional Keywords
  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-300
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-400
  if not keyword_set( NPLT ) then nplt = 15L
  if not keyword_set( FWHM ) then fwhm = 4.

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0
  ydat = x_readspec(yin, INFLG=inflg, head=head, npix=npix, $
                    WAV=xdat, FIL_SIG=ysin, SIG=ysig)

  

  tmp = { velpltstrct }
  tmp1 = { abslinstrct }

; STATE

  state = {             $
            hplt: 5L, $
            nx: 3L, $
            npix: npix, $
            fx: ydat, $
            wave: xdat, $
            sig: fltarr(npix), $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            cursngl: 0L, $
            xpmnx: lonarr(2), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            old_xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            zabs: 0., $
            nlin: 0, $   ; VOIGT LINES
            nset: -1, $   
            lines: replicate(tmp1,500), $ 
            fit: fltarr(n_elements(ydat)) + 1., $
            conti: 1., $  ; ALL STUFF
            curlin: 0L, $
            all_zabs: zin, $
            vmnx: [-1000.,1000.], $
            ntrans: 0L, $  ; PLOTTING LINES
            all_velo: dblarr(5000, 100), $  
            all_pmnx: lonarr(3, 100), $  
            velplt: replicate(tmp, 100), $
            psfile: 0, $
            size: lonarr(2), $
            xpos: 0., $
            ypos: 0., $
            xcurs: 0., $
            ycurs: 0., $
            base_id: 0L, $      ; Widgets
            alldraw_id: 0L, $
            alldrawbase_id: 0L, $
            sngldraw_id: 0L, $
            sngldrawbase_id: 0L, $
            text_id: 0L, $
            vmin_id: 0L, $
            vmax_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            sngdroplist_id: 0L, $
            allzabs_id: 0L, $
            zabs_id: 0L, $
            bval_id: 0L, $
            Ncolm_id: 0L, $
            help_text_id: 0L $
          }

; NORM
  if keyword_set( NRM ) then state.flg_nrm = 1
  if keyword_set( YSIG ) then state.sig = temporary(ysig)

; LINELIST for VELPLT
  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE
  x_llsfit_icmmn
  zabs = state.all_zabs

  ;; VELPLT
  if not keyword_set( tlist ) then $
    tlist = getenv('XIDL_DIR')+'/LLS/Lines/LLSvp_std.lst'
  readcol, tlist, trans, /sil
  state.ntrans = n_elements(trans)
  srt = sort(trans)
  state.velplt[0:state.ntrans-1].wrest = trans[reverse(srt)]
  state.velplt.ymnx = [-0.09, 1.1]
  ;; Names
  for i=0L,state.ntrans-1 do begin
      getfnam, state.velplt[i].wrest, fv, nm
      state.velplt[i].name = nm
  endfor
  
;  state.flg_lines = 2
  

;    WIDGET
  base = WIDGET_BASE( title = 'x_llsfit: LLS fit', /row, $
                    UNAME='BASE', xoffset=100)
  state.base_id = base
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ALL ;;;;

  ;; Toolbar
  alltool = WIDGET_BASE( state.base_id, /column, /frame, /base_align_center,$
                         /align_center)

  ;; Version + Name + trace
  labelbase = widget_base(alltool, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_llsfit', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)

  ;; Vmin, vmax
  vmnx_base = WIDGET_BASE(alltool, /column, /base_align_left,$
                         /align_center)
  state.vmin_id = cw_field(vmnx_base, value=state.vmnx[0], $
                             /float, /return_events, xsize=9,$
                             title='vmin:', UVALUE='VMIN')
  state.vmax_id = cw_field(vmnx_base, value=state.vmnx[1], $
                             /float, /return_events, xsize=9,$
                             title='vmax:', UVALUE='VMAX')

 ;;      BUTTONS
;  butbase = widget_base(alltool, /column, /align_center)
  butbase2 = widget_base(alltool, /column, /align_center)
  print = WIDGET_BUTTON(butbase2, value='PRINT',uvalue='PRINT')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')

  state.allzabs_id = cw_field(alltool, title='zabs', $
                              value=state.all_zabs, /floating, $
                           /column, xsize=10, uvalue='ALLZABS')
  ;; DRAW
  state.alldrawbase_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='ALLDRAW_BASE', frame=2, $
                                   /tracking_events)

  state.alldraw_id = widget_draw(state.alldrawbase_id, xsize=600, ysize=1000L,$
                              /frame, $
                              uvalue='ALLDRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; SNGL ;;;;;;;;;;

  sngl_base = widget_base(base, /column, uvalue='SNGL_BASE', /frame)

  ;; DRAW
  state.size[0] = xsize
  state.size[1] = ysize
  state.sngldrawbase_id = widget_base(sngl_base, /column, /base_align_left, $
                                   uvalue='SNGLDRAW_BASE', frame=2, $
                                   /tracking_events)

  state.sngldraw_id = widget_draw(state.sngldrawbase_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /motion_events, /button_events, uvalue='SNGLDRAW')

  ;; Text
  state.text_id = widget_text(state.sngldrawbase_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')

  ;; Toolbar
  linbar = WIDGET_BASE( sngl_base, /row, /frame, /base_align_center,$
                           /align_center)
  state.zabs_id = cw_field(linbar, title='zabs', value=zabs, /floating, $
                           /column, xsize=10, uvalue='ZABS')
  state.Ncolm_id = cw_field(linbar, title='Ncolm', value=0., /floating, $
                           /column, xsize=10, /return_events, uvalue='NCOLM')
  state.bval_id = cw_field(linbar, title='bval', value=0., /floating, $
                           /column, xsize=10, /return_events, uvalue='BVAL')


  linelist = state.velplt[0:state.ntrans-1].name
  state.sngdroplist_id = widget_droplist(sngl_base, $
                                   frame = 1, $
                                   title = 'Transitions:', $
                                   uvalue = 'SNGLIST', $
                                   value = linelist)

  
; Realize
  WIDGET_CONTROL, base, /realize


  ;; Set SNGL line
  x_llsfit_resetsng, state

  ;; Set VELO
  x_llsfit_setvelo, state

  ;; PLOT
  x_llsfit_allplt, state
  x_llsfit_sngplt, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_llsfit', base

  !P.MULTI= [0,1,1]
  return
end

