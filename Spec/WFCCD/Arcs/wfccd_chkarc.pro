;+ 
; NAME:
; wfccd_chkarc
;    Version 1.0
;
; PURPOSE:
;   Plots a series of arc to allow a quick check
;
; CALLING SEQUENCE:
;   
;   wfccd_chkarc, wfccd, maskid, expsr, XSIZE=, YSIZE=
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
;   wfccd_chkarc, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;
; Events
;;;;

pro wfccd_chkarc_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DRAW' : begin
          case ev.type of
              0 :  ; Button press
              1 :  ; Button release
              2 :  ; Motion
          endcase
      end
      'PRINT' : wfccd_chkarc_print, state
      'NEXT' : begin
          state.curpg = state.curpg + 1
          if state.curpg GT state.npg-1 then state.curpg = 0L
          wfccd_chkarc_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'PREV' : begin
          state.curpg = state.curpg - 1
          if state.curpg LT 0 then state.curpg = state.npg-1
          wfccd_chkarc_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'NSPEC': begin
          widget_control, state.narc_id, get_value=nplt
          state.nplt = nplt
          state.npg =  (state.narc/state.nplt) + (state.narc MOD state.nplt NE 0)
          state.curpg = state.curpg < (state.npg-1)
          wfccd_chkarc_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
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
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkarc_Plot, state
  
; Plot Data

  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  !p.background = clr.white
  !p.color = clr.black
  !P.MULTI= [0,4,state.nplt]
  !x.margin = [0,0]
  !x.style = 1
  !y.style = 1

  pmin = state.curpg*state.nplt
  if state.curpg EQ state.npg - 1 then pmax = state.narc-1 $
  else pmax = (state.curpg+1)*state.nplt-1
  
  spaces = replicate('!17 ',30)
  for i=pmin,pmax do begin

      ; 3600--5550
      ; Plot
      if i NE pmax then begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[3700., 5180.],  $
            ymargin=[0,0], xtickn=spaces, ytickn=spaces, yticks=0
      endif else begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[3700., 5180.],  $
            ymargin=[3,0], ytickn=spaces, yticks=0, charsize=1.5
      endelse
      ; Overplot good lines
      for j=0L,state.nlin-1 do $
        oplot, [state.lines[j], state.lines[j]], [-1e4, 1e5], color=clr.red
      oplot, [4713.14, 4713.14], [-1e4, 1e5], color=clr.red
      xyouts, 3930., 100., 'RMS= '+string(state.rms[i], FORMAT='(f8.3)'), $
        charsize=1.3
      xyouts, 3930., 200., 'N= '+string(i, FORMAT='(i3)'), charsize=1.3

      ; 5650--7400
      ; Plot
      if i NE pmax then begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[5800, 6450.],  $
            ymargin=[0,0], xtickn=spaces, yticks=0, ytickn=spaces
      endif else begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[5800, 6450.],  $
            ymargin=[3,0], ytickn=spaces, yticks=0, charsize=1.5
      endelse
      ; Overplot good lines
      for j=0L,state.nlin-1 do $
        oplot, [state.lines[j], state.lines[j]], [-1e4, 1e5], color=clr.red

      ; 7400--9000
      ; Plot
      if i NE pmax then begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[6480, 7450.],  $
            ymargin=[0,0], xtickn=spaces, ytickn=spaces
      endif else begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[6480, 7450.],  $
            ymargin=[3,0], ytickn=spaces, charsize=1.5
      endelse
      ; Overplot good lines
      for j=0L,state.nlin-1 do $
        oplot, [state.lines[j], state.lines[j]], [-1e4, 1e5], color=clr.red

      ; 7400--9000
      ; Plot
      if i NE pmax then begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[7450, 9000],  $
            ymargin=[0,0], xtickn=spaces, ytickn=spaces
      endif else begin
          plot, state.wv[i,*], state.fx[i,*], xrange=[7450, 9000],  $
            ymargin=[3,0], ytickn=spaces, charsize=1.5
      endelse
      ; Overplot good lines
      for j=0L,state.nlin-1 do $
        oplot, [state.lines[j], state.lines[j]], [-1e4, 1e5], color=clr.red

  endfor

  !P.MULTI= [0,1,1]
  !x.margin = [10,3]

end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkarc_print, state

  widget_control, /hourglass   
; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
;  !p.charthick = 1

  device, decompose=0
  ; Get file name
  ipos = strpos(state.arc_fil, '.fits')
  psfil = strmid(state.arc_fil, 0, ipos)+'.ps'
  ps_open, file=psfil, font=1, /color
  state.psfile = 1
  for qq=0L,state.npg-1 do begin
      state.curpg = qq
      wfccd_chkarc_plot, state
  endfor
  ps_close, /noprint, /noid
  spawn, 'gzip -f '+psfil
  device, decomposed=svdecomp
  state.psfile = 0
;  !p.thick = 1
;  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkarc, wfccd, mask_id, expsr, XSIZE=xsize, YSIZE=ysize, $
                   PRINTONLY=printonly

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_chkarc, wfccd, mask_id, expsr, XSIZE=, YSIZE=,'
    print, '      (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 800

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(expsr) then exp = allexp[expsr] else exp=allexp[0]

; Arc structure
  i = strpos(wfccd[exp].arc_fil, '_')
  arcfil = 'Arcs/ArcS'+strmid(wfccd[exp].arc_fil,i)
  wfccd_readastrct, arcfil, wfarc
  gdarc = where(wfarc.flg_anly NE 0, narc)

; Set min to 0.
  for i=0L,narc-1 do begin
      a = where(wfarc[gdarc[i]].spec LE 0)
      if a[0] NE -1 then wfarc[gdarc[i]].spec[a] = 0.
  endfor

; Lines
  linelist = $
    getenv('XIDL_DIR')+'/Spec/Arcs/Lists/wfccdB_HeNe.lst'
  x_arclist, linelist, lines
  gdlin = where(lines.flg_qual EQ 5, ngd)

; STATE
  ninit = 7L

  state = {             $
            nplt: ninit, $
            npg: (narc/ninit) + (narc MOD ninit NE 0), $
            arc_fil: arcfil, $
            narc: narc, $
            wv: transpose(wfarc[gdarc].wave), $
            fx: transpose(wfarc[gdarc].spec), $
            rms: wfarc[gdarc].fit.rms, $
            nlin: ngd, $
            lines: lines[gdlin].wave, $
            curpg: 0L, $
            pos: [0.04,0.0,1.00,1.0], $ ; Plotting
            tmpxy: fltarr(4), $
            size: fltarr(2), $
            xcurs: 0., $
            xpos: 0., $
            psfile: 0, $
            base_id: 0L, $      ; Widgets
            draw_id: 0L, $
            drawbase_id: 0L, $
            pg_id: 0L, $
            narc_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            help_text_id: 0L $
          }

; PRINTONLY
  if keyword_set( PRINTONLY ) then begin
      state.psfile = 1
      wfccd_chkarc_print, state
      return
  endif
	
;    WIDGET
  base = WIDGET_BASE( title = 'wfccd_chkarc: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='wfccd_chkarc', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
;  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)

;;;;;;;;;;;;;;;;;;
; Pages
  state.pg_id = cw_field(toolbar, value=0L, $
                             /long, /return_events, xsize=3,$
                             title='Page')
  state.narc_id = cw_field(toolbar, value=state.nplt, $
                             /long, /return_events, xsize=3,$
                             title='Nplt', UVALUE='NSPEC')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  state.drawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_id = widget_draw(state.drawbase_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')

; XY position
;  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
;  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
;  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  butbase2 = widget_base(toolbar, /column, /align_center)
  print = WIDGET_BUTTON(butbase2, value='PRINT',uvalue='PRINT')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
;  strhelp = strarr(50)
;  strhelp = ['   Help Menu   ',$
;             'LMB - Truncate/Extend trace', $ 
;             'RMB - Contrast/Brightness', $
;;             'CMB/CMB - Zoom' $ 
;             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

; Realize
  WIDGET_CONTROL, base, /realize

  ; PLOT
  wfccd_chkarc_Plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'wfccd_chkarc', base

  !P.MULTI= [0,1,1]
  return
end

