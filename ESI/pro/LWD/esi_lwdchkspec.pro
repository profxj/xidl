;+ 
; NAME:
; esi_lwdchkspec
;    Version 1.0
;
; PURPOSE:
;   Plots a series of spectra to allow a quick check
;
; CALLING SEQUENCE:
;   
;   esi_lwdchkspec, esi, maskid, expsr, XSIZE=, YSIZE=
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
;   esi_lwdchkspec, esi, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

;pro esi_lwdchkspec_icmmn

;

;common esi_lwdchkspec_cmm, $
;  objstr

;end

;;;;
; Events
;;;;

pro esi_lwdchkspec_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DRAW' : begin
          case ev.type of
              0 : ; Button press
              1 :  ; Button release
              2 : begin ; Motion event
                  state.xcurs = ev.x
;                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state.xcurs,state.pos, $
                                         state.xymnx,$
                                         state.size) 
;                  state.ypos = xgety_plt(state.ycurs,state.pos, $
;                                         state.xymnx,$
;                                         state.size)
                  widget_control, state.xpos_id, set_value=state.xpos
;                  widget_control, state.ypos_id, set_value=state.ypos
              end
          endcase
      end
      'PRINT' : esi_lwdchkspec_print, state
      'NEXT' : begin
          state.curpg = state.curpg + 1
          if state.curpg GT state.npg-1 then state.curpg = 0L
          esi_lwdchkspec_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'PREV' : begin
          state.curpg = state.curpg - 1
          if state.curpg LT 0 then state.curpg = state.npg-1
          esi_lwdchkspec_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'NSPEC': begin
          widget_control, state.nspec_id, get_value=nplt
          state.nplt = nplt
          state.npg =  (state.nobj/state.nplt) + (state.nobj MOD state.nplt NE 0)
          state.curpg = state.curpg < (state.npg-1)
          esi_lwdchkspec_Plot, state
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

pro esi_lwdchkspec_Plot, state
  
; Plot Data

  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  !P.MULTI= [0,1,state.nplt]

  pmin = state.curpg*state.nplt
  if state.curpg EQ state.npg - 1 then pmax = state.nobj-1 $
  else pmax = (state.curpg+1)*state.nplt-1
  
  for i=pmin,pmax do begin
      npix = state.obj[i].npix
      ; Get median flux
      mdwv = where(state.obj[i].wave[0:npix-1] GT state.mdrng[0] AND $
                   state.obj[i].wave[0:npix-1] LT state.mdrng[1] AND $
                   state.obj[i].fx[0:npix-1] GT 0., nwv)
      if nwv NE 0 then mdfx = median(state.obj[i].fx[mdwv]) else mdfx = 1.

      ; Plot
      if i NE pmax then begin
          spaces = replicate('!17 ',30)
          plot, state.obj[i].wave[0:npix-1], $
            state.obj[i].fx[0:npix-1], $
            xrange=[state.xymnx[0],state.xymnx[2]], $
            yrange=[0., 2*mdfx], xtickn=spaces, xmargin=[9,0], ymargin=[0,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.obj[i].wave[0:npix-1], $
            state.obj[i].fx[0:npix-1], $
            xrange=[state.xymnx[0],state.xymnx[2]], $
            yrange=[0., 2*mdfx], xmargin=[9,0], ymargin=[3,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse
      ; Error
      oplot, state.obj[i].wave[0:npix-1], state.sig[i,0:npix-1], $
        color=clr.green
      ; Labels
      xyouts, 0.03*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.75*2*mdfx, $
        strtrim(state.obj[i].slit_id,2)+state.obj[i].obj_id, $
        color=clr.red, charsize=1.5
      xyouts, 0.03*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.45*2*mdfx, $
        'Obj: '+strtrim(i,2), color=clr.blue, charsize=1.5
  endfor

  !P.MULTI= [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro esi_lwdchkspec_print, state

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
      esi_lwdchkspec_plot, state
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

pro esi_lwdchkspec, esi, obj_id, expsr, XSIZE=xsize, YSIZE=ysize, $
                   PRINTONLY=printonly

;common esi_lwdchkspec_cmm

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'esi_lwdchkspec, esi, obj_id, expsr, XSIZE=, YSIZE=,'
    print, '      (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 900
  if not keyword_set( XYMNX ) then xymnx = [3800., 0.0, 8000., 1.0]

; Set exp
  indx = where(esi.flg_anly EQ 1 AND esi.mode EQ 1 AND $
               esi.obj_id EQ obj_id AND esi.type EQ 'OBJ', nindx)
  if not keyword_set(expsr) then expsr = 0L
  exp = indx[expsr]

;  Object Structure
  esiobj = xmrdfits(esi[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)
  objstr = esiobj[where(esiobj.flg_anly GT 1, nobj)]

; STATE

  state = {             $
            nplt: 8L, $
            npg: (nobj/8L) + (nobj MOD 8L NE 0), $
            nobj: nobj, $
            obj: objstr, $
            sig: fltarr(nobj,n_elements(objstr[0].var)), $
            obj_fil: esi[exp].obj_fil, $
            curpg: 0L, $
            mdrng: [5000., 8000.], $
            pos: [0.04,0.0,1.00,1.0], $ ; Plotting
            xymnx: xymnx, $
            tmpxy: fltarr(4), $
            xcurs: 0., $
            xpos: 0., $
            size: lonarr(2), $
            psfile: 0, $
            base_id: 0L, $      ; Widgets
            draw_id: 0L, $
            drawbase_id: 0L, $
            pg_id: 0L, $
            nspec_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            help_text_id: 0L $
          }

; Create Sigma array

  for i=0L,state.nobj-1 do begin
      a = where(state.obj[i].var LE 0., COMPLEMENT=b)
      if a[0] NE -1 then state.sig[i,a] = 0.
      if b[0] NE -1 then state.sig[i,b] = sqrt(state.obj[i].var[b])
  endfor

;    WIDGET
  base = WIDGET_BASE( title = 'esi_lwdchkspec: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='esi_lwdchkspec', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
;  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)

;;;;;;;;;;;;;;;;;;
; Pages
  state.pg_id = cw_field(toolbar, value=0L, $
                             /long, /return_events, xsize=3,$
                             title='Page')
  state.nspec_id = cw_field(toolbar, value=state.nplt, $
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
  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
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
;             'CMB/CMB - Zoom' $ 
;             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

; Realize
  WIDGET_CONTROL, base, /realize

  ; PLOT
  esi_lwdchkspec_Plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'esi_lwdchkspec', base

  !P.MULTI= [0,1,1]
  return
end

