;+ 
; NAME:
; x_msklayout   
;    Version 1.1
;
; PURPOSE:
;    Allows the user to interactively place slitmasks on a field of
;  targets and write out their positions and PA values
;
; CALLING SEQUENCE:
;   x_msklayout, targ, instr
;
; INPUTS:
;   img        - Image(s) for Masking
;   instr      - Name of instrument (LRIS, DEIMOS)
;
; RETURNS:
;
; OUTPUTS:
;   mask       -  Creates a fits table of pixel values to mask
;
; OPTIONAL KEYWORDS:
;  SECTRG  -- 2nd set of targets to plot (in alternate point style)
;  ASTAR   -- File containing the positions of Alignment stars
;  PVOBJ   -- File of objects previously observed
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_msklayout, img, 'SITe1', 'LCO-40'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Events
;;;;

pro x_msklayout_event, ev


  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'NEWBX': begin
          state.nbox = state.nbox + 1
          state.curbox = state.nbox-1
          WIDGET_CONTROL, state.pa_id, get_value=dum
          state.box[2,state.curbox] = dum
          x_msklayout_Plt, state
      end
      'PA': begin
          if state.nbox NE 0 then begin
              WIDGET_CONTROL, ev.id, get_value=dum
              state.box[2,state.curbox] = dum
              x_msklayout_Plt, state
          endif
      end
      'SAVE': begin
          close, /all
          openw, 1, state.outfil
          for qq=0, state.nbox-1 do begin
              ;; RA AND DEC
              x_radec, state.rac, state.decc, rad, decd
              decd = decd + state.box[1,qq]/60.
              rad = rad - state.box[0,qq]/cos(!pi*decd/180.)/60.
              x_radec, ras, decs, rad, decd, /flip
              ;; 
              case state.instr of 
                  'LRIS': $
                    printf, 1, state.box[0,qq], state.box[1,qq], state.box[2,qq],$
                    ras, decs, FORMAT='(2f8.4,1x,f8.3,1x,a12,1x,a12)'
                  'DEIMOS': $
                    printf, 1, state.box[0,qq], state.box[1,qq], state.box[2,qq],$
                    rad/15., decd, ras, decs, $
                    FORMAT='(2f8.4,1x,f8.3,1x,f12.9,1x,f12.8,a12,1x,a12)'
                  else: stop
              endcase
          endfor
          close, /all
      end
      'PRINT': begin
          x_msklayout_print, state
          print, 'x_msklayout: Output in xmask.ps'
      end
      'DONE': begin
          widget_control, ev.top, /destroy
          return
      end
      'DRAW' : begin
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Center
                          if state.nbox NE 0 then begin
                              state.box[0,state.curbox] = xgetx_plt(state, /strct) 
                              state.box[1,state.curbox] = xgety_plt(state, /strct) 
                              x_msklayout_Plt, state
                          endif
                      end 
                      2 : 
                      4 : begin
                          if state.nbox NE 0 then begin
                              state.curbox = state.curbox + 1
                              if state.curbox GE state.nbox then $
                                state.curbox = 0
                              ;; Update widget
                              widget_control, state.pa_id, $
                                set_value=state.box[2,state.curbox]
                              x_msklayout_Plt, state
                          endif
                      end
                  endcase
              end
              1 : ; Button release
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
              end
          endcase
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_msklayout_Plt, state
  

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)
  plot, [0], [0], $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    background=clr.white,  xcharsize=1.5,  ycharsize=1.5, $
    color=clr.black, /nodata, /isotropic

  ;; QSO
  oplot, [0], [0], color=clr.red, psym=7, symsize=1.5

  ;; Alignment stars
  if state.n3 NE 0 then $
    oplot, state.x3, state.y3, color=clr.orange, psym=6

  ;; Secondary
  if state.n2 NE 0 then $
    oplot, state.x2, state.y2, color=clr.blue, psym=2

  ;; Targets
  oplot, state.x1, state.y1, color=clr.green, psym=1

  ;; Previos obj
  if state.n4 NE 0 then begin
      oplot, state.x4, state.y4, color=clr.white, psym=1
      oplot, state.x4, state.y4, color=clr.red, psym=4
  endif

  ;; Boxes
  for qq=0L,state.nbox-1 do begin
      crnx = fltarr(5)
      crny = fltarr(5)
      case strtrim(state.instr,2) of 
          'LRIS': begin
              alph = atan(3.8/2.5)
              radius = sqrt(3.8^2 + 2.5^2)
              ;;
              crnx[0] = cos( alph + !pi*state.box[2,qq]/180.) 
              crny[0] = sin( alph + !pi*state.box[2,qq]/180.)
              crnx[1] = cos( !pi - alph + !pi*state.box[2,qq]/180.)
              crny[1] = sin( !pi - alph + !pi*state.box[2,qq]/180.)
              crnx[2] = cos( !pi + alph + !pi*state.box[2,qq]/180.)
              crny[2] = sin( !pi + alph + !pi*state.box[2,qq]/180.)
              crnx[3] = cos( -alph + !pi*state.box[2,qq]/180.)
              crny[3] = sin( -alph + !pi*state.box[2,qq]/180.)
              crnx[4] = crnx[0]
              crny[4] = crny[0]
          end
          'DEIMOS': begin
              ;; Assuming 16' x 5' field
              alph = atan(8/2.5)
              radius = sqrt(8.^2 + 2.5^2)
              ;;
              crnx[0] = cos( alph + !pi*state.box[2,qq]/180.) 
              crny[0] = sin( alph + !pi*state.box[2,qq]/180.)
              crnx[1] = cos( !pi - alph + !pi*state.box[2,qq]/180.)
              crny[1] = sin( !pi - alph + !pi*state.box[2,qq]/180.)
              crnx[2] = cos( !pi + alph + !pi*state.box[2,qq]/180.)
              crny[2] = sin( !pi + alph + !pi*state.box[2,qq]/180.)
              crnx[3] = cos( -alph + !pi*state.box[2,qq]/180.)
              crny[3] = sin( -alph + !pi*state.box[2,qq]/180.)
              crnx[4] = crnx[0]
              crny[4] = crny[0]
          end
          else: stop
      endcase
      ;; Plot
      oplot, radius*crnx+state.box[0,qq], $
        radius*crny+state.box[1,qq], color=clr.red, linestyle=2
      ;; Curbox
      if qq EQ state.curbox then $
        oplot, radius*crnx+state.box[0,qq], $
        radius*crny+state.box[1,qq], color=clr.black, linestyle=2
      if strtrim(state.instr,2) EQ 'DEIMOS' then begin
          oplot, [radius*crnx[0]+state.box[0,qq]], $
            [radius*crny[0]+state.box[1,qq]], color=clr.purple, psym=6
          oplot, [radius*crnx[3]+state.box[0,qq]], $
            [radius*crny[3]+state.box[1,qq]], color=clr.purple, psym=6
      endif
  endfor
end

pro x_msklayout_print, state
  ;; Open file
  x_psopen, 'xmask.ps', /maxs
  state.psfile = 1
  x_msklayout_Plt, state
  state.psfile = 0
  x_psclose

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_msklayout, targ, instr, SECTRG=sectrg, OUTFIL=outfil, $
                 RAC=rac, DECC=decc, INFIL=infil, ASTAR=astar, $
                 PVOBJ=pvobj


;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'x_msklayout, targ, instr, SECTRG=, OUTFIL= (v1.0)'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( XSIZE ) then    xsize = ysize
  if not keyword_set( RAC ) then rac = '01:01:01.0'
  if not keyword_set( DECC ) then decc = '01:01:01.0'
  if not keyword_set( OUTFIL ) then outfil = 'box.dat'

; Open targets

  case strtrim(instr,2) of 
      'LRIS': readcol, targ, xv, yv, /silent
      'DEIMOS': readcol, targ, xv, yv, id, /silent
      else: stop
  endcase
  if keyword_set( SECTRG ) then begin
      readcol, sectrg, xv2, yv2, id2, /silent
      n2 = n_elements(xv2)
  endif else begin
      xv2 = 0.
      n2 = 0
      yv2 = 0.
  endelse
  if keyword_set( ASTAR ) then begin
      readcol, astar, xv3, yv3
      n3 = n_elements(xv3)
  endif else begin
      xv3 = 0.
      n3 = 0
      yv3 = 0.
  endelse
  if keyword_set( PVOBJ ) then begin
      case strtrim(instr,2) of 
          'LRIS': begin
              readcol, pvobj, xv4, yv4
              n4 = n_elements(xv4)
          end
          'DEIMOS': begin
              readcol, pvobj, pv_id, format='L'
              n4 = n_elements(pv_id)
              xv4 = fltarr(n4)
              yv4 = fltarr(n4)
              for ii=0L,n4-1 do begin
                  a = where(pv_id[ii] EQ ID, na)
                  if na EQ 1 then begin
                      xv4[ii] = xv[a]
                      yv4[ii] = yv[a]
                  endif else begin ;; Secondary
                      if n2 NE 0 then begin
                          a = where(pv_id[ii] EQ ID2, na)
                          if na EQ 1 then begin
                              xv4[ii] = xv2[a]
                              yv4[ii] = yv2[a]
                          endif else stop
                      endif
                  endelse
              endfor
          end
      endcase
  endif else begin
      xv4 = 0.
      n4 = 0
      yv4 = 0.
  endelse

      

;    STATE

  state = { nbox: 0L, $
            box: fltarr(3,100), $
            curbox: 0L, $
            instr: instr, $
            outfil: outfil, $
            rac: rac, $
            decc: decc, $
            x1: xv, $
            y1: yv, $
            n1: n_elements(xv), $
            n2: n2, $
            x2: xv2, $
            y2: yv2, $
            n3: n3, $
            x3: xv3, $
            y3: yv3, $
            n4: n4, $
            x4: xv4, $
            y4: yv4, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plottinxv
            svxymnx: [ min(xv)-0.01*abs(max(xv)-min(xv)), $ ; xmin
                      min(yv)-0.01*abs(max(yv)-min(yv)), $ ; ymin
                      max(xv)+0.01*abs(max(xv)-min(xv)), $ ; xmax
                      max(yv)+0.01*abs(max(yv)-min(yv))], $  ; ymax
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            help: strarr(20), $
            psfile: 0, $
            size: intarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            pa_id: 0L, $
            draw_id: 0L $
          }

  state.xymnx = state.svxymnx

; INFIL
  if keyword_set( INFIL ) then begin
      readcol, infil, xb, yb, pab, /silent, FORMAT='F,F,F'
      state.nbox = n_elements(xb)
      for ii=0L,state.nbox-1 do begin
          state.box[0,ii] = xb[ii]
          state.box[1,ii] = yb[ii]
          state.box[2,ii] = pab[ii]
      endfor
  endif
      

;    WIDGET
  base = WIDGET_BASE( title = 'x_msklayout: ID Stars', /column, $
                      xoffset=100,yoffset=300)
  state.base_id = base
  
;      Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
;        Version + Name
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_msklayout', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Drawing
  drawbase = WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_id = widget_draw(drawbase, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Buttons1
  butbase = widget_base(toolbar, /row, /align_center)
  flipx = WIDGET_BUTTON(butbase, value='NEWBX',uvalue='NEWBX')
  print = WIDGET_BUTTON(butbase, value='PRINT',uvalue='PRINT')
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  state.pa_id = cw_field(butbase, title='PA:', value=0., /floating, $
                           xsize=5, /return_events, uvalue='PA')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
; Realize
  WIDGET_CONTROL, base, /realize

  x_msklayout_Plt, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_msklayout', base

  return
end

