;+ 
; NAME:
; x_editspec   
;   Version 1.1
;
; PURPOSE:
;    Allows user to reset values of a spectrum to reject bad regions
;    using a simple GUI
;
; CALLING SEQUENCE:
;   
;   x_editspec, wv, fx, var, [title], XSIZE=, YSIZE=, ISPEC=, 
;               /BLOCK, NEWVAR=, FLG=, NEWFLG=
;
; INPUTS
;   wv  -- Wavelength array  (expected to be a 2D array [npix, nspec])
;   fx  -- Flux array (expected to be a 2D array [npix, nspec])
;  var  -- Variance array (expected to be a 2D array [npix, nspec])
;  [title] -- List of object ID numbers
;
; RETURNS:
;
; OUTPUTS:
;  NEWVAR=  -- Updated variance array which has bad pixels (regions) masked
;  NEWFLG=  -- Flag used to specify whether to analyse the spectrum
;
; OPTIONAL KEYWORDS:
;   xsize   - Draw window xsize (pixels)
;   ysize   - Draw window ysize (pixels)
;  /BLOCK   - Block the window
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_editspec 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   03-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_editspec_initcommon, wv, fx, var
;
common x_editspec_cmm, $
  all_title, $
  all_wv, $
  all_fx, $
  all_var, $
  all_flg, $
  all_sz

; Set

  all_wv = wv
  all_fx = fx
  all_var = var 
  all_sz = size(all_wv)

end

;;;;
; Events
;;;;

pro x_editspec_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
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
                      4 : begin    ; Set reference line
                          x_editspec_delreg, state 
                          if state.flg_delreg NE 2 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state,  /no_copy
                              return
                          endif else state.flg_delreg = 0
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
;                  widget_control, state.ypos_id, set_value=state.ypos
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
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'N': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              'z': begin  ; Region
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'H': x_helpwidg, state.help
              else:  ; Nothing
          endcase
      end
      'NEXT' : begin
          state.curspec = state.curspec + 1
          if state.curspec GT state.nspec-1 then state.curspec = 0
          x_editspec_Reset, state
      end
      'PREV' : begin
          state.curspec = state.curspec - 1
          if state.curspec LT 0 then state.curspec = state.nspec-1
          x_editspec_Reset, state
      end
      'LIST' : begin
          state.curspec = ev.index
          x_editspec_reset, state
      end
      'NOANLY' : x_editspec_noanly, state
      'UNDO' : x_editspec_Reset, state
      'SAVE' : x_editspec_save, state
      'DONE' : begin
          x_editspec_save, state
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  x_editspec_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_editspec_UpdatePlot, state
  
; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  clr = getcolor(/load)

  plot, state.wave, state.fx, psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=2.0, $
    ycharsize=2.0

  oplot, state.wave, state.sig, psym=state.psym, color=clr.red

end

;;;;;;;;;;;;;;;;;;;;
;;;;;;;
;;  DEL REG
pro x_editspec_delreg, state

  clr = getcolor(/load)
  if state.flg_delreg EQ 0 then begin
      state.tmpxy[0] = xgetx_plt(state, /strct)
      state.flg_delreg = 1
      ; Plot
      oplot, [state.tmpxy[0],state.tmpxy[0]], $
        [state.xymnx[1],state.xymnx[3]], color=clr.blue
  end else begin
      state.tmpxy[2] = xgetx_plt(state, /strct)
      state.flg_delreg = 2
      mn = state.tmpxy[0] < state.tmpxy[2]
      mx = state.tmpxy[0] > state.tmpxy[2]
      pix = where(state.wave GT mn and state.wave LT mx, npix)
      print, 'x_editspec:   ', strtrim(min(pix),2),':', strtrim(max(pix),2)
      if npix NE 0 then state.sig[pix] = 0.
  endelse

  return
end


;;;;;;;;;;;;;;;;;;;;
;  SAVE
;;;;;;;;;;;;;;;;;;;;

pro x_editspec_save, state


common x_editspec_cmm

  ; zero pix
  pix = where(state.sig LE 0., npix)
  if npix NE 0 then all_var[pix,state.curspec] = 0.

end

;;;;;;;;;;;;;;;;;;;;
;  NO ANALYSIS
;;;;;;;;;;;;;;;;;;;;

pro x_editspec_noanly, state


common x_editspec_cmm

  ; zero pix
  all_flg[state.curspec] = 0

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_editspec_Reset, state


common x_editspec_cmm

  ;  Flux, wave, var
  state.fx = all_fx[*,state.curspec]
  state.wave = all_wv[*,state.curspec]
  state.sig = 0.
  a = where(all_var[*,state.curspec] GT 0)
  if a[0] NE -1 then state.sig[a] = sqrt(all_var[a,state.curspec])

  ; Title
  widget_control, state.title_id, set_list_select=state.curspec
  
  ; svxymnx
  mnwv = min(state.wave, max=mxwv)
  mxfx = max(state.fx)
  state.svxymnx = [mnwv, 0., mxwv, mxfx]

; Plotting
  state.xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_editspec, wv, fx, var, title, XSIZE=xsize, YSIZE=ysize, ISPEC=ispec, $
                BLOCK=block, NEWVAR=newvar, FLG=flg, NEWFLG=newflg

common x_editspec_cmm

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_editspec, wv, fx, var, [title], XSIZE=, YSIZE=, ISPEC= [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
;  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
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
;  if not keyword_set( XSIZE ) then    xsize = 1200
;  if not keyword_set( YSIZE ) then    ysize = 800


; Init common

  x_editspec_initcommon, wv, fx, var
  if keyword_set( TITLE ) then all_title = title $
  else all_title = sindgen(all_sz[1])
  if keyword_set( FLG ) then all_flg = flg

;    STATE
  state = { fx: fltarr(all_sz[1]), $
            wave: dblarr(all_sz[1]), $
            sig: fltarr(all_sz[1]), $
            flg_zoom: 0, $
            flg_delreg: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            xpos: 0.d, $
            ypos: 0.d, $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            psym: 10, $
            nspec: 0L, $
            curspec: 0L, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblspec_id: 0L, $
            draw_base_id: 0L, $
            title_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            error_msg_id: 0L $
          }

; Initial spectrum
  if keyword_set( ISPEC ) then state.curspec = ispec

; Nspec
  if all_sz[0] EQ 2 then state.nspec = all_sz[2] else state.nspec = 1

;    WIDGET
  base = WIDGET_BASE( title = 'x_editspec', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
  strlbl = strarr(10)
  strlbl = ['x_editspec', ' ', 'Ver 1.0']
  verslbl = WIDGET_TEXT(toolbar, value=strlbl, xsize=10, ysize=3)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'N -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = '1 -- Al III doublet'
  state.help[19] = '2 -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'C -- AODM'
  state.help[22] = 'q -- Quit '

;;;;;;;;;
;  Toolbar

; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
;  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

  state.title_id = widget_list(toolbar, value=all_title, UVALUE='LIST', $
                         ysize=5,xsize=20)
;;;;;;;;;;;;
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
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
  buttons1 = WIDGET_BASE( toolbar, /column, /base_align_center,$
                           /align_center)
  next = WIDGET_BUTTON(buttons1, value='NEXT',uvalue='NEXT', /align_right)
  prev = WIDGET_BUTTON(buttons1, value='PREV',uvalue='PREV', /align_right)
  buttons2 = WIDGET_BASE( toolbar, /column, /base_align_center,$
                           /align_center)
  undo = WIDGET_BUTTON(buttons2, value='UNDO',uvalue='UNDO', /align_right)
  save = WIDGET_BUTTON(buttons2, value='SAVE',uvalue='SAVE', /align_right)
  noanly = WIDGET_BUTTON(buttons2, value='NOANLY',uvalue='NOANLY', /align_right)
  if not keyword_set( all_flg ) then widget_control, noanly, sensitive=0
  done = WIDGET_BUTTON(buttons2, value='DONE',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Update
  x_editspec_Reset, state
  x_editspec_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set( BLOCK ) then xmanager, 'x_editspec', base, /no_block $
  else xmanager, 'x_editspec', base

; Return
  newvar = temporary(all_var)
  if arg_present(NEWFLG) then newflg = all_flg
  delvarx, all_fx, all_wave, all_sz
  if keyword_set( all_title ) then delvarx, all_title
  if keyword_set( all_flg ) then delvarx, all_flg

return
end
