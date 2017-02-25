;+ 
; NAME:
; x_setapergui   
;       Version 1.1
;
; PURPOSE:
;   Sets objects and apertures interactively using a GUI
;
; CALLING SEQUENCE:
;   x_setapergui, img, XSIZE=, YSIZE=, CLINE=, PEAKFRAC=, OBJSTR=
;
; INPUTS:
;   img   - 1D or 2D image
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   CLINE      - Center line for 2D image to 'extract' 1D [default: sz[1]/2]
;   PEAKFRAC=  - Fraction of peak 
;
; OPTIONAL OUTPUTS:
;  OBJSTR=  -- Object structure with the aperture filled in.  This will
;              include multiple objects if they are identified.
;
; COMMENTS:
;
; EXAMPLES:
;   x_setapergui, 'spec.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;  GETCOLOR (Coyote)
;  DJS_MEDIAN
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_setapergui_initcommon, flg, objstr
;
common x_setapergui_cmmn, sobjstr, nobj

  ;; Set sobjstr
  if flg MOD 4 GE 2 then begin
      sobjstr = objstr
      nobj = n_elements(sobjstr)
  endif else nobj = 0

return
end

;;;;
; Events
;;;;

pro x_setapergui_event, ev

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
                      1 : begin   ; Left Aperture
                          x_setapergui_mvaper, state, 1
                          if state.flg_aper MOD 2 NE 1 then $
                            state.flg_aper = state.flg_aper + 1
                      end 
                      4 : begin ; Right aperture
                          x_setapergui_mvaper, state, 2
                          if state.flg_aper MOD 4 LT 2 then $
                            state.flg_aper = state.flg_aper + 2
                      end 
                      else:
                  endcase
              end
              1 : ; Button release
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
;          if (state.flg_reg EQ 1 AND eventch NE 's') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
;              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
;              return
;          endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              ; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              ;
              'Z': xsetapergui_Zoom, state ; Zoom
              'z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'W': state.xymnx = state.svxymnx ; Reset the screen
              'H': x_helpwidg, state.help
              'A': xsetapergui_Auto, state   ; Automatically define an aperture
              'S': x_setapergui_SetSci, state
              'N': x_setapergui_NewObj, state
              'd': x_setapergui_DelObj, state
              'm': x_setapergui_MoveObj, state
              ; Aperture
              '1': x_setapergui_mvaper, state, 1
              '2': x_setapergui_mvaper, state, 2
              ; Obj
              '=': x_setapergui_nextobj, state, 1
              '-': x_setapergui_nextobj, state, 2
              ; Quit

              'Q': begin
;                  if state.flg_objstr EQ 0 then begin
                  tmp1 = state.aper
                  state.aper[0] = tmp1[0] < tmp1[1]
                  state.aper[1] = tmp1[0] > tmp1[1]
                  *state.pnt_aper = state.aper
;                  endif
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'x_setapergui: Not a valid key!' ; Nothing
          endcase
      end
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
; Update Plot
  xsetapergui_UpdatePlot, ev.top
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xsetapergui_UpdatePlot, base_id
  
common x_setapergui_cmmn

; State

  widget_control, base_id, get_uvalue=state, /no_copy

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  clr = getcolor(/load)

  plot, state.wave, state.fx, psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.5, $
    ycharsize=1.5
      
; Apertures

    if state.flg_aper MOD 2 EQ 1 then begin
        oplot, [state.aper[0], state.aper[0]], $
          [state.xymnx[1],state.xymnx[3]], $
          color=clr.green
    endif
    if state.flg_aper MOD 4 GT 2 then begin
        oplot, [state.aper[1], state.aper[1]], $
          [state.xymnx[1],state.xymnx[3]], $
          color=clr.green
    endif

; Objects
    if state.flg_objstr MOD 4 GE 2 then begin
        ;; Sci obj
        sci = where(strtrim(sobjstr.obj_id,2) EQ 'a')
        oplot, replicate(sobjstr[sci].ycen,2), [-1.e10, 1.e10], $
          color=clr.red
        ;; Current obj
        if sci[0] NE state.curobj then $
          oplot, replicate(sobjstr[state.curobj].ycen,2), [-1.e10, 1.e10], $
          color=clr.blue
        ;; Other obj
        case nobj of
            0: 
            1: 
            else: begin
                for j=0L,nobj-1 do begin
                    if j NE state.curobj AND j NE sci[0] then $
                      oplot, replicate(sobjstr[j].ycen,2), $
                      [-1.e10, 1.e10], color=clr.purple
                endfor
            end
        endcase
    endif
                
; Return State

  widget_control, base_id, set_uvalue=state, /no_copy

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro xsetapergui_Reset, base_id

; State

  widget_control, base_id, get_uvalue=state, /no_copy
;
  state.npix = n_elements(state.fx)

; Plotting
  state.xymnx = state.svxymnx

;
  widget_control, base_id, set_uvalue=state, /no_copy

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Zoom in

pro xsetapergui_Zoom, state

  ; Zoom
  state.xymnx[0] = xgetx_plt(state, /strct) - 40. ; left
  state.xymnx[2] = state.xymnx[0]+80. ; right

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; New Obj

pro x_setapergui_NewObj, state

common x_setapergui_cmmn

  ;; Center
  xval = round(xgetx_plt(state, /strct))
  xmn = 0 > (xval - 5L)
  xmx = (state.npix-1) < (xval + 5L)
  
  cent = x_centspln(state.wave[xmn:xmx], state.fx[xmn:xmx], state.peakfrac, $
                    /FORCE, EDGES=edges)

  ;; Set aper
  if edges[0] NE -1 then state.aper=edges $
    else state.aper = replicate(cent,2) + [-3.,3.]
  state.flg_aper = 3

  ;; Update objstr
  if not keyword_set( sobjstr ) then begin
      tmp = { specobjstrct } 
  endif else begin
      if size(sobjstr[0].wave[0], /type) EQ 5 then tmp = { dblsobjstrct } $
      else tmp = { specobjstrct } 
      tmp.slit_fil = sobjstr[0].slit_fil
      tmp.spec2d_fil = sobjstr[0].spec2d_fil
      tmp.UT = sobjstr[0].UT
      tmp.instr_strct = sobjstr[0].instr_strct
      tmp.field = sobjstr[0].field
  endelse

  tmp.flg_anly = 0

  if nobj GT 0 then sobjstr = [sobjstr,tmp] else begin
      sobjstr = tmp
      state.flg_objstr = state.flg_objstr + 2
  endelse

  state.cent = cent
  sobjstr[nobj].ycen = cent
  sobjstr[nobj].aper[0] = -(state.aper[0]-cent - 5L)
  sobjstr[nobj].aper[1] = state.aper[1]-cent + 5L 
  sobjstr[nobj].obj_id = x_objnumid(nobj)
  state.curobj = nobj
  nobj = nobj + 1

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Del Obj

pro x_setapergui_DelObj, state

common x_setapergui_cmmn

  ;; Center
  xval = round(xgetx_plt(state, /strct))
  mn = min(abs(xval-sobjstr.ycen), imn)
  msk = lonarr(nobj) + 1L
  msk[imn] = 0L

  ;; If SCI, set new sci
  if sobjstr[imn].obj_id EQ 'a' then flg_sci = 1 else flg_sci = 0 

  ;; Keep good ones
  gdobj = where(msk EQ 1)
  sobjstr = sobjstr[gdobj]
  nobj = nobj - 1

  ;; SCI
  if flg_sci EQ 1 then sobjstr[0].obj_id = 'a'
  
  ;; Change current obj
  x_setapergui_nextobj, state, 1
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Move Obj

pro x_setapergui_MoveObj, state

common x_setapergui_cmmn

  ;; Center
  xval = round(xgetx_plt(state, /strct))
  mn = min(abs(xval-sobjstr.ycen), imn)

  state.cent = xval
  sobjstr[imn].ycen = xval
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set Sci -- Sets obj as science obj

pro x_setapergui_SetSci, state

common x_setapergui_cmmn

  ;; Center
  xval = round(xgetx_plt(state, /strct))
  mn = min(abs(xval-sobjstr.ycen), imn)

  ;; Current science
  sci = where(sobjstr.obj_id EQ 'a')

  ;; Swap
  sobjstr[sci].obj_id = sobjstr[imn].obj_id
  sobjstr[imn].obj_id = 'a'
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Move Aper

pro x_setapergui_mvaper, state, flg

common x_setapergui_cmmn

  xval = xgetx_plt(state, /strct)
  case flg of
      1: begin
          state.aper[0] = xval 
          if state.flg_objstr MOD 4 GE 2 then $
            sobjstr[state.curobj].aper[0] = -(state.aper[0]-state.cent)
      end
      2: begin
          state.aper[1] = xval 
          if state.flg_objstr MOD 4 GE 2 then $
            sobjstr[state.curobj].aper[1] = state.aper[1]-state.cent 
      end
      else :
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Move Aper

pro x_setapergui_nextobj, state, flg

common x_setapergui_cmmn

  if state.flg_objstr MOD 4 GE 2 then begin
      case flg of
          1: begin
              state.curobj = state.curobj + 1
              if state.curobj GT (nobj-1) then state.curobj = 0
          end
          2: begin
              state.curobj = state.curobj - 1
              if state.curobj LT 0 then state.curobj = nobj-1
          end
          else :
      endcase
      state.aper = replicate(sobjstr[state.curobj].ycen,2) + $
        sobjstr[state.curobj].aper*[-1.,1.]
      state.cent = sobjstr[state.curobj].ycen
  endif

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Auto Set an aperture

pro xsetapergui_Auto, state

  ; Find the peak and set a region around it
  mx = max(state.fx, imx)
  rgmin = (imx - 8L) > 0
  rgmax = (imx + 8L) < (n_elements(state.fx)-1)

  ; Subtract off the sky
  medlvl = median(state.fx[rgmin:rgmax])

  ; Now Use splines
  mx = x_maxspln(state.wave[rgmin:rgmax], $
                         (state.fx[rgmin:rgmax]-medlvl), $
                         EDGES=edg, EDGVAL=state.peakfrac)
  state.aper = temporary(edg)
  state.flg_aper = 3

  ; Zoom
  rgmin = (imx - 20) > 0
  rgmax = (imx + 20) < (n_elements(state.fx)-1)
  state.xymnx[0] = state.wave[rgmin]
  state.xymnx[2] = state.wave[rgmax]

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

function x_setapergui, img, XSIZE=xsize, YSIZE=ysize, CLINE=cline, $
                    PEAKFRAC=peakfrac, OBJSTR=objstr

common x_setapergui_cmmn
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'aper = x_setapergui(img, XSIZE=,YSIZE=, CLINE=, PEAKFRAC=) [V1.0])'
    return, -1
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-400
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-400
  if not keyword_set( PEAKFRAC ) then    peakfrac = 0.15
  
; Convert to 1D as necessary
  sz = size(img)
  if sz[0] EQ 2 then begin
      if not keyword_set( CLINE ) then cline = sz[1]/2
      ydat = djs_median(img[cline-10:cline+10, *], 1)
  endif else ydat = img

;    Pointer for Aperture
  aper = fltarr(2)
  pnt_aper = PTR_NEW( aper )

; Check for objstr
  flg_objstr = 0
  if keyword_set( OBJSTR ) then flg_objstr = flg_objstr + 2
  if arg_present( OBJSTR ) then flg_objstr = flg_objstr + 1
  x_setapergui_initcommon, flg_objstr, objstr

  
;    STATE
  state = { fx: ydat, $
            wave: findgen(n_elements(ydat)), $
            npix: n_elements(ydat), $
            peakfrac: peakfrac, $  ; Fraction of peak to set aperture automatically
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            svxymnx: [0., $     ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      float(n_elements(ydat)-1), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $ ; ymax
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            pnt_aper: pnt_aper, $
            aper: fltarr(2), $
            flg_objstr: flg_objstr, $
            curobj: 0L, $
            flg_aper: 0, $
            cent: 0., $
            flg_reg: 0, $
            psym: 10, $
            title: '', $
            size: intarr(2), $
            help: strarr(50), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }
  
; Set aper
  if state.flg_objstr MOD 4 GE 2 then begin
      state.aper = replicate(sobjstr[0].ycen,2) + sobjstr[0].aper*[-1.,1.]
      state.flg_aper = 3
      state.cent = sobjstr[0].ycen
  endif
;    Title
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_setapergui', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['x_setaper', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  state.help[0] = '     Help Menu   '
  state.help[1] = 'LMB(1) -- Set Left aperture edge'
  state.help[2] = 'RMB(2) -- Set Right aperture edge'
  state.help[3] = 'Z -- Zoom in'
  state.help[4] = 'lrbt -- Set Left,Right,Bottom,Top'
  state.help[5] = 'z -- Set ymin to 0.'
  state.help[6] = 'W -- Reset screen'
  state.help[7] = '-------------------------'
  state.help[8] = 'N -- New Obj'
  state.help[9] = 'd -- Delete Obj'
  state.help[10] = 'm -- Move Obj'
  state.help[11] = '=/- -- Cycle through obj'
  state.help[12] = 'S -- Set science obj'
  state.help[13] = 'H -- Show this screen'
  state.help[14] = 'Q -- Quit '

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
;  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
  print, 'Press H for help'
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
  loadct, 0, /silent
  
; Update
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  xsetapergui_Reset, base
  xsetapergui_UpdatePlot, base

; Send to the xmanager
  xmanager, 'x_setapergui', base

  ; Pointer
  aper = *pnt_aper
  PTR_FREE, pnt_aper

  ;; Objstr
  if keyword_set( SOBJSTR ) then begin
      srt = sort(sobjstr.obj_id)
      objstr = sobjstr[srt]
      delvarx, sobjstr, nobj
      ;; Reset letters for objstr
      for i=1L,n_elements(objstr)-1 do objstr[i].obj_id = x_objnumid(i)
  endif

  return, aper
end
