;+ 
; NAME:
; x_splot   
;     Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   x_splot, [xdat], ydat, XSIZE=, YSIZE=, TITLE=, XTWO=, YTWO=, PSYM2=,
;      /BLOCK, PSYM1=, YMNX=, XTHR=, YTHR=, PSYM3=
;
; INPUTS:
;   [xdat]     - x values (optional)
;   ydat       - Values 
;
; RETURNS: 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   PSYM1  -- Plot symbol for array1
;   TITLE
;   YTWO   -- 2nd array to plot [If without XTWO, must have same dimension
;              as xdat]
;   XTWO   -- 2nd array of x values to plot [requires YTWO]
;   PSYM2  -- Plot symbol for array 2 [default: 10]
;   YTHR   -- 3nd array to plot [If without XTHR, must have same dimension
;              as xdat]
;   XTHR   -- 3nd array of x values to plot [requires YTHR]
;   PSYM3  -- Plot symbol for array 3 [default: 10]
;   YFOU   -- 4th array to plot [If without XFOU, must have same dimension
;              as xdat]
;   XFOU   -- 4th array of x values to plot [requires YFOU]
;   PSYM4  -- Plot symbol for array 4 [default: 10]
;   YFIV   -- 5th array to plot [If without XFIV, must have same dimension
;              as xdat]
;   XFIV   -- 5th array of x values to plot [requires YFIV]
;   PSYM5  -- Plot symbol for array 5 [default: 10]
;   YSIX   -- 6th array to plot [If without XSIX, must have same dimension
;              as xdat]
;   XSIX   -- 6th array of x values to plot [requires YSIX]
;   PSYM6  -- Plot symbol for array 6 [default: 10]
;
;   YSEV   -- 7th array to plot [If without XSEV, must have same dimension
;              as xdat]
;   XSEV   -- 7th array of x values to plot [requires YSEV]
;   PSYM7  -- Plot symbol for array 7 [default: 10]
;   YEIG   -- 8th array to plot [If without XEIG, must have same dimension
;              as xdat]
;   XEIG   -- 8th array of x values to plot [requires YEIG]
;   PSYM8  -- Plot symbol for array 8 [default: 10]
;   COLOR1 -- Plot color for array1
;   COLOR2 -- Plot color for array2
;   COLOR3 -- Plot color for array3
;   COLOR4 -- Plot color for array4
;   COLOR5 -- Plot color for array5
;   COLOR6 -- Plot color for array6
;   COLOR7 -- Plot color for array7
;   COLOR8 -- Plot color for array8
;   LGND   -- String array of of legend labels.  X_SPLOT will grab
;             colors and symbols itself.
;   /BLOCK
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_splot, y
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;  LEGEND
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP (modified from x1dfit)
;   24-Nov-2001 Added ytwo
;    9-Aug-2007 Added xmnx keyword JFH
;   13-Jun-2011 Added yfiv, ysix, ysev, yeig MMK 
;   19-Jul-2011 Added color and legend options  MMK
;   21-Feb-2014 Handle 2+ portrait monitors, KLC
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_splot_initcomm

common x_splot_common, xdat, ydat, y2, x2, y3, x3, y4, x4,y5, x5,y6, x6,y7,x7,y8,x8,mask, plot_extra

end

;;;;
; Events
;;;;

pro x_splot_event, ev

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
                      1 : state.xymnx[0] = xgetx_plt(state, /strct) ; left
                      2 : state.xymnx[3] = xgety_plt(state, /strct) ; top
                      4 : state.xymnx[2] = xgetx_plt(state, /strct) ; right
                      else :
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
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another s !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'N': state.xymnx[3] = 1.1
              'h': state.psym1 = 10
              'H': x_helpwidg, state.help
              'z': begin  ; Zoom
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, $
                        set_uvalue = state, /no_copy
                      return
                  endif
              end
              'm': begin
                 xmax = xgetx_plt(state, /str)
                 startall = WHERE(state.maskbracket EQ 0.0)
                 startodd = WHERE( (startall MOD 2) EQ 1)
                 IF xmax LE state.maskbracket[startall[startodd[0]]-1] AND $
                    state.maskbracket[startall[startodd[0]]-1] GT 0.0 THEN BEGIN
                    PRINT, "ERROR: max is less than min, max not stored" 
                    printbracket = WHERE(state.maskbracket NE 0.0, nprint)
                    IF nprint NE 0 THEN BEGIN
                       printindex = LINDGEN(printbracket[SIZE(printbracket, /DIM)]+1)
                       PRINT, state.maskbracket[printindex]
                    ENDIF

                 ENDIF ELSE BEGIN
                    state.maskbracket[startall[startodd[0]]] = xmax

                 ENDELSE                 
                 print, "xmax= ",xmax
              end
              'n': begin
                 xmin = xgetx_plt(state, /str)
                 startall = WHERE(state.maskbracket EQ 0.0)
                 starteven = WHERE( (startall MOD 2) EQ 0)
                 IF xmin GE state.maskbracket[startall[starteven[0]]+1] AND $
                    state.maskbracket[startall[starteven[0]]+1] GT 0.0 THEN BEGIN
                    PRINT, "ERROR: min is greater than max, min not stored" 
                    printbracket = WHERE(state.maskbracket NE 0.0, nprint)
                    IF nprint NE 0 THEN BEGIN
                       printindex = LINDGEN(printbracket[SIZE(printbracket, /DIM)]+1)
                       PRINT, state.maskbracket[printindex]
                    ENDIF

                 ENDIF ELSE BEGIN
                    state.maskbracket[startall[starteven[0]]] = xmin
                 ENDELSE
                 print, "xmin= ",xmin
              end
             
              'M': begin  ; Median/mean stats
                  xsplot_stats, state
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              '}': xsplot_pan, state, /noy
              '{': xsplot_pan, state, /left, /noy
              ']': xsplot_pan, state
              '[': xsplot_pan, state, /left
              ' ': state.showval = 1 ; Show value
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'W': state.xymnx = state.svxymnx ; Reset the screen
              'P': xsplot_psfile, state  ; Send Current screen to ps file
              'q': begin
                  if keyword_set( X2 ) then delvarx, x2
                  if keyword_set( Y2 ) then delvarx, y2
                  if keyword_set( XDAT ) then delvarx, xdat
                  if keyword_set( YDAT ) then delvarx, ydat
                  widget_control, ev.top, /destroy
                  return
              end
              else:  ; Nothing
          endcase
      end
      'DONE' : begin
         widget_control, ev.top, /destroy
         return
      end
   endcase

; Update Plot
  xsplot_UpdatePlot, state
;
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xsplot_UpdatePlot, state
  
common x_splot_common

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  color = getcolor(/load)

  ;; Where
  a = where(xdat LE state.xymnx[2] AND xdat GE state.xymnx[0], na)
  plot, [0], [0], psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    title=state.title, background=color.white, $
    xcharsize=1.5, $
    ycharsize=1.5, $
    color=color.black, /nodata, _EXTRA=plot_extra; xlog=tag_exist(plot_extra,'XLOG')

  if na NE 0 then begin
     oplot, xdat[a], ydat[a], psym=state.psym1, color=state.color1
  end

; y2
  if state.flg_y2 then begin
      oplot, x2, y2, color=state.color2, psym=state.psym2
  endif

; y3
  if state.flg_y3 then begin
      oplot, x3, y3, color=state.color3, psym=state.psym3
  endif

; y4
  if state.flg_y4 then begin
      oplot, x4, y4, color=state.color4, psym=state.psym4
   endif

; y5
  if state.flg_y5 then begin
      oplot, x5, y5, color=state.color5, psym=state.psym5
   endif


; y6
  if state.flg_y6 then begin
      oplot, x6, y6, color=state.color6, psym=state.psym6
   endif


; y7
  if state.flg_y7 then begin
      oplot, x7, y7, color=state.color7, psym=state.psym7
  endif

; y8
  if state.flg_y8 then begin
      oplot, x8, y8, color=state.color8, psym=state.psym8
   endif

  if state.flg_leg then begin
     sizeleg = (size(state.legend, /DIM))[0] > 1
     psym = [state.psym1, state.psym2, state.psym3, state.psym4,   $
             state.psym5, state.psym6, state.psym7, state.psym8]
     psym = psym[LINDGEN(sizeleg)]
     replace = WHERE(psym EQ 10, nrep)
     IF nrep NE 0  THEN psym[replace] = -3
     colors = [state.color1, state.color2, state.color3, state.color4,   $
               state.color5, state.color6, state.color7, state.color8]
     colors = colors[LINDGEN(sizeleg)]
     LEGEND, state.legend, psym=psym, colors=colors, outline_color=color.black, textcolors=colors, /top_legend, /center_legend, _extra=plot_extra
  endif

; Show val
  if state.showval EQ 1 then begin
      xyouts, 0.7, 0.97, $
        'Value = '+string(xgetx_plt(state, /strct), $
                          xgety_plt(state, /strct), $
                          format='(2g14.7)'), /NORMAL, charsize=1.5, $
              color=color.black
      state.showval = 0
      print, 'Value = '+string(xgetx_plt(state, /strct), $
                               xgety_plt(state, /strct), $
                               format='(2g14.7)') 
  endif
      
end

;;;;;;;;;;;;;;;;;;;;
;  Pan
;;;;;;;;;;;;;;;;;;;;

pro xsplot_pan, state, LEFT=left, NOY=noy

common x_splot_common

  ; Size of screen
  sz = state.xymnx[2]-state.xymnx[0]

  ; Right or left
  if not keyword_set( LEFT ) then begin ; RIGHT
      state.xymnx[0] = state.xymnx[2] - 0.1*sz
      state.xymnx[2] = state.xymnx[2] + 0.9*sz
  endif else begin
      state.xymnx[2] = state.xymnx[0] + 0.1*sz
      state.xymnx[0] = state.xymnx[0] - 0.9*sz
  endelse

  ; Reset top and bottom
  if not keyword_set( NOY ) then begin
      gdwv = where(xdat GT state.xymnx[0] AND xdat LT state.xymnx[2], ngd)
      if ngd GT 0 then begin
          mn = min(ydat[gdwv], max=mx)
          state.xymnx[1] = mn - (mx-mn)*0.05
          state.xymnx[3] = mx + (mx-mn)*0.05
      endif
  endif

  return
end

      

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;
pro xsplot_Reset, state
common x_splot_common
  state.ntot = n_elements(ydat)
; Plotting
  state.xymnx = state.svxymnx
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro xsplot_psfile, state

; Device
  device, get_decomposed=svdecomp

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  xsplot_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0

end

;;;;;;;;;;;;;;;;;;;;
;  Stats
;;;;;;;;;;;;;;;;;;;;
pro xsplot_stats, state
common x_splot_common
  if state.flg_stats EQ 0 then begin
      state.sv_statx = round(xgetx_plt(state,/strct))
      state.flg_stats = 1
  endif else begin
      x2 = round(xgetx_plt(state,/strct))
      x1 = x2 < state.sv_statx
      x2 = x2 > state.sv_statx
      djs_iterstat, ydat[x1:x2], sigrej=4., mean=mn, median=med, sigma=sig
      print, 'xsplot: Stats (4 sig rejection) -- ;'
      print, ' median = ', med
      print, ' mean = ', mn
      print, ' sigma = ', sig
      state.flg_stats = 0
  endelse

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

pro x_splot, xin, yin, XSIZE=xsize, YSIZE=ysize, TITLE=title, YTWO=ytwo, $
             PSYM2=psym2, BLOCK=block, XTWO=xtwo, PSYM1=psym1, XYOFF=xyoff, $
             XMNX=XMNX, YMNX=ymnx, XTHR=xthr, YTHR=ythr, $
             PSYM3=psym3, XFOU=XFOU, YFOU=YFOU, PSYM4=PSYM4, XFIV=XFIV, $
             YFIV=YFIV, PSYM5=PSYM5, XSIX=XSIX,YSIX=YSIX, $
             PSYM6=PSYM6, XSEV=XSEV, YSEV=YSEV, XEIG=XEIG, YEIG=YEIG,   $
             PSYM7=PSYM7, PSYM8=psym8, BRACKET=bracket, _EXTRA=extra, $
             COLOR1=color1, color2=color2, COLOR3=color3, color4=COLOR4, $
             COLOR5=color5, COLOR6=color6, COLOR7=color7, COLOR8=color8,  $
             LGND=lgnd


common x_splot_common
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_splot, [xdat], ydat, YTWO=ytwo, XSIZE=, YSIZE=, TITLE= '
    print, '            XTWO=, /BLOCK, PSYM_Y2=, PSYM1=, YMNX=, YTHR= [v1.1]'
    return
  endif 

; Set xdat, ydat

  x_splot_initcomm
  if keyword_set(EXTRA) then plot_extra = extra else plot_extra = { dum: 0}

  if keyword_set( yin ) then begin
      if n_elements(xin) NE n_elements(yin) then begin
          print, 'x_splot: Wrong array sizes!'
          return
      endif
      xdat = xin
      ydat = yin
  endif else begin
      xdat = findgen( n_elements(xin) )
      ydat = xin
  endelse

  if keyword_set( YTWO ) then begin
      y2 = ytwo 
      if keyword_set( XTWO ) then x2 = xtwo else x2 = xdat
  endif

  if keyword_set( YTHR ) then begin
      y3 = ythr 
      if keyword_set( XTHR ) then x3 = xthr else x3 = xdat
  endif

  if keyword_set( YFOU ) then begin
      y4 = yfou 
      if keyword_set( XFOU ) then x4 = xfou else x4 = xdat
  endif
  

  if keyword_set( YFIV ) then begin
     y5 = yfiv
     if keyword_set( XFIV ) then x5 = xfiv else x5 = xdat
  endif
  
  if keyword_set( YSIX ) then begin
     y6 = ysix
     if keyword_set( XSIX ) then x6 = xsix else x6 = xdat
  endif
  
  
  if keyword_set( YSEV ) then begin
     y7 = ysev
     if keyword_set( XSEV ) then x7 = xsev else x7 = xdat
  endif
  
  if keyword_set( YEIG ) then begin
     y8 = yeig
     if keyword_set( XEIG ) then x8 = xeig else x8 = xdat
  endif
    
 if keyword_set( bracket ) then begin
    print, "keyword set bracket"
    mask = DBLARR(2,50)
    mask[*] = 1.0e9
 endif
    
    
;  Optional Keywords
 color = getcolor(/load)
  device, get_screen_size=ssz
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
  if not keyword_set( PSYM2 ) then    psym2 = 0
  if not keyword_set( PSYM3 ) then    psym3 = 1
  if not keyword_set( PSYM4 ) then    psym4 = 10
  if not keyword_set( PSYM5 ) then    psym5 = 3
  if not keyword_set( PSYM6 ) then    psym6 = 6
  if not keyword_set( PSYM7 ) then    psym7 = 6
  if not keyword_set( PSYM8 ) then    psym8 = 6
  if not keyword_set( PSYM1 ) then    psym1 = 0
  if not keyword_set( XYOFF ) then    xyoff = [10L,20L]
  if not keyword_set( COLOR1 ) then   color1 = color.black
  if not keyword_set( COLOR2 ) then   color2 = color.red
  if not keyword_set( COLOR3 ) then   color3 = color.blue
  if not keyword_set( COLOR4 ) then   color4 = color.cyan
  if not keyword_set( COLOR5 ) then   color5 = color.limegreen
  if not keyword_set( COLOR6 ) then   color6 = color.magenta
  if not keyword_set( COLOR7 ) then   color7 = color.orange
  if not keyword_set( COLOR8 ) then   color8 = color.dodgerblue
  if not keyword_set( LGND )   then   lgnd = [' ']

;    STATE
  state = { ntot: n_elements(ydat), $
            reg: fltarr(100,2), $
            flg_y2: 0, $
            flg_y3: 0, $
            flg_y4: 0, $
            flg_y5: 0, $
            flg_y6: 0, $
            flg_y7: 0, $
            flg_y8: 0, $
            flg_leg: 0, $
            psym1: psym1, $
            psym2: psym2, $
            psym3: psym3, $
            psym4: psym4, $
            psym5: psym5, $
            psym6: psym6, $
            psym7: psym7, $
            psym8: psym8, $
            color1: color1, $
            color2: color2, $
            color3: color3, $
            color4: color4, $
            color5: color5, $
            color6: color6, $
            color7: color7, $            
            color8: color8, $
            flg_zoom: 0, $
            flg_flip: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            svxymnx: [ min(xdat)-0.01*abs(max(xdat)-min(xdat)), $ ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      max(xdat)+0.01*abs(max(xdat)-min(xdat)), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $  ; ymax
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            flg_stats: 0, $
            sv_statx: 0L, $
            help: strarr(20), $
            showval: 0, $
            psym: psym1, $
            psfile: 0, $
            title: '', $
            size: intarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L, $
            maskbracket: DBLARR(2,50),  $
            legend: lgnd   $
          }
  
; YMNX
  if keyword_set( YMNX ) then begin
      state.svxymnx[1] = ymnx[0]
      state.svxymnx[3] = ymnx[1]
  endif

  state.svxymnx[0] = min(xdat)
  state.svxymnx[2] = max(xdat)

; XMNX
  if keyword_set( XMNX ) then begin
      state.svxymnx[0] = xmnx[0]
      state.svxymnx[2] = xmnx[1]
  endif

      
; Flag

  if keyword_set( YTWO ) then state.flg_y2 = 1
  if keyword_set( YTHR ) then state.flg_y3 = 1
  if keyword_set( YFOU ) then state.flg_y4 = 1
  if keyword_set( YFIV ) then state.flg_y5 = 1
  if keyword_set( YSIX ) then state.flg_y6 = 1
  if keyword_set( YSEV ) then state.flg_y7 = 1
  if keyword_set( YEIG ) then state.flg_y8 = 1
  if lgnd[0] NE ' ' then state.flg_leg = 1
  
; Set xvxymnx[0,2]


;    Title
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_splot', /column, xoffset=xyoff[0], $
                    yoffset=xyoff[1])
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)

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
;  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

;        Help
  state.help[0] = '     Help Menu   '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'N -- Set ymax to 1.1'
  state.help[9] = 'h -- Switch to histogram mode'
  state.help[10] = 'H -- Show this screen'
  state.help[11] = 'W -- Reset screen'
  state.help[11] = 'M -- stats'
  state.help[12] = 'q -- Quit '

  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
  
; Update
  xsplot_Reset, state
  xsplot_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set( BLOCK ) then xmanager, 'x_splot', base, /no_block $
    else xmanager, 'x_splot', base


return
end
