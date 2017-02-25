;+ 
; NAME:
; cldy_plot   
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a CLOUDY grid and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   cldy_plot, ingrid, /INFIL
;
; INPUTS:
;   flux  - Flux array (or FITS file)
;   [inflg]  - Sigma array (or FITS file)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   INFIL       -   set to indicate that the 'ingrid' is a file, and should be 
;                   read in
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_plot, 'grid.fits', /INFIL
;   
;   cldy_plot, grid
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;  -- Built by GEP
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro cldy_plot_initcommon
;
common cldy_plot_lines, $
  flg_lines, $
  lines, $
  zabs

  flg_lines = 0
  zabs = 0.

end

;;;;
; Events
;;;;

pro cldy_plot_event, ev

common cldy_plot_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval
  case uval of
      'PICVARS': begin
          state.usevar = state.picvars[long(ev.index)]
      end
      'U': begin
          state.ou = state.u
          state.u = float(state.urange[long(ev.index)])
      end
      'FeH': begin
          state.ofeh = state.feh
          state.feh = float(state.fehrange[long(ev.index)])
      end
      'nH': begin
          state.onh = state.nh
          state.nh = float(state.nhrange[long(ev.index)])
      end
      'NHI': begin
          state.onhi = state.nhi
          state.nhi = float(state.nhirange[long(ev.index)])
      end
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
              end
              1 : begin ; Button Release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state,  /no_copy
                  return
              end
              2 : begin ; Motion event
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'IONWIN': begin
          didx = widget_info(state.ionwin_id,/list_select)
          idx = state.ionwin[didx]
          if idx[0] NE -1 then begin
              first = lindgen(n_elements(didx))
              second = lindgen(n_elements(didx))
              state.nions=n_elements(didx)
              for ii=0L,n_elements(didx)-1 do begin
                  case idx[ii] of 
                      'HI': begin
                          first[ii] = 1
                          second[ii] = 1
                      end
                      'SiII': begin
                          first[ii] = 14
                          second[ii] = 2
                      end
                      'SiIV': begin
                          first[ii] = 14
                          second[ii] = 4
                      end
                      'CII': begin
                          first[ii] = 6
                          second[ii] = 2
                      end
                      'CIV': begin
                          first[ii] = 6
                          second[ii] = 4
                      end
                      'AlIII': begin
                          first[ii] = 13
                          second[ii] = 3
                      end
                      'OI': begin
                          first[ii] = 8
                          second[ii] = 1
                      end
                      'FeII': begin
                          first[ii] = 26
                          second[ii] = 2
                      end
                      else: stop
                  endcase
              endfor
              for jj=0, n_elements(didx)-1 do begin
                  state.ions[0,jj] = first[jj]
                  state.ions[1,jj] = second[jj]
              endfor
          endif
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          case eventch of
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'cldy_plot: Not a valid key!' ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  cldyplot_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro cldyplot_UpdatePlot, state
  
common cldy_plot_lines

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  case state.usevar of
      'U': begin
          ii = where(state.grid.nh EQ state.nh and state.grid.nhi EQ state.nhi $
                     and state.grid.feh EQ state.feh,nii)
          if nii EQ 0 then begin
              state.nh = state.onh
              state.nhi = state.onhi
              state.feh = state.ofeh
              state.u = state.ou
              ii = where(state.grid.nh EQ state.nh and state.grid.nhi EQ state.nhi $
                     and state.grid.feh EQ state.feh,nii)
          endif
          xmin = min(state.grid[ii].u,max=xmax)
          xplt = state.grid[ii].u
          state.xtitle = 'U'
      end
      'nH': begin
          ii = where(state.grid.u EQ state.u and state.grid.nhi EQ state.nhi $
                     and state.grid.feh EQ state.feh,nii)
          if nii EQ 0 then begin
              state.nh = state.onh
              state.nhi = state.onhi
              state.feh = state.ofeh
              state.u = state.ou
              ii = where(state.grid.u EQ state.u and state.grid.nhi EQ state.nhi $
                     and state.grid.feh EQ state.feh,nii)
          endif
          xmin = min(state.grid[ii].nh,max=xmax)
          xplt = state.grid[ii].nh
          state.xtitle = 'NH'
      end
      'NHI': begin
          ii = where(state.grid.nh EQ state.nh and state.grid.u EQ state.u $
                     and state.grid.feh EQ state.feh,nii)
          if nii EQ 0 then begin
              state.nh = state.onh
              state.nhi = state.onhi
              state.feh = state.ofeh
              state.u = state.ou
              ii = where(state.grid.nh EQ state.nh and state.grid.u EQ state.u $
                     and state.grid.feh EQ state.feh,nii)
          endif
          xmin = min(state.grid[ii].nhi,max=xmax)
          xplt = state.grid[ii].nhi
          state.xtitle = 'NHI'
      end
      'FeH': begin
          ii = where(state.grid.nh EQ state.nh and state.grid.u EQ state.u $
                     and state.grid.nhi EQ state.nhi,nii)
          if nii EQ 0 then begin
              state.nh = state.onh
              state.nhi = state.onhi
              state.feh = state.ofeh
              state.u = state.ou
              ii = where(state.grid.nh EQ state.nh and state.grid.u EQ state.u $
                     and state.grid.nhi EQ state.nhi,nii)
          endif
          xmin = min(state.grid[ii].feh,max=xmax)
          xplt = state.grid[ii].feh
          state.xtitle = 'FEH'
      end
      else: stop
  endcase


  xcolors = x_setclrs()
  clr = getcolor(/load)

  dumx = fltarr(2)
  dumy = fltarr(2)
  dumx[0] = -99.
  dumy[0] = -99.

  for q = 0, state.nions-1 do begin
      ;   Abund
      getabnd, elm, state.ions[0,q], abnd, flag=1
      yplt = state.nhi - state.grid[ii].X[1,1]-12. + abnd + $
          state.grid[ii].X[state.ions[0,q],state.ions[1,q]]
      ;  Metallicity
      if state.ions[0,q] NE 1 then yplt = yplt + state.FeH
      ;
      if q EQ 2 then psym = 1 else psym=q+1
      psym = psym<7

      if q EQ 0 then plot, xplt, yplt, psym=psym, xrange=[xmin,xmax], $
          yrange=[state.ymnx[0],state.ymnx[1]], color=clr.black, $
          background=clr.white, ystyle=1, xstyle=1, $
          xtitle=state.xtitle, ytitle='N(X)' $
      else oplot, xplt, yplt, psym=psym, color=xcolors[q]

      ;  Label
      case state.ions[1,q] of
          1: nm = elm+'I'
          2: nm = elm+'II'
          3: nm = elm+'III'
          4: nm = elm+'IV'
          5: nm = elm+'V'
          6: nm = elm+'VI'
          7: nm = elm+'VII'
          else: stop
      endcase

      dumx[1] = xmax-0.5
      dumy[1] = state.ymnx[0]+3-q*0.3
      oplot, dumx, dumy, psym=psym, color=xcolors[q]
      xyouts, xmax-0.4, state.ymnx[0]+2.9-q*0.3, nm, charsize=1.5, $
        color=xcolors[q]
  endfor


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;  MAIN  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro cldy_plot, ingrid, INFIL=infil, NHI=nhi, U=u, HDEN=hden, FEH=feh

;common x_specplot_lines

;
  if  N_params() LT 1  then begin 
    print,'Syntax - cldy_plot, ingrid, /infil'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200

; Read in the Data
  if keyword_set( INFLG ) then begin
      grid = xmrdfits(ingrid,1,strctyp='strctcldy')
  endif else grid = ingrid

; Values
  s = sort(grid.nhi)
  r = uniq(grid[s].nhi)
  nhirange = strtrim(string(grid[s[r]].nhi),2)
  s = sort(grid.nh)
  r = uniq(grid[s].nh)
  nhrange = strtrim(string(grid[s[r]].nh),2)
  s = sort(grid.feh)
  r = uniq(grid[s].feh)
  fehrange = strtrim(string(grid[s[r]].feh),2)
  s = sort(grid.u)
  r = uniq(grid[s].u)
  urange = strtrim(string(grid[s[r]].u),2)


; Init common

  cldy_plot_initcommon

;    STATE
  state = { grid: grid, $
            nhi: 17., $
            nh: -2., $
            onhi: 17., $
            onh: -2., $
            nions: 8L, $
            U: -3., $
            oU: -3., $
            FeH: 0., $
            oFeH: 0., $
            ions: [[1,1],[14,2],[14,4],[6,2],[6,4],[13,3],[8,1],[26,2]], $
            psfile: 0, $ ; Postscript
            ymnx: fltarr(2), $
            use_var: 'U', $
            nhirange: nhirange, $
            nhrange: nhrange, $
            fehrange: fehrange, $
            urange: urange, $
            psym: 10, $
            xtitle: 'U', $
            base_id: 0L, $      ; Widgets
            vars_id: 0L, $   ;  vars
            u_id: 0L, $ 
            nh_id: 0L, $
            nhi_id: 0L, $
            feh_id: 0L, $
            size: lonarr(2), $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            picvars_id: 0L, $
            usevar: 'U',$
            ionwin_id: 0L, $
            ionwin: ['HI','SiII','SiIV','CII','CIV','AlIII','OI','FeII'], $
            picvars: ['U','NHI','nH','FeH'] $
          }

  
; Set ymnx

  if not keyword_set(YMNX) then state.ymnx = [8.,18.]

  if keyword_set(NHI) then state.nhi = nhi else nhi = state.nhi
  if keyword_set(HDEN) then state.nh = hden else hden = state.nh
  if keyword_set(U) then state.u = u else u = state.u
  if keyword_set(FEH) then state.feh = feh else feh = state.feh

  state.onhi = state.nhi
  state.onh = state.nh
  state.ou = state.u
  state.ofeh = state.feh

;    WIDGET
  base = WIDGET_BASE( title = 'cldy_plot', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)


;;;;;;;;;
;  Toolbar

  state.picvars_id = WIDGET_LIST(toolbar, VALUE=state.picvars, $
                                     uvalue='PICVARS', ysize=4)
  state.ionwin_id = WIDGET_LIST(toolbar, VALUE=state.ionwin, $
                                     uvalue='IONWIN', ysize=4, /MULTIPLE)

  


; Values

  vars_id = widget_base(toolbar, /row, /frame,/base_align_center, /align_center)
  state.u_id = WIDGET_LIST(vars_id, VALUE=state.urange, uvalue='U',ysize=4)
  state.feh_id = WIDGET_LIST(vars_id, VALUE=state.fehrange, uvalue='FeH',ysize=4)
  state.nh_id = WIDGET_LIST(vars_id, VALUE=state.nhrange, uvalue='nH',ysize=4)
  state.nhi_id = WIDGET_LIST(vars_id, VALUE=state.nhirange, uvalue='NHI',ysize=4)

  
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
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Update
  cldyplot_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'cldy_plot', base, /no_block $
  else xmanager, 'cldy_plot', base

return
end
