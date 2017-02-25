;+ 
; NAME:
; ximgd_setzoom
;    Version 1.1
;
; PURPOSE:
;  Sets up zooming for either Image viewing or a standard Plot
;  window
;
; CALLING SEQUENCE:
;   ximgd_setzoom, state, [flg], /PLOT
;
; INPUTS:
;   state   -- Structure with tv dependent info or Plot info
;   [FLG]    --  If 1, set the intital values, otherwise do the zoom
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /PLOT -- Zooming in plot window not image window
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ximgd_setzoom, state, flg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_setzoom, state, flg, PLOT=plot

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ximgd_setzoom, state, [flg], /plot  [v1.1]'
    return
  endif 

  if keyword_set( PLOT ) then begin
      if(state.flg_zoom EQ 0) then begin
          state.tmpxy[0] = xgetx_plt(state, /strct) ; left
          state.tmpxy[1] = xgety_plt(state, /strct) ; bottom
          state.flg_zoom = 1
      endif else begin
          state.tmpxy[2] = xgetx_plt(state, /strct) ; right
          state.tmpxy[3] = xgety_plt(state, /strct) ; top
          state.flg_zoom = 0
          ; Allow for flipping
          idx = tag_exist(state, 'flg_flip')
          if idx NE 0 then begin
              if state.flg_flip MOD 2 EQ 0 then begin
                  state.xymnx[0] = state.tmpxy[0] < state.tmpxy[2]
                  state.xymnx[2] = state.tmpxy[0] > state.tmpxy[2]
              endif else begin
                  state.xymnx[2] = state.tmpxy[0] < state.tmpxy[2]
                  state.xymnx[0] = state.tmpxy[0] > state.tmpxy[2]
              endelse
          endif else begin ; NO FLIP
              state.xymnx[0] = state.tmpxy[0] < state.tmpxy[2]
              state.xymnx[2] = state.tmpxy[0] > state.tmpxy[2]
          endelse
          ; Standard y values
          state.xymnx[1] = state.tmpxy[1] < state.tmpxy[3]
          state.xymnx[3] = state.tmpxy[1] > state.tmpxy[3]
      endelse
  endif else begin

      if flg EQ 1 then begin
          state.tmpreg[0] = x_tvx(state.tv, /intg)
          state.tmpreg[1] = x_tvy(state.tv, /intg)
          state.zoom.flg = 1
      endif else begin
          state.tmpreg[2] = x_tvx(state.tv, /intg)
          state.tmpreg[3] = x_tvy(state.tv, /intg)
          
          if (state.tmpreg[0] NE state.tmpreg[2] AND $
              state.tmpreg[1] NE state.tmpreg[3]) then begin
              state.zoomreg[0] = state.tmpreg[0] < state.tmpreg[2]
              state.zoomreg[2] = state.tmpreg[0] > state.tmpreg[2]
              state.zoomreg[1] = state.tmpreg[1] < state.tmpreg[3]
              state.zoomreg[3] = state.tmpreg[1] > state.tmpreg[3]
              
              state.zoomreg[0] = state.zoomreg[0] > 0
              state.zoomreg[2] = state.zoomreg[2] > 0
              state.zoomreg[1] = state.zoomreg[1] < state.sz_img[0]-1
              state.zoomreg[3] = state.zoomreg[3] < state.sz_img[1]-1
          endif
          x_setgridxy, state, state.zoomreg, /FILL
      endelse
      state.tv.xymnx = float(state.zoomreg)
  endelse
  return

end
