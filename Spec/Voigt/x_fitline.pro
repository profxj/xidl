;+ 
; NAME:
; x_fitline   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   x_fitline, ydat, [head], XSIZE=, YSIZE=, TITLE=, WAVE=
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array
;   ERR        - Error array (fits or image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fitline, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xfitline_UpdatePlot, state
  
; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

  plot, state.wave, state.fx, psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.7, $
    ycharsize=1.7

  ; Plot Error array
  if state.flg_sig EQ 1 then $
    oplot, state.wave, state.sig, psym=state.psym, color=clr.red

  ;; FIT
  if state.nlin NE 0 then begin
      oplot, state.wave, state.fit, color=clr.green
      ;; Mark current line
      oplot, replicate( (state.lines[state.curlin].zabs+1.)*$
                        state.lines[state.curlin].wrest, 2), $
        [state.xymnx[1], state.xymnx[3]], color=clr.blue, linestyle=1
  endif
  ;; Continuum
  oplot, state.wave, state.conti, color=clr.purple, linestyle=1


end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro xfitline_Reset, state


; Plotting
  state.xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_psfile, state

; Device
  device, get_decomposed=svdecomp

  !p.thick = 3
  !p.charthick = 3

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  xfitline_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0
  !p.thick = 1
  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;
;  New Lya line
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_newlin, state

  ;; Grab x,y pos
  xpt = state.xpos

  ;; Setup HI
  tmp = x_setline(1215.670d)

  ;; SET
  state.nset = state.nset + 1

  ;; Get z
  state.lines[state.nlin] = tmp
  state.lines[state.nlin].zabs = (xpt / 1215.6701) - 1.
  state.lines[state.nlin].N = 20.3
  state.lines[state.nlin].b = 30.0
  state.lines[state.nlin].set = state.nset

  ;; Update text windows
  widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
  widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
  widget_control, state.bval_id, set_value=state.lines[state.nlin].b

  ;; Current line and increment
  state.curlin = state.nset
  state.nlin = state.nlin + 1

  ;; Update fit
  x_fitline_updfit, state
  return
end

;;;;;;;;;;;;;;;;;;;;
;  New LLS
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_newlls, state

  widget_control, /hourglass   
  ;; Get z
  xpt = state.xpos
  ztmp = (xpt / 1215.6701) - 1.
  ;; NSET
  state.nset = state.nset + 1

  ;; Get LLS lines
  tmplin = x_mknewlls(ztmp, set=state.nset)
  nlin = n_elements(tmplin)

  ;; Add to state
  state.lines[state.nlin:state.nlin+nlin-1] = tmplin

  ;; Update window
  widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
  widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
  widget_control, state.bval_id, set_value=state.lines[state.nlin].b

  ;; Current line and increment
  state.curlin = state.nset
  state.nlin = state.nlin + nlin

  ;; Update fit
  x_fitline_updfit, state
  return
end

;;;;;;;;;;;;;;;;;;;;
;  Delete line
;     Deletes a set of lines not just one
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_dellin, state

  ;; No lines
  if state.nset LT 0 then return

  ;; 1 line
  if state.nset EQ 0 then begin
      state.nset = -1
      state.nlin = 0
      state.curlin = state.nset
      return
  endif

  ;; Multiple lines
  gdlin = where(state.lines.set NE state.curlin, ngd)
  state.lines[0:ngd-1] = state.lines[gdlin]

  ;; Reset set
  if state.curlin NE state.nset then begin
      big = where(state.lines.set GT state.curlin)
      state.lines[big].set = state.lines[big].set - 1
  endif
  state.nset = state.nset - 1
  state.curlin = state.curlin < state.nset

  ;; Update fit
  x_fitline_updfit, state
  return
end

;;;;;;;;;;;;;;;;;;;;
;  Update fit
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_updfit, state, FLG_PLT=flg_plt

  ;; flg_plt
  flg_plt = 1

;  xpt = state.xpos
  ;; Calculate
;  xmin = state.wave[0] > (xpt - 40.*(1+min(state.lines[0:state.nlin-1].zabs)))
;  xmax = state.wave[state.npix-1] < (xpt + 40.* $
;                                     (1+max(state.lines[0:state.nlin-1].zabs)))

;  mn = min(abs(state.wave-xmin), imn)
;  mx = min(abs(state.wave-xmax), imx)

  ;; Voigt
  state.fit = 1.
  state.fit = x_allvoigt(state.wave, state.lines[0:state.nlin-1], $
                                  SIGMA=state.FWHM) 

  ;; Continuum
  if state.conti[0] NE 1. then state.fit = state.fit * state.conti

  return
end

;;;;;;;;;;;;;;;;;;;;
;  ASCII Output
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_output, state

  close, 11
  openw, 11, 'fort.13'
  for q=0L, state.nlin-1 do begin
      printf, 11, q, state.lines[q].ion, state.lines[q].wrest, $
        state.lines[q].zabs, state.lines[q].N, state.lines[q].b, $
        FORMAT='(i2,1x,a11,1x,f9.3,1x,f9.7,1x,f5.2,1x,f5.2)'
  endfor
  close, 11

  return
end

;;;;;;;;;;;;;;;;;;;;
;  IDL Output
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_idlout, state

  lines = state.lines[0:state.nlin-1]
  conti = state.conti
  if tag_exist(state, 'cstr') EQ 1 then cstr = state.cstr else cstr = 0.
  save, cstr, conti, lines, filename='fort.idl'
  delvarx, lines

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Update window
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_updwin, state

  ;; Update text windows
  widget_control, state.zabs_id, set_value=state.lines[state.curlin].zabs
  widget_control, state.Ncolm_id, set_value=state.lines[state.curlin].N
  widget_control, state.bval_id, set_value=state.lines[state.curlin].b

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Continuum
;;;;;;;;;;;;;;;;;;;;

pro x_fitline_continuum, state, flg


; CASE

  case flg of 
      0: begin ; Set point
          i = state.cstr.npts 
          state.cstr.xval[i] = state.xpos
          state.cstr.yval[i] = state.ypos
          state.cstr.msk[i] = 1L
          state.cstr.npts =  state.cstr.npts  + 1
      end
      1: begin ; Move point
          if state.cstr.npts EQ 0 then return
          gd = where(state.cstr.msk EQ 1)
          mn = min(abs(state.cstr.xval[gd]-state.xpos), imn)
          state.cstr.yval[gd[imn]] = state.ypos
      end
      else: stop
  endcase

  gdmsk = where(state.cstr.msk EQ 1, nmsk)
  if nmsk LT 3 then state.conti = mean(state.cstr.yval[gdmsk]) $
  else begin ;; SPLINE
      mn = min(state.cstr.xval[gdmsk], imn)
      ;; Low
      low = where(state.wave LE mn, nlow)
      if nlow NE 0 then state.conti[low] = state.cstr.yval[gdmsk[imn]]
      ;; High
      mx = max(state.cstr.xval[gdmsk], imx)
      high = where(state.wave GE mx, nhigh)
      if nhigh NE 0 then state.conti[high] = state.cstr.yval[gdmsk[imx]]
      ;; Interp
      gd = where(state.wave GT mn AND state.wave LT mx, ngd)
      if ngd NE 0 then begin
          ;; SORT
          srt = sort(state.cstr.xval[gdmsk])
          state.conti[gd] = spline(state.cstr.xval[gdmsk[srt]], $
                                   state.cstr.yval[gdmsk[srt]], $
                                   state.wave[gd])
          ;; Update Fit
          if state.nlin NE 0 then x_fitline_updfit, state
      endif
  endelse   
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fitline, state, eventch, FLG_PLT=flg_plt

  ;; FLG_PLT
  flg_plt = 1

  if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
      print, 'Set the other region!'
      flg_plt = 0
      return
  endif

  case eventch of
      'b': state.xymnx[1] = state.ypos
      'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
      'l': state.xymnx[0] = state.xpos
      'r': state.xymnx[2] = state.xpos
      't': state.xymnx[3] = state.ypos
      'T': state.xymnx[3] = 1.1 ; Set ymax to 1.1
      ;; ZOOMING
      'i': x_speczoom, state, 0 ; Zoom in
      'o': x_speczoom, state, 1 ; Zoom out
      'z': begin
          ximgd_setzoom, state, /plot
          if state.flg_zoom EQ 1 then begin
              flg_plt = 0
              return
          endif
      end
      ;; Reset
      'w': state.xymnx = state.svxymnx ; Reset the screen
      ;; PANNING
      '}': x_specpan, state, /noy
      '{': x_specpan, state, /left, /noy
      ']': x_specpan, state
      '[': x_specpan, state, /left
      ;; HEL	P
      'H': begin 
          x_helpwidg, state.help
          flg_plt = 0
      end
      ;; Add new line
      'L': x_fitline_newLLS, state
      'c': x_fitline_newlin, state
      'd': x_fitline_dellin, state
      ;; Continuum
      'C': begin ;; Set all to a constant
          state.fit = state.fit * state.ypos / state.conti
          state.conti = state.ypos
          if tag_exist(state, 'CONTI_ID') then $
            widget_control, state.conti_id, set_value=state.ypos
          ;; Zero out continuum old points
          if tag_exist(state, 'cstr') then begin
              state.cstr.msk[*] = 0L
              state.cstr.npts = 0L
          endif
      end
      '3': begin ;; Add point
          if n_elements(state.conti) NE 1 then $
            x_fitline_continuum, state, 0L
      end
      '4': begin ;; Move a point
          if n_elements(state.conti) NE 1 then $
            x_fitline_continuum, state, 1L
      end
      ;; Colm
      'n': begin
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = state.lines[gdset[0]].N - 0.05
          x_fitline_updfit, state
          widget_control, state.Ncolm_id, $
            set_value=state.lines[gdset[0]].N
      end
      'N': begin
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = state.lines[gdset[0]].N + 0.05
          widget_control, state.Ncolm_id, $
            set_value=state.lines[gdset[0]].N
          x_fitline_updfit, state
      end
      ;; b-value
      'v': begin
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = state.lines[gdset[0]].b - 1.0
          x_fitline_updfit, state
          widget_control, state.bval_id, $
            set_value=state.lines[state.curlin].b
      end
      'V': begin
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = state.lines[gdset[0]].b + 1.0
          widget_control, state.bval_id, $
            set_value=state.lines[state.curlin].b
          x_fitline_updfit, state
      end
      ;; Switch lines
      '=': begin
          state.curlin = (state.curlin+1) < state.nset 
          x_fitline_updwin, state
      end
      '-': begin
          state.curlin = (state.curlin-1) > 0
          x_fitline_updwin, state
      end
      ;; Output
      'O': x_fitline_output, state
      'I': x_fitline_idlout, state
                                ; Postscript
      'P': x_fitline_psfile, state  
                                ; QUIT
      else:  print, 'x_fitline: Not a valid key!' ; Nothing
  endcase

  return
end
  
