;+ 
; NAME:
; x_identify   
;   Version 1.1
;
; PURPOSE:
;    Interactively identifies and fits reference lines
;
; CALLING SEQUENCE:
;   
;   x_identify, spec, calib, LINELIST=, XSIZE=, YSIZE=, REJ=, AUTO=,
;                 DISP=, LINEROOT=
;
; INPUTS:
;   spec       - 1D Spectrum
;
; RETURNS:
;
; OUTPUTS:
;   calib      - Structure describing the wavelength solution
;
; OPTIONAL KEYWORDS:
;   LINELIST   - Reference list
;   DISP       - Guess at dispersion (A or km/s per pix [pos/neg
;                value])
;   FLUX       - Assumes wide slit and therefore wide arc lines
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_identify, arc, calib
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Dec-2001 Written by JXP
;   19-Dec-2001 Added fitting
;   01-Feb-2003 Added input
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro x_identify_initcommon

;

common x_identify_fit, svffit, fin_fit, tcalib, flg_calib, $
  xid_pdmenu, spec_msk, gd_lines

  ; Wavelength solution
  tcalib = { fitstrct }
  tcalib.func = 'POLY'
  tcalib.nord = 2
  tcalib.niter = 3
  flg_calib=0
  xid_pdmenu = 0
  gd_lines = 0

end

;;;;
; Events
;;;;

pro x_identify_event, ev

common x_identify_fit

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
                      1 : x_identify_SetLine, state  ; Mark a line
                      2 : state.xymnx = state.svxymnx ; Reset the screen
                      4 : begin ; Zoom region
                          x_identify_SetZoom, state
                          if state.flg_reg EQ 1 then begin
                              WIDGET_CONTROL, state.base_id, set_uvalue = $
                                state, /no_copy
                              return
                          endif
                      end
                  endcase
              end
              1 : begin ; Button release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
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
          if (state.flg_reg EQ 1 AND eventch NE 's') then begin
              widget_control, state.error_msg_id, $
                set_value='Expecting another s !!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          case eventch of
              ;; Set edges
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'l': state.xymnx[0] = xgetx_plt(state, /STRCT) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              's': begin  ; Zoom Region
                  x_identify_SetZoom, state
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = $
                        state, /no_copy
                      return
                  endif
              end
              ;; Panning
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'W': x_identify_Scr, state ; Reset the screen
              'z': x_identify_ZoomLine, state, xgetx_plt(state, /strct) ; Zoom
              'f': x_identify_Flip, state ; Flip the x-axis
              ; LINES
              'm': x_identify_SetLine, state        ; Mark one line
              'd': x_identify_DelLine, state        ; Delete one line 
              'D': x_identify_DelLine, state, /all  ; Delete all lines
              'F': x_identify_FitLines, state       ; Fit current lines
              'L': x_identify_AutoID, state         ; Auto id more lines 
              'A': x_identify_Auto, state           ; Auto id the template lines
              ; Postscript
              'P': x_identify_psfile, state
              ; QUIT
              'q': begin
                  gdlin = where(state.lines.flg_plt MOD 2 EQ 1, nlin)
                  gd_lines = state.lines[gdlin]
                  widget_control, ev.top, /destroy
                  return
              end
              else: begin ; Nothing
                  widget_control, state.error_msg_id, $
                    set_value='Keystroke not defined: '+eventch
              end
          endcase
      end
      'FLIP': x_identify_Flip, state
      'PIXWAV': begin  ; Set Pixel/Wavelength
          WIDGET_CONTROL, ev.id, get_value = val
          case val of 
              0: begin
                  state.xymnx[0] = min(state.xdat, max=mx)
                  state.xymnx[2] = temporary(mx)
                  if state.flg_plot MOD 2 EQ 1 then $
                    state.flg_plot = state.flg_plot-1
              end
              1: begin ; Set to wavelength
                  if flg_calib NE 1 then $
                    WIDGET_CONTROL, ev.id, set_value = 0 $
                  else begin
                      if state.flg_plot MOD 2 NE 1 then $
                        state.flg_plot=state.flg_plot+1
                      state.xymnx[0] = min(state.wave, max=mx)
                      state.xymnx[2] = temporary(mx)
                  endelse
              end
              else:
          endcase
      end
      'RESDISPV': begin  ; Set dispersion/resolution
          widget_control, state.resdispv_id, get_value=tmp
          state.disp = temporary(tmp)
      end
      'LINELIST': begin
          state.listid = ev.index
          x_identify_PrsList, state
      end
      'DONE' : begin
          ;; Save good lines
          gdlin = where(state.lines.flg_plt MOD 2 EQ 1, nlin)
          gd_lines = state.lines[gdlin]
          widget_control, ev.top, /destroy
          return
      end
      else:
  endcase

; Update Plot
  x_identify_UpdatePlot, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_identify_UpdatePlot, state
  
common x_identify_fit

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  color = getcolor(/load)

  if state.flg_plot MOD 2 EQ 0 then begin  ; PIXELS

      ;  Arc spectrum
      plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
        /nodata, xstyle=1, ystyle=1, $
        xrange=[state.xymnx[0], state.xymnx[2]], $
        position=state.pos,  background=color.white, color=color.black 
      oplot, state.xdat, state.fx, psym=10, color=color.black 

;  LINES

      if state.nlin NE 0 then begin
          gdlin = where(state.lines.flg_plt MOD 2 EQ 1, nlin)
          if nlin NE 0 then begin
              peaks = state.fx[round(state.lines[gdlin].pix)]
              top = where(peaks LT 0.8*state.xymnx[3], ntop, $
                       COMPLEMENT=bot, NCOMPLEMENT=nbot)
              ; Most lines
              if ntop GT 0 then begin
                  oplot, [state.lines[gdlin[top]].pix], $
                    [peaks[top] + 0.1*(state.xymnx[3]-state.xymnx[1])], $
                    psym=7, color=color.red
                  xyouts, [state.lines[gdlin[top]].pix], $
                    [peaks[top] + 0.13*(state.xymnx[3]-state.xymnx[1])], $
                    [strtrim(state.lines[gdlin[top]].name,2)+' '+$
                     strmid(strtrim(state.lines[gdlin[top]].wave,2),0,7)], $
                    orientation=90., color=color.red, charsize=1.3
              endif
              ; High lines
              if nbot GT 0 then begin
                  oplot, [state.lines[gdlin[bot]].pix], $
                    [replicate(state.xymnx[1] + $
                               0.13*(state.xymnx[3]-state.xymnx[1]),nbot)], $
                    psym=7, color=color.red
                  xyouts, [state.lines[gdlin[bot]].pix-3], $
                    [replicate(state.xymnx[1]+$
                               0.10*(state.xymnx[3]-state.xymnx[1]),nbot)], $
                    [strtrim(state.lines[gdlin[bot]].name,2)+' '+$
                     strmid(strtrim(state.lines[gdlin[bot]].wave,2),0,7)], $
                    orientation=90., color=color.red, charsize=1.3, $
                    alignment=1.0
              endif
                
          endif
      endif

  endif else begin  ; WAVELENGTH

;  Arc spectrum
      plot, [state.xymnx[0],state.xymnx[2]], [state.xymnx[1],state.xymnx[3]], $
        /nodata, xstyle=1, ystyle=1, $
        xrange=[state.xymnx[0], state.xymnx[2]], $
        position=state.pos,  background=color.white, color=color.black 
      oplot, state.wave, state.fx, psym=10, color=color.black 

;  LINES

      if state.nlin NE 0 then begin
          gdlin = where(state.lines[0:state.nlin-1].flg_plt MOD 2 EQ 1, nlin)
          if nlin NE 0 then begin
              peaks = state.fx[round(state.lines[gdlin].pix)]
              oplot, [state.lines[gdlin].wave_fit], $
                [peaks + 0.1*(state.xymnx[3]-state.xymnx[1])], $
                psym=7, color=color.red
              xyouts, [state.lines[gdlin].wave_fit], $
                [peaks + 0.13*(state.xymnx[3]-state.xymnx[1])], $
                [strtrim(state.lines[gdlin].name,2)+' '+$
                 strmid(strtrim(state.lines[gdlin].wave,2),0,7)], $
                orientation=90., color=color.red, charsize=1.3
                
          endif
      endif
  endelse


end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_identify_Reset, state

; gdpix Flag
;   1 = Good
;   2 = In region
;   4 = Rejected by fit
  state.gdpix[0:n_elements(state.xdat)-1] = 1

; Plotting
  state.xymnx = state.svxymnx

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;
;  Set List
;;;;;;;;;
   
pro x_identify_PrsList, state

; Read the List 
  linfil = strjoin([state.lineroot, state.lists[state.listid]])
  readcol, linfil, FORMAT='D,I,A', wav, flg, nam
  state.nlin = n_elements(wav)
  state.lines[0:state.nlin-1].wave = temporary(wav)
  state.lines[0:state.nlin-1].flg_qual = temporary(flg)
  state.lines[0:state.nlin-1].name = temporary(nam)

; Reset to No Plotting
  state.lines.flg_plt = 0
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;
;  Mark Line
;;;;;;;;;
   
pro x_identify_SetLine, state

common x_identify_fit

  ; Make sure a list has been loaded in
  if state.nlin EQ 0 then begin
      print, 'x_identify_SetLine: Need to read a list in first'
      widget_control, state.error_msg_id, $
        set_value='x_identify_SetLine: Need to read a list in first'
      return
  endif
      
  ; Wave or Pixel?

  widget_control, state.text_id, sensitive = 0
  if state.flg_plot EQ 0 then begin ; PIXEL
      pix = xgetx_plt(state.xcurs,state.pos,state.xymnx, state.size)  
  endif else begin  ; WAVE
      wv = xgetx_plt(state.xcurs,state.pos,state.xymnx, state.size)  
      pix = x_fndfitval(wv, tcalib, [0., state.norg-1.])
  endelse

  ; CENTROID
  center = x_centspln(state.xdat[round(pix)-state.linpix:round(pix)+state.linpix], $
                      state.fx[round(pix)-state.linpix:round(pix)+state.linpix])
  ; CHOOSE ARC LINE
  line = x_slctline(state.lines, NLIN=state.nlin, ILIN=ilin, $
                    PDMENU=xid_pdmenu)
      
  ; Set the values
  state.lines[ilin].flg_plt = 1
  state.lines[ilin].pix = center

  ; Set wave as applicable
  if flg_calib EQ 1 then state.lines[ilin].wave_fit = $
    x_calcfit(center, FITSTR=tcalib)

  widget_control, state.text_id, /sensitive, /input_focus
  
  return
end

;;;;;;;;;;;;;;;;;;;;
;  Delete
;;;;;;;;;;;;;;;;;;;;

pro x_identify_DelLine, state, ALL=all

  ; All
  if keyword_set( ALL ) then begin
      state.lines.flg_plt = 0
      return
  end
  
  ; Find mouse
  
  x = xgetx_plt(state, /strct)

  ; Find existing lines
  gd = where(state.lines.flg_plt MOD 2 EQ 1, ngd)
  if ngd EQ 0 then return

  ; PIX OR WAVE?
  if state.flg_plot MOD 2 EQ 0 then begin ; PIXELS
      sep = abs(state.lines[gd].pix - x)
      minsep = min(sep, imin)
      if minsep LT 5. then state.lines[gd[imin]].flg_plt = 0 $
      else widget_control, state.error_msg_id, $
        set_value='DelLine: No nearby line'
  endif else begin ; WAVELENGTH
      sep = abs(state.lines[gd].wave - x)
      minsep = min(sep, imin)
      if minsep LT 3. then state.lines[gd[imin]].flg_plt = 0 $
      else widget_control, state.error_msg_id, $
        set_value='DelLine: No nearby line'
  endelse

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;
;  Fit Lines
;;;;;;;;;

pro x_identify_FitLines, state

common x_identify_fit

  ; ID lines
  gdlin = where(state.lines.flg_plt MOD 2 EQ 1, count)
  if count LT 3 then begin
      print, 'Need to ID at least 3 lines!'
      return
  endif

  ; FIT
  fit = x1dfit(state.lines[gdlin].pix, state.lines[gdlin].wave, /inter, $
              DELPTS=delpts, FITSTR=tcalib)

  ; Set wavelengths
  state.wave = x_calcfit(state.xdat, FITSTR=tcalib)
  state.svwvxymnx[0] = min(state.wave, MAX=mxwv)
  state.svwvxymnx[2] = mxwv

  flg_calib = 1

  WIDGET_CONTROL, state.pixwav_id, set_value = 1 ; Set switch to wavelength

  ; DELETE LINES
  if delpts[0] NE -1 then begin
      state.lines[gdlin[delpts]].flg_plt = 0
      for i=0L,n_elements(delpts)-1 do $
        print, 'Deleted ', state.lines[gdlin[delpts[i]]].name+' '+$
        strmid(strtrim(state.lines[gdlin[delpts[i]]].wave,2),0,7)
      ; Reset gdlin
      gdlin = where(state.lines.flg_plt MOD 2 EQ  1, count)
  endif

  ; PLOT AS WAVE
  if state.flg_plot MOD 2 EQ 0 then state.flg_plot = state.flg_plot + 1
  state.xymnx[0] = min(state.wave, max=mx)
  state.xymnx[2] = temporary(mx)

  ; RESET FLIP
  state.flg_flip = 0

  ; ID LINES
  state.lines[gdlin].wave_fit = x_calcfit(state.lines[gdlin].pix,$
                                          FITSTR=tcalib)

  ; DISP/RES
  x_identify_SetDispRes, state

  return
end

;;;;;;;;;
;  Auto ID :: Finds new lines given a first solution
;;;;;;;;;

pro x_identify_AutoID, state

common x_identify_fit

  widget_control, /hourglass

  ;; Find Peaks
  x_fndpeaks, state.fx, peak, NSIG=5., MSK=spec_msk
  npk = n_elements(peak)

  ; Add to lines
  for i=0L,npk-1 do begin
      ;; Require that a calibration line lies within +/- 5 pixels
      gdln = where( (state.lines[0:state.nlin-1].wave - $
                     state.wave[(peak[i]-5)>0])* $
                    (state.lines[0:state.nlin-1].wave - $
                     state.wave[(peak[i]+5)<(state.norg-1)]) LE 0 AND $
                    state.lines[0:state.nlin-1].flg_qual NE 0, nmatch)
      if nmatch NE 0 then begin
          ;; Find the closest calibration line
          sep = abs(state.lines[gdln].wave - state.wave[peak[i]])
          minsep = min(sep, imin)
          state.lines[gdln[imin]].flg_plt = 1 
          ;; Center up on the peak
          center = peak[i]
          state.lines[gdln[imin]].pix = center 
          state.lines[gdln[imin]].wave_fit = x_calcfit(center, FITSTR=tcalib)
      endif
  endfor

  return
end

;;;;;;;;;
;  Auto :: Finds a small set of lines given A/pix or km/s/pix
;     Keys on the lines in the line list with flag=5
;;;;;;;;;

pro x_identify_Auto, state

common x_identify_fit

  widget_control, /hourglass

; Error Checking
  ; Check dispersion
  if state.disp LE 0. then begin
      print, 'x_identify_Auto: Need to Set dispersion first!'
      widget_control, state.error_msg_id, $
        set_value='x_identify_Auto: Need to Set dispersion first!'
      return
  endif
  ; Not ready for resolution
  if state.flg_disp EQ 1 then begin
      print, 'x_identify_Auto: Not ready for resolution!'
      return
  endif

  ; Find the template lines
  gdlin = where(state.lines.flg_qual EQ 5, ngd)
  if ngd EQ 0 then begin
      print, 'x_identify_Auto: No lines for template!'
      widget_control, state.error_msg_id, $
        set_value='x_identify_Auto: No lines for template!'
      return
  endif
  
  ; Fit the spectrum with a crude BSPLIN
  autofit = x_fitrej(state.xdat, state.fx, 'BSPLIN', 31, $
                     hsigma=2., lsigma=5., rms=rms)

  ; Find all 100 sigma peaks and avoid edges
  gdpix = where( state.fx GT autofit+100.*rms AND $
                 state.xdat GT 5 AND $
                 state.xdat LT state.norg-6, ngpix)

  ; Include only the inflection points (max of 5 pts)
  npk = 0L
  peak = lonarr(state.norg)
  for i=0L,ngpix-1 do begin
      if ( state.fx[gdpix[i]-1] GT state.fx[gdpix[i]-2] AND $
           state.fx[gdpix[i]] GT state.fx[gdpix[i]-1] AND $
           state.fx[gdpix[i]] GT state.fx[gdpix[i]+1] AND $
           state.fx[gdpix[i]+1] GT state.fx[gdpix[i]+2] ) then begin
          peak[npk] = gdpix[i]
          npk = npk + 1
      endif
  endfor
  peak = temporary(peak[0:npk-1])

  ; Center up
  center = fltarr(npk)
  for i=0L,npk-1 do begin
      ; Center up on the peak
      center[i] = x_centspln(state.xdat[peak[i]-3:peak[i]+3], $
                          state.fx[peak[i]-3:peak[i]+3])
      ; Catch bad centroids (Not applicable anymore)
      if center[i] EQ -1 OR abs(center[i]-state.xdat[peak[i]]) GT 1. then begin
          center[i] = state.xdat[peak[i]]
      endif 
  endfor

  ; Flip now!
  if state.flg_bluered EQ 1 then center = (state.norg-1.d)-center

  ; Loop parameters
  nzro = (state.norg-1)/2L
  ndlmb = 20
  nnonl = 50  ; Strect (non-linear term)
  all_chi = fltarr(nzro,ndlmb,nnonl)


  ; Search limits
  min_zro = 3000.d  ; Atmosphere
  max_zro = min_zro + nzro*state.disp/2.
  min_dlmb = state.disp*0.90
  max_dlmb = state.disp*1.10
  min_nonl = -(300.d*state.disp)/(state.norg-1)^2  ; Less than 300 pix stretch
  max_nonl = (300.d*state.disp)/(state.norg-1)^2  ; Less than 300 pix stretch

  ; Begin the loop
  for q1=0L,nzro-1 do begin
      Zro = min_zro + q1*(max_zro-min_zro)/float(nzro)
      for q2=0L,ndlmb-1 do begin
          dlmb = min_dlmb + q2*(max_dlmb-min_dlmb)/float(ndlmb)
          for q3=0L,nnonl-1 do begin
              nonl = min_nonl + q3*(max_nonl-min_nonl)/float(nnonl)

              ; Calculate lambda for the peaks 
              lambda = Zro + center*dlmb + (center^2)*nonl 

              ; Find closest for each template line and calulate chisq
              chisq = 0.d
              for i=0L,ngd-1 do begin
                  min_sep = min(abs(lambda-state.lines[gdlin[i]].wave))
                  chisq = chisq + min_sep^2/(state.disp/2.)^2
              endfor
              ; Save chisq
              all_chi[q1,q2,q3] = chisq
          endfor
      endfor
  endfor

  ; Take min
  min_chi = min(all_chi, imin)
  
  ; Find the parameters
  q1 = imin MOD nzro
  q2 = (imin MOD (nzro*ndlmb))/nzro
  q3 = imin/(ndlmb*nzro)
  Zro = min_zro + q1*(max_zro-min_zro)/float(nzro)
  dlmb = min_dlmb + q2*(max_dlmb-min_dlmb)/float(ndlmb)
  nonl = min_nonl + q3*(max_nonl-min_nonl)/float(nnonl)

  ; Find the best lines
  lambda = Zro + center*dlmb + center^2*nonl
  for i=0L,ngd-1 do begin
      min_sep = min(abs(lambda-state.lines[gdlin[i]].wave), imin)
      ; Keep it? (5 pix separation)
      if min_sep LT 5*state.disp then begin
          ; Deal with the dreaded flip
          if state.flg_bluered EQ 0 then state.lines[gdlin[i]].pix = center[imin] $
          else state.lines[gdlin[i]].pix = (state.norg-1.d)-center[imin]
          ; Save the line!
          state.lines[gdlin[i]].flg_plt = 1
      endif
  endfor

  return
end

;;;;;;;;;;;;;;;;;
; Set Dispersion/Resolution automatically
pro x_identify_SetDispRes, state


  ;  Avoids the ends of the spectrum where 
     ;     wavelength error may be high
  nstrt = long(state.norg/10)
  nend = long(9*state.norg/10)

  ; Dispersion
  d1 = state.wave[nstrt+1]-state.wave[nstrt]
  d2 = state.wave[nend] - state.wave[nend-1]
  diff_disp = abs(d1-d2)/abs(d1+d2)

  ; Resolution (dlambda/lambda)
  r1 = (state.wave[nstrt+1]-state.wave[nstrt])/state.wave[nstrt]
  r2 = (state.wave[nend] - state.wave[nend-1])/state.wave[nend]
  diff_res = abs(r1-r2)/abs(r1+r2)

  if diff_disp LT diff_res then begin
      state.flg_disp = 0
      state.disp = abs(d1)  ; A/pix [zero point term]
      widget_control, state.resdispb_id, set_value=0
      widget_control, state.resdispv_id, set_value=state.disp
  endif else begin
      state.flg_disp = 1
      state.disp = abs(r1)*2.99E5  ; km/s/pix [zero point term]
      widget_control, state.resdispb_id, set_value=1
  endelse

end

;;;;;;;;;;;;;;;;;;;;;
;WINDOW STUFF
;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;
; Set Zoom region
pro x_identify_SetZoom, state

  if(state.flg_reg EQ 0) then begin
      state.tmpxy[0] = xgetx_plt(state, /strct) ; left
      state.tmpxy[1] = xgety_plt(state, /strct) ; bottom
      state.flg_reg = 1
  end else begin
      state.tmpxy[2] = xgetx_plt(state, /strct) ; right
      state.tmpxy[3] = xgety_plt(state, /strct) ; top
      state.flg_reg = 0
      if state.flg_flip MOD 2 EQ 0 then begin
          state.xymnx[0] = state.tmpxy[0] < state.tmpxy[2]
          state.xymnx[2] = state.tmpxy[0] > state.tmpxy[2]
      endif else begin
          state.xymnx[2] = state.tmpxy[0] < state.tmpxy[2]
          state.xymnx[0] = state.tmpxy[0] > state.tmpxy[2]
      endelse
      state.xymnx[1] = state.tmpxy[1] < state.tmpxy[3]
      state.xymnx[3] = state.tmpxy[1] > state.tmpxy[3]
  endelse
end

;;;;;;;;;;;;;;;;;
; Zoom on a line
pro x_identify_ZoomLine, state, x 

  if state.flg_plot MOD 2 EQ 0 then begin ; PIXELS
      xmin = (round(x) - 7) > 0
      xmax = (round(x) + 7) < (state.norg-1)
      ymin = min(state.fx[xmin:xmax], MAX=ymax)

      state.xymnx[0] = xmin-0.5
      state.xymnx[2] = xmax+0.5
      state.xymnx[1] = ymin-(ymax-ymin)*0.05
      state.xymnx[3] = ymax+(ymax-ymin)*0.35
  endif else begin
      wavsep = abs(state.wave - x)
      minsep = min(wavsep, iwav)
      xmin = (iwav - 7) > 0
      xmax = (iwav + 7) < (state.norg-1)
      ymin = min(state.fx[xmin:xmax], MAX=ymax)

      state.xymnx[0] = state.wave[xmin] < state.wave[xmax]
      state.xymnx[2] = state.wave[xmin] > state.wave[xmax]
      state.xymnx[1] = ymin-(ymax-ymin)*0.05
      state.xymnx[3] = ymax+(ymax-ymin)*0.35
      delvarx, wavsep
  endelse

  if state.flg_flip MOD 2 EQ 1 then x_identify_Flip, state, /NOFLG
      
end

;;;;;;;;;;;;;;;;;
; Flip
pro x_identify_Flip, state, NOFLG=noflg

  ; Flag
  if not keyword_set(NOFLG) then state.flg_flip = state.flg_flip+1

  ; Flip
  tmp = state.xymnx[0]
  state.xymnx[0] = state.xymnx[2]
  state.xymnx[2] = tmp

end

;;;;;;;;;;;;;;;;;
; Reset Window
pro x_identify_Scr, state 

  ; PIX or WAVE?
  if state.flg_plot MOD 2 EQ 0 then state.xymnx = state.svxymnx $
  else state.xymnx=state.svwvxymnx 

  ; Flip
  if state.flg_flip MOD 2 EQ 1 then x_identify_Flip, state, /NOFLG

end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro x_identify_psfile, state

; Device
  device, get_decomposed=svdecomp

  !p.thick = 3
  !p.charthick = 3

  device, decompose=0
  ps_open, file='idl.ps', font=1, /color
  state.psfile = 1
  x_identify_UpdatePlot, state
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
  state.psfile = 0
  !p.thick = 1
  !p.charthick = 1

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

pro x_identify, spec, calib, XSIZE=xsize, YSIZE=ysize, LINELIST=linelist, $
                REJ=rej, AUTO=auto, DISP=disp, LINEROOT=lineroot, $
                DEBUG=debug, MSK=msk, OUTLIN=outlin, $
                REDBLUE=redblue, WAVE=wave, FLUX=flux, INLIN=inlin

common x_identify_fit


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_identify, spec, calib, LINELIST=, XSIZE=, YSIZE= '
    print, '     /REJ, /AUTO, DISP=, LINEROOT=, /FLUX, INLIN=, WAVE=, /REDBLUE [v1.1]'
    return
  endif 

; Initialize the common block
  x_identify_initcommon
  if keyword_set( MSK ) then spec_msk = msk

;  Optional Keywords

  if not keyword_set( LINEROOT ) then $
    lineroot = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/'
  if keyword_set( FLUX ) then linpix = 10L else linpix=5L


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Non-interactive

  if keyword_set (AUTO) then begin
      print, 'Not yet!'
      return
  endif else begin
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; INTERACTIVE

  if not keyword_set( XSIZE ) then    xsize = 1200
  if not keyword_set( YSIZE ) then    ysize = 800

  tmp = { arclinstrct }

;    STATE
  state = { $
            fx: spec, $
            xdat: findgen(n_elements(spec)), $
            wave: dblarr(n_elements(spec)), $
            norg: n_elements(spec), $
            lines: replicate(tmp, 20000), $  ; LINES
            nlin: 0L, $
            linpix: linpix, $
            nlist: 0, $
            listid: 0, $
            lists: strarr(100), $
            lineroot: lineroot, $  ; Path to lists
            cwlabel: strarr(10000), $
            flg_plot: 0, $  ; 0 = Pixel, 1 = Wavelength
            flg_disp: 0, $  ; 0= Dispersion, 1=Resolution
            disp: 0.d, $
            flg_bluered: 0, $  ; 0= Blue-Red; 1=Red-Blue
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            svxymnx: fltarr(4), $
            svwvxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            psfile: 0, $
            flg_flip: 0, $
            flg_reg: 0, $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            linelist_id: 0L, $
            pixwav_id: 0L, $
            redblue_id: 0L, $
            guess_id: 0L, $
            guessbut_id: 0L, $
            resdispb_id: 0L, $
            resdispv_id: 0L, $
            error_msg_id: 0L, $
            gdpix: intarr(200000) $
          }


  ; svxymnx
  state.svxymnx = [min(state.xdat)-0.01*abs(max(state.xdat)-min(state.xdat)), $
                  min(state.fx)-0.01*abs(max(state.fx)-min(state.fx)), $
                  max(state.xdat)+0.01*abs(max(state.xdat)-min(state.xdat)), $
                  max(state.fx)+0.01*abs(max(state.fx)-min(state.fx))]
  state.svwvxymnx = state.svxymnx
;    WIDGET
  base = WIDGET_BASE( title = 'x_identify: Interactive Mode', /column, $
                    xoffset=300, yoffset=200)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
  strlbl = strarr(10)
  strlbl = ['x_identify', ' ', 'Ver 1.0']
  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;;; LINELISTS

  a = findfile(strjoin([state.lineroot,'*.lst']), count=nlist)
  state.nlist = nlist

  if state.nlist NE 0 then $
    for i = 0,state.nlist-1 do state.lists[i] = fileandpath(a[i]) $
  else begin
      print, 'x_identify: No line lists!'
      return
  endelse

  state.linelist_id = WIDGET_LIST(toolbar, VALUE=state.lists[0:state.nlist-1],$
                                  uvalue='LINELIST', ysize=6)

  ; Set initial linelist
  if keyword_set( LINELIST ) then begin
      init = where(strtrim(linelist,2) EQ a, count)
      if count NE 0 then begin
          widget_control, state.linelist_id, set_list_select=init[0] 
          state.listid = init[0]
          x_identify_PrsList, state
      endif else print, linelist+' not found!'
  endif
      
;;; Pixel/Wavelength Toggle

  xaxis_id = WIDGET_BASE(toolbar, /column, /align_center, /frame)
  state.pixwav_id = CW_BGROUP(xaxis_id, ['Pixel', 'Wave'], $
                              row=2, /exclusive, $
                              set_value=0,  uvalue='PIXWAV', $
                              LABEL_TOP="x-axis")
  state.redblue_id = CW_BGROUP(xaxis_id, ['Blue-Red', 'Red-Blue'], $
                              row=2, /exclusive, $
                              set_value=0,  uvalue='BLUERED')

;;;; GUESS
  state.guess_id = WIDGET_BASE(toolbar, /column, /base_align_center, /frame, $
                               /align_center)
  state.resdispb_id = CW_BGROUP(state.guess_id, ['Dispersion', 'Resolution'], $
                               /exclusive, $
                               row=2, /frame, set_value=0, uvalue='RESDISPB')
  state.resdispv_id = CW_FIELD(state.guess_id, value=0., xsize=10, ysize=1, $
                               title='RES/DISP', /floating, $
                               /return_events, uvalue='RESDISPV')
  state.guessbut_id = WIDGET_BUTTON(state.guess_id, value='Guess',$
                                    uvalue='GUESS', sensitive=0)
  if state.disp GT 0 then widget_control, state.guessbut_id, /sensitive
  
  
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'Mouse:::::::::::::::::',$
             'LMB -- Mark a line',$
             'RMB+RMB -- Zoom in region',$
             'CMB -- Reset screen',$
             'Window Stuff:::::::::: ',$
             'l -- Set Left Window',$
             'r -- Set Right Window',$
             'b -- Set Bottom Window',$
             't -- Set Top Window',$
             'W -- Reset screen',$
             'f -- Flip x-axis',$
             's,s -- Zoom',$
             'Lines::::::::::::::::: ',$
             'm -- mark a line (LMB)',$
             'F -- Fit current lines',$
             'L -- Auto ID more lines (given a fit)',$
             'A -- Auto ID all lines',$
             'd -- Delete a line',$
             'D -- Delete all lines',$
             'Saving:::::::::::::::: ',$
             'q -- Quit and save']
  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
                             /scroll)
;      Error base
  error_base_id  = widget_base(toolbar, /column)
  error_btn_id = widget_button(error_base_id, value='Erase Error Msg', $
                               uvalue='ERRORB')
;          Error Messages
  state.error_msg_id = widget_text(error_base_id, value='', xsize=30, $
                                   ysize=2, /scroll)
  
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)
  
  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /button_events, /motion_events, $
                              uvalue='DRAW') ;, retain=2)
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
  
; Include all elements at first
  state.gdpix[0:n_elements(xdat)] = 1
  
; Color Table 

; Update
  x_identify_Reset, state

  ; Blue/Red
  if keyword_set(REDBLUE) then begin
      x_identify_Flip, state
      state.flg_bluered = 1
      widget_control, state.redblue_id, set_value=1
      widget_control, state.error_msg_id, $
        set_value='We are flipped! '
  endif

  ; Dispersion
  if keyword_set(DISP) then begin
      if disp LT 0 then state.flg_disp = 1 ; Resolution (km/s/pix)
      state.disp = abs(disp)
      widget_control, state.resdispv_id, set_value=state.disp 
  endif

; INLIN
  if keyword_set(INLIN) then begin
      state.nlin = n_elements(inlin)
      state.lines[0:state.nlin-1] = inlin
  endif

  x_identify_UpdatePlot, state

; Debug
  if keyword_set( DEBUG ) then begin
      state.listid = 0
      x_identify_PrsList, state
      state.lines[6].flg_plt = 1
      state.lines[6].pix = 868.456
      state.lines[26].flg_plt = 1
      state.lines[26].pix = 458.463
      state.lines[66].flg_plt = 1
      state.lines[66].pix = 92.0117
  endif
      
  ; 
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
; Send to the xmanager
  xmanager, 'x_identify', base
  
; Passing back
  
  endelse

  calib = tcalib
;  Wavelength array
  if arg_present(WAVE) AND flg_calib NE 0 then $
      wave = x_calcfit(findgen(n_elements(spec)), FITSTR=tcalib)
  delvarx, tcalib

  ;; OUTLIN (good lines)
  outlin = temporary(gd_lines)

  return
end
