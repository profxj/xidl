;+ 
; NAME:
; x_fitlls   
;   Version 1.11
;
; PURPOSE:
;    GUI used to fit DLA profiles interactively.  The user determines
;    the continuum at the same time.
;
; CALLING SEQUENCE:
;  x_fitlls, flux_fil, [err_fil], XSIZE=, YSIZE=, TITLE=, 
;               WAVE=, /BLOCK, FWHM=, INIFIT=, INFLG=
;
; INPUTS:
;   flux_fil   - FITS file (or array) containing flux
;   [err_fil]  - FITS error array
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array 
;   FWHM=      - Resoltuion (FWHM in pixels)  [default: 3]
;   INFLG=     - Usual flag for reading the FITS file (see x_readspec)
;   INIFIT=    - IDL file containing a saved version of the fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fitlls, 'spec.fits'
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

;;;;
; Common
;;;;

pro x_specplot_initcommon
;
common x_specplot_lines, $
   flg_lines, $
   lines, $
   linid, $
   zabs

   flg_lines = 0
   zabs = 0.

end


;;;;
; Events
;;;;

pro x_fitlls_event, ev

common x_specplot_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval
  flg_plt = 1

  case uval of
      'CRUDEHI': begin
          state.crude_val[1] = ev.value
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDELO': begin
          state.crude_val[0] = ev.value
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDE': begin
          state.crude = ev.value
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      ;; 
      'MAINWIN': begin
          state.all_xymnx[*,state.curwin] = state.xymnx
          state.all_xpmnx[*,state.curwin] = state.xpmnx
          x_fitlls_idx, state, ev.index, 0, id
          state.curwin = id
          state.xymnx[*] = state.all_xymnx[*,state.curwin]
          state.xpmnx[*] = state.all_xpmnx[*,state.curwin]
          flg_plt = 3
      end
      'OTHWIN': begin
          ;; Set flag fiwn
          idx = widget_info(state.othwin_id,/list_select)
          state.flgwin[*] = 0B
          if idx[0] NE -1 then begin
              for ii=0L,n_elements(idx)-1 do begin
                  x_fitlls_idx, state, idx[ii], 1, id
                  state.flgwin[id] = 1B 
              endfor
          endif 
          ;; Open (as necessary)
          x_fitlls_openwin, state
          if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
      end
      'EXACT': begin
          state.exact = ev.value
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'GCONT': x_fitlls_autocont, state
      'UREGIONS': begin
          state.uregions = ev.value
          x_fitlls_updfit, state
          x_fitlls_UpdatePlot, state
          x_fitlls_PlotOth, state
      end
      'INDIV': begin
          state.indiv = ev.value
          x_fitlls_UpdatePlot, state
          x_fitlls_PlotOth, state
      end
      'ZABS' : begin
         widget_control, /hourglass
         widget_control, state.zabs_id, get_value=tmp
         gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
         state.lines[gdset].zabs = tmp
         state.vpstr.z[state.curlin] = state.lines[gdset[0]].zabs
         x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
         if state.curlin EQ 0L then x_fitlls_setxy, state
         if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
      end
      'CONTI' : begin
          widget_control, state.conti_id, get_value=tmp
          state.conti = tmp
      end
      'BVAL' : begin
          widget_control, state.bval_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = tmp
          state.vpstr.b[state.curlin] = state.lines[gdset[0]].b
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
          if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
      end
      'NCOLM' : begin
          widget_control, state.Ncolm_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = tmp
          state.vpstr.n[state.curlin] = state.lines[gdset[0]].N
          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
          if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
      end
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
                      4 : begin
                          mn = min(abs(state.xpos-(state.lines.zabs+1)* $
                                       state.lines.wrest),imn)
                          idx = where(state.lines.set EQ $
                                      state.lines[imn].set)
                          state.lines[idx].zabs = $
                            (state.xpos / state.lines[imn].wrest) - 1.
                          x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
                          state.vpstr.z[state.lines[imn].set] = state.lines[idx[0]].zabs
                          widget_control, state.zabs_id, $
                            set_value=state.lines[state.curlin].zabs
                          if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
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
                  widget_control, state.ypos_id, set_value=state.ypos
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
          if eventch EQ  'q' then begin
              widget_control, ev.top, /destroy
              return
          end
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
              's': begin ; Regions
                  x_fitlls_GetReg, state
                  if state.flg_reg EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
              'S': x_fitlls_DelReg, state ; Delete a region
              'F': x_fitlls_DelReg, state, /all  ; Delete all regions
              ;; Add new LLS
              'L': begin
                  x_fitlls_newLLS, state
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
               end
              ;'Y': begin
              ;    x_fitlls_newsubLLS, state
              ;    if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              ;end
              'c': x_fitlls_newlin, state
              '2': x_fitlls_newlin, state, 12.
              '3': x_fitlls_newlin, state, 13.
              '4': x_fitlls_newlin, state, 14.
              '5': x_fitlls_newlin, state, 15.
              ;;
              'd': x_fitlls_dellin, state
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
                  ;; Plot
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
               end
              ;'3': begin ;; Add point
              ;    if n_elements(state.conti) NE 1 then $
              ;      x_fitline_continuum, state, 1L
              ;end
              ;'4': begin ;; Move a point
              ;    if n_elements(state.conti) NE 1 then $
              ;      x_fitline_continuum, state, 2L
              ;end
              ;; Exact
              'X': begin
                  x_fitlls_updfit, state, /EXACT
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              end
              ;; Colm
              'n': begin
                 widget_control, /hourglass
                  if state.nlin NE 0 then begin
                      gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
                      state.lines[gdset].N = state.lines[gdset[0]].N - 0.05
                      state.vpstr.n[state.curlin] = state.lines[gdset[0]].N
                      x_fitlls_updfit, state, EXACT=state.exact
                      widget_control, state.Ncolm_id, $
                        set_value=state.lines[gdset[0]].N
                      if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
                  endif
                  x_fitlls_updwin, state
              end
              'N': begin
                 widget_control, /hourglass
                  if state.nlin NE 0 then begin
                      gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
                      state.lines[gdset].N = state.lines[gdset[0]].N + 0.05
                      state.vpstr.n[state.curlin] = state.lines[gdset[0]].N
                      widget_control, state.Ncolm_id, $
                        set_value=state.lines[gdset[0]].N
                      x_fitlls_updfit, state, EXACT=state.exact
                      if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
                  endif
                  x_fitlls_updwin, state
              end
              ;; b-value
              'v': begin
                  gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
                  state.lines[gdset].b = state.lines[gdset[0]].b - 1.0
                  state.vpstr.b[state.curlin] = state.lines[gdset[0]].b
                  x_fitlls_updfit, state, EXACT=state.exact
                  widget_control, state.bval_id, $
                    set_value=state.lines[state.curlin].b
                  x_fitlls_updwin, state
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              end
              'V': begin
                  gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
                  state.lines[gdset].b = state.lines[gdset[0]].b + 1.0
                  state.vpstr.b[state.curlin] = state.lines[gdset[0]].b
                  widget_control, state.bval_id, $
                    set_value=state.lines[state.curlin].b
                  x_fitlls_updfit, state, EXACT=state.exact
                  x_fitlls_updwin, state
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              end
              ;; Switch lines
              '=': begin
                  if state.curlin EQ state.nset then state.curlin = 0 $
                  else state.curlin = (state.curlin+1) < state.nset 
                  x_fitlls_updwin, state
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              end
              '-': begin
                  if state.curlin EQ 0 then state.curlin = state.nset $
                  else state.curlin = (state.curlin-1) > 0
                  x_fitlls_updwin, state
                  if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
              end
              ;; non-HI lines
              'a': x_fitlls_setline, state
              ;; Output
              'O': x_fitline_output, state
              'I': x_fitline_idlout, state
                                ; Postscript
              'P': x_fitline_psfile, state  
              ;; Plot all
              'p': if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
                                ; QUIT
              else:  print, 'x_fitline: Not a valid key!' ; Nothing
          endcase
          ;; New!
          if (eventch EQ 'L' or eventch EQ 'c') AND state.nset EQ 0 then begin
              ;; Set xymnx
              x_fitlls_setxy, state
              ;; Plot
              if (flg_plt MOD 4) LT 2 then flg_plt = flg_plt + 2
          endif
      end
      ;; BUTTONS
      'IDLOUT' : x_fitline_idlout, state, /ERR
      'VPOUT' : x_fitlls_vpout, state
      'SPECPLT' : begin
          if state.nlin NE 0 then begin 
              gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
              if ngd NE 0 then zin = state.lines[gdlin[0]].zabs else zin = 0
          endif else zin = 0
          x_specplot, state.fx, state.sig, wave=state.wave, inflg=4, $
            zin=zin, /lls, /block
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

  ;; PXMNX
  if state.xymnx[0] NE state.old_xymnx[0] OR $
    state.xymnx[2] NE state.old_xymnx[2] then begin
      state.old_xymnx = state.xymnx
      state.xpmnx = x_getxpmnx(state)
  endif

; Update Plot
  if (flg_plt MOD 2) EQ 1 then x_fitlls_UpdatePlot, state
  if (flg_plt MOD 4) GE 2 then x_fitlls_PlotOth, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  updfit
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_updfit, state, FLG_PLT=flg_plt, EXACT=exact, WVOFF=wvoff

common x_specplot_lines
  ;; flg_plt
  flg_plt = 1

  ;wvoff = 10. ;; Ang
  if not keyword_set( wvoff ) then wvoff = 30. ;; Ang

  ;; Voigt
  state.fit = 1.
  ;; FIRST CHECK IF THERE ARE REGIONS BOUNDS TO APPLY
  if state.uregions EQ 0 then begin
      if keyword_set(EXACT) then begin
         ;; Determine region to generate lines in
         ;stop
         npix =  n_elements(state.wave)
         flg = replicate(0B, npix)
         for kk=0L,state.nlin-1 do begin
            wvobs = state.lines[kk].wrest*(1.+state.lines[kk].zabs) 
;            mnwv = min(state.lines[0:state.nlin-1].wrest $
;                       *(1.+state.lines[0:state.nlin-1].zabs), max=mxwv)
            mn = min(abs(state.wave - (wvobs - wvoff)), mnpx) 
            mx = min(abs(state.wave - (wvobs + wvoff)), mxpx)
            mnpx = mnpx > 0
            mxpx = mxpx < (npix-1)
            flg[mnpx:mxpx] = 1B
         endfor
         ;; Generate the lines
;         state.fit[mnpx:mxpx] = x_voigt(state.wave[mnpx:mxpx], $
;                                        state.lines[0:state.nlin-1], $
;                                        FWHM=state.FWHM)
         gdp = where(flg)
         state.fit[gdp] = x_voigt(state.wave[gdp], $
                                        state.lines[0:state.nlin-1], $
                                        FWHM=state.FWHM)
         ;; Individual
         if state.flg_indiv then begin
            for i=0, state.nset do begin
               gdlin = where(i EQ state.lines[0:state.nlin-1].set, ngd)
               if ngd EQ 0 then stop
               mnwv = min(state.lines[gdlin].wrest $
                          *(1.+state.lines[gdlin].zabs), max=mxwv)
               mn = min(abs(state.wave - mnwv + wvoff), mnpx)
               mx = min(abs(state.wave - mxwv - wvoff), mxpx)
               state.ifit[i,mnpx:mxpx] = x_voigt(state.wave[mnpx:mxpx], $
                                                 state.lines[gdlin], $
                                                 FWHM=state.FWHM)
            endfor
         endif
      endif else begin
         stop ;; I do not trust x_allvoigt.  Reconsider JXP 18 September 2013
         state.fit = x_allvoigt(state.wave, state.lines[0:state.nlin-1], $
                                SIGMA=state.FWHM)
         for i=0, state.nset do begin
            gdlin = where(i EQ state.lines[0:state.nlin-1].set, ngd)
            if ngd EQ 0 then stop
            state.ifit[i,0:n_elements(state.fit)-1] = x_allvoigt(state.wave, $
                                                                 state.lines[gdlin], $
                                                                 SIGMA=state.FWHM)
         endfor
      endelse
  endif else begin
      if keyword_set(EXACT) then begin
          for i=0,state.nreg-1 do begin
              gd = where(state.wave GE state.regions[0,i] AND $
                         state.wave LE state.regions[1,i],ngd)
              if ngd LE 33 then continue
              state.fit[gd] = x_voigt(state.wave[gd], state.lines[0:state.nlin-1], $
                                      FWHM=state.FWHM)
              for j=0, state.nset do begin
                  gdlin = where(j EQ state.lines[0:state.nlin-1].set, ngdlin)
                  if ngdlin EQ 0 then stop
                  state.ifit[j,gd] = x_voigt(state.wave[gd], state.lines[gdlin], $
                                             FWHM=state.FWHM)
              endfor
          endfor
      endif else begin ;; Use an ok approximation
          for i=0,state.nreg-1 do begin
              gd = where(state.wave GE state.regions[0,i] AND $
                         state.wave LE state.regions[1,i],ngd)
              if ngd LE 33 then continue
              stop ;; I do not trust x_allvoigt.  Reconsider JXP 18 September 2013
              state.fit[gd] = x_allvoigt(state.wave[gd], state.lines[0:state.nlin-1], $
                                         SIGMA=state.FWHM)
              for j=0,state.nset do begin
                  gdlin = where(j EQ state.lines[0:state.nlin-1].set, ngdlin)
                  if ngdlin EQ 0 then stop
                  state.ifit[j,gd] = x_allvoigt(state.wave[gd], state.lines[gdlin], $
                                                SIGMA=state.FWHM)
              endfor
          endfor
      endelse
  endelse

  ;; Add in Continuum opacity
  alin = state.lines[0:state.nlin-1]
  lya = where( abs(alin.wrest-1215.6701) LT 1e-3 and $ 
               alin.N GT 16.5, nlls)
  constc = 2.9979246d+10
  consth = 6.6260690d-27
  constev = 1.6021772d-12
  for jj=0L,nlls-1 do begin
     low_wv = where(state.wave LT (1+alin[lya[jj]].zabs)*911.7641, npix)
     if npix GT 10 then begin
         tau_eng = consth*constc/((state.wave[low_wv]/(1+alin[lya[jj]].zabs))*1d-8) / constev ; eV
         tau_LL = 10.^(alin[lya[jj]].N) * x_photocross(1,1, tau_eng) 
;     tau_LL = 10.^(16+tau_id*delta_N) * (tau_wv/911.7641)^3 $
;              / (1.5748d17)
         state.fit[low_wv] = state.fit[low_wv] * exp(-tau_LL)
      endif
  endfor

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  SetLine, ala JMO
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_autocont, state

common x_specplot_lines
  nbin = floor((state.wave[state.npix-1]-state.wave[0])/25.0)

  cxpos = fltarr(nbin+2)
  cypos = fltarr(nbin+2)

  cxpos[0] = state.wave[0]
  cxpos[nbin+1] = state.wave[state.npix-1]

  cypos[0] = state.fx[0]
  cypos[nbin+1] = state.fx[state.npix-1]

  for qq=1,nbin do begin
      cxpos[qq] = state.wave[0] + 25.0*(qq)
  endfor
 
  for qq=1,nbin do begin
      a = where(state.wave GT cxpos[qq-1] AND state.wave LT cxpos[qq+1])
      tmpf = state.fx[a]
      s = sort(tmpf)
      tmpf2 = tmpf[s]
      val = floor(0.75*n_elements(tmpf))

      cypos[qq] = tmpf2[val]
  endfor

  for qq=0, n_elements(cxpos)-1 do begin
    state.xpos = cxpos[qq]
    state.ypos = cypos[qq]
    x_fitline_continuum, state, 1
  endfor

  x_fitlls_updateplot, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  SetLine
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_setline, state

common x_specplot_lines

  xpt = state.xpos
  flg_lines = 1
  llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
  lines = x_setllst(llist,0)
  ; Set Line with gui
  setwave = x_slctline(lines, /ISM)
  ;; Get lines
  gd = where(abs(state.all_lin.wrest-setwave) LE 0.01,ngd)
  if ngd EQ 0 then stop
  tmp = state.all_lin[gd]
  ion = state.all_lin[gd].ion
  crd = strsplit(ion,' ',/extract)

  use = where(strmatch(strtrim(state.all_lin.ion,2),'*'+crd[0]+' *') EQ 1, nuse)
  if nuse EQ 0 then stop

  tmplin = state.all_lin[use]
  tmplin.N = 13.
  tmplin.b = 15.0
  tmplin.zabs = (xpt / setwave) - 1.

  ;; NSET
  state.nset = state.nset + 1
  nlin = n_elements(tmplin)

  ;; Add to state
  state.lines[state.nlin:state.nlin+nlin-1] = tmplin
  state.lines[state.nlin:state.nlin+nlin-1].set = state.nset

  ;; Add to vpstr
  tmp = crd[0]
  bg = strpos(tmp,'I')
  nm = strmid(tmp,0,bg)
  nm2 = strmid(tmp,bg)
  
  if bg EQ 1 then $
      state.vpstr.ion[state.vpstr.nion] = nm+' '+nm2 else $
      state.vpstr.ion[state.vpstr.nion] = tmp
  state.vpstr.z[state.vpstr.nion] = state.lines[state.nlin].zabs
  state.vpstr.n[state.vpstr.nion] = state.lines[state.nlin].N
  state.vpstr.b[state.vpstr.nion] = state.lines[state.nlin].b
  state.vpstr.nion = state.vpstr.nion + 1

  ;; Update window
  widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
  widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
  widget_control, state.bval_id, set_value=state.lines[state.nlin].b

  ;; Current line and increment
  state.curlin = state.nset
  state.nlin = state.nlin + nlin

  ;; Update fit
  x_fitlls_updfit, state, EXACT=state.exact

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  dellin
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_dellin, state

common x_specplot_lines
  ;; No lines
  if state.nset LT 0 then return

  ;; 1 line
  if state.nset EQ 0 then begin
      state.nset = -1
      state.nlin = 0
      state.curlin = state.nset
      state.vpstr.nion = 0
      state.vpstr.ion[0] = ''
      state.vpstr.z[0] = 0.
      state.vpstr.n[0] = 0.
      state.vpstr.b[0] = 0.
      return
  endif

  ;; Multiple lines
  gdlin = where(state.lines.set NE state.curlin, ngd)
  state.lines[0:ngd-1] = state.lines[gdlin]

  ;; VPSTR
  for i=state.curlin, state.vpstr.nion-2 do begin
      state.vpstr.ion[i] = state.vpstr.ion[i+1]
      state.vpstr.z[i] = state.vpstr.z[i+1]
      state.vpstr.n[i] = state.vpstr.n[i+1]
      state.vpstr.b[i] = state.vpstr.b[i+1]
  endfor
  state.vpstr.ion[state.vpstr.nion-1] = ''
  state.vpstr.z[state.vpstr.nion-1] = 0.
  state.vpstr.n[state.vpstr.nion-1] = 0.
  state.vpstr.b[state.vpstr.nion-1] = 0.
  state.vpstr.nion = state.vpstr.nion - 1

  ;; Reset set
  if state.curlin NE state.nset then begin
      big = where(state.lines.set GT state.curlin)
      state.lines[big].set = state.lines[big].set - 1
  endif
  state.nset = state.nset - 1
  state.curlin = state.curlin < state.nset

  ;; Update fit
  x_fitline_updfit, state, EXACT=state.exact

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  newlls
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_newlls, state
  
common x_specplot_lines
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
    state.lines[state.nlin:state.nlin+nlin-1].set = state.nset

    ;; Update window
    widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
    widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
    widget_control, state.bval_id, set_value=state.lines[state.nlin].b

    ;; Add to vpstr
    state.vpstr.ion[state.vpstr.nion] = 'H I'
    state.vpstr.z[state.vpstr.nion] = ztmp
    state.vpstr.n[state.vpstr.nion] = state.lines[state.nlin].N
    state.vpstr.b[state.vpstr.nion] = state.lines[state.nlin].b
    state.vpstr.nion = state.vpstr.nion + 1

    ;; Current line and increment
    state.curlin = state.nset
    state.nlin = state.nlin + nlin
   
    ;; Update fit
    x_fitlls_updfit, state, EXACT=state.exact
    return
 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  newlls  -- Only through Lyd
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_newsublls, state
  
common x_specplot_lines
    widget_control, /hourglass
    ;; Get z
    xpt = state.xpos
    ztmp = (xpt / 1215.6701) - 1.
    ;; NSET
    state.nset = state.nset + 1
  
    ;; Get LLS lines
    tmplin = x_mknewlls(ztmp, set=state.nset)
    nlin = n_elements(tmplin)

    ;; Cut
    tmplin = tmplin[nlin-4:*]
    nlin = 4L

    ;; Add to state
    state.lines[state.nlin:state.nlin+nlin-1] = tmplin
    state.lines[state.nlin:state.nlin+nlin-1].set = state.nset

    ;; Update window
    widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
    widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
    widget_control, state.bval_id, set_value=state.lines[state.nlin].b

    ;; Add to vpstr
    state.vpstr.ion[state.vpstr.nion] = 'H I'
    state.vpstr.z[state.vpstr.nion] = ztmp
    state.vpstr.n[state.vpstr.nion] = state.lines[state.nlin].N
    state.vpstr.b[state.vpstr.nion] = state.lines[state.nlin].b
    state.vpstr.nion = state.vpstr.nion + 1

    ;; Current line and increment
    state.curlin = state.nset
    state.nlin = state.nlin + nlin
   
    ;; Update fit
    x_fitlls_updfit, state, EXACT=state.exact
    return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  fitlls_newlin
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_newlin, state, NHI

common x_specplot_lines
  ;; Grab x,y pos
  xpt = state.xpos

  if not keyword_set(NHI) then NHI = 13.

  ;; Setup HI
  tmp = x_setline(1215.670d)
  tmp2 = x_setline(1025.7223d)
  tmp3 = x_setline(972.5368d)
  tmp4 = x_setline(949.7431d)

  ;; SET
  state.nset = state.nset + 1
  state.curlin = state.nset
  idx = state.nlin + lindgen(4)
 
  ;; Get z
  state.lines[idx] = [tmp,tmp2,tmp3,tmp4]
  state.lines[idx].zabs = (xpt / 1215.6701) - 1.
  state.lines[idx].N = NHI
  state.lines[idx].b = 25.0
  state.lines[idx].set = state.nset

  ;; Add to vpstr
  state.vpstr.ion[state.vpstr.nion+lindgen(4)] = 'H I'
  state.vpstr.z[state.vpstr.nion+lindgen(4)] = state.lines[state.nlin].zabs
  state.vpstr.n[state.vpstr.nion+lindgen(4)] = state.lines[state.nlin].N
  state.vpstr.b[state.vpstr.nion+lindgen(4)] = state.lines[state.nlin].b
  state.vpstr.nion = state.vpstr.nion + 4

  ;; Update text windows
  widget_control, state.zabs_id, set_value=state.lines[state.nlin].zabs
  widget_control, state.Ncolm_id, set_value=state.lines[state.nlin].N
  widget_control, state.bval_id, set_value=state.lines[state.nlin].b

  ;; Current line and increment
  state.nlin = state.nlin + 4

  ;; Update fit
  x_fitlls_updfit, state, EXACT=state.exact
  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  vpout
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_vpout, state

common x_specplot_lines
    g_mkvpin, state.vpstr, FIL=state.vpfil

    return 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  GetReg
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_GetReg, state
common x_specplot_lines
    ;;
    if state.flg_reg EQ 0 then begin
        state.flg_reg = 1
        state.regions[0,state.nreg] = xgetx_plt(state,/strct)
        state.vpstr.reg_beg[state.vpstr.nreg] = state.regions[0,state.nreg]
    endif else begin
        state.flg_reg = 0
        state.regions[1,state.nreg] = xgetx_plt(state,/strct)
        state.vpstr.reg_end[state.vpstr.nreg] = state.regions[1,state.nreg]
        if state.regions[1,state.nreg] LT state.regions[0,state.nreg] then begin
            tmp = state.regions[1,state.nreg]
            state.regions[1,state.nreg] = state.regions[0,state.nreg]
            state.regions[0,state.nreg] = tmp
            state.vpstr.reg_end[state.vpstr.nreg] = state.regions[1,state.nreg]
            state.vpstr.reg_beg[state.vpstr.nreg] = state.regions[0,state.nreg]
        endif
        if state.regions[1,state.nreg] LT state.midwv then $
            state.vpstr.fluxfil[state.vpstr.nreg] = state.regfil1 else $
            state.vpstr.fluxfil[state.vpstr.nreg] = state.regfil2
        state.nreg = state.nreg + 1
        state.vpstr.nreg = state.vpstr.nreg + 1
    endelse

    return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Delete Region
;;;;;;;;;;;;;;;;;;;;
pro x_fitlls_DelReg, state, ALL=all

common x_specplot_lines
  if state.nreg EQ 0 then return

  ;; Deletes 1 or All regions

   if keyword_set(ALL) OR state.nreg EQ 1 then begin
       for i=0, state.nreg-1 do begin
           state.regions[0,i] = 0.
           state.regions[1,i] = 0.
           state.vpstr.reg_beg[i] = 0.
           state.vpstr.reg_end[i] = 0.
           state.vpstr.fluxfil[i] = ''
       endfor
       state.nreg = 0
       state.vpstr.nreg = 0
   endif else begin
       xpos = xgetx_plt(state, /strct)
       fndreg = where( xpos GE state.regions[0,0:state.nreg-1] AND $
                       xpos LE state.regions[1,0:state.nreg-1], count)
       if count EQ 0 then begin
           widget_control, state.error_msg_id, set_value='No regions found!'
           return
       endif
       bad = fndreg[0]
       for i=bad,state.nreg-2 do begin
           state.regions[0,i] = state.regions[0,i+1]
           state.regions[1,i] = state.regions[1,i+1]
           state.vpstr.reg_beg[i] = state.vpstr.reg_beg[i+1]
           state.vpstr.reg_end[i] = state.vpstr.reg_end[i+1]
           state.vpstr.fluxfil[i] = state.vpstr.fluxfil[i+1]
       endfor
       state.regions[0,state.nreg-1] = 0.
       state.regions[1,state.nreg-1] = 0.
       state.vpstr.reg_beg[state.nreg-1] = 0.
       state.vpstr.reg_end[state.nreg-1] = 0.
       state.vpstr.fluxfil[state.nreg-1] = ''
       state.nreg = state.nreg-1
       state.vpstr.nreg = state.vpstr.nreg-1
   endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  updwin
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_updwin, state

common x_specplot_lines
   gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
   if ngd NE 0 then begin
       widget_control, state.zabs_id, set_value=state.lines[gdset[0]].zabs
       widget_control, state.Ncolm_id, set_value=state.lines[gdset[0]].N
       widget_control, state.bval_id, set_value=state.lines[gdset[0]].b
   endif

end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_UpdatePlot, state
  
common x_specplot_lines
; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

;  if state.flg_dum EQ 1 then stop
;  state.flg_dum = 0
;  stop
  plot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.fx[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, $
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
    oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.sig[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, color=clr.red

; REGIONS
  if state.nreg GT 0 then begin
      for i=0,state.nreg-1 do begin
          gd = where(state.wave GE state.regions[0,i] AND state.wave LE state.regions[1,i],ngd)
          if ngd EQ 0 then continue
          oplot, state.wave[gd], state.fx[gd], color=clr.red, psym=state.psym
      endfor
  endif

  ;; FIT
  if state.nlin NE 0 then begin
      if state.indiv EQ 0 then begin
          oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
            state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
            state.fit[state.xpmnx[0]:state.xpmnx[1]], color=clr.blue, thick=3
      endif else begin
          for i=0L, state.nset do begin
              oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
                state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
                state.ifit[i,state.xpmnx[0]:state.xpmnx[1]], color=clr.blue, thick=3
          endfor
          oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
            state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
            state.ifit[state.curlin,state.xpmnx[0]:state.xpmnx[1]], color=clr.red
          oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
            state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
            state.fit[state.xpmnx[0]:state.xpmnx[1]], color=clr.purple
      endelse
      ;; Mark all lines
      for i=0L,state.nlin-1 do begin
          oplot, replicate( (state.lines[i].zabs+1.)*$
                            state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.blue, linestyle=1
      endfor
      ;; Mark current set
      gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin, ngd)
      for i=0L,ngd-1 do begin
          oplot, replicate( (state.lines[gdlin[i]].zabs+1.)*$
                            state.lines[gdlin[i]].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.red, linestyle=1
      endfor
  endif

; CONTINUUM
  ;; Line
  oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.conti[state.xpmnx[0]:state.xpmnx[1]], color=clr.purple, $
    linestyle=1
  ;; Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.orange, symsize=5
  endif

  
end

;;;;
;; Other windows
;;;;

pro x_fitlls_PlotOth, state
  
common x_specplot_lines
  if state.psfile EQ 1 then return
  clr = getcolor(/load)

  ;; Loop on Windows
  gd = where(state.flgwin NE 0, ngd)
  for ii=0L,ngd-1 do begin
      jj=gd[ii]

      ;; Set Window
      wset, jj

      ;; Reset y-range by continuum
      contiv = state.conti[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]]
      mxc = max(contiv)
      state.all_xymnx[1,jj] = -1. * 0.05 * mxc
      state.all_xymnx[3,jj] = 1.1 * mxc

      ;if jj EQ 1 then stop

      ;; Plot Data
      plot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
        state.fx[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
        psym=state.psym, position=state.pos, $
        xrange=[state.all_xymnx[0,jj], $
                state.all_xymnx[2,jj]], $
        yrange=[state.all_xymnx[1,jj], $
                state.all_xymnx[3,jj]], xstyle=1, ystyle=1, $
        xtitle='!17Wavelength', ytitle='Flux', $
        title=state.title, $
        background=clr.white, $
        color=clr.black, $
        xcharsize=1.7, $
        ycharsize=1.7

      ;; Plot Error array
      if state.flg_sig EQ 1 then $
        oplot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
        state.sig[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
        psym=state.psym, color=clr.red
      
      ;; FIT
      if state.nlin NE 0 then begin
          if state.indiv EQ 0 then begin
              oplot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                state.conti[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]]* $
                state.fit[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                color=clr.green
          endif else begin
              for i=0L, state.nset do begin
                  oplot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                    state.conti[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]]* $
                    state.ifit[i,state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                    color=clr.green
              endfor
              oplot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                state.conti[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]]* $
                state.ifit[state.curlin,state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                color=clr.red
              oplot, state.wave[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                state.conti[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]]* $
                state.fit[state.all_xpmnx[0,jj]:state.all_xpmnx[1,jj]], $
                color=clr.blue
          endelse
          ;; Mark all lines
          for i=0L,state.nlin-1 do begin
              oplot, replicate( (state.lines[i].zabs+1.)*$
                                state.lines[i].wrest, 2), $
                [state.all_xymnx[1,jj], state.all_xymnx[3,jj]], $
                color=clr.blue, linestyle=1
          endfor
          ;; Mark current set
          gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin, ngd)
          for i=0L,ngd-1 do begin
              oplot, replicate( (state.lines[gdlin[i]].zabs+1.)*$
                                state.lines[gdlin[i]].wrest, 2), $
                [state.all_xymnx[1,jj], state.all_xymnx[3,jj]], $
                color=clr.red, linestyle=1
          endfor
      endif
      
; CONTINUUM
      ;; Line
      oplot, state.wave, state.conti, color=clr.purple, linestyle=1
      ;; Points
      if state.cstr.npts NE 0 then begin
          gdc = where(state.cstr.msk EQ 1)
          oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
            color=clr.orange, symsize=5
      endif
  endfor
  
  return

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_Reset, state

common x_specplot_lines

; Plotting
  state.xymnx = state.svxymnx
  state.old_xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;
;  VP IN
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_vpin, state, vpin

common x_specplot_lines
  a = findfile(vpin, count=na)
  if na EQ 0 then begin
      print, 'x_fitlls: FILE ', vpin, ' does not exist!'
      stop
  endif

  g_vpparse, vpin, VPSTR=invpstr

  state.vpstr = invpstr

  name = invpstr.ion[0:invpstr.nion-1]

  bd = where(strmid(name,1,1) EQ ' ', nbd)
  if nbd NE 0 then $     ;  fixing names to match 'all_lin'
      name[bd] = strmid(name[bd],0,1)+strmid(name[bd],2)

  lls = where(name EQ 'HI',nlls)
  metals = where(name NE 'HI', nmetals)
  nmet = 0
  for i=0,nmetals-1 do begin
    use = where(strmatch(strtrim(state.all_lin.ion,2),'*'+name[metals[i]]+' *') EQ 1,nuse)
    nmet = nmet+nuse
  endfor

  nlines = nlls*19+nmet
  nset = n_elements(name)

  ;; SET LINES

  nlin=0
  for i=0,nset-1 do begin
      if name[i] EQ 'HI' then begin
          state.nset = state.nset + 1
          tmplin = x_mknewlls(invpstr.z[i],set=state.nset)
          num = n_elements(tmplin)
          state.lines[state.nlin:state.nlin+num-1] = tmplin
          state.lines[state.nlin:state.nlin+num-1].set = state.nset
          state.lines[state.nlin:state.nlin+num-1].N = invpstr.n[i]
          state.lines[state.nlin:state.nlin+num-1].b = invpstr.b[i]
          state.curlin = state.nset
          state.nlin = state.nlin + num
          x_fitlls_setxy, state
      endif else begin
          use = where(strmatch(strtrim(state.all_lin.ion,2),'*'+name[i]+' *') EQ 1,nuse)
          state.nset = state.nset+1
          tmplin = state.all_lin[use]
          num = n_elements(tmplin)
          state.lines[state.nlin:state.nlin+num-1] = tmplin
          state.lines[state.nlin:state.nlin+num-1].set = state.nset
          state.lines[state.nlin:state.nlin+num-1].N = invpstr.n[i]
          state.lines[state.nlin:state.nlin+num-1].zabs = invpstr.z[i]
          state.lines[state.nlin:state.nlin+num-1].b = invpstr.b[i]
          state.curlin = state.nset
          state.nlin = state.nlin + num
      endelse
  endfor

  ;; SET REGIONS
  for i=0, invpstr.nreg-1 do begin
      state.regions[0,i] = invpstr.reg_beg[i]
      state.regions[1,i] = invpstr.reg_end[i]
      state.nreg = state.nreg+1
  endfor

  state.uregions=1

  x_fitlls_updfit, state, EXACT=state.exact
  x_fitlls_UpdatePlot, state

  return

end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_inifit, state, inifit

common x_specplot_lines

  ;; File?
  a = findfile(inifit, count=na)
  if na EQ 0 then begin
      print, 'x_fitlls: FILE ', inifit, ' does not exist!'
      stop
  endif

  ;; Restore
  restore, inifit

  ;; Error
;  if keyword_set( SIG ) then begin
;      if sig[0] NE 0. then begin
;          state.crude_val = sig
;          state.crude = 1
;          widget_control, state.crude_id, set_value=1
;      endif
;  endif

  ;; Continuum
  if keyword_set(conti) then state.conti[*] = conti else state.conti[*] = 1.

  ;; conti_str
  if keyword_set( CSTR ) then state.cstr = cstr

  ;; Lines
  a = where(lines.zabs GT 0., nlin)
  state.nlin = nlin
  state.lines[0:state.nlin-1] = lines[0:nlin-1]
  state.curlin = 0L

  ;; Zoom in 
  cen = (1.+state.lines[state.curlin].zabs)*1215.6701
  gd = where(abs(state.wave-cen) LT 200.)
  mx = max(state.conti[gd])
  state.xymnx = [cen-200., -0.1*mx, cen+200., mx*2.]

  ;; Zabs, N
  widget_control, state.Ncolm_id, set_value=state.lines[state.curlin].N
  widget_control, state.bval_id, set_value=state.lines[state.curlin].b
  widget_control, state.zabs_id, set_value=state.lines[state.curlin].zabs

  ;; nset
  state.nset = state.lines[state.nlin-1].set

  ;; 
  x_fitlls_setxy, state
  x_fitlls_updfit, state, EXACT=state.exact

end

;;;;;;;
;; Open Windows (as necessary)
;;;;;;;

pro x_fitlls_openwin, state

common x_specplot_lines
  ;; Check
  gd = where(state.flgwin, ngd)
  for ii=0L,ngd-1 do begin
      ;; Open
      case gd[ii] of
          0: lbl='Lya'
          1: lbl='Lyb'
          2: lbl='Lyg'
          3: lbl='Lyd'
          4: lbl='Lye'
          5: lbl='Lyf'
          19:lbl= 'LL'
          else: stop
      endcase
      if x_chkwindow(gd[ii]) NE 1 then window, gd[ii], title=lbl
  endfor
  return
end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_setxy, state

common x_specplot_lines
  ;; Grab redshift
  id = where(state.lines.set EQ 0L)
  z = state.lines[id[0]].zabs
  for qq=0L,19 do begin
      case qq of
          0: lmb = 1215.6701 ;;  Lya
          1: lmb = 1025.722
          2: lmb = 972.5368
          3: lmb = 949.7431
          4: lmb = 937.8035
          5: lmb = 930.7483
          19: lmb = 914.039
          else:
      endcase
      ;; Velocity range
      case qq of
          19: begin
              vmin = -400.
              vmax = 800.
          end
          else: begin
              ;vmin = -200.
              ;vmax = 200.
              vmin = state.lym_vmnx[0]
              vmax = state.lym_vmnx[1]
          end
      endcase
      ;; Setup
      state.all_xymnx[0,qq] = lmb*(1+z) * (1. + vmin/3e5)
      state.all_xymnx[2,qq] = lmb*(1+z) * (1. + vmax/3e5)
      state.xymnx = state.all_xymnx[*,qq]
      state.all_xpmnx[*,qq] = x_getxpmnx(state)
;      id = where(state.wave GT state.all_xymnx[0,qq] AND $
;                 state.wave LT state.all_xymnx[2,qq],npix)
;      srt = sort(state.fx[id])
;      med = state.fx[id[srt[0.9*npix]]]

      contiv = state.conti[state.all_xpmnx[0,qq]:state.all_xpmnx[1,qq]]
      mxc = max(contiv)
      state.all_xymnx[1,qq] = -1. * 0.05 * mxc
      state.all_xymnx[3,qq] = 1.1 * mxc
      ;state.xymnx = state.all_xymnx[*,qq]
      ;; Kludge?
      ;state.all_xpmnx[*,qq] = x_getxpmnx(state)
      ;state.all_xymnx[1,qq] = state.ypos
      ;state.all_xymnx[3,qq] = state.ypos
      ;state.xymnx = state.all_xymnx[*,qq]
  endfor
  state.xymnx = state.all_xymnx[*,state.curwin]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fitlls_idx, state, evid, flg, id, LBL=lbl

common x_specplot_lines
  case flg of
      0: lbl = state.tranwin[evid]
      1: lbl = state.othwin[evid]
      else: stop
  endcase
 
  ;; Index
  case lbl of
      'Lya': id= 0 
      'Lyb': id= 1
      'Lyg': id= 2
      'Lyd': id= 3
      'Lye': id= 4
      'Lyf': id= 5
      'LL':  id= 19
      else: stop
  endcase

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; MAIN ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; /INDIV -- Calculate for the individual screens
pro x_fitlls, yin, ysin, XSIZE=xsize, YSIZE=ysize, TITLE=title, $
              WAVE=wave, BLOCK=block, FWHM=fwhm, INIFIT=inifit, $
              INFLG=inflg, ZIN=zin, EXFIL=exfil, FX=fx, VPIN=vpin, $
              VPOUTFIL=vpoutfil, YOUT=yout, EXOUT=exout, $
              CONTI_FIL=conti_fil, exact=exact, WVOFF=wvoff, $
              LYM_VMNX=LYM_VMNX, INDIV=indiv

common x_specplot_lines
  x_specplot_initcommon
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_fitlls, flux, [error], XSIZE=,YSIZE=, TITLE=, WAVE=, '
    print, '            INIFIT=, INFLG=, FWHM=, /BLOCK, /exact, /INDIV, WVOFF=, LYM_VMNX=) [v1.2]'
    return
  endif 
;  Optional Keywords

  if not keyword_set(EXACT) then exact = 0

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = round(ssz[0]*2/3.)
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200L
  if not keyword_set( FWHM ) then    fwhm = 3.
  if not keyword_set( VPOUTFIL ) then vpoutfil = 'fort.13'

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0
  ydat = x_readspec(yin, INFLG=inflg, head=head, $
                  NPIX=npix, WAV=xdat, FIL_SIG=ysin, SIG=ysig)
  if not keyword_set(YSIG) then ysig = replicate(1., npix)

  if keyword_set( EXFIL ) then begin
      ydat2 = x_readspec(exfil, INFLG=inflg, head=head, $
                        NPIX=npix2, WAV=xdat2)
      gd = where(ydat NE 0,ngd)
      if ngd GT 0 then begin 
          ydat = ydat[gd]
          xdat = xdat[gd]
      endif
      mx1 = max(xdat)
      mn2 = min(xdat2)
      wv1 = (mx1+mn2)/2.
      ;gd1 = where(xdat GE mn2,ngd1)
      ;gd2 = where(xdat2 LE mx1 AND xdat2 GT 0.,ngd2)
      ;wv1 = xdat[gd1[long((ngd1+ngd2)/2.)]]
      use1 = where(xdat LT wv1)
      use2 = where(xdat2 GE wv1)
      ydat = [ydat[use1],ydat2[use2]]
      xdat = [xdat[use1],xdat2[use2]]
      srt = sort(xdat)
      xdat = xdat[srt]
      ydat = ydat[srt]
      npix = n_elements(ydat)
      ysig = replicate(0., npix)
  endif else wv1 = 0.


  tmp1 = { newabslinstrct }

;  tmp2 = { conti_str, $
  tmp2 = { npts: 0L, $
           xval: fltarr(1000), $
           yval: fltarr(1000), $
           msk: lonarr(1000) }

  tmp3 = {vpstrct}

  tmp4 = xmrdfits(getenv('XIDL_DIR')+'/Spec/Lines/all_lin.fits', 1,/silent)

  if not keyword_set(YOUT) then yout = yin
  if not keyword_set(EXOUT) AND keyword_set(EXFIL) then exout = exfil $
  else exout = 'unk.dat'


;JO DID THIS
wv1=0.

;    STATE
  state = { fx: ydat, $
            wave: xdat, $
            sig: ysig, $
            regions: fltarr(2,500), $  ;  array of regions for vpfit output
            nreg: 0L, $
            flg_reg: 0L, $
            regfil1: yout, $
            regfil2: exout, $ 
            all_lin: tmp4, $
            npix: npix, $
            midwv: wv1, $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            nlin: 0, $ ; DLA flag
            vpstr: tmp3, $ ;  VPSTRCT for inputing into vpfit
            vpfil : vpoutfil, $
            EXACT: exact,$
            wvoff: 30.,$ ;; Ang
            lym_vmnx: [-200, 200.], $ ;; Plot range for Lyman series
            INDIV: 0,$
            UREGIONS: 0,$
            GCONT: 0,$
            lines: replicate(tmp1,1600), $ ; DLA flag
            fit: fltarr(n_elements(ydat)) + 1., $
            ifit: fltarr(100,n_elements(ydat)) + 1., $
            flg_indiv: keyword_set(INDIV), $
            conti: replicate(1.,npix), $
            cstr: tmp2, $
            curlin: 0L, $
            curwin: 0L, $  ;; Transition for Main Window  0=Lya; 1=Lyb
            othwin: ['Lya', 'Lyb', 'Lyg', 'Lyd','Lye','Lyf', 'LL'], $
            tranwin: ['Lya', 'Lyb','Lyg', 'Lyd','Lye','Lyf','LL'], $
            flgwin: bytarr(20), $
            nset: -1, $
            crude: 0, $  ; Crude Error
            crude_val: [-0.1, 0.1], $
            xpos: 0.d, $
            ypos: 0.d, $
            flg_dum: 1, $
            psfile: 0, $ ; Postscript
            svxymnx: [0., $     ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      float(n_elements(ydat)-1), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $ ; ymax
            xymnx: fltarr(4), $
            old_xymnx:fltarr(4), $
            all_xymnx:fltarr(4,20), $
            tmpxy: fltarr(4), $
            xpmnx: lonarr(2), $
            all_xpmnx: lonarr(2,20), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            bval_id: 0L, $
            Ncolm_id: 0L, $
            lines_id: 0L, $
            crude_id: 0L, $
            exact_id: 0L, $
            indiv_id: 0L, $
            reg_id: 0L, $
            crudehi_id: 0L, $
            crudelo_id: 0L, $
            mainwin_id: 0L, $
            othwin_id: 0L, $
            error_msg_id: 0L $
          }

; WAVE
  if keyword_set(WAVE) then state.wave = wave
  if keyword_set(LYM_VMNX) then state.lym_vmnx = LYM_VMNX
  if keyword_set(WVOFF) then state.wvoff = wvoff
  if keyword_set( YSIN ) then state.flg_sig = 1

  if keyword_set(CONTI_FIL) then begin
      readcol, conti_fil, x_spl, y_spl, /sile
      nspl = n_elements(x_spl)
      state.cstr.npts = nspl
      state.cstr.xval[0:nspl-1] = x_spl
      state.cstr.yval[0:nspl-1] = y_spl
      state.cstr.msk[0:nspl-1] = 1L
  endif

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE
; Set svxymnx[0,2]

  state.svxymnx[0] = min(state.wave)
  state.svxymnx[2] = max(state.wave)

;    Title
  if size(yin, /type) EQ 7 then state.title = yin
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'x_fitlls', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['x_fitlls', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'T -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = 'A -- Al III doublet'
  state.help[19] = 'C -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'N -- AODM'
  state.help[22] = 'q -- Quit '


;;;;;;;;;
;  Toolbar

; zabs

  linbar = WIDGET_BASE( toolbar, /row, /frame, /base_align_center,$
                           /align_center)
  state.zabs_id = cw_field(linbar, title='zabs', value=0., /floating, $
                           /column, /return_events, xsize=10, uvalue='ZABS')
  state.Ncolm_id = cw_field(linbar, title='Ncolm', value=zabs, /floating, $
                           /column, xsize=7, /return_events, uvalue='NCOLM')
  state.bval_id = cw_field(linbar, title='bval', value=zabs, /floating, $
                           /column, xsize=6, /return_events, uvalue='BVAL')



; continuum
;  state.conti_id = cw_field(toolbar, title='conti',
;  value=state.conti, /floating, $
;                           /column, xsize=12, /return_events, uvalue='CONTI')



; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Main Window
  state.mainwin_id = WIDGET_LIST(toolbar, VALUE=state.tranwin, $
                                  uvalue='MAINWIN',  ysize=4)
  widget_control, state.mainwin_id, set_list_select=0

; Other windows
  state.othwin_id = WIDGET_LIST(toolbar, VALUE=state.othwin, $
                                  uvalue='OTHWIN',  ysize=4, /MULTIPLE)

  ;; Initialize
  widget_control, state.othwin_id, set_list_select=[0,1]
  if x_chkwindow(0) NE 1 then window, 0, title='Lya'
  if x_chkwindow(1,get_all_active=gaa) NE 1 then window, 1, title='Lyb'
  state.flgwin[0:1] = 1B
  if gaa[0] NE -1 then for qq=0L,n_elements(gaa)-1 do wshow, gaa[qq]

; Exact
  ;exbase_id = widget_base(toolbar, /column, /align_center, frame=2)
  exbase_id = widget_base(toolbar, /row, /align_center, frame=2)
  state.exact_id = CW_BGROUP(exbase_id, ['N', 'Y'], label_top='Ex?', $
                              row=2, /exclusive, /no_release, $
                              set_value=state.exact,  uvalue='EXACT')
  state.indiv_id = CW_BGROUP(exbase_id, ['N', 'Y'], label_top='Idividual Fits?', $
                              row=2, /exclusive, /no_release, $
                              set_value=state.indiv,  uvalue='INDIV')
  state.reg_id = CW_BGROUP(exbase_id, ['N', 'Y'], label_top='Fit Regions?', $
                              row=2, /exclusive, /no_release, $
                              set_value=state.uregions,  uvalue='UREGIONS')

  
  
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
  buttbase = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  specplt = WIDGET_BUTTON(buttbase, value='Specplt',uvalue='SPECPLT', $
                          /align_right)
  idlout = WIDGET_BUTTON(buttbase, value='IDL Out',uvalue='IDLOUT', /align_right)
  vpout = WIDGET_BUTTON(buttbase, value='VP Out', uvalue='VPOUT', /align_right)
  gcont = WIDGET_BUTTON(buttbase, value='Auto Cont', uvalue='GCONT', /align_right)
  done = WIDGET_BUTTON(buttbase, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
;  if not keyword_set( NOCTB ) then loadct, 2, /silent

; HELP
  state.help[0] = '  :::Help Menu for x_fitlls::: '
  state.help[1] = 'RMB -- Recenter line'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'lrbt -- Set Left Right Bottom Top'
  state.help[4] = 'T -- Set ymax to 1.1'
  state.help[5] = 'zz -- Zoom corners'
  state.help[6] = 'io -- Zoom in/out'
  state.help[7] = '{}[] -- Pan'
  state.help[8] = 'w -- Reset the screen'
  state.help[9] = 'H -- This widget'
  state.help[11] = 'c -- Add new Lya line'
  state.help[12] = 'L -- Add new LLS'
  state.help[13] = 'd -- Delete current line'
  state.help[14] = 'C -- Set continnuum' 
  state.help[15] = 'nN -- Adjust colm'
  state.help[16] = 'vV -- Adjust bvalue'
  state.help[17] = '=- -- Loop through lines'
  state.help[19] = 'P -- Print to postrscipt'
  state.help[20] = 'I -- IDL output'
  state.help[21] = '3,4 -- Add/move continuum point'

; Update
  x_fitlls_Reset, state

  ;; Zin
  if keyword_set( ZIN ) then begin
      state.xpos = (1.+zin)*1215.6701
      x_fitlls_newLLS, state
  endif

  ;; Init Fit
  if keyword_set(INIFIT) then x_fitlls_inifit, state, inifit
  ;; Set pmnx
  state.xpmnx = x_getxpmnx(state)
;  if keyword_set(INIFIT) then begin
;    for jj=0,19 do begin
;      state.all_xpmnx[0,jj] = state.xpmnx[0]
;;      state.all_xpmnx[1,jj] = state.xpmnx[1]
;    endfor
;  endif
  ;; Plot
  if keyword_set(VPIN) then x_fitlls_vpin, state, vpin

  x_fitlls_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'x_fitlls', base, /no_block $
  else xmanager, 'x_fitlls', base

return
end
