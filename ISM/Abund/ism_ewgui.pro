;+ 
; NAME:
; ism_ewgui
;    Version 1.0
;
; PURPOSE:
;   Visually check SDSS CIV absorption detected with sdss_fndciv
;   Edit EW and continuum as desired
;
; CALLING SEQUENCE:
;   
;   ism_ewgui, 
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
;   ism_ewgui 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   1-Dec-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro ism_ewgui_icmmn, specfil, outfil, IFLG=iflg, clobber=clobber

  common ism_ewgui_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    fit_px, $
    fit_fx, $
    all_fit, $
    sng_str, $
    fitstr


  sng_str = { $
            parms: dblarr(5,10), $
            sig_parms: dblarr(5,10), $
            ew: fltarr(10), $
            sigew: fltarr(10), $
            boxew: 0., $
            boxsigew: 0., $
            flg_ew: 0, $
            region: dblarr(2), $
            fitnum: 0L $
            }
  if not keyword_set(IFLG) then auto = 1 else auto = 1
  fx = x_readspec(specfil, INFLG = iflg, NPIX = npix $
                  , WAV = wv, FIL_SIG = ysin, SIG = sig  $
                  , auto=auto)
  conti = replicate(1., npix)
  all_fit = replicate(0., npix)
  
  if x_chkfil(outfil+'*', /silent) and not keyword_set(CLOBBER) then begin
      print, 'ism_ewguil:  Reading (and will overwrite) ', outfil
      fitstr = xmrdfits(outfil,1)
  endif else begin
      fitstr = [sng_str]
  endelse

  return
end
  

;;;;
; Events
;;;;

pro ism_ewgui_event, ev

  common ism_ewgui_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'IVALU' : begin
          widget_control, state.ivalu_id, get_value=tmp
          mt = where(fitstr.fitnum EQ tmp, nmt)
          if nmt EQ 1 then begin
              state.curfit = mt[0]
              state.xymnx[0] = fitstr[state.curfit].region[0] - 20
              state.xymnx[2] = fitstr[state.curfit].region[1] + 20
              ism_ewgui_update, state, FLG=1
          endif else begin
              print, 'ism_ewgui: Bad fitnum'
              print, 'ism_ewgui: Resetting to previous value'
              widget_control, state.ivalu_id, get_value=fitstr[state.curfit].fitnm
          endelse
      end
      'SPLT': x_specplot, fx, sig, wav=wv, inf=4, /bloc
      'BACK': ism_ewgui_back, state
      'DONE' : begin
          ism_ewgui_save, state
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          print, 'ism_ewgui: Saving..'
          ism_ewgui_save, state
          print, 'ism_ewgui: Done!'
      end
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'SI2DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  mn = min(abs(wv- xgetx_plt(state, /strct)),imn)
                  case ev.press of
                      1 : begin
                          fitstr[state.curfit].region[0] = $
                            (wv[imn] < fitstr[state.curfit].region[1])
                      end
                      4 : begin
                          fitstr[state.curfit].region[1] = $
                            (wv[imn] > fitstr[state.curfit].region[0])
                      end
                      else: 
                  endcase
                  ism_ewgui_update, state, flg=1
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
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'T': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              ;; Continuum
              '3': ism_ewgui_continuum, state, 1L ;; Add point
              '4': ism_ewgui_continuum, state, 2L ;; Move a point
              ; ZOOMING
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'w': state.xymnx = state.svxymnx ; Reset the screen
              ; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              ;; Continuum
              'C': ism_ewgui_conti, state
              ; Crude Analysis
              'g': ism_ewgui_next, state, /GAUSS
              'B': ism_ewgui_next, state, /BOX
              'N': ism_ewgui_newreg, state
              'D': ism_ewgui_delreg, state
              ;;
              else:  print, 'ism_ewgui: Not a valid key!' ; Nothing
          endcase
          ism_ewgui_update, state, flg=1
      end


      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro ism_ewgui_plot, state
  
  common ism_ewgui_cmm

  ;; Set plot window
  clr = getcolor(/load)
  

  ;; MgII Plot
  !p.multi = [0,1,1]
  widget_control, state.si2_id, get_value=wind
  wset, wind
  plot, wv, fx, color=clr.black, background=clr.white, psym=10, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1

  ;; Continuum and error
  oplot, wv, (1.-all_fit)*conti, color=clr.red, thick=1
  oplot, wv, conti, color=clr.gray, linesty=1, thick=3
  oplot, wv, sig, color=clr.orange, thick=1, psym=10


  if fitstr[state.curfit].fitnum GT 0 then begin
      ;; Region
      oplot, replicate(fitstr[state.curfit].region[0],2), $
             [-9e9,9e9], color=clr.green, linestyle=2, thick=3
      oplot, replicate(fitstr[state.curfit].region[1],2), $
             [-9e9,9e9], color=clr.green, linestyle=2, thick=3
      ;; Fit
      if not keyword_set(fit_px) or not keyword_set(fit_fx) then begin 
          fit_px = where(wv GE fitstr[state.curfit].region[0] AND $
                         wv LE fitstr[state.curfit].region[1] , np)
          acoeff = fitstr[state.curfit].parms[0:2,0]
          fit_fx = acoeff[0]*exp(-0.5*((wv[fit_px]-acoeff[1])/acoeff[2])^2)
          all_fit[fit_px] = fit_fx
      endif
      oplot, wv[fit_px], (1.-fit_fx)*conti[fit_px], color=clr.blue
  endif

  ;; Continuum
  ;; Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.cyan, symsize=5, thick=3
  endif

end

;;;;;;;;;;;;;;;;;;;;
;  Continuum
;;;;;;;;;;;;;;;;;;;;

pro ism_ewgui_continuum, state, flg

  common ism_ewgui_cmm
  if not keyword_set( flg ) then return

; CASE

  case flg of 
      0: stop ;  RESET to 1/2
      1: begin ; Set point
          if state.cstr.npts EQ 9 then return
          i = state.cstr.npts 
          state.cstr.xval[i] = state.xpos
          state.cstr.yval[i] = state.ypos
          state.cstr.msk[i] = 1L
          state.cstr.npts =  state.cstr.npts  + 1
      end
      2: begin ; Move point
          if state.cstr.npts EQ 0 then return
          gd = where(state.cstr.msk EQ 1)
          mn = min(abs(state.cstr.xval[gd]-state.xpos), imn)
          state.cstr.yval[gd[imn]] = state.ypos
          state.cstr.xval[gd[imn]] = state.xpos
      end
      3: ; Do nothing just respline
      else: stop
  endcase

  ;; Save
  mgiistr[state.cursi2].EW[40+lindgen(state.cstr.npts)] = $
    state.cstr.xval[lindgen(state.cstr.npts)]
  mgiistr[state.cursi2].sigEW[40+lindgen(state.cstr.npts)] = $
    state.cstr.yval[lindgen(state.cstr.npts)]

  ;; Redo the continuum
  if state.cstr.npts NE 0 then begin
      gdmsk = where(state.cstr.msk EQ 1, nmsk)
      if nmsk LT 3 then conti[*] = mean(state.cstr.yval[gdmsk]) $
      else begin ;; SPLINE
          mn = min(state.cstr.xval[gdmsk], imn)
          ;; Low
          low = where(wv LE mn, nlow)
          if nlow NE 0 then conti[low] = state.cstr.yval[gdmsk[imn]]
          ;; High
          mx = max(state.cstr.xval[gdmsk], imx)
          high = where(wv GE mx, nhigh)
          if nhigh NE 0 then conti[high] = state.cstr.yval[gdmsk[imx]]
          ;; Interp
          gd = where(wv GT mn AND wv LT mx, ngd)
          if ngd NE 0 then begin
              ;; SORT
              srt = sort(state.cstr.xval[gdmsk])
              swv = sort(wv[gd])
              twv = wv[gd[swv]]
              conti[gd[swv]] = xspline(state.cstr.xval[gdmsk[srt]], $
                                             state.cstr.yval[gdmsk[srt]], $
                                             twv)
              ;; Update Fit
;         if state.nlin NE 0 then x_fitline_updfit, state, EXACT=state.exact
          endif
      endelse   
  endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro ism_ewgui_update, state, FLG=flg
  common ism_ewgui_cmm
 
  ;; Update flag
  widget_control, state.bgroup_id, get_value=v_lls
  fitstr[state.curfit].flg_ew = v_lls

  if flg EQ 1 then begin
      ism_ewgui_fit, state
      ism_ewgui_boxcar, state
      ism_ewgui_updinfo, state
      ism_ewgui_plot, state
  endif
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro ism_ewgui_setup, state, flg

  common ism_ewgui_cmm
    
  ;; xymnx
  state.xymnx[0] = min(wv, max=mx)
  state.xymnx[2] = mx
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  ;; Flush/fix continuum points
  state.cstr.npts = 0
  state.cstr.xval = 0.
  state.cstr.yval = 0.
  state.cstr.msk = 0

  ;; Recover fit limits
  srt = sort(fx[gd])
  ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
  state.xymnx[1] = 0.
  state.xymnx[3] = ymd*1.5
  state.svxymnx = state.xymnx

  flg = 1
  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro ism_ewgui_next, state, BOX=box, GAUSS=gauss, HAND=hand, NG=ng

  common ism_ewgui_cmm

  widget_control, /hourglass

  ;; Boxcar
  if keyword_set(BOX) then begin
      fitstr[state.curfit].flg_ew = 1
  endif

  ;; Gaussian
  if keyword_set(GAUSS) then begin
      fitstr[state.curfit].flg_ew = 0
  endif

  ;; Proceed
  ism_ewgui_update, state, FLG=flg

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Back 
pro ism_ewgui_back, state
  common ism_ewgui_cmm

  ;; Match?
  state.curfit = (state.curfit-1) > 0
  ism_ewgui_update, state, FLG=1

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Back 
pro ism_ewgui_conti, state
  common ism_ewgui_cmm

  ;; Match?
  if fitstr[state.curfit].fitnum GT 0 then begin
      fitstr[state.curfit].parms[3] = xgety_plt(state, /strct) 
      px = where(wv GE fitstr[state.curfit].region[0] AND $
                 wv LE fitstr[state.curfit].region[1] , np)
      conti[px] = fitstr[state.curfit].parms[3]
  endif
      
  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro ism_ewgui_updinfo, state
  common ism_ewgui_cmm

  ;; Name
;  widget_control, state.name_id, $
;    set_value=strtrim(dlastr[state.curdla].qso,2)

  ;; Fit
  widget_control, state.ivalu_id, set_value=fitstr[state.curfit].fitnum
  
  ;; EW
  widget_control, state.ew_id, $
                  set_value=string(fitstr[state.curfit].ew[0],format='(f7.3)')+ $
                  ' '+string(fitstr[state.curfit].sigew[0],format='(f5.3)')
  ;; Box EW
  widget_control, state.boxew_id, $
                  set_value=string(fitstr[state.curfit].boxew,format='(f7.3)')+' '+$
                  string(fitstr[state.curfit].boxsigEW,format='(f5.3)')

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Boxcar
pro ism_ewgui_boxcar, state
  common ism_ewgui_cmm

  ;; Simple sum

  sumpix = where(wv GE fitstr[state.curfit].region[0] AND $
             wv LE fitstr[state.curfit].region[1] , np)
  if np NE 0 then begin
      dwv = wv[sumpix[np-1]+1] - wv[sumpix[np-1]]
      dwv = abs(dwv)
      ;; Value
      state.boxEW = total(1.-(fx[sumpix]>0.)/conti[sumpix])*dwv 

      ;; ERROR
      sumvar = total((sig[sumpix]/conti[sumpix])^2)
      state.sigboxEW = sqrt(sumvar)*dwv
  endif
  
  ;; Save
  ism_ewgui_updinfo, state
  fitstr[state.curfit].boxew = state.boxEW
  fitstr[state.curfit].boxsigew = state.sigboxEW

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New Region
pro ism_ewgui_newreg, state
  common ism_ewgui_cmm

  ;; First one?
  if fitstr[state.curfit].fitnum EQ 0 then begin
      fitstr[state.curfit].fitnum = 1 
  endif else begin
      a = where(fitstr.fitnum GT 0, na)
      state.curfit = na
      fitstr = [fitstr, sng_str]
      fitstr[state.curfit].fitnum = max(fitstr.fitnum)+1
  endelse

  mn = min(abs(wv- xgetx_plt(state, /strct)),imn)
  fitstr[state.curfit].region[0] = wv[(imn-5) > 0]
  fitstr[state.curfit].region[1] = wv[(imn+5) < (npix-1)]

  state.xymnx[0] = fitstr[state.curfit].region[0] - 20
  state.xymnx[2] = fitstr[state.curfit].region[1] + 20

  ;; Fit flag (Gaussian = default)
  widget_control, state.bgroup_id, set_value=0

  ;; Do the boxcar
  ism_ewgui_boxcar, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New Region
pro ism_ewgui_delreg, state
  common ism_ewgui_cmm

  ;; First one?
  if fitstr[state.curfit].fitnum EQ 0 then return
  a = where(fitstr.fitnum GT 0, na)
  high = where(fitstr.fitnum GT fitstr[state.curfit].fitnum,nhigh)

  ;; Shift
  fitstr[state.curfit:na-2] = fitstr[state.curfit+1:na-1]

  ;; Delete last one
  fitstr[na-1].fitnum = 0
  
  state.curfit = 0 > state.curfit < (na-2)

  ;; Update
  state.xymnx[0] = fitstr[state.curfit].region[0] - 20
  state.xymnx[2] = fitstr[state.curfit].region[1] + 20

  ;; Fit flag (Gaussian = default)
  widget_control, state.bgroup_id, set_value=fitstr[state.curfit].flg_ew

  ;; Do the boxcar
  ism_ewgui_boxcar, state

  return

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_ewgui_save, state

  common ism_ewgui_cmm

  gd = where(fitstr.fitnum GT 0,ngd)
  if ngd NE 0 then mwrfits, fitstr[gd], state.outfil, /create
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_ewgui_fit, state

  common ism_ewgui_cmm

  if fitstr[state.curfit].fitnum EQ 0 then return

  px = where(wv GE fitstr[state.curfit].region[0] AND $
             wv LE fitstr[state.curfit].region[1] , np)
  wcen = round(mean(fitstr[state.curfit].region))
  mn = min(abs(wv-wcen),cpix)
  subpx = px
  fit_px = px

  ;; Profile
  gprof = -1.*(fx[px] / conti[px] - 1.)

  dwv = abs(wv[cpix]-wv[cpix+1])
  gsssig = 2.

  ;; FIT
  fit_wv = wv[px]
  fit_fx = x_gaussfit(fit_wv, gprof, acoeff, $
                    estimates=[max(gprof), wcen, gsssig], $
                    sigma=sigma, nterms=3, COVAR=covar, $
                    measure_errors=sig[px]/conti[px])
  all_fit[px] = fit_fx

  ;; Save EW
  ewval = acoeff[0] * acoeff[2] * sqrt(!pi*2.) ; A
  sigew1 = sqrt(total( ((sig[subpx]/conti[subpx])*dwv)^2)) ; A
  sigew2 = sqrt(!pi*2.) * sqrt( (acoeff[0]*sigma[2])^2 + $
                                (acoeff[2]*sigma[0])^2 )
;      sigew = sigew1 > sigew2
  sigew = sigew1 < sigew2

  ;; Error checking (for a bad fit)
  ew_chk = total((1. - fx[px]/conti[px])*abs(dwv[px])) ;; Obs Ang
  if abs(ew_chk-ewval) GT 5.*sigew then begin
      print, 'Bad EW value!', ew_chk, ewval
      print, 'Probably a bad Gaussian fit'
      print, 'Setting EW to the lower value', ew_chk, ewval
      ewval = (ewval < ew_chk) > 0.
  endif

  state.gaussEW = ewval 
  state.siggaussEW = sigew 

  ;; Save params
  fitstr[state.curfit].parms[0:2] = acoeff
  fitstr[state.curfit].ew[0] = ewval
  fitstr[state.curfit].sigew[0] = sigew

      
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; MAIN PROGRAM ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_ewgui, specfil, outfil, $
               XSIZE=xsize, YSIZE=ysize, IFLG=iflg, CLOBBER=clobber

  common ism_ewgui_cmm
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'ism_ewgui, specfil, outfil, IFLG=, /CLOBBER [v1.0]'
    return
  endif 


  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
;  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = 0L

; Initialize the common blcok
  ism_ewgui_icmmn, specfil, outfil, IFLG=iflg, CLOBBER=clobber

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }
  tmp3 = { contistrct }

; STATE

  state = {             $
          curfit: 0L, $
          outfil: outfil, $
          model: fltarr(npix), $
          gaussEW: 0., $
          siggaussEW: 0., $
          boxEW: 0., $
          sigboxEW: 0., $
          cstr: tmp3, $
          flg_redo: keyword_set(redo), $
          xpos: 0.0, $
          ypos: 0.0, $
          ipress: 0L, $
          pos: [0.1,0.1,0.95,0.95], $ ; Plotting
          flg_zoom: 0, $
          psfile: 0, $
          help: strarr(50), $
          svxymnx: fltarr(4), $
          xymnx: fltarr(4), $
          tmpxy: fltarr(4), $
          xcurs: 0., $
          ycurs: 0., $
          size: lonarr(2), $
          base_id: 0L, $        ; Widgets
          ldraw_id: 0L, $       ; Lya
          text_id: 0L, $       ; 
          draw_base_id: 0L, $
          fxval_id: 0L, $
          iwvval_id: 0L, $
          mdraw_id: 0L, $       ; Spec Window
          mdrawbase_id: 0L, $
          swvval_id: 0L, $
          xmax_id: 0L, $
          name_id: 0L, $
          nspec_id: 0L, $
          lines_id: 0L, $
          si2_id: 0L, $
          tdraw_id: 0L, $
          left_id: 0L, $
          right_id: 0L, $
          rhs_id: 0L, $
          info_id: 0L, $
          quality_id: 0L, $
          ew_id: 0L, $
          ivalu_id: 0L, $
          boxew_id: 0L, $
          hew_id: 0L, $
          bgroup_id: 0L, $
          hsew_id: 0L, $
          help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES

;    WIDGET
  base = WIDGET_BASE( title = 'ism_ewgui: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.info_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, ysize=ysize/5.)
  state.draw_base_id = widget_base(state.base_id, /column, /base_align_left, $
                                   xsize=xsize, $
                                   ysize=round(4*ysize/5), $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)
  state.si2_id = widget_draw(state.draw_base_id, xsize=xsize, $
                              ysize=round(4*ysize/5), /frame, retain=2, $
                              uvalue='SI2DRAW', /button_even, /motion_events)
  state.text_id = widget_text(state.draw_base_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
  state.size[0] = xsize
  state.size[1] = round(4*ysize/5)

;;;;;; Info window ;;;;;;;;;;;
;               ysize=round(ysize/3.))
  ;; Info
  civinf = widget_base(state.info_id, /column, /align_center, frame=2)
;  state.name_id = cw_field(civinf, title='Obj ', value=' ', xsize=18)
  zinf = widget_base(civinf, /row, /align_center, frame=2)
  state.ew_id = cw_field(civinf, title='EW: ', value=' ', xsize=13)
  state.boxew_id = cw_field(civinf, title='Box EW: ', value=' ', xsize=13)
  handinf = widget_base(state.info_id, /column, /align_center, frame=2)

  ;; BUTTONS

  ;; Indexing
  index = widget_base(state.info_id, /column, /align_center, frame=2)
  state.ivalu_id = cw_field(index, title='Indx: ', value=state.curfit, $
                            xsize=6, /return_events, uvalue='IVALU')
  back  = WIDGET_BUTTON(index, value='BACK', uvalue='BACK')

  state.bgroup_id = cw_bgroup(state.info_id, ['G', 'B'], $
                           /exclusive, row=2, UVALUE='BGROUP')
  widget_control, state.bgroup_id, set_value=0

  ;; Saving
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  splt = WIDGET_BUTTON(butbase2, value='SPLT',uvalue='SPLT')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB - Truncate/Extend trace', $ 
             'RMB - Contrast/Brightness', $
             'CMB/CMB - Zoom' $ 
             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ;; PLOT
  ism_ewgui_setup, state
  ism_ewgui_updinfo, state
  ism_ewgui_plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'ism_ewgui', base
  delvarx, fx, wv, npix, sig, fitstr

  return
end
	
