;+ 
; NAME:
; galx_fit_emiss
;    Version 1.0
;
; PURPOSE:
;   Visually check SDSS CIV absorption detected with sdss_fndciv
;   Edit EW and continuum as desired
;
; CALLING SEQUENCE:
;   
;   galx_fit_emiss, 
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
;   galx_fit_emiss 
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

pro galx_fit_emiss_icmmn, specfil, outfil, IFLG=iflg, clobber=clobber

  common galx_fit_emiss_cmm, $
    npix, $
     fx, $
     wv, $
     sig, $
     conti, $
     fit_px, $
     fit_fx, $
     all_fit, $
     sng_str, $
     flg_lines, $
     lines, $
     fitstr

  flg_lines = 0

  sng_str = { emissstrct }

  if not keyword_set(IFLG) then auto = 1 else auto = 1
  fx = x_readspec(specfil, INFLG = iflg, NPIX = npix $
                  , WAV = wv, FIL_SIG = ysin, SIG = sig  $
                  , auto=auto)
  conti = replicate(1., npix)
  all_fit = replicate(0., npix)
  
  if x_chkfil(outfil+'*', /silent) and not keyword_set(CLOBBER) then begin
      print, 'galx_fit_emissl:  Reading (and will overwrite) ', outfil
      fitstr = xmrdfits(outfil,1)
  endif else begin
      fitstr = [sng_str]
  endelse

  return
end
  

;;;;
; Events
;;;;

pro galx_fit_emiss_event, ev

  common galx_fit_emiss_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'FGROUP' : begin
         widget_control, state.fgroup_id, get_value=tmp
         state.flg_fit = tmp
      end
      'BGROUP' : begin  ; Continuum
         widget_control, state.bgroup_id, get_value=v_lls
         fitstr[state.curline].flg_c = v_lls
         galx_fit_emiss_getconti, state
      end
      'ZGAL' : begin
          widget_control, state.zgal_id, get_value=tmp
          state.zgal = tmp
          galx_fit_emiss_update, state, flg=1
      end
      'WVALU' : begin
         widget_control, state.wvalu_id, get_value=tmp
         mt = where(abs(fitstr.wrest-tmp) LT 1, nmt)
         case nmt of 
            0: begin
              print, 'galx_fit_emiss: No line with that wavelength'
              print, 'galx_fit_emiss: Try again'
              widget_control, state.wvalu_id, set_value=fitstr[state.curline].wrest
           end
            1: begin ; Match
               state.curline = mt[0]
               state.xymnx[0] = fitstr[state.curline].f_reg[0] - 20
               state.xymnx[2] = fitstr[state.curline].f_reg[1] + 20
               if fitstr[state.curline].flg_flux EQ 1 then state.flg_fit = 1
               widget_control, state.fgroup_id, set_value=1
            end
            else: stop
         endcase
      end
      'SPLT': x_specplot, fx, sig, wav=wv, inf=4, /bloc
      'BACK': galx_fit_emiss_back, state
      'DONE' : begin
          galx_fit_emiss_save, state
          print, 'galx_fit_emiss: Exiting..'
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          galx_fit_emiss_save, state
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
                          fitstr[state.curline].f_reg[0] = $
                            (wv[imn] < fitstr[state.curline].f_reg[1])
                      end
                      4 : begin
                          fitstr[state.curline].f_reg[1] = $
                            (wv[imn] > fitstr[state.curline].f_reg[0])
                      end
                      else: 
                  endcase
                  galx_fit_emiss_update, state, flg=1
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
              'T': state.xymnx[3] = 2 * state.xymnx[3]
              ;; Continuum
              '1': galx_fit_emiss_continuum, state, 1L ;; Modify left point 1
              '2': galx_fit_emiss_continuum, state, 2L ;; Modify left point 2
              '3': galx_fit_emiss_continuum, state, 3L ;; Modify right point 1
              '4': galx_fit_emiss_continuum, state, 4L ;; Modify right point 2
              ;; Help
              'H': x_helpwidg, state.help
              ;; Fitting
              'F': begin
                 state.flg_fit = 1
                 widget_control, state.fgroup_id, set_value=1
              end
              'f': begin 
                 state.flg_fit = 0
                 widget_control, state.fgroup_id, set_value=0
              end
              ;; Upper limit
              'U': begin 
                 fitstr[state.curline].flg_flux = -1*fitstr[state.curline].flg_flux  
                 print, 'galx_fit_emiss: Changing flag on this line', $
                        fitstr[state.curline].flg_flux
                 ;; Turn off fitting too
                 state.flg_fit = 0
                 widget_control, state.fgroup_id, set_value=0
              end
              ;; Boxcar
              'B': begin 
                 print, 'galx_fit_emiss: Using boxcar EW for the Flux and Error'
                 print, 'Gauss = ', fitstr[state.curline].flux, 'and box = ', state.boxEW
                 fitstr[state.curline].flg_flux = 2 ; Default is Gaussian
                 fitstr[state.curline].flux = state.boxEW
                 fitstr[state.curline].fsig = state.sigboxEW 
                 ;; Turn off fitting too
                 state.flg_fit = 0
                 widget_control, state.fgroup_id, set_value=0
              end
              ; ZOOMING
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'w': state.xymnx = state.svxymnx ; Reset the screen
              ; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              ; Crude Analysis
              'N': galx_fit_emiss_newline, state
              'D': galx_fit_emiss_delreg, state
              'M': galx_fit_emiss_modify, state ; Switch to nearby line
              else:
           end
          galx_fit_emiss_update, state, flg=1
       end
      else :
   endcase
      
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro galx_fit_emiss_plot, state
  
  common galx_fit_emiss_cmm

  ;; Set plot window
  clr = getcolor(/load)
  

  ;; MgII Plot
  !p.multi = [0,1,1]
  widget_control, state.si2_id, get_value=wind
  wset, wind
  ymnx = [state.xymnx[1], state.xymnx[3]]
  plot, wv, fx, color=clr.black, background=clr.white, psym=10, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=ymnx, xstyle=1, ystyle=1

  ;; Data
;  oplot, wv, (1.-all_fit)*conti, color=clr.red, thick=1
  oplot, wv, conti, color=clr.gray, linesty=1, thick=3
  oplot, wv, sig, color=clr.orange, thick=1, psym=10

  ;; Fitted Lines
  gdl = where(fitstr.wcen LT state.xymnx[2] and $
              fitstr.wcen GT state.xymnx[0], nlin)
  for ii=0L,nlin-1 do begin
     idx = gdl[ii]
     if idx EQ state.curline then begin
        ;; Continuum regions
        for jj=0,3 do  $
           oplot, replicate(fitstr[idx].c_reg[jj],2), ymnx, color=clr.red, lines=3
        
        ;; Fit regions
        for jj=0,1 do $
           oplot, replicate(fitstr[idx].f_reg[jj],2), ymnx, color=clr.darkgreen, $
                  lines=2
     endif
     
     ;; Generate the fit
     wvp = where(wv GT fitstr[idx].c_reg[0,1] and $ 
                 wv LT fitstr[idx].c_reg[1,0] )
     wvv = fitstr[idx].c_reg[0,1] + $
           (fitstr[idx].c_reg[1,0]-fitstr[idx].c_reg[0,1])*findgen(200)/199
     fit_fx = fitstr[idx].A_gauss* $
              exp(-0.5*((wvv-fitstr[idx].wcen)/fitstr[idx].vsigma)^2)
     conti_i = interpol(conti[wvp], wv[wvp], wvv)
     oplot, wvv, fit_fx+conti_i, color=clr.blue

  endfor

  ;; Continuum
  oplot, wv, conti, color=clr.cyan, thick=2

  ;; 
  galx_fit_emiss_pltlines, state

end

pro galx_fit_emiss_pltlines, state, IMG=img

  common galx_fit_emiss_cmm

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.blue
  endif else pclr = 3

  ; Find Min and Max in region
  wvmin = state.xymnx[0]/(state.zgal+1.)
  wvmax = state.xymnx[2]/(state.zgal+1.)

  ; Parse list
  allwv = where(lines.wave GE wvmin AND lines.wave LE wvmax, cntwv)

  ; Plot
  ymax = state.xymnx[1] + 0.02*(state.xymnx[3]-state.xymnx[1])
  ymax2 = state.xymnx[1] + 0.8*(state.xymnx[3]-state.xymnx[1])

  for q=0L,cntwv-1 do begin
      xplt = lines[allwv[q]].wave*(1.+state.zgal)
      ;; Name
      xyouts, xplt, ymax2, $ ; GAL
              strtrim(lines[allwv[q]].name,2), color=pclr, $
              orientation=90., charsize=1.5 
      ;; Dotted line
      oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
        color=pclr, linestyle=1
  endfor

  return
end


;;;;;;;;;;;;;;;;;;;;
;  Continuum
;;;;;;;;;;;;;;;;;;;;

pro galx_fit_emiss_continuum, state, flg

  common galx_fit_emiss_cmm
  if not keyword_set( flg ) then return

; CASE
;  mn = min(abs(wv-state.xpos),imn)

  case flg of 
      1: begin 
         fitstr[state.curline].c_reg[0,0] = (state.xpos < $
                                 (fitstr[state.curline].c_reg[0,1]-state.wv_off)) > $
                                            wv[0]
      end
      2: begin 
         fitstr[state.curline].c_reg[0,1] = (double(state.xpos) > $
                                 (fitstr[state.curline].c_reg[0,0]+state.wv_off)) < $
                                 (fitstr[state.curline].f_reg[0]-state.wv_off)
      end
      3: begin 
         fitstr[state.curline].c_reg[1,0] = (state.xpos < $
                                 (fitstr[state.curline].c_reg[1,1]-state.wv_off)) > $
                                 (fitstr[state.curline].f_reg[1]+state.wv_off)
      end
      4: begin 
         fitstr[state.curline].c_reg[1,1] = (state.xpos > $
                                 (fitstr[state.curline].c_reg[1,0]+state.wv_off)) < $
                                            max(wv)
      end
      else: stop
   endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro galx_fit_emiss_update, state, FLG=flg
  common galx_fit_emiss_cmm
 

  if flg EQ 1 and fitstr[state.curline].wrest GT 0. then begin
      galx_fit_emiss_getconti, state
      if state.flg_fit then galx_fit_emiss_fit, state
;      galx_fit_emiss_boxcar, state
      galx_fit_emiss_updinfo, state
  endif
  galx_fit_emiss_plot, state

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro galx_fit_emiss_setup, state, flg

  common galx_fit_emiss_cmm
    
  ;; xymnx
  state.xymnx[0] = min(wv, max=mx)
  state.xymnx[2] = mx
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  ;; Recover fit limits
  srt = sort(fx[gd])
  ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
  state.xymnx[1] = 0.
  state.xymnx[3] = ymd*1.5
  state.svxymnx = state.xymnx

  if fitstr[state.curline].wrest GT 0 then begin
     if fitstr[state.curline].flg_flux EQ 1 then state.flg_fit = 1 $
     else state.flg_fit = 0
     widget_control, state.fgroup_id, set_value=state.flg_fit
     widget_control, state.bgroup_id, set_value=fitstr[state.curline].flg_c
  endif

  flg = 1
  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro galx_fit_emiss_next, state, BOX=box, GAUSS=gauss, HAND=hand, NG=ng

  common galx_fit_emiss_cmm

  widget_control, /hourglass

  ;; Boxcar
  if keyword_set(BOX) then begin
      fitstr[state.curline].flg_ew = 1
  endif

  ;; Gaussian
  if keyword_set(GAUSS) then begin
      fitstr[state.curline].flg_ew = 0
  endif

  ;; Proceed
  galx_fit_emiss_update, state, FLG=flg

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Back 
pro galx_fit_emiss_back, state
  common galx_fit_emiss_cmm

  ;; Match?
  state.curline = (state.curline-1) > 0
  galx_fit_emiss_update, state, FLG=1

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Back 
pro galx_fit_emiss_getconti, state
  common galx_fit_emiss_cmm

  pix1 = where(wv GE fitstr[state.curline].c_reg[0,0] AND $
               wv LE fitstr[state.curline].c_reg[0,1] , np1)
  pix2 = where(wv GE fitstr[state.curline].c_reg[1,0] AND $
               wv LE fitstr[state.curline].c_reg[1,1] , np2)
  allp = where(wv GE fitstr[state.curline].c_reg[0,0] AND $
               wv LE fitstr[state.curline].c_reg[1,1] , nallp)
  if np1*np2*nallp EQ 0 then stop

  case fitstr[state.curline].flg_c of
     0: begin  ;; Take median of all values
        cval = median( [fx[pix1], fx[pix2]] )
        conti[allp] = cval
     end
     1: begin  ;; Linear fit
        fitp = [pix1,pix2]
        line = linfit( wv[fitp], fx[fitp])
        conti[allp] = line[0] + line[1]*wv[allp]
     end
     else: stop
  endcase
      
  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro galx_fit_emiss_updinfo, state
  common galx_fit_emiss_cmm

  ;; Name
;  widget_control, state.name_id, $
;    set_value=strtrim(dlastr[state.curdla].qso,2)

  ;; Fit
;  widget_control, state.ivalu_id, set_value=fitstr[state.curline].fitnum
  
  ;; EW
  widget_control, state.flux_id, $
                  set_value=string(fitstr[state.curline].flux,format='(f9.3)')+ $
                  ' '+string(fitstr[state.curline].fsig,format='(f7.3)')

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New Region
pro galx_fit_emiss_newline, state
  common galx_fit_emiss_cmm

  sv_cur = state.curline
  state.curline = n_elements(fitstr) 
  fitstr = [fitstr, sng_str]

  ;;  Identify the line
  mn = min(abs(wv-xgetx_plt(state, /strct)),imn)
  wclick = wv[imn]
  mn = min(abs(lines.wave*(1+state.zgal) - wclick), mnlin) 
  fitstr[state.curline].wcen = wclick
  if mn LT 2. then begin
     mt = where(abs(fitstr.wrest-lines[mnlin].wave) LT 5., nmt)
     if nmt GT 0 then begin
        print, 'You have already fit this line!'
        print, 'Delete it first and then refit as a new line if you wish.'
        print, 'Or just modify it by typing an "M"'
        fitstr = fitstr[0:state.curline-1]
        state.curline = sv_cur
        return
     endif
     fitstr[state.curline].wrest = lines[mnlin].wave 
     fitstr[state.curline].ion = lines[mnlin].name 
  endif else begin
     val = x_guinum(1,TITLE='Enter rest wavelength:')
     fitstr[state.curline].wrest = val
  endelse
  widget_control, state.wvalu_id, set_value=fitstr[state.curline].wrest

  ;; Fit region
  fitstr[state.curline].f_reg[0] = wv[(imn-5) > 0]
  fitstr[state.curline].f_reg[1] = wv[(imn+5) < (npix-1)]

  ;; Continuum
  fitstr[state.curline].c_reg[0,0] = wv[(imn-25) > 0]
  fitstr[state.curline].c_reg[0,1] = wv[(imn-15) > 0]
  fitstr[state.curline].c_reg[1,0] = wv[(imn+15) < (npix-1)]
  fitstr[state.curline].c_reg[1,1] = wv[(imn+25) < (npix-1)]

  ;; Plotting
  state.xymnx[0] = fitstr[state.curline].c_reg[0,0] - 30
  state.xymnx[2] = fitstr[state.curline].c_reg[1,1] + 30

  ;; Continuum flag (Gaussian = default)
  widget_control, state.bgroup_id, set_value=0
  ;; Fit
  state.flg_fit = 0
  widget_control, state.fgroup_id, set_value=1

  ;; Do the boxcar
;  galx_fit_emiss_boxcar, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New Region
pro galx_fit_emiss_modify, state
  common galx_fit_emiss_cmm

  mt = where(abs(fitstr.wrest-state.xpos/(1+state.zgal)) LT 5, nmt)
  case nmt of 
     0: begin
        print, 'galx_fit_emiss: Get closer to the line'
        print, 'galx_fit_emiss: And try it again'
     end
     1: begin
        state.curline = mt[0]
        state.xymnx[0] = fitstr[state.curline].f_reg[0] - 20
        state.xymnx[2] = fitstr[state.curline].f_reg[1] + 20
        if fitstr[state.curline].flg_flux EQ 1 then state.flg_fit = 1 $
        else state.flg_fit = 0
        ;; Reset continuum
        galx_fit_emiss_getconti, state
        ;; GUI
        widget_control, state.fgroup_id, set_value=state.flg_fit
        widget_control, state.bgroup_id, set_value=fitstr[state.curline].flg_c
     end
                    ;;
     else:  print, 'galx_fit_emiss: Not a valid key!' ; Nothing
  endcase
  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New Region
pro galx_fit_emiss_delreg, state
  common galx_fit_emiss_cmm

  ;; First one?
  if fitstr[state.curline].wrest LE 0 then return

  ;; Throw out
  msk = replicate(1B,n_elements(fitstr))
  msk[state.curline] = 0B
  keep = where(msk, nlin)
  if nlin GT 0 then fitstr = fitstr[keep] else begin
   state.curline = 0L
   state.flg_fit = 0L
   fitstr[0].wrest = 0.
   fitstr[0].wcen = 0.
   return
  endelse
  state.curline = n_elements(fitstr)-1

  state.flg_fit = 0 
;  if fitstr[state.curline].wrest LE 0 then state.flg_fit = 0 else begin
;     if fitstr[state.curline].flg_flux EQ 1 then state.flg_fit = 1 $
;     else state.flg_fit = 0
;  endelse
  widget_control, state.fgroup_id, set_value=state.flg_fit

  ;; Update
  if fitstr[state.curline].wrest GT 0 then begin
     state.xymnx[0] = fitstr[state.curline].f_reg[0] - 20
     state.xymnx[2] = fitstr[state.curline].f_reg[1] + 20
  endif
  widget_control, state.wvalu_id, set_value=fitstr[state.curline].wrest

  ;; Fit flag (Gaussian = default)
  widget_control, state.bgroup_id, set_value=fitstr[state.curline].flg_c

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro galx_fit_emiss_save, state

  common galx_fit_emiss_cmm

  ;; Sort
  srt = sort(fitstr.wrest)
  fitstr = fitstr[srt]
  ;; Parse
  gd = where(fitstr.wrest GT 0,ngd)
  print, 'galx_fit_emiss: Saving..', state.outfil
  if ngd NE 0 then mwrfits, fitstr[gd], state.outfil, /create
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro galx_fit_emiss_fit, state

  common galx_fit_emiss_cmm

  if fitstr[state.curline].wrest LE 0 then return

  px = where(wv GE fitstr[state.curline].f_reg[0] AND $
             wv LE fitstr[state.curline].f_reg[1] , np)
  wcen = round(mean(fitstr[state.curline].f_reg))
  mn = min(abs(wv-wcen),cpix)
  subpx = px
  fit_px = px

  ;; Profile
  gprof = (fx[px] - conti[px]) 

  dwv = abs(wv[cpix]-wv[cpix+1])
  gsssig = 2.

  ;; FIT
  fit_wv = wv[px]
  fit_fx = x_gaussfit(fit_wv, gprof, acoeff, $
                    estimates=[max(gprof), wcen, gsssig], $
                    sigma=sigma, nterms=3, COVAR=covar, $
                    measure_errors=sig[px])
  all_fit[px] = fit_fx

  ;; Boxcar 
  state.boxEW = total(gprof*dwv)  ; Ang

  ;; Flux
  flux = acoeff[0] * acoeff[2] * sqrt(!pi*2.) ; A
  sig_f1 = sqrt(total( ((sig[subpx]/conti[subpx])*dwv)^2)) ; A
 ; sig_f1 = sqrt(mean( ((sig[subpx]/conti[subpx])*dwv)^2)) ; A
  sig_f1 = sqrt(total( ((sig[subpx]*dwv)^2))) ; A  [Boxcar]
  state.sigboxEW = sig_f1
  sig_f2 = sqrt(!pi*2.) * sqrt( (acoeff[0]*sigma[2])^2 + $
                                (acoeff[2]*sigma[0])^2 )
;      sigew = sigew1 > sigew2
;  sigew = sigew1 < sigew2
  sig_f = sig_f1 < sig_f2


  ;; EW

  ;; Error checking (for a bad fit)
;  ew_chk = total((1. - fx[px]/conti[px])*abs(dwv[px])) ;; Obs Ang
;  if abs(ew_chk-ewval) GT 5.*sigew then begin
;      print, 'Bad EW value!', ew_chk, ewval
;      print, 'Probably a bad Gaussian fit'
;      print, 'Setting EW to the lower value', ew_chk, ewval
;      ewval = (ewval < ew_chk) > 0.
;  endif

;  state.gaussEW = ewval 
;  state.siggaussEW = sigew 

  ;; Save params
  fitstr[state.curline].wcen = acoeff[1]
  fitstr[state.curline].vsigma = acoeff[2]
  fitstr[state.curline].A_gauss = acoeff[0]

  fitstr[state.curline].sig_wc = sigma[1]
  fitstr[state.curline].sig_vsig = sigma[2]
  fitstr[state.curline].sig_A = sigma[0]

  fitstr[state.curline].flg_flux = 1  ; Default is Gaussian
  fitstr[state.curline].flux = double(flux)  
  fitstr[state.curline].fsig = double(sig_f)

  fitstr[state.curline].zem = acoeff[1]/fitstr[state.curline].wrest - 1
  fitstr[state.curline].zsig = sigma[1]/fitstr[state.curline].wrest 
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
pro galx_fit_emiss, specfil, outfil, ZGAL=zgal, $
               XSIZE=xsize, YSIZE=ysize, IFLG=iflg, CLOBBER=clobber, $
                    UPDATE=update

  common galx_fit_emiss_cmm
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'galx_fit_emiss, specfil, outfil, IFLG=, /CLOBBER, /UPDATE [v1.0]'
    return
  endif 


  device, get_screen_size=ssz
  if not keyword_set( ZGAL ) then    zgal = 0.
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
;  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = 0L

; Initialize the common blcok
  galx_fit_emiss_icmmn, specfil, outfil, IFLG=iflg, CLOBBER=clobber

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }

; STATE

  state = {             $
          curline: 0L, $
          outfil: outfil, $
          model: fltarr(npix), $
          zgal: double(zgal), $
          gaussEW: 0., $
          siggaussEW: 0., $
          boxEW: 0.d, $
          sigboxEW: 0.d, $
          flg_redo: keyword_set(redo), $
          flg_fit: 0, $
          wv_off: 2., $ ; Ang
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
          zgal_id: 0L, $
          lines_id: 0L, $
          si2_id: 0L, $
          tdraw_id: 0L, $
          left_id: 0L, $
          right_id: 0L, $
          rhs_id: 0L, $
          info_id: 0L, $
          quality_id: 0L, $
          flux_id: 0L, $
          wvalu_id: 0L, $
          boxew_id: 0L, $
          hew_id: 0L, $
          bgroup_id: 0L, $
          fgroup_id: 0L, $
          hsew_id: 0L, $
          help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES

;    WIDGET
  base = WIDGET_BASE( title = 'galx_fit_emiss: Check spectra', /column, $
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
  state.zgal_id = cw_field(civinf, title='z', value=state.zgal, /floating, $
                           /column, xsize=8, /return_events, uvalue='ZGAL')
  state.flux_id = cw_field(civinf, title='Flux: ', value=' ', xsize=17)
;  state.boxew_id = cw_field(civinf, title='Box EW: ', value=' ', xsize=13)
  handinf = widget_base(state.info_id, /column, /align_center, frame=2)

  ;; BUTTONS

  ;; Indexing
  index = widget_base(state.info_id, /column, /align_center, frame=2)
  state.wvalu_id = cw_field(index, title='wrest: ', value=state.curline, /floating, $
                            xsize=7, /return_events, uvalue='WVALU')
;  back  = WIDGET_BUTTON(index, value='BACK', uvalue='BACK')

  state.bgroup_id = cw_bgroup(state.info_id, ['C', 'L'], $
                           /exclusive, row=2, label_top='Conti', UVALUE='BGROUP')
  widget_control, state.bgroup_id, set_value=0
  state.fgroup_id = cw_bgroup(state.info_id, ['N', 'Y'], $
                           /exclusive, row=2, label_top='Fit?', UVALUE='FGROUP')
  widget_control, state.fgroup_id, set_value=0

  ;; Saving
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  splt = WIDGET_BUTTON(butbase2, value='SPLT',uvalue='SPLT')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = ['   Help Menu   ',$
             'l -- Set Left ', $
             'r -- Set Right ', $
             'b -- Set Bottom ', $
             't -- Set Top ', $
             'Z -- Set ymin to 0.', $
             'T -- Set ymax to 1.1', $
             'U -- Increase ymax by 2x', $
             'w -- Reset the screen', $
             'L -- Set redshift with a line', $
             'i -- Zoom in', $
             'o -- Zoom out', $
             'z -- Zoom region', $
             '[ -- Pan left', $
             '] -- Pan right', $
             '{ -- Pan left, no y-rescale', $
             '} -- Pan right, no y-rescale', $
             'H -- Show this screen', $
             'N -- New line to analyze', $
             'M -- Change line to analyze', $
             'D -- Delete current line', $
             'F/f -- Start/stop fitting', $
             '1-4 -- Modify continuum region', $
             'LMB/RMB -- Modify fit region' $
            ]
  nhelp = n_elements(strhelp)
  for ii=0l,nhelp-1 do state.help[ii] = strhelp[ii]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  if not keyword_set(UPDATE) then WIDGET_CONTROL, base, /realize

  ;; Galaxy line list
  llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/gal.lst'
  lines = x_setllst(llist, 1)
  flg_lines = 1

  ;; PLOT
  galx_fit_emiss_setup, state
  if not keyword_set(UPDATE) then begin
     galx_fit_emiss_updinfo, state
     galx_fit_emiss_plot, state
  endif else begin ;; Just evaluate all the lines
     gd = where(fitstr.wrest GT 0,ngd)
     if ngd EQ 0 then begin
        print, 'galx_fit_emiss: No lines to update!  Returning...'
        return
     endif
     print, 'galx_fit_emiss: Updating ', ngd, ' lines'
     for ii=0L,ngd-1 do begin
        state.curline = gd[ii]
        sv_flg = fitstr[state.curline].flg_flux
        galx_fit_emiss_getconti, state
        galx_fit_emiss_fit, state
        fitstr[state.curline].flg_flux = sv_flg
        if sv_flg EQ 2 then begin
           fitstr[state.curline].flux = state.boxEW
           fitstr[state.curline].fsig = state.sigboxEW 
        endif
     endfor
     print, 'galx_fit_emiss: Writing ', state.outfil
     mwrfits, fitstr[gd], state.outfil, /create
     return
  endelse
     
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'galx_fit_emiss', base
  delvarx, fx, wv, npix, sig, fitstr

  return
end
	
