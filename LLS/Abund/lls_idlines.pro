;+ 
; NAME:
;  lls_idlines
;    GUI based on ism_idlines to vette QPQ metal transitions
; REVISION HISTORY:
;   11-Jun-2013 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro lls_idlines_icmmn, specfil, outfil, IFLG=iflg, DATA=data, CFIL=cfil

  common lls_idlines_cmm, $
    npix, $
     fx, $
     FSCALE, $
     wv, $
     sig, $
     conti, $
     OUT_CFLG, $
     sng_str, $
     flg_lines, $
     lines, $
     all_str

  flg_lines = 0
  out_cflg = 0


  ;; Read specfil
  if not keyword_set(DATA) then begin
     if not keyword_set(IFLG) then auto = 1 else begin
        auto = 1
        if (iflg EQ 2) or (iflg EQ 11) then  auto = 0 else auto = 1
     endelse
     fx = x_readspec(specfil, INFLG = iflg, NPIX = npix $
                     , WAV = wv, FIL_SIG = ysin, SIG = sig  $
                     , auto=auto)
     if keyword_set(CFIL) then begin
        stop
     endif else conti = replicate(1., npix)
  endif else begin ;; QPQ format
     fx = data.flux_bg
     wv = data.wave
     sig = data.sig_bg
     conti = data.conti
     npix = n_elements(fx)
  endelse

  if median(fx) LT 1e-5 then begin
     FSCALE = median(fx)
     fx = fx / FSCALE
     sig = sig / FSCALE
     ;; Deal with continuum if there is one
  endif else FSCALE = 1.

  return
end
  

;;;;
; Events
;;;;

pro lls_idlines_event, ev

  common lls_idlines_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
     'LIM_FLAG': begin
        if ev.select EQ 1 then begin
           widget_control, state.limit_id, get_value=limit_val
           all_str[state.curline].flg_limit = limit_val+1
           lls_idlines_update, state, flg=1
        endif
     end
     'ANLY_FLAG': begin
        if ev.select EQ 1 then begin
           widget_control, state.anly_id, get_value=anly_val
           all_str[state.curline].flg_anly = anly_val
           lls_idlines_update, state, flg=1
        endif
     end
     'EYE_FLAG': begin
        widget_control, state.eye_id, get_value=eye_val
        all_str[state.curline].flg_eye = 1*eye_val[0] + 2*eye_val[1] + 4*eye_val[2]
        lls_idlines_update, state, flg=1
     end
     'LLIST': begin
        state.lastline = state.curline
        state.curline = where(abs(all_str.wrest-lines[ev.index].wave) LT 1e-3,na)
        if na NE 1 then stop
        lls_idlines_setwin, state
        lls_idlines_update, state, flg=1
     end
      'SPLT' : x_specplot, fx, sig, wave=wv, zin=state.zgal, /lls, /block, infl=4
      'VPLT' : x_velplt, fx, state.zgal, wave=wv, infl=4, /lls, /un_norm
      'ZGAL' : begin
          widget_control, state.zgal_id, get_value=tmp
          state.zgal = tmp
          lls_idlines_update, state, flg=1
      end
      'WVALU' : begin
         widget_control, state.wvalu_id, get_value=tmp
         mt = where(abs(fitstr.wrest-tmp) LT 1, nmt)
         case nmt of 
            0: begin
              print, 'lls_idlines: No line with that wavelength'
              print, 'lls_idlines: Try again'
              widget_control, state.wvalu_id, set_value=fitstr[state.curline].wrest
           end
            1: begin ; Match
               state.lastline = state.curline
               state.curline = mt[0]
               state.xymnx[0] = fitstr[state.curline].f_reg[0] - 20
               state.xymnx[2] = fitstr[state.curline].f_reg[1] + 20
               if fitstr[state.curline].flg_flux EQ 1 then state.flg_fit = 1
               widget_control, state.fgroup_id, set_value=1
            end
            else: stop
         endcase
      end
      'DONE' : begin
          lls_idlines_save, state
          print, 'lls_idlines: Exiting..'
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          lls_idlines_save, state
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
                  xwv=state.xpos/(1+state.zgal)
                  case ev.press of
                      1 : begin ; LMB
                         v0 = 3e5 * (xwv-all_str[state.curline].wrest) / $
                              all_str[state.curline].wrest
                         all_str[state.curline].dv[0] = v0 < $
                                                        (all_str[state.curline].dv[1] - 1.)
                         widget_control, state.vmin_id, set_value=v0
                      end
                      4 : begin ; RMB
                         v1 = 3e5 * (xwv-all_str[state.curline].wrest) / $
                              all_str[state.curline].wrest
                         all_str[state.curline].dv[1] = v1 > $
                                                        (all_str[state.curline].dv[0] + 1.)
                         widget_control, state.vmax_id, set_value=v1
                      end
                      else: 
                  endcase
                  lls_idlines_update, state, flg=1
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
              ;; Switch
              '-': lls_idlines_prev, state
              '=': lls_idlines_next, state
              'm': begin ;; Work on the nearest line
                 rwv=state.xpos/(1+state.zgal)
                 mn = min(abs(all_str.wrest-rwv), imn)
                 state.lastline = state.curline
                 state.curline = imn
                 lls_idlines_setup, state, /NOWIN
              end
              ;; Help
              'H': x_helpwidg, state.help
              ;; Blend
              'B': begin 
                 ;; Set flag
                 if (all_str[state.curline].flg_eye MOD 4) LT 2 then $
                    all_str[state.curline].flg_eye += 2
                 ;; Button
                 widget_control, state.eye_id, get_value=eye_val
                 eye_val[1] = 1
                 widget_control, state.eye_id, set_value=eye_val
              end
              ;; Saturated
              'S': begin 
                 ;; Set flag
                 if (all_str[state.curline].flg_eye MOD 8) LT 4 then $
                    all_str[state.curline].flg_eye += 4
                 ;; Button
                 widget_control, state.eye_id, get_value=eye_val
                 eye_val[2] = 1
                 widget_control, state.eye_id, set_value=eye_val
                 ;; Lower limit
                 all_str[state.curline].flg_limit = 2
                 widget_control, state.limit_id, set_value=1
              end
              ;; Continuum
              '3': lls_idlines_continuum, state, 1L
              '4': lls_idlines_continuum, state, 2L
              ;; Detection
              'D': begin 
                 ;; Set flag
                 if (all_str[state.curline].flg_eye MOD 2) LT 1 then $
                    all_str[state.curline].flg_eye += 1
                 ;; Button
                 widget_control, state.eye_id, get_value=eye_val
                 eye_val[0] = 1
                 widget_control, state.eye_id, set_value=eye_val
              end
              'N': begin ;; No Good
                 all_str[state.curline].flg_anly = 0
                 widget_control, state.anly_id, set_value=0
              end
              'U': begin ;; Upper limit
                 if all_str[state.curline].flg_limit NE 3 then begin
                    all_str[state.curline].flg_limit = 3
                    widget_control, state.limit_id, set_value=2
                 endif else begin
                    all_str[state.curline].flg_limit = 1
                    widget_control, state.limit_id, set_value=0
                 endelse
              end
              'L': begin ;; Lower Limit
                 if all_str[state.curline].flg_limit NE 2 then begin
                    all_str[state.curline].flg_limit = 2
                    widget_control, state.limit_id, set_value=1
                 endif else begin
                    all_str[state.curline].flg_limit = 1
                    widget_control, state.limit_id, set_value=0
                 endelse
              end
              ; ZOOMING
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'w': state.xymnx = state.svxymnx ; Reset the screen
              'W': begin
                 state.xymnx[0] = min(wv, max=mxw)
                 state.xymnx[2] = mxw
                 gdw = where(wv GT state.xymnx[0] and wv LT state.xymnx[2], npx)
                 srt = sort(fx[gdw])
                 md = fx[gdw[srt[round(0.9*npx)]]]
                 state.xymnx[3] = md*1.5
                 state.xymnx[1] = -0.1 * md
              end
              ;; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              ;; Regions
              '0': all_str[state.curline].dv = all_str[state.lastline].dv
              'T': begin ;; Set all lines to this value
                 old_vmnx = state.def_vmnx
                 state.def_vmnx = all_str[state.curline].dv
                 ;; Tied
                 tied = where(abs(all_str.dv[0]-old_vmnx[0]) LT 1e-2 and $
                              abs(all_str.dv[1]-old_vmnx[1]) LT 1e-2, ntied)
                 if ntied GT 0 then all_str[tied].dv = state.def_vmnx
              end
              'F': begin ;; Force all lines to this value
                 state.def_vmnx = all_str[state.curline].dv
                 all_str.dv = state.def_vmnx
              end
              'M': begin ;; Set all metal-lines to this value
                 old_vmnx = state.def_vmnx
                 state.def_vmnx = all_str[state.curline].dv
                 ;; Tied
                 tied = where(abs(all_str.dv[0]-old_vmnx[0]) LT 1e-2 and $
                              abs(all_str.dv[1]-old_vmnx[1]) LT 1e-2 and $
                              strmid(all_str.ionnm,0,2) NE 'HI', ntied)
                 if ntied GT 0 then all_str[tied].dv = state.def_vmnx
              end
              'f': all_str[state.curline].dv = state.def_vmnx  ;; Set current line to the saved dv
              'q': begin
                 lls_idlines_save, state
                 print, 'lls_idlines: Exiting..'
                 widget_control, ev.top, /destroy
                 return
              end
              else:
           end
          lls_idlines_update, state, flg=1
       end
      else :
   endcase
      
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro lls_idlines_plot, state
  
  common lls_idlines_cmm

  ;; Set plot window
  clr = getcolor(/load)
  

  ;; MgII Plot
  !p.multi = [0,1,1]
  widget_control, state.si2_id, get_value=wind
  wset, wind
  ymnx = [state.xymnx[1], state.xymnx[3]]
  plot, wv, fx, color=clr.black, background=clr.white, psym=10, $
        position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=ymnx, xstyle=1, ystyle=1, xtit='Wavelength (Ang)', $
        ytit='Relative Flux'

  ;; Data
;  oplot, wv, (1.-all_fit)*conti, color=clr.red, thick=1
;  oplot, wv, conti, color=clr.gray, linesty=1, thick=3
  oplot, wv, sig, color=clr.orange, thick=1, psym=10
  oplot, wv, conti, color=clr.cyan, thick=2, psym=10

  ; Find Min and Max in region
  wvmin = state.xymnx[0]/(state.zgal+1.)
  wvmax = state.xymnx[2]/(state.zgal+1.)

  ; Parse list
  allwv = where(all_str.wrest GE wvmin AND all_str.wrest LE wvmax, cntwv)

  ; Plot

  for q=0L,cntwv-1 do begin
     lsty=0
     ymax = state.xymnx[1] + 0.5*(state.xymnx[3]-state.xymnx[1])
     if allwv[q] EQ state.curline then pclr=clr.blue $ 
     else begin
        lsty=2
        pclr = clr.darkgreen
        ymax = state.xymnx[1] + 0.2*(state.xymnx[3]-state.xymnx[1])
     endelse

     ;; Region
     wvmnx = (all_str[allwv[q]].dv/3e5 + 1) * all_str[allwv[q]].wrest
;     if cntwv GT 1 then stop
     gdr = where(wv GT wvmnx[0]*(1+state.zgal) and wv LT wvmnx[1]*(1+state.zgal), ngdr)
     if ngdr GT 0 then begin
        rclr = clr.blue
        ;; Color
        if (all_str[allwv[q]].flg_eye mod 2) GE 1 then rclr = clr.green ; detection
        if (all_str[allwv[q]].flg_eye mod 4) GE 2 then rclr = clr.orange ; blend
         if (all_str[allwv[q]].flg_eye mod 8) GE 4 then rclr = clr.purple ; saturated
        if all_str[allwv[q]].flg_anly EQ 0 then rclr = clr.red ; NG
        ;; Paint
        oplot, wv[gdr], fx[gdr], color=rclr, psym=10, linesty=lsty, thick=2
     endif

     ;; Line
     xplt = all_str[allwv[q]].wrest*(1.+state.zgal)
     xyouts, xplt, ymax, $ ; QAL
             strtrim(all_str[allwv[q]].ionnm,2), $;+' '+string(lines[allwv[q]].fval,format='(f6.4)'), $
             color=pclr, $
             orientation=90., charsize=1.5, align=0.5
     oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
            color=pclr, linestyle=lsty

  endfor

  ;; Continuum Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.green, symsize=5
   endif

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro lls_idlines_update, state, FLG=flg
  common lls_idlines_cmm
 

  lls_idlines_plot, state
  lls_idlines_velplt, state

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup 
pro lls_idlines_setup, state, NOWIN=nowin

  common lls_idlines_cmm
  if not keyword_set(NOWIN) then lls_idlines_setwin, state
  lls_idlines_update, state

  ;; Velocity
  widget_control, state.vmin_id, set_value=all_str[state.curline].dv[0]
  widget_control, state.vmax_id, set_value=all_str[state.curline].dv[1]

  ;; Droplist
  a = where(abs(lines.wave-all_str[state.curline].wrest) LT 1e-3, na)
  if na NE 1 then stop
  widget_control, state.droplist_id, set_droplist_select=a[0]

  ;; Buttons
  widget_control, state.limit_id, set_value=all_str[state.curline].flg_limit-1 
  widget_control, state.anly_id, set_value=all_str[state.curline].flg_anly 

  eye_val = lonarr(3)
  eye_val[0] = (all_str[state.curline].flg_eye MOD 2) EQ 1
  eye_val[1] = (all_str[state.curline].flg_eye MOD 4) GE 2
  eye_val[2] = (all_str[state.curline].flg_eye MOD 8) GE 4
  widget_control, state.eye_id, set_value=eye_val

  ;; Line list
  llist = 'lls_plt_lines.lst'
  plt_lines = x_setllst(llist, 0)
  state.ntrans = n_elements(plt_lines)
  state.velplt[0:state.ntrans-1].wrest = plt_lines.wave
  state.velplt[0:state.ntrans-1].name = plt_lines.name
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.39]

  ;; Set flg
  gd = where(state.velplt.wrest*(state.zgal+1) GT min(wv) AND $
             state.velplt.wrest*(state.zgal+1) GT (state.zgal+1.)*1215.6701 AND $
             state.velplt.wrest*(state.zgal+1) LT max(wv), na)
  if na GT 0 then begin
     state.velplt[*].flg = 0
     state.nplt = na
     state.velplt[gd].flg = 1
     
     ;; Setup the plots
     state.all_velo[*,gd] = x_allvelo(wv, state.zgal, $
                                      state.velplt[gd].wrest,$
                                      state.velp_vmnx, $
                                      all_pmnx=all_pmnx, NPIX=5000L)
     ;; Set view
     for jj=0L,na-1 do begin
        ;; Continuum exists?
        if abs(median(conti[all_pmnx[0,jj]:all_pmnx[1,jj]])-1.) GT 1.e-2 then continue
        ;;
        all_fx = fx[all_pmnx[0,jj]:all_pmnx[1,jj]]
        np = n_elements(all_fx)
        if np GT 1 then begin  ;; Some fall in artificial gaps
           srt = sort(all_fx)
           ymx = 1.1 * all_fx[srt[round(0.9*np)]]
           ymn = -0.1 * all_fx[srt[round(0.9*np)]]
           state.velplt[gd[jj]].ymnx = [ymn, ymx]
           ;print, state.velplt[gd[jj]].wrest, ymn, ymx
           ;stop
        endif
     endfor
     ;; Pixels
     state.all_pmnx[*,gd] = all_pmnx
  endif
    
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup 
pro lls_idlines_setwin, state

  common lls_idlines_cmm

  wvobs = all_str[state.curline].wrest*(1+state.zgal)
  dv = all_str[state.curline].dv
  dv[0] = state.scale_win*dv[0] < ((-50.)*state.scale_win)
  dv[1] = state.scale_win*dv[1] > dv[1]+50.*state.scale_win
  dlam = dv/3e5 * wvobs

  state.xymnx[0] = wvobs + dlam[0]
  state.xymnx[2] = wvobs + dlam[1]
  gdw = where(wv GT state.xymnx[0] and wv LT state.xymnx[2], npx)
  srt = sort(fx[gdw])
  md = fx[gdw[srt[round(0.9*npx) < (npx-1)]]]
  state.xymnx[3] = md*1.5
  state.xymnx[1] = -0.1 * md

  state.svxymnx = state.xymnx
    
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro lls_idlines_next, state

  common lls_idlines_cmm

  state.lastline = state.curline
  state.curline = (state.curline+1) < (n_elements(all_str) - 1)
  lls_idlines_setup, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Back 
pro lls_idlines_prev, state
  common lls_idlines_cmm

  ;; Match?
  state.lastline = state.curline
  state.curline = (state.curline-1) > 0
  lls_idlines_setup, state

  return

end

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Velocity Plot
pro lls_idlines_velplt, state
  
  common lls_idlines_cmm

  ; Set plot window
  widget_control, state.velplt_id, get_value=wind
  wset, wind

  clr = getcolor(/load)
  
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)

  !p.multi = [0,1,ngd+1,0,1]

  for j=0L,ngd-1 do begin
      i = gdlin[j]
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]

      ;; Plot
      if  (j NE ngd-1 ) then begin
          spaces = replicate('!17 ',30)
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fx[pixmin:pixmax]/conti[pixmin:pixmax], xrange=state.velp_vmnx, $
            yrange=state.velplt[i].ymnx, xtickn=spaces, xmargin=[9,3], $
            ymargin=[0,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fx[pixmin:pixmax]/conti[pixmin:pixmax], xrange=state.velp_vmnx, $
            yrange=state.velplt[i].ymnx, xmargin=[9,3], ymargin=[3,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse

      ;; Labels
      xyouts, 0.07*(state.velp_vmnx[1]-state.velp_vmnx[0])+state.velp_vmnx[0], $
        state.velplt[i].ymnx[0]+ $
        (state.velplt[i].ymnx[1]-state.velplt[i].ymnx[0])*0.05, $
        strtrim(state.velplt[i].name,2), $
        color=clr.black, charsize=1.5

      ;; Color region
      wrest = state.velplt[i].wrest
      rclr = clr.blue
      if abs(wrest-1215.6701) LT 1e-3 then begin
         dv = HI_str.dv 
      endif else begin
         aidx = where(abs(all_str.wrest-wrest) LT 1e-3, nmt)
         if nmt NE 1 then stop
         if (all_str[aidx].flg_eye mod 2) GE 1 then rclr = clr.green ; detection
         if (all_str[aidx].flg_eye mod 4) GE 2 then rclr = clr.orange ; blend
         if (all_str[aidx].flg_eye mod 8) GE 4 then rclr = clr.purple ; saturated
         if all_str[aidx].flg_anly EQ 0 then rclr = clr.red           ; NG
         dv = all_str[aidx].dv
      endelse

      gdr = where(state.all_velo[0:state.all_pmnx[2,i],i] GT dv[0] AND $ 
                  state.all_velo[0:state.all_pmnx[2,i],i] LT dv[1], ngdr) 
      if ngdr GT 0 then $
         oplot, state.all_velo[gdr,i], fx[pixmin+gdr]/conti[pixmin+gdr], $ 
                color=rclr, psym=10, thick=2
      
      ;; Lines
      oplot, [0., 0.], state.velplt[i].ymnx, color=clr.blue, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.green, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  endfor

  !p.multi = [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;
;  Continuum
;;;;;;;;;;;;;;;;;;;;

pro lls_idlines_continuum, state, flg

  common lls_idlines_cmm
  if not keyword_set( flg ) then return

; CASE

  case flg of 
      0: stop ;  RESET to 1/2
      1: begin ; Set point
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

  if state.cstr.npts NE 0 then begin
      gdmsk = where(state.cstr.msk EQ 1, nmsk)
      if nmsk GE 3 then begin
         OUT_CFLG = 1
         mn = min(state.cstr.xval[gdmsk], imn)
         mx = max(state.cstr.xval[gdmsk], imx)
         ;; Interp
         gd = where(state.wave GT mn AND state.wave LT mx, ngd)
         if ngd NE 0 then begin
            ;; SORT
            srt = sort(state.cstr.xval[gdmsk])
            swv = sort(state.wave[gd])
            twv = state.wave[gd[swv]]
            state.conti[gd[swv]] = xspline(state.cstr.xval[gdmsk[srt]], $
                                           state.cstr.yval[gdmsk[srt]], $
                                           twv)
         endif
      endif
   endif

  ;; Write to common
  conti = state.conti

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lls_idlines_save, state

  common lls_idlines_cmm

  ;; Sort
  srt = sort(all_str.wrest)
  all_str = all_str[srt]
  print, 'lls_idlines: Saving..', state.outfil
  mwrfits, all_str, state.outfil, /create
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
;; lls_idlines, 'J0042-1037_nbin3_coadd.fits', '~/tmp.fits', 0.095d, infl=3
pro lls_idlines, specfil, outfil, zgal, LLIST=llist, $
                 CFIL=cfil, CLOBBER=clobber, $
                 XSIZE=xsize, YSIZE=ysize, INFLG=inflg, $
                 INSTR=instr, FLG_LLST=flg_llst, ILINES=ilines, $
                 FIRST_LINE=first_line, DATA=data, DEF_VMNX=def_vmnx, $
                 STRCT=strct, SCALE_WIN=scale_win, SKIP_NSTOP=skip_nstop, $
                 IHI_STR=ihi_str, CONTI_STRCT=conti_strct, EDIT=edit

  common lls_idlines_cmm
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'lls_idlines, specfil, outfil, zgal, LLIST=, INFLG=, INSTR= [v1.0]'
    return
  endif 

  ;; Check to see if outfil exists
  a = file_search(outfil+'*', count=nfil)
  if nfil GT 0 and not keyword_set(CLOBBER) and not keyword_set(EDIT) then begin
     print, 'lls_idlines: Outfil exists.  Use /CLOBBER to over-write!'
     return
  endif
  ;; Screen size
  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100

; Initialize the common blcok
  if not keyword_set( DEF_VMNX ) then  def_vmnx=[-600., 600]
  if not keyword_set( SCALE_WIN ) then  scale_win = 2. ;; Viewing parameter

  lls_idlines_icmmn, specfil, outfil, IFLG=inflg, DATA=data, CFIL=cfil

  ;; COS/Halos line list
  if not keyword_set(ILINES) then begin
     if not keyword_set(LLIST) then begin
        if not keyword_set(FLG_LLST) then flg_llst = 0L
        case flg_llst of
           0: llist = 'lls_lines.lst'
           1: stop ; llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/mstrong_dla.lst'
           else: stop
        endcase
     endif
     lines = x_setllst(llist, 0)
  endif else lines = ilines
  flg_lines = 1

  gdl = where(lines.wave GT min(wv,max=mxw)/(1+zgal) and $
              lines.wave LT mxw/(1+zgal), nlines)
  lines = lines[gdl]

  ;; Initialize the structure
  if not keyword_set(INSTR) and not keyword_set(EDIT) then begin
     all_str = replicate({ismabndstrct},nlines)
     all_str.datfil = specfil
     all_str.wrest = lines.wave
     all_str.ionnm = lines.name
     all_str.dv = def_vmnx
     all_str.flg_anly = 2
     all_str.flg_eye = 0
     all_str.flg_limit = 1
     all_str.zabs = zgal
  endif else begin
     if keyword_set(EDIT) and not keyword_set(INSTR) then instr=outfil
     all_str = xmrdfits(instr, 1, /silen)
     if n_elements(all_str) NE nlines then begin
        print, 'lls_idlines: Your linelist is smaller/bigger than mine!'
        print, 'lls_idlines: Am going to create a new one if you choose to proceed'
        if not keyword_set(SKIP_NSTOP) then stop
        ;; 
        tmp_str = all_str
        all_str = replicate({ismabndstrct},nlines)
        all_str.wrest = lines.wave
        all_str.ionnm = lines.name
        all_str.dv = def_vmnx
        all_str.flg_anly = 2
        all_str.flg_eye = 0
        all_str.flg_limit = 1
        all_str.zabs = zgal
        ;; Fill in
        for ii=0L,n_elements(tmp_str)-1 do begin
           mt = where(abs(tmp_str[ii].wrest - all_str.wrest) LT 1d-3, nmt)
           case nmt of 
              0: 
              1: all_str[mt] = tmp_str[ii]
              else: stop
           endcase
        endfor
     endif
  endelse

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }
  tmp3 = { contistrct }

  ;; STATE
  state = {             $
          curline: 0L, $
          lastline: 0L, $
          outfil: outfil, $
          specfil: specfil, $
          inflg: 0, $
          model: fltarr(npix), $
          zgal: double(zgal), $
          def_vmnx: def_vmnx, $
          wave: wv, $
          fx: fx, $
          conti: conti, $
          cstr: tmp3, $
          velp_vmnx: [-500., 500], $
          xpos: 0.0, $
          ypos: 0.0, $
          ipress: 0L, $
          pos: [0.1,0.1,0.95,0.95], $ ; Plotting
          flg_zoom: 0, $
          psfile: 0, $
          ntrans: 0L, $               ; PLOTTING LINES
          nplt: 0L, $               ; PLOTTING LINES
          velplt: replicate(tmp, 300), $
          all_velo: dblarr(5000, 300), $  
          all_pmnx: lonarr(3, 300), $  
          help: strarr(50), $
          scale_win: scale_win, $ ;; Automated viewing
          mxxymnx: fltarr(4), $
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
          droplist_id: 0L, $
          tdraw_id: 0L, $
          anly_id: 0L, $
          eye_id: 0L, $
          info_id: 0L, $
          si2_id: 0L, $
          limit_id: 0L, $
          wvalu_id: 0L, $
          vmin_id: 0L, $
          vmax_id: 0L, $
          velbase_id: 0L, $
          velplt_id: 0L, $
          main_id: 0L, $
          hsew_id: 0L, $
          help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES

  if keyword_set(INFLG) then state.inflg = inflg

;    WIDGET
  base = WIDGET_BASE( title = 'lls_idlines: Check spectra', /row, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  ;; Info + Main draw window
  xmain = round(2*xsize/3.)
  state.main_id = $
    WIDGET_BASE( state.base_id, /column, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, ysize=ysize, xsize=xmain)
  state.info_id = $
    WIDGET_BASE( state.main_id, /row, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, ysize=ysize/5.)
  state.draw_base_id = widget_base(state.main_id, /column, /base_align_left, $
                                   xsize=xmain, $
                                   ysize=round(4*ysize/5), $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)
  state.si2_id = widget_draw(state.draw_base_id, xsize=xmain, $
                              ysize=round(4*ysize/5), /frame, retain=2, $
                              uvalue='SI2DRAW', /button_even, /motion_events)
  state.text_id = widget_text(state.draw_base_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
  state.size[0] = xmain
  state.size[1] = round(4*ysize/5)

  ;; Velocity plot
  state.velbase_id = widget_base(state.base_id, /column, /base_align_right, $
                                   xsize=xsize-xmain, $
                                   ysize=ysize, $
                                   uvalue='VEL_BASE', frame=2)
  state.velplt_id = widget_draw(state.velbase_id, xsize=xsize-xmain, $
                                ysize=ysize, /frame, retain=2, $
                                uvalue='VEL_DRAW')

;;;;;; Info window ;;;;;;;;;;;
;               ysize=round(ysize/3.))
  ;; Info
  civinf = widget_base(state.info_id, /column, /align_center, frame=2)
;  state.name_id = cw_field(civinf, title='Obj ', value=' ', xsize=18)
  state.zgal_id = cw_field(civinf, title='z_galaxy', value=state.zgal, /floating, $
                           /column, xsize=8, /return_events, uvalue='ZGAL')
  handinf = widget_base(state.info_id, /column, /align_center, frame=2)


  ;; ;;;;;;;;
  ;; Line
  wave = widget_base(state.info_id, /column, /align_center, frame=2)
  linelist = lines.name
  state.droplist_id = widget_droplist(wave, $
                                      frame = 1, $
                                      title = 'Lines:', $
                                      uvalue = 'LLIST', $
                                      value = linelist)
  if not keyword_set(FIRST_LINE) then first_line = 1215.6701d
  ;; Initialize to Lya if possible
  a = where(abs(lines.wave-first_line) LT 1e-3, na)
  if na NE 1 then a = 0L ;; Just take the really first line
  widget_control, state.droplist_id, set_droplist_select=a[0]
  b = where(abs(all_str.wrest-lines[a].wave) LT 1e-3, nb)
  if nb NE 1 then stop
  state.curline = b[0]

  ;; ;;;;;;;;
  ;; Flags
  flags = widget_base(state.info_id, /row, /align_center, frame=2)
  state.limit_id = cw_bgroup(flags, ['Normal','Lower', 'Upper'], $
                             /exclusive, row=3, UVALUE='LIM_FLAG', /frame)
  state.anly_id = cw_bgroup(flags, ['NG','Measure','Use'], $
                            /exclusive, row=3, UVALUE='ANLY_FLAG', /frame)
  state.eye_id = cw_bgroup(flags, ['Detected','Blended','Saturated'], /frame, /nonexclu, $
                           row=2, UVALUE='EYE_FLAG', SET_VALUE=[0,0,0])
  widget_control, state.limit_id, set_value=0
  widget_control, state.anly_id, set_value=2

  ;; ;;;
  ;; Region
  region = widget_base(state.info_id, /column, /align_center, frame=2)
  state.vmin_id = cw_field(region, title='vmin: ', value=state.def_vmnx[0], /floating, $
                           xsize=7, /return_events, uvalue='VMIN')
  state.vmax_id = cw_field(region, title='vmax: ', value=state.def_vmnx[1], /floating, $
                           xsize=7, /return_events, uvalue='VMAX')



  ;; Saving
  butbase2 = widget_base(state.info_id, /column, /align_center, frame=2)
  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  vplt = WIDGET_BUTTON(butbase2, value='VPLT',uvalue='VPLT')
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
             'LMB/RMB -- Modify velocity region' $
            ]
  nhelp = n_elements(strhelp)
  for ii=0l,nhelp-1 do state.help[ii] = strhelp[ii]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  if not keyword_set(UPDATE) then WIDGET_CONTROL, base, /realize


  ;; PLOT
  lls_idlines_setup, state
  if not keyword_set(UPDATE) then begin
     lls_idlines_update, state
  endif else begin ;; Just evaluate all the lines
     stop ;; Not functional
     gd = where(fitstr.wrest GT 0,ngd)
     if ngd EQ 0 then begin
        print, 'lls_idlines: No lines to update!  Returning...'
        return
     endif
     print, 'lls_idlines: Updating ', ngd, ' lines'
     for ii=0L,ngd-1 do begin
        state.curline = gd[ii]
        sv_flg = fitstr[state.curline].flg_flux
        lls_idlines_getconti, state
        lls_idlines_fit, state
        fitstr[state.curline].flg_flux = sv_flg
        if sv_flg EQ 2 then begin
           fitstr[state.curline].flux = state.boxEW
           fitstr[state.curline].fsig = state.sigboxEW 
        endif
     endfor
     print, 'lls_idlines: Writing ', state.outfil
     mwrfits, fitstr[gd], state.outfil, /create
     return
  endelse
     
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'lls_idlines', base
  delvarx, fx, wv, npix, sig, fitstr
  strct = all_str

  conti_strct = { $
              conti: CONTI, $
              flg_conti: OUT_CFLG $
              }

  return
end
	
