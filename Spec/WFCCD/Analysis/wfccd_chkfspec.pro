;+ 
; NAME:
; wfccd_chkfspec
;    Version 1.0
;
; PURPOSE:
;   Plots a series of spectra to allow a quick check
;
; CALLING SEQUENCE:
;   
;   wfccd_chkfspec, wfccd, maskid, expsr, XSIZE=, YSIZE=
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
;   wfccd_chkfspec, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Events
;;;;

pro wfccd_chkfspec_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DRAW' : begin
          case ev.type of
              0 : ; Button press
              1 :  ; Button release
              2 : begin ; Motion event
                  state.xcurs = ev.x
;                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state.xcurs,state.pos, $
                                         state.xymnx,$
                                         state.size) 
;                  state.ypos = xgety_plt(state.ycurs,state.pos, $
;                                         state.xymnx,$
;                                         state.size)
                  widget_control, state.xpos_id, set_value=state.xpos
;                  widget_control, state.ypos_id, set_value=state.ypos
              end
              else:
          endcase
      end
      'PRINT' : wfccd_chkfspec_print, state
      'NEXT' : begin
          state.curpg = state.curpg + 1
          if state.curpg GT state.npg-1 then state.curpg = 0L
          wfccd_chkfspec_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'PREV' : begin
          state.curpg = state.curpg - 1
          if state.curpg LT 0 then state.curpg = state.npg-1
          wfccd_chkfspec_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'PAGE': begin
          widget_control, state.pg_id, get_value=pg
          state.curpg = 0 > pg < (state.npg-1)
          wfccd_chkfspec_Plot, state
      end
      'NSPEC': begin
          widget_control, state.nspec_id, get_value=nplt
          state.nplt = nplt
          state.npg =  (state.nobj/state.nplt) + (state.nobj MOD state.nplt NE 0)
          state.curpg = state.curpg < (state.npg-1)
          wfccd_chkfspec_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkfspec_Plot, state
  
; Plot Data

  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  !P.MULTI= [0,1,state.nplt]

  pmin = state.curpg*state.nplt
  if state.curpg EQ state.npg - 1 then pmax = state.nobj-1 $
  else pmax = (state.curpg+1)*state.nplt-1
  
  for i=pmin,pmax do begin
      npix = state.obj[i].npix
      ; Get median flux
      mdwv = where(state.obj[i].wave[0:npix-1] GT state.mdrng[0] AND $
                   state.obj[i].wave[0:npix-1] LT state.mdrng[1] AND $
                   state.obj[i].fx[0:npix-1] GT 0., nwv)
      if nwv NE 0 then mdfx = median(state.obj[i].fx[mdwv]) else mdfx = 1.

      ; Plot
      if i NE pmax then begin
          spaces = replicate('!17 ',30)
          plot, state.obj[i].wave[0:npix-1], $
            state.obj[i].fx[0:npix-1], $
            xrange=[state.xymnx[0],state.xymnx[2]], $
            yrange=[0., 2*mdfx], xtickn=spaces, xmargin=[14,0], ymargin=[0,0], $
            charsize=1.3, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.obj[i].wave[0:npix-1], $
            state.obj[i].fx[0:npix-1], $
            xrange=[state.xymnx[0],state.xymnx[2]], $
            yrange=[0., 2*mdfx], xmargin=[14,0], ymargin=[3,0], $
            charsize=1.3, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse
      ; Error
      oplot, state.obj[i].wave[0:npix-1], state.sig[0:npix-1,i], $
        color=clr.green

      ; Zfind
      if state.flg_zfind EQ 1 AND state.obj[i].zans.z_err NE 0. then begin
          synth = x_synthspec(state.obj[i].zans, $
                            loglam=alog10(state.obj[i].wave[0:npix-1]))
          oplot, state.obj[i].wave[0:npix-1], synth, color=clr.red
          xyouts, 0.13*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], $
            0.85*2*mdfx, 'z= '+string(state.obj[i].zans.z,'(f7.4)'), $
            color=clr.red, charsize=1.5
          xyouts, 0.13*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], $
            0.70*2*mdfx, 'chi= '+string(state.obj[i].zans.rchi2,'(f5.2)'), $
            color=clr.red, charsize=1.5
      ; Key lines
      ; Key lines
          for jj=0L,state.nkey-1 do begin
              wv=state.keylines[jj]*(state.obj[i].zans.z+1)
              oplot, [wv,wv], [-1e20,1e20], color=clr.gray, linestyle=2
          endfor
      endif

      ; Labels
      xyouts, 0.03*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.75*2*mdfx, $
        strtrim(state.obj[i].slit_id,2)+state.obj[i].obj_id, $
        color=clr.green, charsize=1.5
      xyouts, 0.03*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.45*2*mdfx, $
        'Obj: '+strtrim(i,2), color=clr.blue, charsize=1.5
  endfor

  !P.MULTI= [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkfspec_print, state

  widget_control, /hourglass   
; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
;  !p.charthick = 1

  device, decompose=0
  ; Get file name
;  ipos = strpos(state.obj_fil, '.fits')
;  psfil = strmid(state.obj_fil, 0, ipos)+'.ps'
  psfil = 'fspec_chk.ps'
  ps_open, file=psfil, font=1, /color
  state.psfile = 1
  for qq=0L,state.npg-1 do begin
      state.curpg = qq
      wfccd_chkfspec_plot, state
  endfor
  ps_close, /noprint, /noid
  spawn, 'gzip -f '+psfil
  device, decomposed=svdecomp
  state.psfile = 0
;  !p.thick = 1
;  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_chkfspec, fspec_fil, ZFIND=zfind, XSIZE=xsize, YSIZE=ysize, $
                   PRINTONLY=printonly, NPLT=nplt

;common wfccd_chkfspec_cmm

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_chkfspec, fspec_fil, /zfind, XSIZE=, YSIZE=, /PRINTONLY (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 900
  if not keyword_set( XYMNX ) then xymnx = [3800., 0.0, 8700., 1.0]
  if not keyword_set( NPLT ) then nplt = 6L

; Check for file
;  flg = x_chkfil(fspec_fil)
;  if flg NE 1 then return
      
  
;  Object Structure
  wfccd_wrfspec, wffspec, fspec_fil, /read

  gdfspec = where(wffspec.flg_anly NE 0, nobj)

; STATE

  state = {             $
            nplt: nplt, $
            npg: (nobj/nplt) + (nobj MOD nplt NE 0), $
            nobj: nobj, $
            obj: wffspec[gdfspec], $
            sig: fltarr(5000,nobj), $
            curpg: 0L, $
            mdrng: [5000., 8000.], $
            pos: [0.04,0.0,1.00,1.0], $ ; Plotting
            flg_zfind: 0, $
            keylines: fltarr(100), $
            nkey: 0L, $
            xymnx: xymnx, $
            tmpxy: fltarr(4), $
            xcurs: 0., $
            xpos: 0., $
            size: lonarr(2), $
            psfile: 0, $
            base_id: 0L, $      ; Widgets
            draw_id: 0L, $
            drawbase_id: 0L, $
            pg_id: 0L, $
            nspec_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            help_text_id: 0L $
          }

; Create Sigma array

  for i=0L,state.nobj-1 do begin
      a = where(state.obj[i].var LE 0., COMPLEMENT=b)
      if a[0] NE -1 then state.sig[a,i] = 0.
      if b[0] NE -1 then state.sig[b,i] = sqrt(state.obj[i].var[b])
  endfor

; Zfind
  if keyword_set( ZFIND ) then begin
      state.flg_zfind = 1
      state.nkey = 9L
      state.keylines = [ 3727., $
                         4861.3, $
                         5007., $
                         6562.8, $
                         3933.7, $
                         3968.5, $
                         4000.0, $
                         4304.4, $
                         5892.5]
  endif
                         

;    WIDGET
  base = WIDGET_BASE( title = 'wfccd_chkfspec: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L, yoffset=50L)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='wfccd_chkfspec', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
;  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)

;;;;;;;;;;;;;;;;;;
; Pages
  state.pg_id = cw_field(toolbar, value=0L, $
                             /long, /return_events, xsize=3,$
                             title='Page', UVALUE='PAGE')
  state.nspec_id = cw_field(toolbar, value=state.nplt, $
                             /long, /return_events, xsize=3,$
                             title='Nplt', UVALUE='NSPEC')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  state.drawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_id = widget_draw(state.drawbase_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')

; XY position
  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
;  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  butbase2 = widget_base(toolbar, /column, /align_center)
  print = WIDGET_BUTTON(butbase2, value='PRINT',uvalue='PRINT')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
;  strhelp = strarr(50)
;  strhelp = ['   Help Menu   ',$
;             'LMB - Truncate/Extend trace', $ 
;             'RMB - Contrast/Brightness', $
;             'CMB/CMB - Zoom' $ 
;             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

; Realize
  WIDGET_CONTROL, base, /realize

  ; PLOT
  wfccd_chkfspec_Plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'wfccd_chkfspec', base, /no_block

  !P.MULTI= [0,1,1]
  return
end

