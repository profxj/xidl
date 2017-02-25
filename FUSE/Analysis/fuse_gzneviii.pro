;+ 
; NAME:
; fuse_gzneviii
;    Version 1.1
;
; PURPOSE:
;   GUI which plots the spectra in two strips corresponding to the
;   doublet of OVI.  User can indicate the regions where OVI
;   could be detected (visually).
;
; CALLING SEQUENCE:
;   
;   fuse_gzneviii, datfil, XSIZE=, YSIZE=, OUTFIL=, GZFIL=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XSIZE  -- Size of gui x-pixels [default: 80% of screen] 
;  YSIZE  -- Size of gui y-pixels [default: 80% of screen]
;  GZFIL  -- File containing saved regions
;
; OPTIONAL OUTPUTS:
;  SVSTATE  -- Save the state to return to search later
;  OUTFIL   -- File to contain saved regions 
;              [default: 'find_neviii.fits']
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_gzneviii, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro fuse_gzneviii_icmmn, datfil, gzfil

  common fuse_gzneviii_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    gzreg

  ;; Read FUSE data
  fx = x_readspec(datfil, inflg=3, sig=sig, wav=wv, npix=npix)

  if keyword_set( GZFIL ) then $
    gzreg = xmrdfits(gzfil, 1, /silent) else begin
      tmp2 = { $
               zmin: 0.d, $
               zmax: 0.d, $
               flg: 0 $
             }
      gzreg = replicate(tmp2, 1000L)
  endelse

  ;; Sort
  srt = sort(wv)
  wv = wv[srt]
  fx = fx[srt]
  sig = sig[srt]

  return
end
  

;;;;
; Events
;;;;

pro fuse_gzneviii_event, ev

  common fuse_gzneviii_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SAVE': fuse_gzneviii_svgz, state
      'NEXT': begin
          fuse_gzneviii_next, state
          fuse_gzneviii_update, state
      end
      'PREV': begin
          fuse_gzneviii_prev, state
          fuse_gzneviii_update, state
      end
      'LDRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.ltext_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'LDRAW' : begin
          widget_control, state.ltext_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin         ; Button press
                  state.ipress = ev.press
                  case state.ipress of
                      1 : begin ;; Set lhs
                          gzreg[state.ngz].zmin = xgetx_plt(state, /strct) 
                      end
                      2 : begin ;; New region
                          state.ngz = state.ngz + 1
                          gzreg[state.ngz].zmin = $
                            xgetx_plt(state, /strct) - 0.0005
                          gzreg[state.ngz].zmax = $
                            xgetx_plt(state, /strct) + 0.0005
                          gzreg[state.ngz].flg = 1
                      end
                      4 : begin ;; Set rhs
                          gzreg[state.ngz].zmax = xgetx_plt(state, /strct) 
                      end
                  endcase
                  fuse_gzneviii_update, state
              end
              1 :
              2 : begin         ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
              end
          endcase
      end
      'LTEXT' : begin
          eventch = string(ev.ch)
          case eventch of
              ; ZOOM
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; top
              'W': state.xymnx = state.svxymnx
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              'd': fuse_gzneviii_rmreg, state
              else: 
          endcase
          fuse_gzneviii_update, state
      end
;;;;;;;;; Metals ;;;;;;;;;;
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          zabs = tmp
          flg_lines = 1
      end
      'SDRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.stext_id, sensitive = 0
          endif
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
;  Lya Plot
;;;;;;;;;;;;;;;;;;;;


pro fuse_gzneviii_spectra, state
  
  common fuse_gzneviii_cmm

  ;; Top window
  widget_control, state.ldraw_id, get_value=wind
  wset, wind

  clr = getcolor(/load)
  
  ;; Plot top
  plot, state.z1arr, fx, $
    xrange=[state.xymnx[0], state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
    charsize=1.2, psym=10, background=clr.white, color=clr.black, $
    xmargin=[0,0], ymargin=[0,0], xstyle=1,$
    ystyle=1, xtickn=spaces, /nodata

  ;; Regions
  gd = where(gzreg.flg NE 0 AND gzreg.zmin GT state.xymnx[0] AND $
             gzreg.zmax LT state.xymnx[2], ngd)
  for jj=0L,ngd-1 do begin
      id = gd[jj]
      if id EQ state.ngz then fclr = clr.orange else fclr=clr.yellow
      polyfill, [gzreg[id].zmin, gzreg[id].zmin, $
                 gzreg[id].zmax, gzreg[id].zmax], $
        [-9e9, 9e9, 9e9, -9e9], color=fclr
  endfor

  ; SIGMA
  oplot, state.z1arr, sig, psym=10, color=clr.red

  ;; H2
  oplot, (state.h2plotx/770.409)-1., state.h2ploty, color=clr.green, psym=10, $
    linestyle=1

  ;; Data
  oplot, state.z1arr, fx, color=clr.black, psym=10


  ;; Plot bottom
  
  widget_control, state.ldraw2_id, get_value=wind
  wset, wind
  plot, state.z2arr, fx, $
    xrange=[state.xymnx[0], state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos2, $
    charsize=1.2, psym=10, background=clr.white, color=clr.black, $
    xtitle='!17zabs', xmargin=[0,0], ymargin=[4,0], xstyle=1,$
    ystyle=1, /nodata

  ;; Regions
  gd = where(gzreg.flg NE 0 AND gzreg.zmin GT state.xymnx[0] AND $
             gzreg.zmax LT state.xymnx[2], ngd)
  for jj=0L,ngd-1 do begin
      id = gd[jj]
      if id EQ state.ngz then fclr = clr.orange else fclr=clr.yellow
      polyfill, [gzreg[id].zmin, gzreg[id].zmin, $
                 gzreg[id].zmax, gzreg[id].zmax], $
        [-9e9, 9e9, 9e9, -9e9], color=fclr
  endfor
  ;; SIGMA
  oplot, state.z2arr, sig, psym=10, color=clr.red
  
  ;; H2
  oplot, (state.h2plotx/780.324)-1., state.h2ploty, color=clr.green, psym=10, $
    linestyle=1

  ;; Data
  oplot, state.z2arr, fx, color=clr.black, psym=10
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro fuse_gzneviii_Metals, state
  
  common fuse_gzneviii_cmm

  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.mdraw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)
  
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)

  ny = ngd / 2 + (ngd MOD 2 EQ 1)

  !p.multi = [0,2,ny,0,1]

  for j=0L,ngd-1 do begin
      i = gdlin[j]
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]

      ;; Plot
      if (j NE ny-1) AND (j NE ngd-1 ) then begin
          spaces = replicate('!17 ',30)
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fx[pixmin:pixmax]/conti[pixmin:pixmax], xrange=state.vmnx, $
            yrange=state.velplt[i].ymnx, xtickn=spaces, xmargin=[9,3], $
            ymargin=[0,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fx[pixmin:pixmax]/conti[pixmin:pixmax], xrange=state.vmnx, $
            yrange=state.velplt[i].ymnx, xmargin=[9,3], ymargin=[3,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse

      ;; Labels
      xyouts, 0.07*(state.vmnx[1]-state.vmnx[0])+state.vmnx[0], $
        state.velplt[i].ymnx[0]+ $
        (state.velplt[i].ymnx[1]-state.velplt[i].ymnx[0])*0.05, $
        strtrim(state.velplt[i].name,2), $
        color=clr.black, charsize=1.5
      
      ;; Lines
      oplot, [0., 0.], state.velplt[i].ymnx, color=clr.blue, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.green, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  endfor

  !p.multi = [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro fuse_gzneviii_update, state
;  fuse_gzneviii_updfit, state
  fuse_gzneviii_setxym, state
  fuse_gzneviii_spectra, state
;  fuse_gzneviii_Metals, state
  return
end

;;;;;;;;;;;;;;;
;;; SET XY LIM
pro fuse_gzneviii_setxym, state

  state.xymnx[0] = state.zabs-state.dz
  state.xymnx[2] = state.zabs+state.dz

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro fuse_gzneviii_next, state
  common fuse_gzneviii_cmm

  state.zabs = state.zabs+state.dz
  fuse_gzneviii_setxym, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro fuse_gzneviii_prev, state
  common fuse_gzneviii_cmm

  state.zabs = state.zabs - state.dz
  fuse_gzneviii_setxym, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro fuse_gzneviii_updinfo, state
  common fuse_gzneviii_cmm

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs

  return

end

;;;;;;;;;;;;;;;
;;; Remove region
pro fuse_gzneviii_rmreg, state
  common fuse_gzneviii_cmm

  zv = xgetx_plt(state, /strct) 
  a = where(gzreg.zmin LT zv and gzreg.zmax GT zv, na)
  if na NE 0 then gzreg[a].flg = 0

  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro fuse_gzneviii_svgz, state
  common fuse_gzneviii_cmm

  ;; zabs
  mwrfits, gzreg, state.outfil, /create

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

pro fuse_gzneviii, datfil, XSIZE=xsize, YSIZE=ysize, $
                   OUTFIL=outfil, GZFIL=gzfil

  common fuse_gzneviii_cmm
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'fuse_gznevii, datfil, OUTFIL=, XSIZE=, YSIZE=, GZFIL=, [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( YSIZE ) then    ysize = round(0.8*ssz[1])
  if not keyword_set( XSIZE ) then xsize = round(0.8*ssz[0])
  if not keyword_set( DELTAZ ) then deltaz = 0.0025
  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( OUTFIL ) then begin
      if keyword_set( GZFIL ) then outfil= gzfil else $
        outfil = 'fin_neviii.fits'
  endif

; Initialize the common blcok
  fuse_gzneviii_icmmn, datfil, gzfil

  if not keyword_set( STRTZ) then strtz = (wv[0]+1.)/770.409 - 1.

  ;; H2
  h2lines = fuse_h2lin()

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }

; STATE
  a = where(gzreg.flg NE 0,na)
  if na NE 0 then ngz = max(a) else ngz = 0L

  state = {             $
            nqal: 0L, $
            zabs: strtz, $
            dz: deltaz, $
            z1arr: (wv/770.409)-1., $
            z2arr: (wv/780.324)-1., $
            ngz: ngz, $
            outfil: outfil, $
            patt: bytarr(10,10), $
            h2lines: h2lines, $
            h2plotx: fltarr(2*n_elements(h2lines)), $
            h2ploty: fltarr(2*n_elements(h2lines)), $
            ntrans: 0L, $       ; PLOTTING LINES
            vmnx: [-800., 800.], $
            nplt: 0, $
            dla_lin: tmp2, $
            dla_conti: 0., $
            all_velo: dblarr(5000, 300), $  
            all_pmnx: lonarr(3, 300), $  
            velplt: replicate(tmp, 300), $
            xpos: 0.0, $
            ypos: 0.0, $
            ipress: 0L, $
            pos: [0.1,0.0,0.95,1.00], $ ; Plotting
            pos2: [0.1,0.1,0.95,1.00], $ ; Plotting
            flg_zoom: 0, $
            psfile: 0, $
            help: strarr(50), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            xcurs: 0., $
            ycurs: 0., $
            size: lonarr(2), $
            base_id: 0L, $      ; Widgets
            ldraw_id: 0L, $    ; Spectra
            ldraw2_id: 0L, $    
            ltext_id: 0L, $    ; 
            ldrawbase_id: 0L, $
            fxval_id: 0L, $
            iwvval_id: 0L, $
            mdraw_id: 0L, $       ; Spec Window
            mdrawbase_id: 0L, $
            swvval_id: 0L, $
            zabs_id: 0L, $
            xmax_id: 0L, $
            name_id: 0L, $
            nspec_id: 0L, $
            pmin_id: 0L, $
            pmax_id: 0L, $
            lines_id: 0L, $
            lhs_id: 0L, $
            rhs_id: 0L, $
            info_id: 0L, $
            quality_id: 0L, $
            scr1_id: 0L, $
            scr2_id: 0L, $
            hits_id: 0L, $
            NHI_id: 0L, $
            NHIb_id: 0L, $
            mtl_id: 0L, $
            stat_id: 0L, $
            ra_id: 0L, $
            dec_id: 0L, $
            mag_id: 0L, $
            help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES

  state.patt[*] = 255
  state.patt[5,5] = 0
  state.xymnx[1] = -0.1
  state.xymnx[3] = 1.2
  state.svxymnx = state.xymnx

;; H2
  for ii=0L,n_elements(h2lines)-1 do begin
      state.h2plotx[ii*2:ii*2+1] = h2lines[ii].wrest
      state.h2ploty[ii*2] = -9.e9
      state.h2ploty[ii*2+1] = 9.e9
  endfor


;    WIDGET
  base = WIDGET_BASE( title = 'fuse_gzneviii: Find NeVIII', /row, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.lhs_id = WIDGET_BASE( state.base_id, /column, $
                              /base_align_center,/align_center, $
                              uvalue='RHS_BASE', frame=2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Lya DRAW
  state.ldrawbase_id = $
    WIDGET_BASE( state.lhs_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='LDRAW_BASE')

  state.ldraw_id = widget_draw(state.ldrawbase_id, xsize=round(xsize*3./5), $
                              ysize=round(ysize/2.), $
                              /button_events, /motion_events, uvalue='LDRAW')

  state.size[0]=round(xsize*3/5.)
  state.size[1]=2*ysize/3.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Lya TEXT
  state.ltext_id = widget_text(state.ldrawbase_id, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'LTEXT', $
                              value = '')

;;;;;; Info window ;;;;;;;;;;;
  ldraw2base_id = $
    WIDGET_BASE( state.lhs_id, /row, /base_align_center,/align_center, $
               uvalue='LDRAW2_BASE')
  state.ldraw2_id = widget_draw(ldraw2base_id, xsize=round(xsize*3./5), $
                              ysize=round(ysize/2.),  uvalue='LDRAW2')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Metals DRAW
  state.mdrawbase_id = $
    WIDGET_BASE( state.base_id, /column, /base_align_center,/align_center, $
               uvalue='SDRAW_BASE', frame=2)

  state.mdraw_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./5), $
                              ysize=ysize*4/5, /frame, retain=2, $
                              uvalue='SDRAW')

;      BUTTONS
  butbase = widget_base(state.mdrawbase_id, /row, /align_center, frame=2)
  state.zabs_id = cw_field(butbase, title='zabs: ', value=0., xsize=3)
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'LMB - Truncate/Extend trace', $ 
             'RMB - Contrast/Brightness', $
             'CMB/CMB - Zoom' $ 
             ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ; Load data
;  fuse_gzneviii_setup, state

  ; PLOT
  fuse_gzneviii_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'fuse_gzneviii', base
  delvarx, fx, wv, npix, sig

  return
end
	
