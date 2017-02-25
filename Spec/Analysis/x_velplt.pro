;+ 
; NAME:
; x_velplt
;    Version 1.1
;
; PURPOSE:
;  GUI used to display velocity profiles for absorption line system.
;  The user can control a number of options for the plot, output ps
;  files and PG plot input files, etc.
;
; CALLING SEQUENCE:
; x_velplt, flux_fil, zin, VMNX=vmnx, INFLG=inflg, $
;             TITLE=title, NPLT=nplt, NRM=nrm, XSIZE=xsize, YSIZE=ysize, $
;             SIG_FIL=, /ESIDLA, SVSTATE=, /LLS
;
; INPUTS:
;  flux_fil -- Input FITS file for the flux (header includes
;              wavelength info)
;  zin -- Redshift of absorption lines
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  VMNX= -- Velocity range to plot [default [-300, 300]]
;  NPLT= -- Number of plots to show [default: 15L]
; TITLE= -- Title to give GUI
;  XSIZE      - Size of gui in screen x-pixels [default = 700]
;  YSIZE      - Size of gui in screen y-pixels [default = ssz -200]
;  /LLS -- Use LLS line list
;  /SUBLLS -- Use smaller LLS line list
;  /ESIDLA -- Use ESI DLA line list
;  /CII - subset of the ESIDLA line list for measuring CII
;  /SHRT - subset of ESIDLA line list for a shorter list (and hiz)
;  SIG_FIL= -- Sigma file
;  INFLG=  -- Flag for flux_fil data (see x_readspec)
;  SVSTATE= -- Input IDL file that saved the state of the GUI
;              previously
;  WAVE -- wavelength array (flux_fil should be array)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_velplt, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   29-Oct-2002 Written by JXP
;   13-Jan-2006 Add wave keyword and enable passing of arrays, KLC
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro x_velplt_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  flg_allplt = 0

  case uval of
      'SPPLOT': x_specplot, state.fx, state.sig, wave=state.wave, inflg=4, $
        /block, zin=state.zabs, /lls
      'VELDROP': begin
          gd = where(state.velplt.flg EQ 1)
          state.curlin = gd[ev.index]
          widget_control, state.ymin_id, $
            set_value=state.velplt[state.curlin].ymnx[0]
          widget_control, state.ymax_id, $
            set_value=state.velplt[state.curlin].ymnx[1]
      end
      'YMIN' : begin
          widget_control, state.ymin_id, get_value=tmp
          state.velplt[state.curlin].ymnx[0] = tmp
          flg_allplt = 1
      end
      'YMAX' : begin
          widget_control, state.ymax_id, get_value=tmp
          state.velplt[state.curlin].ymnx[1] = tmp
          flg_allplt = 1
      end
      'NCLM' : begin  ;; Ncolm
          widget_control, state.nclm_id, get_value=tmp
          state.nx = tmp
          flg_allplt = 1
      end
      'NHPLT' : begin  ;; Ncolm
          widget_control, state.nhplt_id, get_value=tmp
          state.hplt = tmp
          flg_allplt = 1
      end
      'LLIST': begin
          state.llist = ev.index
          x_velplt_llist, state
      end
      'LIST': begin
          state.velplt.flg = 0
          indx = widget_info(state.list_id, /list_select)
          state.velplt[indx].flg = 1
          ;; VELDROP
          widget_control, state.veldroplist_id, $
            set_value=state.velplt[indx].name
          ;; NPG
          state.nplt = n_elements(indx)
          x_velplt_setnpg, state
          flg_allplt = 1
      end
      'NEXT' : begin
          state.curpg = state.curpg + 1
          if state.curpg GT state.npg-1 then state.curpg = 0L
          widget_control, state.pg_id, set_value=state.curpg
          flg_allplt = 1
      end
      'PREV' : begin
          state.curpg = state.curpg - 1
          if state.curpg LT 0 then state.curpg = state.npg-1
          widget_control, state.pg_id, set_value=state.curpg
          flg_allplt = 1
      end
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          state.zabs = tmp
          x_velplt_setvelo, state
          flg_allplt = 1
      end
      'VMIN' : begin
          widget_control, state.vmin_id, get_value=tmp
          state.vmnx[0] = tmp
          x_velplt_setvelo, state
          flg_allplt = 1
      end
      'VMAX' : begin
          widget_control, state.vmax_id, get_value=tmp
          state.vmnx[1] = tmp
          x_velplt_setvelo, state
          flg_allplt = 1
      end
      'PRINT' : x_velplt_print, state
      'SVPGPLOT' : begin
          x_velplt_svpgplot, state
          print, 'x_velplt: PGPLOT info written to vel_plt.par'
      end
      'SVSTATE' : begin
          save, state, filename='xvp_state.idl', /compress
          print, 'x_velplt: State written to xvp_state.idl'
      end
      'DONE' : begin
          widget_control, state.base_id, /destroy
          return
      end
      else :
  endcase

  ;; PLOT
  if flg_allplt EQ 1 then x_velplt_allplt, state

  ;; RETURN
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SET VELO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_velplt_setvelo, state

  ;; Just do the plots
  gd = where(state.velplt.flg EQ 1)

  ;; Call x_allvelo
  state.all_velo[*,gd] = x_allvelo(state.wave, state.zabs, $
                                                 state.velplt[gd].wrest,$
                                                 state.vmnx, $
                                                 all_pmnx=all_pmnx, NPIX=5000L)

  state.all_pmnx[*,gd] = all_pmnx

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;

pro x_velplt_allplt, state
  
  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.alldraw_id, get_value=wind
      wset, wind
  endif

  ;; GDLIN
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)
  if ngd NE state.nplt then stop

  ; PLOT
  clr = getcolor(/load)

  !P.MULTI= [0,state.nx,state.hplt,0,1]

  pmin = state.curpg*(state.nx*state.hplt)
  if state.curpg EQ state.npg - 1 then pmax = state.nplt-1 $
  else pmax = (state.curpg+1)*(state.nx*state.hplt)-1
;  pmax = state.nx*state.hplt - 1

  ;; LOOP
  for j=pmin, pmax do begin
      i = gdlin[j]
      ;; Checks
      if i EQ state.ntrans then break
      if state.all_pmnx[2,i] EQ 0 then continue
      if state.velplt[i].flg MOD 2 NE 1 then continue

      ;; Get pxmin, pxmax + velo array
;      x_pixminmax, state.wave, lines[gdlin[i]].wave, state.zabs, $
;        state.svxymnx[0], state.svxymnx[2], PIXMIN=pixmin, PIXMAX=pixmax, $
;        VELO=velo
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]
      
      ;; Plot
      if (j NE pmax) AND ((j-pmin) NE state.hplt-1 ) then begin
          spaces = replicate('!17 ',30)
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            state.fx[pixmin:pixmax], xrange=state.vmnx, $
            yrange=state.velplt[i].ymnx, xtickn=spaces, xmargin=[9,3], $
            ymargin=[0,0], $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            state.fx[pixmin:pixmax], xrange=state.vmnx, $
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
  
  !P.MULTI= [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro x_velplt_print, state

  widget_control, /hourglass   
; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
;  !p.charthick = 1

  device, decompose=0
  ; Get file name
  psfil = 'xvp.ps'
  ps_open, file=psfil, font=1, /color, /portrait
  state.psfile = 1
  for qq=0L,state.npg-1 do begin
      state.curpg = qq
      x_velplt_allplt, state
  endfor
  ps_close, /noprint, /noid
;  spawn, 'gzip -f '+psfil
  device, decomposed=svdecomp
  state.psfile = 0
;  !p.thick = 1
;  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;
;  Set Line list
;;;;;;;;;;;;;;;;;;;;

pro x_velplt_llist, state

  case state.llist of
      0: begin                  ; All DLA
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qal.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end
      1: begin                  ; ESI DLA
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/esi_dla.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end
      3: begin                  ; LLS
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end
      4: begin                  ; H2
          h2list = fuse_h2lin()
          tmp = { lliststrct }
          lines = replicate(tmp, n_elements(h2list))
          lines.wave = h2list.wrest
          lines.name = h2list.label
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = $
            lines.name+' '+strtrim(round(lines.wave),2)
          delvarx, lines, h2list
      end
      5: begin                  ; Subset of LLS.LST
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls_sub.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end 
      6: begin                  ; Subset of esi_dla_sub for CII
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/esi_dla_sub.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end 
      7: begin                  ; Subset of esi_dla_sub for HIZ (short)
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/dla_shrt.lst'
          lines = x_setllst(llist, 0)
          state.ntrans = n_elements(lines)
          state.velplt[0:state.ntrans-1].wrest = lines.wave
          state.velplt[0:state.ntrans-1].name = lines.name
          delvarx, lines
      end
      else: stop
  endcase

  ;; Default is Normalized data
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.09]

  
  ;; Set flg
  a = where(state.velplt.wrest*(state.zabs+1.) GT min(state.wave) AND $
            state.velplt.wrest*(state.zabs+1.) LT max(state.wave[0:state.npix-1]), $
            na )
  if na EQ 0 then stop
  state.nplt = na
  state.curlin = a[0]
  
  state.velplt[a].flg = 1

  ;; Un-normalized
  if state.flg_norm EQ 2 then begin
     tmp_velo = x_allvelo(state.wave, state.zabs, $
                          state.velplt[0:state.ntrans-1].wrest,$
                          [-400, 400], all_pmnx=tmp_pmnx, NPIX=5000L)
     for jj=0L,na-1 do begin
        ii = a[jj]
        all_fx = state.fx[tmp_pmnx[0,ii]:tmp_pmnx[1,ii]]
        npix = n_elements(all_fx)
        srt = sort(all_fx)
        ymx = 1.1 * all_fx[srt[round(0.9*npix)]]
        ymn = -0.1 * all_fx[srt[round(0.9*npix)]]
        state.velplt[ii].ymnx = [ymn, ymx]
     endfor
  endif

end

;;;;;;;;;;;;;;;;;;;;
;  Set npg
;;;;;;;;;;;;;;;;;;;;

pro x_velplt_setnpg, state

  ;; Value
  state.npg =  state.nplt/(state.nx*state.hplt) + $
    (state.nplt MOD (state.nx*state.hplt) NE 0)
  ;; 
  widget_control, state.npg_id, set_value=state.npg-1

  ;; Curpg
  state.curpg = state.curpg < (state.npg-1)

end

;;;;;;;;;;;;;;;;;;;;
;  PGPLOT -- create vel_plt.par
;;;;;;;;;;;;;;;;;;;;

pro x_velplt_svpgplot, state

  ;; Open file
  close, 11
  openw, 11, 'vel_plt.par'
  printf, 11, '0.5, 1.1'
  ;; File
  if size(state.fil, /type) EQ 7 then printf, 11, state.fil else $
    printf, 11, 'filename'
  ;; Label
  printf, 11, 'Label!'

  ;; Nplt
  printf, 11, strtrim(state.nplt,2)

  ;; z
  printf, 11, strtrim(state.zabs,2)

  ;; vmin,vmax
  printf, 11, FORMAT='(f9.2,1a,f9.2)', state.vmnx[0], ',', state.vmnx[1]

  ;; LOOP
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)
  if ngd NE state.nplt then stop

  for q=0L,state.nplt-1 do begin
      gd = gdlin[q]
      ;; wrest
      printf, 11, strtrim(state.velplt[gd].wrest,2)
      ;; ymin, ymax
      printf, 11, FORMAT='(f5.2,1a,f5.2)', state.velplt[gd].ymnx[0], ',', $
        state.velplt[gd].ymnx[1]
      ;; Label
      printf, 11, '0.19, -0.8'
      printf, 11, '0'
  endfor
  close, 11
      
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_velplt, flux_fil, zin, VMNX=vmnx, INFLG=inflg, $
              TITLE=title, NPLT=nplt, XSIZE=xsize, YSIZE=ysize, $
              SIG_FIL=ysin, ESIDLA=esidla, SVSTATE=svstate, LLS=lls, $
              WAVE=wave, AIR=air, H2=h2, SUBLLS=sublls, CII=CII, SHRT=SHRT, $
              UN_NORM=un_norm

;
  if  N_params() LT 2  and not keyword_set(SVSTATE) then begin 
    print,'Syntax - ' + $
      'x_velplt, spec, zin,  VMNX=, INFLG=, TITLE=, SVSTATE='
    print, '      /ESIDLA, /LLS, /SUBLLS, NPLT=, /AIR, /H2, /CII, /SHRT, /UN_NORM [v1.2]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( ZIN ) then zin = 3.0
  if not keyword_set( XSIZE ) then xsize = 700L
  device, get_screen_size=ssz
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( VMNX ) then vmnx = [-300., 300.]
  if not keyword_set( NPLT ) then nplt = 15L
  if not keyword_set( FWHM ) then fwhm = 4.
  if not keyword_set( LSTFONT ) then lstfont = '6x10'


; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0
  if size(flux_fil,/type) eq 7 then $
    ydat = x_readspec(flux_fil, INFLG=inflg, head=head, NPIX=npix, $
                    WAV=xdat, FIL_SIG=ysin, SIG=ysig) $
  else begin 
      ydat = flux_fil
      flux_fil = 'array'
      npix = n_elements(ydat)
  endelse 
  
  tmp = { velpltstrct }

  ;; WAVE
  if keyword_set(WAVE) then xdat = wave
      
; STATE
      
  nstate = {  $
            fil: flux_fil, $
            hplt: 8L, $
            nx: 2L, $
            npg: 0L, $
            curpg:0L, $
            nplt: 0L, $
            npix: npix, $
            fx: ydat, $
            wave: xdat, $
            sig: fltarr(npix), $
            flg_sig: 0, $
            flg_zoom: 0, $
           flg_norm: 1, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            xpmnx: lonarr(2), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            old_xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            curlin: 0L, $
            zabs: zin, $
            vmnx: vmnx, $
            llist: 0L, 	$
            ntrans: 0L, $       ; PLOTTING LINES
            all_velo: dblarr(5000, 1000), $  
            all_pmnx: lonarr(3, 1000), $  
            velplt: replicate(tmp, 1000), $
            psfile: 0, $
            size: lonarr(2), $
            xpos: 0., $
            ypos: 0., $
            xcurs: 0., $
            ycurs: 0., $
            base_id: 0L, $      ; Widgets
            alldraw_id: 0L, $
            alldrawbase_id: 0L, $
            list_id: 0L, $
            vmin_id: 0L, $
            vmax_id: 0L, $
            droplist_id: 0L, $
            npg_id: 0L, $
            pg_id: 0L, $
            zabs_id: 0L, $
            ymin_id: 0L, $
            ymax_id: 0L, $
            veldroplist_id: 0L, $
            bval_id: 0L, $
            nclm_id: 0L, $
            nhplt_id: 0L, $
            help_text_id: 0L $
          }

  ;; SVSTATE
  if keyword_set(SVSTATE) then begin
      a = findfile(svstate,count=na)
      if na EQ 0 then begin
          print, 'x_velplt: SVSTATE ', svstate, ' does not exist! Returning..'
          return
      endif
      restore, svstate
      copy_struct, state, nstate
      state = temporary(nstate)
  endif else begin

      ;; SORT
      srt = sort(nstate.wave)
      nstate.wave = nstate.wave[srt]
      nstate.fx = nstate.fx[srt]
      nstate.sig = nstate.sig[srt]

      if keyword_set(AIR) then begin
          tmp = nstate.wave
          airtovac, tmp
          nstate.wave = tmp
      endif

      ;;
      state = temporary(nstate)
      if keyword_set(UN_NORM) then state.flg_norm = 2
; NORM
      if keyword_set( YSIG ) then state.sig = temporary(ysig)

; LINELIST for VELPLT

      resolve_routine, 'x_specplot', /NO_RECOMPILE
      ;; Line list
      state.llist = 0L
      if keyword_set(ESIDLA) then state.llist = 1
      if keyword_set(LLS) then state.llist = 3
      if keyword_set(H2) then state.llist = 4
      if keyword_set(SUBLLS) then state.llist = 5
      if keyword_set(CII) then state.llist = 6
      if keyword_set(SHRT) then state.llist = 7
      x_velplt_llist, state
  endelse
  
;    WIDGET
  base = WIDGET_BASE( title = 'x_velplt: LLS fit', /row, $
                      UNAME='BASE', xoffset=100)
  state.base_id = base
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ALL ;;;;

  ;; DRAW
  state.alldrawbase_id = widget_base(base, /column, /base_align_left, $
                                     uvalue='ALLDRAW_BASE', frame=2, $
                                     /tracking_events)
  
  state.alldraw_id = widget_draw(state.alldrawbase_id, xsize=xsize, $
                                 ysize=ysize, /frame, retain=2, $
                                 uvalue='ALLDRAW')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TRANS LIST
  list = state.velplt[0:state.ntrans-1].name
  
  ysz = 30 < n_elements(list)
  state.list_id = widget_list(base, value=list, xsize=20L, ysize=ysz, $
                              FONT=lstfont, uvalue = 'LIST', /MULTIPLE)
  ;; Select them all
  gd = where(state.velplt.flg EQ 1)
  widget_control, state.list_id, set_list_select=gd

;;;;;;;;;;;FIDDLE ;;;;
  fiddle_base = WIDGET_BASE( state.base_id, /column, /frame, /base_align_center,$
                         /align_center)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; VELPLT DROPLIST
  veldrop = WIDGET_BASE( fiddle_base, /column, /frame, /base_align_center,$
                         /align_center)
  vellst = state.velplt[gd].name
  state.veldroplist_id = widget_droplist(veldrop, $
                                      frame = 1, $
                                      title = 'Transition:', $
                                      uvalue = 'VELDROP', $
                                      value = vellst)

  state.ymin_id = cw_field(veldrop, title='ymin', /return_events, $
                           value=state.velplt[gd[0]].ymnx[0], /float, $
                           /column, xsize=10, uvalue='YMIN')
  state.ymax_id = cw_field(veldrop, title='ymax', /return_events, $
                           value=state.velplt[gd[0]].ymnx[1], /float, $
                           /column, xsize=10, uvalue='YMAX')
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Toolbar
  alltool = WIDGET_BASE( fiddle_base, /column, /frame, /base_align_center,$
                         /align_center)
  
  ;; Version + Name + trace
  labelbase = widget_base(alltool, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_velplt', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
  
  ;; ZABS
  state.zabs_id = cw_field(alltool, title='zabs', /return_events, $
                           value=state.zabs, /float, $
                           /column, xsize=10, uvalue='ZABS')
  
  ;; LINE LIST
  linelist = ['All DLA', 'ESI DLA', 'Low Z DLA', 'ALL LLS', 'H2', 'CII']
  state.droplist_id = widget_droplist(alltool, $
                                      frame = 1, $
                                      title = 'Linelist:', $
                                      uvalue = 'LLIST', $
                                      value = linelist)
  widget_control, state.droplist_id, set_droplist_select=state.llist
  
;;;;;;;;;;;;;;;;;;
; Pages
  pgpair = WIDGET_BASE(alltool, /row, /base_align_left,$
                        /align_center, /frame)

  npgpair = WIDGET_BASE(pgpair, /row, /base_align_left,$
                        /align_center)
  state.npg_id = cw_field(pgpair, value=state.npg, $
                          /long, /return_events, xsize=3,$
                          title='Npg:')
  state.pg_id = cw_field(pgpair, value=state.curpg, $
                         /long, /return_events, xsize=3,$
                             title='Page')
  npltpair = WIDGET_BASE(alltool, /row, /base_align_left,$
                        /align_center)
  state.nclm_id = cw_field(npltpair, value=state.nx, $
                         /long, /return_events, xsize=3,$
                             title='Nx', uvalue='NCLM')
  state.nhplt_id = cw_field(npltpair, value=state.hplt, $
                         /long, /return_events, xsize=3,$
                             title='Ny', uvalue='NHPLT')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Vmin, vmax
  vmnx_base = WIDGET_BASE(alltool, /column, /base_align_left,$
                          /align_center)
  state.vmin_id = cw_field(vmnx_base, value=state.vmnx[0], $
                           /float, /return_events, xsize=9,$
                           title='vmin:', UVALUE='VMIN')
  state.vmax_id = cw_field(vmnx_base, value=state.vmnx[1], $
                           /float, /return_events, xsize=9,$
                           title='vmax:', UVALUE='VMAX')
  
  ;;      BUTTONS
  butbase = widget_base(alltool, /row, /align_center)
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  butbase2 = widget_base(alltool, /row, /align_center)
  print = WIDGET_BUTTON(butbase2, value='PRINT',uvalue='PRINT')
  plot = WIDGET_BUTTON(butbase2, value='SPPLOT',uvalue='SPPLOT')
  dum_svstate = WIDGET_BUTTON(butbase2, value='SVSTATE',uvalue='SVSTATE')
  dum_svpgplot = WIDGET_BUTTON(butbase2, value='SVPGPLOT',uvalue='SVPGPLOT')
  done = WIDGET_BUTTON(alltool, value='DONE',uvalue='DONE')

; Realize
  WIDGET_CONTROL, base, /realize


; INIT


  if not keyword_set(SVSTATE) then begin
      ;; npg
      x_velplt_setnpg, state
      
      ;; Set VELO
      x_velplt_setvelo, state
  endif
      
  ;; PLOT
  x_velplt_allplt, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_velplt', base

  !P.MULTI= [0,1,1]

; Reset
  if flux_fil eq 'array' then flux_fil = ydat 
  return
end

