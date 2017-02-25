;+ 
; NAME:
; sdss_finchk
;    Version 1.1
;
; PURPOSE:
;   GUI used to make a final check of metal-strong candidates
;
; CALLING SEQUENCE:
;  sdss_finchk, dlafil, con_dir, MAG=mag, XSIZE=, YSIZE=, 
;                 OUTFIL=outfil, RA=ra
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  MAG=  -- Largest magnitude to include in the check
;  RA=   -- 2-element array giving constraint on RA
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_finchk, 'dla_dr1.fits', 'ABSLIN/'
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

pro sdss_finchk_icmmn, dlafil, STRNG=strng, MAG=mag, RA=ra, DEC=dec

  common sdss_finchk_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    fit, $
    dlastr, $
    gddla

  ;; QALSTR
;  qalstr = xmrdfits(qalfil, 1, /silent)

  if not keyword_set( STRNG ) then strng = 4
  if not keyword_set( MAG ) then mag = 99.
;  if not keyword_set( DRA ) then dra = 2.

  ;; DLA str
  if x_chkfil(dlafil+'*', /silent) EQ 1 then $
    dlastr = xmrdfits(dlafil, 1, /silent) $
  else stop


  gddla = where(dlastr.flg_mtl GE strng AND dlastr.rmag LT mag)

  ;; RA
  if keyword_set(RA) then begin
 ;     tmp =  where(abs(dlastr[gddla].ra/15.  - RA) LT dra, ntmp)
      tmp=where(dlastr[gddla].ra/15. ge RA[0]/15. and dlastr[gddla].ra/15. le RA[1]/15., ntmp)
    if ntmp EQ 0 then stop
      gddla = gddla[tmp]
  endif

  ;; Dec
  if keyword_set( DEC ) then begin
      tmp = where(dlastr[gddla].dec LT DEC, ntmp)
      if ntmp EQ 0 then stop
      gddla = gddla[tmp]
  endif

  ;; Sort on Mag
  srt = sort(dlastr[gddla].rmag)
  gddla = gddla[srt]

  print, 'ntarg = ', n_elements(gddla)
  return
end
  

;;;;
; Events
;;;;

pro sdss_finchk_event, ev

  common sdss_finchk_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'PRINT': sdss_finchk_print, state
      'SPLT': x_specplot, fx, sig, wave=wv, zin=state.zabs, /qal, /block, $
        inflg=4
      'SAVE': sdss_finchk_svdla, state
      'NEXT': begin
          sdss_finchk_next, state
          sdss_finchk_setup, state
          sdss_finchk_update, state
      end
      'PREV': begin
          sdss_finchk_prev, state
          sdss_finchk_setup, state
          sdss_finchk_update, state
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
                      1 : state.dla_lin.N = state.dla_lin.N + 0.05
                      2 : begin
                          state.zabs = state.xpos / 1215.6701 - 1.
                          state.dla_lin.zabs = state.zabs
                          state.dla_conti = state.ypos
                          sdss_finchk_updfit, state
                          widget_control, state.zabs_id, $
                            set_value=strtrim(state.zabs,2)
                      end
                      4 : state.dla_lin.N = state.dla_lin.N - 0.05
                  endcase
                  sdss_finchk_updinfo, state
                  sdss_finchk_updfit, state
                  sdss_finchk_Lya, state
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
              else: 
          endcase
          sdss_finchk_Lya, state
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
          sdss_finchk_svdla, state
          widget_control, ev.top, /destroy
          return
      end
      'DNSV' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  sdss_finchk_upddstr, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Lya Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_finchk_Lya, state
  
  common sdss_finchk_cmm

  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.ldraw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)
  
  ; Plot
  plot, wv, fx, $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], pos=state.pos, $
    charsize=1.2, psym=10, background=clr.white, color=clr.black, $
    xtitle='!17Wavelength', xmargin=[0,0], ymargin=[0,0], xstyle=1,$
    ystyle=1 

  ;; Fit
  oplot, wv, fit, color=clr.blue

  ; Plot Lines as required
;  if flg_lines EQ 1 then x_specplot_PltLines, state, /IMG

  ;; center
  oplot, 1215.6701*[state.zabs+1.,state.zabs+1], [-1e5, 1e5], color=clr.green,$ 
    linestyle=2

  ; SIGMA
  oplot, wv, sig, psym=10, color=clr.red

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_finchk_Metals, state
  
  common sdss_finchk_cmm

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
pro sdss_finchk_update, state
  sdss_finchk_updfit, state
  sdss_finchk_updinfo, state
  sdss_finchk_Lya, state
  sdss_finchk_Metals, state
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_finchk_updfit, state

  common sdss_finchk_cmm

  fit[*] = 1.
  fit = x_allvoigt(wv, state.dla_lin, SIGMA=3.)
;  fit = x_voigt(wv, state.dla_lin, FWHM=3.)


  fit = fit*state.dla_conti

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro sdss_finchk_setup, state
  common sdss_finchk_cmm

  ;; Kludge for continuum
  id = long(strsplit(dlastr[state.curdla].qso_name,'-',/extract))
  nid = n_elements(id)
  sdss_objinf, [id[0],id[nid-1]], /dr1, zem=zem

  ;; Continuum
  if keyword_set( ZEM ) then cdir = getenv('SDSSPATH')+'/DR1_QSO/ABSLIN/' $
  else cdir = state.con_dir

  ;; Read data
  parse_sdss, getenv('SDSSPATH')+strtrim(dlastr[state.curdla].sdss_obs[0],2), $
    fx, wv, conti, SIG=sig, NPIX=npix, CDIR=cdir

  if  n_elements(conti) LT n_elements(fx) then begin
      tmpc = fltarr(n_elements(fx))
      tmpc[0:n_elements(conti)-1] = conti
      conti = temporary(tmpc)
  endif

  ;; Set zabs
  state.zabs = dlastr[state.curdla].zabs
  state.zqso = dlastr[state.curdla].z_qso

  ;; xymnx
  state.xymnx[0] = 1215.6701*(1.+state.zabs) - 200.
  state.xymnx[2] = 1215.6701*(1.+state.zabs) + 200.
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  if ngd GT 1 then begin
      srt = sort(fx[gd])
      ymd = fx[gd[srt[round(0.9*ngd)]]]
  endif else ymd = 0.
  state.xymnx[1] = -1.
  state.xymnx[3] = ymd*1.5
  state.svxymnx = state.xymnx

  ;; DLA
  state.dla_lin.N = 20.3
  state.dla_lin.b = 30.
  state.dla_lin.zabs = state.zabs

  ;; Fit
  state.dla_conti = ymd
  fit = replicate(state.dla_conti, npix)

  ;; Set flg
  gd = where(state.velplt.wrest*(state.zabs+1) GT min(wv) AND $
             state.velplt.wrest*(state.zabs+1) GT (state.zqso+1.)*1215.6701 AND $
             state.velplt.wrest*(state.zabs+1) LT max(wv), na)
  if na EQ 0 then stop
  state.velplt[*].flg = 0
  state.nplt = na
  state.velplt[gd].flg = 1

  ;; Just do the plots
  state.all_velo[*,gd] = x_allvelo(wv, state.zabs, $
                                   state.velplt[gd].wrest,$
                                   state.vmnx, $
                                   all_pmnx=all_pmnx, NPIX=5000L)
  state.all_pmnx[*,gd] = all_pmnx

  ;; Set DLA structure
;  sdss_finchk_setdla, state

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set Lines
pro sdss_finchk_llist, state

  llist = getenv('XIDL_DIR')+'/SDSS/sdss_dla.lst'
  lines = x_setllst(llist, 0)
  state.ntrans = n_elements(lines)
  state.velplt[0:state.ntrans-1].wrest = lines.wave
  state.velplt[0:state.ntrans-1].name = lines.name
  delvarx, lines
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.39]

  weak = where(abs(state.velplt.wrest - 1808.0130d) LT 0.01 OR $
               abs(state.velplt.wrest - 1611.2005d) LT 0.01 OR $
               abs(state.velplt.wrest - 2026.136d) LT 0.01 OR $
               abs(state.velplt.wrest - 2260.7805d) LT 0.01 )
  
  state.velplt[weak].ymnx = [0.7, 1.1]
               
  

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_finchk_next, state
  common sdss_finchk_cmm

  idx = where(state.curdla EQ gddla)
  idx = idx + 1
  ndla = n_elements(gddla)
  if idx[0] EQ ndla then idx[0] = 0L else idx = idx[0]
  state.curdla = gddla[idx]

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro sdss_finchk_prev, state
  common sdss_finchk_cmm

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_finchk_updinfo, state
  common sdss_finchk_cmm

  ;; Namej
  widget_control, state.name_id, $
    set_value=strtrim(dlastr[state.curdla].qso_name,2)

  ;; RA, DEC
  widget_control, state.mag_id, set_value=dlastr[state.curdla].rmag
  widget_control, state.ra_id, set_value=dlastr[state.curdla].ra/15.
  widget_control, state.dec_id, set_value=dlastr[state.curdla].dec

  ;; Quality
  widget_control, state.quality_id, set_value=dlastr[state.curdla].quality

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs

  ;; Metals and NHI
  widget_control, state.mtl_id, set_value=dlastr[state.curdla].flg_mtl
  widget_control, state.NHI_id, set_value=dlastr[state.curdla].NHI
  widget_control, state.NHIb_id, set_value=dlastr[state.curdla].flg_NHI

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Update DLA structure
pro sdss_finchk_upddstr, state

  common sdss_finchk_cmm

  ;; zabs, conti
  dlastr[state.curdla].zabs = state.zabs
  dlastr[state.curdla].conti = state.dla_conti

  ;; Metals
  widget_control, state.mtl_id, get_value=v_mtl
  dlastr[state.curdla].flg_mtl = v_mtl

  ;; NHI
  widget_control, state.NHIb_id, get_value=v_NHI
  dlastr[state.curdla].flg_NHI = v_NHI
  dlastr[state.curdla].NHI = state.dla_lin.N

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Print
pro sdss_finchk_print, state

  common sdss_finchk_cmm

  x_radec, ras, decs, dlastr[state.curdla].ra, dlastr[state.curdla].dec, /FLIP
  if (dlastr[state.curdla].zabs+1)*1215.67 GT 3850. then $
    NHI = strtrim(dlastr[state.curdla].NHI,2) else NHI = '-99.99'
  ;; Create line
  lin = strtrim(dlastr[state.curdla].qso_name,2)+' '+$
    strtrim(ras,2)+' '+strtrim(decs,2)+' '+$
    strtrim(dlastr[state.curdla].Rmag,2)+' '+$
    strtrim(dlastr[state.curdla].z_qso,2)+' '+$
    strtrim(dlastr[state.curdla].zabs,2)+' '+NHI

  openw, 55, state.outfil, /append
  printf, state.filnum, lin
  close, 55

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
;sdss_finchk, 'dr3_qal_chkd.fits', 'ABSLIN/', MAG=19., DEC=20., RA=[0., 8*15]

pro sdss_finchk, dlafil, con_dir, MAG=mag, DEC=dec, $
                 XSIZE=xsize, YSIZE=ysize, OUTFIL=outfil, RA=ra

  common sdss_finchk_cmm
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_finchk, dlafil, con_dir, MAG=, RA=, OUTFIL= [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( OUTFIL ) then outfil='finchk.lst'
  if not keyword_set( FILNUM ) then filnum = 55L

; Initialize the common blcok
  sdss_finchk_icmmn, dlafil, MAG=mag, RA=ra, DEC=dec

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }

; STATE

  state = {             $
            ndla: 0L, $
            curdla: gddla[0], $
            dlafil: dlafil, $
            flg_new: 0, $
            zabs: 0., $
            zqso: 0., $
            con_dir: '', $
            minval: minval, $
            ntrans: 0L, $       ; PLOTTING LINES
            vmnx: [-800., 800.], $
            nplt: 0, $
            dla_lin: tmp2, $
            dla_conti: 0., $
            outfil: outfil, $
            filnum: filnum, $
            all_velo: dblarr(5000, 300), $  
            all_pmnx: lonarr(3, 300), $  
            velplt: replicate(tmp, 300), $
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
            base_id: 0L, $      ; Widgets
            ldraw_id: 0L, $    ; Lya
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

  sdss_finchk_llist, state
  state.dla_lin = x_setline(1215.6701d)

; Other setup
  state.ndla = n_elements(dlastr)
  state.con_dir = con_dir

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_finchk: Check spectra', /row, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.lhs_id = WIDGET_BASE( state.base_id, /column, $
                              /base_align_center,/align_center, $
                              uvalue='RHS_BASE', frame=2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Lya DRAW
  state.ldrawbase_id = $
    WIDGET_BASE( state.lhs_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='LDRAW_BASE', frame=2)

  state.ldraw_id = widget_draw(state.ldrawbase_id, xsize=round(xsize*3./5), $
                              ysize=round(2*ysize/3.), /frame, $
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
  state.info_id = $
    WIDGET_BASE( state.lhs_id, /column, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, xsize=round(xsize*3./5), $
               ysize=round(ysize/3.))
  ;; Info
  dlainf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(dlainf, title='Obj ', value=' ', xsize=18)
  state.zabs_id = cw_field(dlainf, title='zabs: ', value=state.zabs, xsize=7)
  state.quality_id = cw_field(dlainf, title='Quality: ', value=0., xsize=5)
  state.scr1_id = cw_field(dlainf, title='NHI Scr: ', value=0., xsize=4)
  state.scr2_id = cw_field(dlainf, title='Mtl Scr: ', value=0., xsize=4)
  state.hits_id = cw_field(dlainf, title='Mtl Hit: ', value=0, xsize=3)

  ;; RA, DEC
  radeci = widget_base(state.info_id, /row, /align_center, frame=2)
  state.mag_id = cw_field(radeci, title='MAG: ', value=0., xsize=10)
  state.ra_id = cw_field(radeci, title='RA: ', value=0., xsize=10)
  state.dec_id = cw_field(radeci, title='DEC: ', value=0., xsize=10)

  ;; Lya
  lyainf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.NHI_id = cw_field(lyainf, title='NHI: ', value=20.3, xsize=7)
  state.NHIb_id = cw_bgroup(lyainf, ['NG', 'S', 'M', 'GD'], /exclusive, column=4, $
                           UVALUE='NHIBGROUP')
  ;; Metals
  mtlinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.mtl_id = cw_bgroup(mtlinf, ['??', 'None','Weak', 'Med', 'Str', 'VStr'], $
                           /exclusive, column=6, UVALUE='MTLBGROUP')

;      BUTTONS
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  state.stat_id = cw_field(butbase, title='Status: ', value='New', xsize=3)
  print = WIDGET_BUTTON(butbase, value='PRINT',uvalue='PRINT')
  splt = WIDGET_BUTTON(butbase, value='SPLT',uvalue='SPLT')
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  dnsv = WIDGET_BUTTON(butbase, value='DNSV',uvalue='DNSV')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Metals DRAW
  state.mdrawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               uvalue='SDRAW_BASE', frame=2)

  state.mdraw_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./5), $
                              ysize=ysize, /frame, retain=2, $
                              uvalue='SDRAW')


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

  ; Set qso and qal
;  if qalstr[0].DLA_z[0] EQ 0. then sdss_finchk_next, state

  ; Load data
  sdss_finchk_setup, state

  ; PLOT
  sdss_finchk_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_finchk', base
  delvarx, fx, wv, npix, sig

  return
end
	
