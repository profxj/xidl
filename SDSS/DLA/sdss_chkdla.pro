;+ 
; NAME:
; sdss_chkdla
;    Version 1.1
;
; PURPOSE:
;   Launches a GUI that is used to check the DLA candidates 
;  (and metal-line systems) identified by the automated algorithm
;
; CALLING SEQUENCE:
; sdss_chkdla, qalfil, dlafil, con_dir, IQSO=, /FNEW, XSIZE=, YSIZE=
;
; INPUTS:
; qalfil  -- FITS file containing the QAL structure
; dlafil  -- FITS file containing the SDSS DLA structure (may be old) 
; con_dir -- Path to the directory containting the continuum files
;
; RETURNS:
;
; OUTPUTS:
; dlafil -- FITS file containing the SDSS DLA structure
;
; OPTIONAL KEYWORDS:
;  /FNEW -- Start with first quasar that hasnt been checked
;  OBJNM= --
;  XSIZE -- Size of gui in screen x-pixels [default: screen-200]
;  YSIZE -- Size of gui in screen y-pixels [default: screen-200]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_chkdla, 'sdss_dr3_QAL.fits',  'dr3_qal_chkd.fits', 'ABSLIN/',
;   /FNEW
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

pro sdss_chkdla_icmmn, qalfil, dlafil

  common sdss_chkdla_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    fit, $
    qalstr, $
    dlastr

  ;; QALSTR
  qalstr = xmrdfits(qalfil, 1, /silent)
  print, 'n total = ', n_elements(qalstr)

  ;; DLA str
  if x_chkfil(dlafil+'*', /silent) EQ 1 then $
    dlastr = xmrdfits(dlafil, 1, /silent) $
  else begin
      dlastr = { sdssdlastrct }
  endelse

  return
end
  

;;;;
; Events
;;;;

pro sdss_chkdla_event, ev

  common sdss_chkdla_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  if not keyword_set(UVAL) then stop
  case uval of
      'SPLT': x_specplot, fx, sig, wave=wv, zin=state.zabs, /lls, /block, $
        inflg=4
      
      'SAVE': sdss_chkdla_svdla, state
      
      'NEXT': begin
          sdss_chkdla_next, state
          sdss_chkdla_setup, state
          sdss_chkdla_update, state
      end
      'NGNEW': begin
         ;; Set the values
         dlastr[state.curdla].flg_mtl = 1
         dlastr[state.curdla].flg_NHI = 0
         dlastr[state.curdla].NHI = state.dla_lin.N
         ;; Continue
         sdss_chkdla_next, state
         sdss_chkdla_new, state
         sdss_chkdla_update, state
      end
      'NEW': begin
          sdss_chkdla_next, state
          sdss_chkdla_new, state
          sdss_chkdla_update, state
      end

      'PREV': begin
          sdss_chkdla_prev, state
          sdss_chkdla_setup, state
          sdss_chkdla_update, state
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
                          sdss_chkdla_updfit, state
                          widget_control, state.zabs_id, $
                            set_value=strtrim(state.zabs,2)
                      end
                      4 : state.dla_lin.N = state.dla_lin.N - 0.05
                  endcase
                  sdss_chkdla_updinfo, state
                  sdss_chkdla_updfit, state
                  sdss_chkdla_Lya, state
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
              'x': begin
                  sdss_chkdla_updfit, state, /exact
                  sdss_chkdla_Lya, state
              end
              else: 
          endcase
          sdss_chkdla_Lya, state
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
          sdss_chkdla_svdla, state
          widget_control, ev.top, /destroy
          return
      end
      'DNSV' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  sdss_chkdla_upddstr, state

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Lya Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_chkdla_Lya, state
  
  common sdss_chkdla_cmm

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


pro sdss_chkdla_Metals, state
  
  common sdss_chkdla_cmm

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
pro sdss_chkdla_update, state
  sdss_chkdla_updfit, state
  sdss_chkdla_updinfo, state
  sdss_chkdla_Lya, state
  sdss_chkdla_Metals, state
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkdla_updfit, state, EXACT=exact

  common sdss_chkdla_cmm

  if not keyword_set( wvoff ) then wvoff = 100.

  ;; FIT
  fit[*] = 1.
  if keyword_set( EXACT ) then begin
      mnwv = state.dla_lin[0].wrest*(1.+state.dla_lin[0].zabs)
      mn = min(abs(wv - mnwv + wvoff), mnpx)
      mx = min(abs(wv - mnwv - wvoff), mxpx)
      fit[mnpx:mxpx] = x_voigt(wv[mnpx:mxpx], $
                               state.dla_lin[0], FWHM=3.)
  endif else fit = x_allvoigt(wv, state.dla_lin, SIGMA=3.)

  fit = fit*state.dla_conti

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro sdss_chkdla_setup, state
  common sdss_chkdla_cmm

  ;; Read data
  parse_sdss, getenv(state.sdsspath) $
    +strtrim(qalstr[state.curqso].file_name,2), fx, wv, conti, $
    SIG=sig, NPIX=npix, CDIR=state.con_dir

  ;; Set zabs
  state.zabs = qalstr[state.curqso].dla_z[state.curqal]
  state.zqso = qalstr[state.curqso].z_qso

  ;; xymnx
  state.xymnx[0] = 1215.6701*(1.+state.zabs) - 200.
  state.xymnx[2] = 1215.6701*(1.+state.zabs) + 200.
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  if ngd GT 1 then begin
      srt = sort(fx[gd])
      ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
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
  sdss_chkdla_setdla, state

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set Lines
pro sdss_chkdla_llist, state

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
pro sdss_chkdla_next, state
  common sdss_chkdla_cmm

  if qalstr[state.curqso].DLA_z[state.curqal+1] GT 0. then begin
      state.curqal = state.curqal + 1 
  endif else begin
      state.curqso = state.curqso + 1
      ;; Check for end
      nqso = n_elements(qalstr)
      if state.curqso EQ nqso then state.curqso = 0L
      state.curqal = 0L
  endelse

  if qalstr[state.curqso].DLA_z[state.curqal] LT 0.001 then sdss_chkdla_next, state

  ;; Quality assessment
  if qalstr[state.curqso].dla_quality[state.curqal] LT state.minval then begin
      if state.flg_lya EQ 0 then sdss_chkdla_next, state $ 
      else begin
          a = where(abs(qalstr[state.curqso].dla_z[state.curqal] - $
                        qalstr[state.curqso].dla_zabs1) LT 0.02, na)
          if na EQ 0 then sdss_chkdla_next, state
      endelse
  endif

  widget_control, state.nqal_id, set_value=state.curqso
  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New
pro sdss_chkdla_new, state
  common sdss_chkdla_cmm

  na = 1
  while( na NE 0 ) do begin
      a = where(abs(dlastr.ra - qalstr[state.curqso].ra) LT 0.0005 AND $
                abs(dlastr.dec - qalstr[state.curqso].dec) LT 0.0005 AND $
                abs(dlastr.zabs - qalstr[state.curqso].dla_z[state.curqal]) $
                LT 0.01,  na)
      if na NE 0 then sdss_chkdla_next, state
  endwhile
  sdss_chkdla_setup, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro sdss_chkdla_prev, state
  common sdss_chkdla_cmm

  if state.curqal EQ 0 then begin
      state.curqso = state.curqso - 1
      if state.curqso EQ -1 then stop
      gd = where(qalstr[state.curqso].dla_z GT 0. AND $
                 qalstr[state.curqso].dla_quality GT state.minval, ngd)
      ;; Good?
      if ngd EQ 0 then sdss_chkdla_prev, state $
      else state.curqal = ngd-1
  endif else begin
      state.curqal = state.curqal - 1
      if qalstr[state.curqso].dla_quality[state.curqal] LT state.minval then $
        sdss_chkdla_prev, state
  endelse

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chkdla_updinfo, state
  common sdss_chkdla_cmm

  ;; Namej
  widget_control, state.name_id, $
    set_value=strtrim(qalstr[state.curqso].qso_name,2)

  ;; RA, DEC
  widget_control, state.mag_id, set_value=qalstr[state.curqso].qso_mag
  widget_control, state.ra_id, set_value=qalstr[state.curqso].ra/15.
  widget_control, state.dec_id, set_value=qalstr[state.curqso].dec

  ;; Quality
  widget_control, state.quality_id, $
    set_value=qalstr[state.curqso].dla_quality[state.curqal]

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs

  ;; Metals
  mn = min(abs(qalstr[state.curqso].dla_zabs2-state.zabs), imn)
  if mn LT 0.01 then begin
      widget_control, state.hits_id, $
        set_value=round(qalstr[state.curqso].dla_hits[imn])
      widget_control, state.scr2_id, set_value=qalstr[state.curqso].dla_score2[imn]
  endif else begin
      widget_control, state.hits_id, set_value=0.
      widget_control, state.scr2_id, set_value=0.
  endelse

  ;; NHI
  mn = min(abs(qalstr[state.curqso].dla_zabs1-state.zabs), imn)
  if mn LT 0.01 then $
    widget_control, state.scr1_id, set_value=qalstr[state.curqso].dla_score1[imn] $
  else $
    widget_control, state.scr1_id, set_value=0.

  widget_control, state.NHI_id, set_value=state.dla_lin.N

  if state.flg_new NE 0 then begin
      if state.dla_lin.N GE 20.8 then $
        widget_control, state.NHIb_id, set_value=3
      if state.dla_lin.N LE 20.2 then $
        widget_control, state.NHIb_id, set_value=1
      if state.dla_lin.N GT 20.2 AND state.dla_lin.N LT 20.8 then $
        widget_control, state.NHIb_id, set_value=2
  endif 

  ;; Metals
  if state.flg_new NE 0 then $
      widget_control, state.mtl_id, set_value=0

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set DLA structure
pro sdss_chkdla_setdla, state

  common sdss_chkdla_cmm

  ;; New?
  a = where(abs(dlastr.ra - qalstr[state.curqso].ra) LT 0.00001 AND $
            abs(dlastr.dec - qalstr[state.curqso].dec) LT 0.0001 AND $
            abs(dlastr.zabs - state.zabs) LT 0.01,  na)
  case na of 
      0: begin  ;; New
;stop          
tmp = { sdssdlastrct }
          dlastr = [dlastr, tmp]
          ndla = n_elements(dlastr)
          state.curdla = ndla-1
          state.flg_new = 1
          ;; Name, etc
          dlastr[state.curdla].qso_name = qalstr[state.curqso].qso_name
          dlastr[state.curdla].ra = qalstr[state.curqso].ra
          dlastr[state.curdla].Rmag = qalstr[state.curqso].qso_mag
          dlastr[state.curdla].dec = qalstr[state.curqso].dec
          dlastr[state.curdla].sdss_obs[0] = qalstr[state.curqso].file_name
          dlastr[state.curdla].zabs = state.zabs
          dlastr[state.curdla].quality = $
            qalstr[state.curqso].dla_quality[state.curqal]
          dlastr[state.curdla].z_qso = qalstr[state.curqso].z_qso
          dlastr[state.curdla].conti = state.dla_conti
          ;; Widget
          widget_control, state.stat_id, set_value='New'
          widget_control, state.mtl_id, set_value=0
      end
      1: begin  ;; Old
;stop    
          state.curdla = a[0]
          state.dla_lin.N = dlastr[state.curdla].NHI
          state.dla_conti = dlastr[state.curdla].conti
          state.flg_new = 0
          ;; Widgets
          widget_control, state.NHIb_id, set_value=dlastr[state.curdla].flg_NHI
          widget_control, state.mtl_id, set_value=dlastr[state.curdla].flg_mtl
          ;; Name
          b = where(dlastr[state.curdla].sdss_obs EQ $
                    qalstr[state.curqso].file_name, nb)
          if nb EQ 0 then begin
              bb = where(strlen(strtrim(dlastr[state.curdla].sdss_obs,2)) NE 0)
              dlastr[state.curdla].sdss_obs[bb[0]] = $
                qalstr[state.curqso].file_name
          endif $
          else dlastr[state.curdla].sdss_obs[0] = qalstr[state.curqso].file_name
          ;; Widget
          widget_control, state.stat_id, set_value='Old'
      end
      else: stop
  endcase
          
  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Update DLA structure
pro sdss_chkdla_upddstr, state

  common sdss_chkdla_cmm

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkdla_svdla, state

  common sdss_chkdla_cmm

  sdss_chkdla_upddstr, state
  mwrfits, dlastr, state.dlafil, /create
  ;; Compress
  print, 'sdss_chkdla:  Saving and compressing '+state.dlafil
  spawn, 'gzip -f '+state.dlafil

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

pro sdss_chkdla, qalfil, dlafil, con_dir, IQSO=iqso, FNEW=fnew, $
              XSIZE=xsize, YSIZE=ysize, REDLA=redla, LYA=lya

  common sdss_chkdla_cmm
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'sdss_chkdla, qalfil, dlafil, con_dir, /FNEW, IQSO=, XSIZE=,YSIZE= ' + $
      '/REDLA, /LYA [v1.2]'
    return
  endif 

;  stop  ;; FIX BUG ON PLATE/FIB/MJD

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then xsize = ssz[0]-100
  if not keyword_set( YSIZE ) then ysize = ssz[1]-100
  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = 0L
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'   

; Initialize the common blcok
  sdss_chkdla_icmmn, qalfil, dlafil

  tmp = { velpltstrct }
  tmp2 = { newabslinstrct }

; STATE

  state = {             $
            nqal: 0L, $
            curqso: iqso, $
            curqal: 0, $
            curdla: 0L, $
            dlafil: dlafil, $
            flg_lya: keyword_set(LYA), $
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
            sdsspath: sdsspath, $
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
            nqal_id: 0L, $
            help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES
  sdss_chkdla_llist, state
  state.dla_lin = x_setline(1215.6701d)

; Other setup
  state.nqal = n_elements(qalstr)
  state.con_dir = con_dir

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_chkdla: Check spectra', /row, $
                    UNAME='BASE')
;                    UNAME='BASE', /tlb_size_events, xoffset=200L)
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
                              ysize=round(ysize/2.), /frame, $
                              /button_events, /motion_events, uvalue='LDRAW')

  state.size[0]=round(xsize*3/5.)
  state.size[1]=ysize/2.

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
               ysize=round(ysize/2.))
  ;; Info
  dlainf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(dlainf, title='O', value=' ', xsize=18)
  state.zabs_id = cw_field(dlainf, title='za:', value=state.zabs, xsize=7)
  state.quality_id = cw_field(dlainf, title='Qual:', value=0., xsize=5)
  state.scr1_id = cw_field(dlainf, title='NHI S:', value=0., xsize=4)
  state.scr2_id = cw_field(dlainf, title='M S:', value=0., xsize=4)
  state.hits_id = cw_field(dlainf, title='M H:', value=0, xsize=3)

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
  ;; Counter
  qalinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.nqal_id = cw_field(qalinf, title='QAL #:', value=0., xsize=5)

;      BUTTONS
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  state.stat_id = cw_field(butbase, title='Status: ', value='New', xsize=3)
  splt = WIDGET_BUTTON(butbase, value='SPLT',uvalue='SPLT')
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  ngnew  = WIDGET_BUTTON(butbase, value='NGNEW' ,uvalue='NGNEW')
  new  = WIDGET_BUTTON(butbase, value='NEW' ,uvalue='NEW')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  dnsv = WIDGET_BUTTON(butbase, value='DONOTSV',uvalue='DNSV')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Metals DRAW
  state.mdrawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               uvalue='SDRAW_BASE', frame=2)

  state.mdraw_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./5), $
                              ysize=round(ysize*0.9), /frame, retain=2, $
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
  if qalstr[0].DLA_z[0] EQ 0. then sdss_chkdla_next, state

  ; Load data
  sdss_chkdla_setup, state

  ; New
  if keyword_set( FNEW ) then sdss_chkdla_new, state
  
  ; PLOT
  sdss_chkdla_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chkdla', base
  delvarx, fx, wv, npix, sig, qalstr

  return
end
	
