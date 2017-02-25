;+ 
; NAME:
; sdss_chkmgii
;    Version 1.0
;
; PURPOSE:
;   Visually check SDSS Mg II absorption with a GUI
;
; CALLING SEQUENCE:
;   
;   sdss_chkmgii, x, maskid, expsr, XSIZE=, YSIZE=
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
;   sdss_chkmgii, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Dec-2003 Written by GEP/SHF
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_chkmgii_icmmn, qalfil, mgiifil, badmgiifil

  common sdss_chkmgii_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    fit, $
    qalstr, $
    mgiistr, $
    badmgii

  ;; QALSTR
  qalstr = xmrdfits(qalfil, 2, /silent)

  ;; MGII str
  if x_chkfil(mgiifil, /silent) EQ 1 then $
    mgiistr = xmrdfits(mgiifil, 1, /silent) $
  else begin
    mgiistr = { sdssmgiistrct }
  endelse
  if x_chkfil(badmgiifil, /silent) EQ 1 then $
    badmgii = xmrdfits(badmgiifil, 1, /silent) $
  else begin
    badmgii = { sdssmgiistrct }
  endelse

  return
end
  

;;;;
; Events
;;;;

pro sdss_chkmgii_event, ev

  common sdss_chkmgii_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SPLT': x_specplot, fx, sig, wave=wv, zin=state.zabs, /qal, /block, $
        inflg=4
      'GOOD': begin
          sdss_chkmgii_setmgii, state
          sdss_chkmgii_next, state
          sdss_chkmgii_setup, state
          sdss_chkmgii_update, state
      end
      'BAD': begin
          sdss_chkmgii_setbad, state
          sdss_chkmgii_next, state
          sdss_chkmgii_setup, state
          sdss_chkmgii_update, state
      end
      'DONE' : begin
          sdss_chkmgii_svmgii, state
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          sdss_chkmgii_svmgii, state
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_chkmgii_Metals, state
  
  common sdss_chkmgii_cmm

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
pro sdss_chkmgii_update, state
  sdss_chkmgii_updinfo, state
  sdss_chkmgii_Metals, state
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro sdss_chkmgii_setup, state
  common sdss_chkmgii_cmm
    
  ingood = where(qalstr[state.curqso].qso_name EQ mgiistr.qso_name AND $
                 qalstr[state.curqso].zabs EQ mgiistr.zabs, nig)
  inbad = where(qalstr[state.curqso].qso_name EQ badmgii.qso_name AND $
                 qalstr[state.curqso].zabs EQ badmgii.zabs, nib)

  while(nig NE 0 OR nib NE 0) do begin
    state.curqso = state.curqso + 1
    if state.curqso EQ state.nqal then begin
      sdss_chkmgii_svmgii, state
      widget_control, ev.top, /destroy
      return
    endif
    ingood = where(qalstr[state.curqso].qso_name EQ mgiistr.qso_name AND $
                   qalstr[state.curqso].zabs EQ mgiistr.zabs, nig)
    inbad = where(qalstr[state.curqso].qso_name EQ badmgii.qso_name AND $
                   qalstr[state.curqso].zabs EQ badmgii.zabs, nib)
  endwhile

  ;; Read data
  file_name = '      ' + strtrim(qalstr[state.curqso].sdss_obs[0],2)
  strput, file_name, '/Users/prochter'
  parse_sdss, file_name, fx, wv, conti, $
    SIG=sig, NPIX=npix, CDIR=state.con_dir

  ;; Set zabs
  state.zabs = qalstr[state.curqso].zabs
  state.zqso = qalstr[state.curqso].z_qso

  ;; Set EW
  state.ew = qalstr[state.curqso].ew[1]

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

  ;; MGII 
  state.mgii_lin.N = 20.3
  state.mgii_lin.b = 30.
  state.mgii_lin.zabs = state.zabs

  ;; Fit
  state.mgii_conti = ymd
  fit = replicate(state.mgii_conti, npix)

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


  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set Lines
pro sdss_chkmgii_llist, state

  llist = getenv('XIDL_DIR')+'/SDSS/sdss_mgii.lst'
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
pro sdss_chkmgii_next, state
  common sdss_chkmgii_cmm

  state.curqso = state.curqso + 1 

  if state.curqso EQ state.nqal then begin
    sdss_chkmgii_svmgii, state
    widget_control, ev.top, /destroy
    return
  endif

  ingood = where(qalstr[state.curqso].qso_name EQ mgiistr.qso_name AND $
                 qalstr[state.curqso].zabs EQ mgiistr.zabs, nig) 
  inbad = where(qalstr[state.curqso].qso_name EQ badmgii.qso_name AND $
                 qalstr[state.curqso].zabs EQ badmgii.zabs, nib) 

  while(nig NE 0 OR nib NE 0) do begin
    state.curqso = state.curqso + 1 
    if state.curqso EQ state.nqal then begin
      sdss_chkmgii_svmgii, state
      widget_control, ev.top, /destroy
      return
    endif
    ingood = where(qalstr[state.curqso].qso_name EQ mgiistr.qso_name AND $
                   qalstr[state.curqso].zabs EQ mgiistr.zabs, nig) 
    inbad = where(qalstr[state.curqso].qso_name EQ badmgii.qso_name AND $
                   qalstr[state.curqso].zabs EQ badmgii.zabs, nib) 
  endwhile

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chkmgii_updinfo, state
  common sdss_chkmgii_cmm

  ;; Namej
  widget_control, state.name_id, $
    set_value=strtrim(qalstr[state.curqso].qso_name,2)

  ;; RA, DEC
  widget_control, state.mag_id, set_value=qalstr[state.curqso].rmag
  widget_control, state.ra_id, set_value=qalstr[state.curqso].ra
  widget_control, state.dec_id, set_value=qalstr[state.curqso].dec
  widget_control, state.ew_id, set_value=qalstr[state.curqso].ew[1]

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set BAD MGII structure
pro sdss_chkmgii_setbad, state

  common sdss_chkmgii_cmm

  if badmgii[0].zabs NE 0 then begin
    tmp = { sdssmgiistrct }
    badmgii = [badmgii, tmp]
  endif
  nmgii = n_elements(badmgii)
  state.curmgii = nmgii-1
  state.flg_new = 1
  ;; Name, etc
  badmgii[state.curmgii].qso_name = qalstr[state.curqso].qso_name
  badmgii[state.curmgii].ra = qalstr[state.curqso].ra
  badmgii[state.curmgii].dec = qalstr[state.curqso].dec
  badmgii[state.curmgii].rmag = qalstr[state.curqso].rmag
  badmgii[state.curmgii].sdss_obs = qalstr[state.curqso].sdss_obs
  badmgii[state.curmgii].z_qso = qalstr[state.curqso].z_qso
  badmgii[state.curmgii].zabs = qalstr[state.curqso].zabs
  badmgii[state.curmgii].ew = qalstr[state.curqso].ew
  badmgii[state.curmgii].sigew = qalstr[state.curqso].sigew
  badmgii[state.curmgii].wrest = qalstr[state.curqso].wrest

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set MGII structure
pro sdss_chkmgii_setmgii, state

  common sdss_chkmgii_cmm

  if mgiistr[0].zabs NE 0 then begin
    tmp = { sdssmgiistrct }
    mgiistr = [mgiistr, tmp]
  endif
  nmgii = n_elements(mgiistr)
  state.curmgii = nmgii-1
  state.flg_new = 1
  ;; Name, etc
  mgiistr[state.curmgii].qso_name = qalstr[state.curqso].qso_name
  mgiistr[state.curmgii].ra = qalstr[state.curqso].ra
  mgiistr[state.curmgii].dec = qalstr[state.curqso].dec
  mgiistr[state.curmgii].rmag = qalstr[state.curqso].rmag
  mgiistr[state.curmgii].sdss_obs = qalstr[state.curqso].sdss_obs
  mgiistr[state.curmgii].z_qso = qalstr[state.curqso].z_qso
  mgiistr[state.curmgii].zabs = qalstr[state.curqso].zabs
  mgiistr[state.curmgii].ew = qalstr[state.curqso].ew
  mgiistr[state.curmgii].sigew = qalstr[state.curqso].sigew
  mgiistr[state.curmgii].wrest = qalstr[state.curqso].wrest

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkmgii_svmgii, state

  common sdss_chkmgii_cmm

  mwrfits, mgiistr, state.mgiifil, /create
  mwrfits, badmgii, state.badmgiifil, /create
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

pro sdss_chkmgii, qalfil, mgiifil, badmgiifil, con_dir, IQSO=iqso, FNEW=fnew, $
              XSIZE=xsize, L_YSIZE=i_ysize, M_YSIZE=s_ysize, $
              YSIZE=ysize, OBJNM=objnm, ZIN=zin, XMAX=xmax, LLIST=llist

  common sdss_chkmgii_cmm
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'sdss_chkmgii, qalfil, mgiifil, badmgiifil, con_dir, /FNEW'
    print, '        I_YSIZE=, S_YSIZE= [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 800
  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = 0L

; Initialize the common blcok
  sdss_chkmgii_icmmn, qalfil, mgiifil, badmgiifil

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }

; STATE

  state = {             $
            nqal: 0L, $
            curqso: iqso, $
            curqal: 0, $
            curmgii: 0L, $
            mgiifil: mgiifil, $
            badmgiifil: badmgiifil, $
            flg_new: 0, $
            zabs: 0., $
            zqso: 0., $
            con_dir: '', $
            minval: minval, $
            ntrans: 0L, $       ; PLOTTING LINES
            vmnx: [-800., 800.], $
            nplt: 0, $
            mgii_lin: tmp2, $
            mgii_conti: 0., $
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
            top_id: 0L, $
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
            ew_id: 0L, $
            help_text_id: 0L, $
	    ew: 0L$
          }

;;;;;;;;;;;;;;
; SETUP LINES
  sdss_chkmgii_llist, state
  state.mgii_lin = x_setline(1215.6701d)

; Other setup
  state.nqal = n_elements(qalstr)
  state.con_dir = con_dir

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_chkmgii: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.top_id = WIDGET_BASE( state.base_id, /column, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(2*ysize/3.), $
                              uvalue='TOP_BASE', frame=2)

;;;;;; Info window ;;;;;;;;;;;
  state.info_id = $
    WIDGET_BASE( state.base_id, /column, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, xsize=xsize, $
               ysize=round(ysize/3.))
  ;; Info
  mgiiinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(mgiiinf, title='Obj ', value=' ', xsize=18)
  state.zabs_id = cw_field(mgiiinf, title='zabs: ', value=state.zabs, xsize=7)
  state.ew_id = cw_field(mgiiinf, title='EW: ', value=state.ew, xsize=7)

  ;; RA, DEC
  radeci = widget_base(state.info_id, /row, /align_center, frame=2)
  state.mag_id = cw_field(radeci, title='MAG: ', value=0., xsize=10)
  state.ra_id = cw_field(radeci, title='RA: ', value=0., xsize=10)
  state.dec_id = cw_field(radeci, title='DEC: ', value=0., xsize=10)

  ;; Lya
  lyainf = widget_base(state.info_id, /row, /align_center, frame=2)

;      BUTTONS
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  good = WIDGET_BUTTON(butbase, value='GOOD',uvalue='GOOD')
  bad  = WIDGET_BUTTON(butbase, value='BAD', uvalue='BAD')
  save = WIDGET_BUTTON(butbase, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  splt = WIDGET_BUTTON(butbase, value='SPLT',uvalue='SPLT')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Metals DRAW
  state.mdrawbase_id = $
    WIDGET_BASE( state.top_id, /row, /base_align_center,/align_center, $
               uvalue='SDRAW_BASE', frame=2)

  state.mdraw_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./3), $
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
  if qalstr[0].zabs[0] EQ 0. then sdss_chkmgii_next, state

  ; Load data
  sdss_chkmgii_setup, state

  ; PLOT
  sdss_chkmgii_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chkmgii', base
  delvarx, fx, wv, npix, sig, qalstr

  return
end
	
