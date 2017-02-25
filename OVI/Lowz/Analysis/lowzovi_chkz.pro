;+ 
; NAME:
; lowzovi_chkz
;    Version 1.0
;
; PURPOSE:
;   GUI used to check the redshift for a galaxy in the LCOOVI survey.
;
; CALLING SEQUENCE:
;   lowzovi_chkz, lowzovi, XSIZE=, YSIZE=, /PRINTONLY, NPLT=, SOBJ=
;     BADFIL=
;
; INPUTS:
;  lowzovi -- LCO galaxy survey structure
;
; RETURNS:
;
; OUTPUTS:
;  BADFIL= -- ASCII file containing the list of bad redshifts
;
; OPTIONAL KEYWORDS:
; SOBJ=  -- Starting object [default: the first]
; NPLT=  -- Number of galaxies to plot at a time (with and without
;           fit shown)  [default: 1L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lowzovi_chkz, lowzovi, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lowzovi_chkz_initlzovin

common lowzovi_chkz_common, lzovi_wff, lzovi_zans, lzovi_svbad

if keyword_set( lzovi_svbad ) then delvarx, lzovi_svbad

return
end

;;;;
; Events
;;;;

pro lowzovi_chkz_event, ev

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
      'PRINT' : 
;      'PRINT' : lowzovi_chkz_print, state
      'EXP_Y' : begin ;; Expand the y-axis scale
          state.scl_y = state.scl_y * 2
          lowzovi_chkz_Plot, state
      end
      'LOW_Y' : begin ;; Expand the y-axis scale
          state.scl_y = state.scl_y / 2
          lowzovi_chkz_Plot, state
      end
      'GOOD' : begin
          state.curobj =  state.curobj + 1
          if state.curobj EQ state.nobj then state.curobj = 0L
          lowzovi_chkz_Setup, state
          lowzovi_chkz_Plot, state
      end
      'BAD' : begin
          lowzovi_chkz_Svbad, state
          state.curobj =  state.curobj + 1
          if state.curobj EQ state.nobj then state.curobj = 0L
          lowzovi_chkz_Setup, state
          lowzovi_chkz_Plot, state
      end
      'TEMPL' : begin
	  state.flg_zfind = 1 - state.flg_zfind
          lowzovi_chkz_Plot, state
      end
      'PLOT2D' : begin
          widget_control, /hourglass   
          wfccd_pltobj, state.sv_fil, state.mod_obj, 0L, /fspec
      end
      'SMOOTH' : begin
          widget_control, state.smooth_id, get_value=smth
          state.smooth = smth
          lowzovi_chkz_Plot, state
      end
      'PREV' : begin
          state.curobj =  state.curobj - 1
          if state.curobj EQ -1 then state.curobj = state.nobj-1
          lowzovi_chkz_Setup, state
          lowzovi_chkz_Plot, state
      end
      'OBJ': begin
          widget_control, state.obj_id, get_value=obj
          all_obj = strtrim(state.galstr.id,2)+strtrim(state.galstr.obj_id,2)
          indx = where(all_obj EQ strtrim(obj[0],2), nindx)
          if nindx EQ 0 then begin
              print, 'lowzovi_chkz: Object not identified!  Ignoring...'
              stop
          endif else begin
              state.curobj = indx[0]
              lowzovi_chkz_Setup, state
              lowzovi_chkz_Plot, state
          endelse
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

pro lowzovi_chkz_Plot, state
  
; Plot Data

  common lowzovi_chkz_common

  ; PSFILE
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  ; PLOT
  clr = getcolor(/load)

  !P.MULTI= [0,1,state.nplt*2]

  
;  for i=pmin,pmax do begin
  npix = state.npix
                                ; Get median flux
  mdwv = where(state.wave[0:npix-1] GT state.mdrng[0] AND $
               state.wave[0:npix-1] LT state.mdrng[1] AND $
               state.fx[0:npix-1] GT 0., nwv)
  if nwv NE 0 then mdfx = median(state.fx[mdwv]) else mdfx = 1.
  mdfx = mdfx * state.scl_y

  ;; Plot smooth
  spaces = replicate('!17 ',30)
  plot, state.wave[0:npix-1], $
    median(state.fx[0:npix-1], state.smooth), $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[0., 2*mdfx], xtickn=spaces, xmargin=[14,0], ymargin=[0,0], $
    charsize=1.5, psym=10, background=clr.white, color=clr.black, $
    xstyle=1, ystyle=1

  ;; Zfind
  if state.flg_zfind EQ 1 AND lzovi_zans.z_err NE 0. then begin
      synth = synthspec(lzovi_zans, $
                        loglam=alog10(state.wave[0:npix-1]))
      oplot, state.wave[0:npix-1], synth, color=clr.blue
  endif

  ;; Key lines
  for jj=0L,state.nkey-1 do begin
      wv=state.keylines[jj]*(lzovi_zans.z+1)
      oplot, [wv,wv], [-1e20,1e20], color=clr.gray, linestyle=2
  endfor

  ;; Plot Normal
  plot, state.wave[0:npix-1], $
    state.fx[0:npix-1], $
    xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[0., 2*mdfx], xmargin=[14,0], ymargin=[3,0], $
    charsize=1.5, psym=10, background=clr.white, color=clr.black, $
    xstyle=1, ystyle=1
  ;; Error
  oplot, state.wave[0:npix-1], state.sig[0:npix-1], $
    color=clr.red

  ;; Zfind
  if state.flg_zfind EQ 1 AND lzovi_zans.z_err NE 0. then begin
      synth = synthspec(lzovi_zans, $
                        loglam=alog10(state.wave[0:npix-1]))
      oplot, state.wave[0:npix-1], synth, color=clr.blue
      xyouts, 0.83*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], $
        0.85*2*mdfx, 'z= '+string(lzovi_zans.z,'(f7.4)'), $
        color=clr.blue, charsize=1.5
      xyouts, 0.83*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], $
        0.70*2*mdfx, 'chi= '+string(lzovi_zans.rchi2,'(f5.2)'), $
        color=clr.red, charsize=1.5
  endif
  ;; Key lines
  for jj=0L,state.nkey-1 do begin
      wv=state.keylines[jj]*(lzovi_zans.z+1)
      oplot, [wv,wv], [-1e20,1e20], color=clr.gray, linestyle=2
  endfor

  ;; Labels
  xyouts, 0.90*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.75*2*mdfx, $
    state.obj, color=clr.green, charsize=2.5
;  xyouts, 0.03*(state.xymnx[2]-state.xymnx[0])+state.xymnx[0], 0.45*2*mdfx, $
;    'Obj: '+strtrim(i,2), color=clr.blue, charsize=1.5

  !P.MULTI= [0,1,1]

end

;;;;;;;;;;;;;;;;;;;;
;  Setup
;;;;;;;;;;;;;;;;;;;;

pro lowzovi_chkz_Setup, state

  common lowzovi_chkz_common
  ;; Check fspec_fil
  ffil = where(strlen(strtrim(state.galstr[state.curobj].fspec_fil,2)) GT 0, nfil)
  case nfil of
      0: stop
      1: 
      else: print, 'lowzovi_chkz_Setup: Multiple fspec_fil.  Taking first! '
  endcase

  ;; File
  findx = ffil[0]
  fspec_fil = strtrim(state.galstr[state.curobj].fspec_fil[findx],2)

  if state.sv_fil NE fspec_fil then begin
      ;; save
      state.sv_fil = fspec_fil
      ;; Read in 
      wfccd_wrfspec, wffspec, fspec_fil, /read
      lzovi_wff = temporary(wffspec)
  endif

  ;; Find obj
  state.mod_obj = strtrim(state.galstr[state.curobj].id MOD 10000L,2)+ $
    strtrim(state.galstr[state.curobj].obj_id,2)
  state.obj = strtrim(state.galstr[state.curobj].id,2)+ $
    strtrim(state.galstr[state.curobj].obj_id,2)
  indx = x_getobjnm(lzovi_wff, state.mod_obj)
  if indx LT 0 then stop
  widget_control, state.obj_id, set_value=state.obj

  ;; Load it up
  np = lzovi_wff[indx].npix

  state.npix = np
  state.wave[0:np-1] = lzovi_wff[indx].wave[0:np-1]
  state.fx[0:np-1] = lzovi_wff[indx].fx[0:np-1]

  ;; Sig
  state.sig[0:np-1] = fltarr(np)
  a = where(lzovi_wff[indx].var[0:np-1] GT 0.)
  state.sig[a] = float(sqrt(lzovi_wff[indx].var[a]))

  ;; zans
  lzovi_zans = lzovi_wff[indx].zans

  ;; Redshift
  widget_control, state.z_id, set_value=lzovi_zans.z

  ;; Redshift
  widget_control, state.indx_id, set_value=state.curobj

end

;;;;;;;;;;;;;;;;;;;;
;  Print
;;;;;;;;;;;;;;;;;;;;

pro lowzovi_chkz_print, state

  widget_control, /hourglass   
; Device
  device, get_decomposed=svdecomp

;  !p.thick = 1
;  !p.charthick = 1

  device, decompose=0
  ; Get file name
;  ipos = strpos(state.obj_fil, '.fits')
;  psfil = strmid(state.obj_fil, 0, ipos)+'.ps'
  psfil = 'lowzovi_chkz.ps'
  ps_open, file=psfil, font=1, /color
  state.psfile = 1
  lowzovi_chkz_plot, state
  ps_close, /noprint, /noid
  spawn, 'gzip -f '+psfil
  device, decomposed=svdecomp
  state.psfile = 0
;  !p.thick = 1
;  !p.charthick = 1

end

;;;;;;;;;;;;;;;;;;;;
;  Save bad
;;;;;;;;;;;;;;;;;;;;

pro lowzovi_chkz_svbad, state

  common lowzovi_chkz_common

  ;; 
  if not keyword_set( lzovi_svbad ) then lzovi_svbad = state.obj $
  else lzovi_svbad = [lzovi_svbad, state.obj]

  close, 77
  openw, 77, state.bad_fil
  for i=0L,n_elements(lzovi_svbad)-1 do printf, 77, lzovi_svbad[i]

  close, 77
  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lowzovi_chkz, lowzovi_fil, XSIZE=xsize, YSIZE=ysize, $
                   PRINTONLY=printonly, NPLT=nplt, SOBJ=sobj, BADFIL=badfil

common lowzovi_chkz_common
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'lowzovi_chkz, lowzovi_fil, XSIZE=, YSIZE=, /PRINTONLY '+$
	'BADFIL=, NPLT=, SOBJ= [v1.1]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( BADFIL ) then badfil='badz.lst'

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then xsize = ssz[0]-200 ;1200
  if not keyword_set( YSIZE ) then ysize = ssz[1]-200 ; 900
  if not keyword_set( XYMNX ) then xymnx = [3800., 0.0, 8700., 1.0]
  if not keyword_set( NPLT ) then nplt = 1L

; Check for file
;  flg = x_chkfil(fspec_fil)
;  if flg NE 1 then return
      
  
  lowzovi_chkz_initlzovin

  if size(lowzovi_fil, /type) EQ 7 then $
    lzovi_gal = xmrdfits(lowzovi_fil, 1, structyp='galsurveystrct',/silent) $
  else lzovi_gal = lowzovi_fil

  ;; PARSE
  gdgal = where(lzovi_gal.obj_id EQ 'a' AND $
                lzovi_gal.flg_anly MOD 4 GT 1, ngd) 
  lzovi_gal = lzovi_gal[gdgal]

;                lzovi_gal.flg_anly MOD 8 GT 3 AND $
  

; STATE

  state = {             $
            nplt: nplt, $
            npg: (ngd/nplt) + (ngd MOD nplt NE 0), $
            nobj: ngd, $
            galstr: lzovi_gal, $
            wave: dblarr(5000), $
            fx: fltarr(5000), $
            sig: fltarr(5000), $
            smooth: 7L, $
            npix: 0L, $
            obj: ' ', $
            mod_obj: ' ', $
            curobj: 0L, $
            mdrng: [5000., 8000.], $
            pos: [0.04,0.0,1.00,1.0], $ ; Plotting
            flg_zfind: 1, $
            sv_fil: ' ', $
            bad_fil: badfil, $
            keylines: fltarr(100), $
            nkey: 0L, $
          scl_y: 1., $
            xymnx: xymnx, $
            tmpxy: fltarr(4), $
            xcurs: 0., $
            xpos: 0., $
            size: lonarr(2), $
            psfile: 0, $
            base_id: 0L, $      ; Widgets
            draw_id: 0L, $
            drawbase_id: 0L, $
            obj_id: 0L, $
            nspec_id: 0L, $
            xpos_id: 0L, $
            indx_id: 0L, $
            ntot_id: 0L, $
            z_id: 0L, $
            smooth_id: 0L, $
            help_text_id: 0L $
          }

; Zfind
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

;    WIDGET
  base = WIDGET_BASE( title = 'lowzovi_chkz: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L, yoffset=50L)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='lowzovi_chkz', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)
;  state.name_id = WIDGET_LABEL(labelbase, value=state.title, /align_center)

;;;;;;;;;;;;;;;;;;
; Pages
  objbase = widget_base(toolbar, /column, /base_align_center, /align_center)
  state.obj_id = cw_field(objbase, value='', /return_events, xsize=8,$
                          title='Obj', UVALUE='OBJ', $
                          /string, font='Courier*Bold')
  state.indx_id = cw_field(objbase, title='Index', value=state.curobj, xsize=7)
  state.ntot_id = cw_field(objbase, title='Ntot', value=state.nobj, xsize=7)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  state.drawbase_id = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_id = widget_draw(state.drawbase_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')

; XY position + redshift
  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.z_id = cw_field(xy_id, title='z:', value=0., xsize=13)
  state.smooth_id = cw_field(xy_id, title='Smooth', value=state.smooth, xsize=5,$
                            /return_events, /long, UVALUE='SMOOTH')
;  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

;      BUTTONS
  butbase = widget_base(toolbar, /column, /align_center)
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  good = WIDGET_BUTTON(butbase, value='GOOD',uvalue='GOOD')
  bad = WIDGET_BUTTON(butbase, value='BAD',uvalue='BAD')
  butbase2 = widget_base(toolbar, /column, /align_center)
  plot2d = WIDGET_BUTTON(butbase2, value='PLOT2D',uvalue='PLOT2D')
  templt = WIDGET_BUTTON(butbase2, value='TEMPL',uvalue='TEMPL')
  print = WIDGET_BUTTON(butbase2, value='PRINT',uvalue='PRINT')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  butbase3 = widget_base(toolbar, /column, /align_center)
  expandy = WIDGET_BUTTON(butbase3, value='EXP_Y',uvalue='EXP_Y')
  lowandy = WIDGET_BUTTON(butbase3, value='LOW_Y',uvalue='LOW_Y')

  
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

  ; Setup
  lowzovi_chkz_Setup, state

  ; PLOT
  lowzovi_chkz_Plot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'lowzovi_chkz', base, /no_block

  !P.MULTI= [0,1,1]
  return
end

