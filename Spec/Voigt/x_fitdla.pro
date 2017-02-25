;+ 
; NAME:
; x_fitdla   
;   Version 1.12
;
; PURPOSE:
;    GUI used to fit DLA profiles interactively.  The user determines
;    the continuum at the same time.
;
; CALLING SEQUENCE:
;  x_fitdla, flux_fil, [err_fil], XSIZE=, YSIZE=, TITLE=, 
;               WAVE=, /BLOCK, FWHM=, INIFIT=, INFLG=
;
; INPUTS:
;   flux_fil   - FITS file (or array) containing flux
;   [err_fil]  - FITS error array
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array 
;   FWHM=      - Resoltuion (FWHM in pixels)  [default: 2]
;   INFLG=     - Usual flag for reading the FITS file (see x_readspec)
;   INIFIT=    - IDL file containing a saved version of the fit
;   BIN=       - bin by # pixels specified
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;  A qso template can be overplotted to help defining the continuum.
;  To do this, first create a file called 'template.txt' in the
;  directory where x_fitdla will be run.  It must have two columns -
;  rest wavelength in Angstroms and flux.  Then create another file
;  called 'redshifts.txt' in the same directory.  It must also have
;  two columns - the first is the spectrum filename (just the
;  filename, no path) and the second is the redshift of the spectrum
;  in that file.
;
; EXAMPLES:
;   x_fitdla, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;   02-Aug-2009 Mendez altered the matrix
;   07-Sep-2011 added keyword smooth=# which smooths the data by # pixels (MN)
;   03-Apr-2014 added keyword bin=# which bins the data by # pixels (MR & MN)
;-
;------------------------------------------------------------------------------

;;;;
; Events
;;;;

pro x_fitdla_initcommon, nspec, INIFIL=inifil

  common x_fitdla_common, $
     fd_nspec, $
     fd_all_cstr, $
     fd_all_npix, $
     fd_all_conti, $
     fd_all_fx, $
     fd_all_sig, $
     fd_all_wave, $
     fd_svxymnx, $
     fd_zqso, $
     OUT_CONTI

  if not keyword_set(INIFIL) then begin
     fd_nspec = nspec
     fd_all_npix = lonarr(nspec)
     fd_all_fx = fltarr(300000L, nspec)
     fd_all_sig = fltarr(300000L, nspec)
     fd_all_conti = fltarr(300000L, nspec)
     fd_all_wave = dblarr(300000L, nspec)
     fd_zqso = fltarr(nspec)
     fd_svxymnx = fltarr(4,nspec)
     fd_all_cstr = replicate({contistrct}, nspec)
  endif else begin
     restore, inifil
  endelse

  return
end

pro x_fitdla_event, ev

common x_fitdla_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval
  flg_plt = 0

  case uval of
      'LYB': begin
         ;print, 'x_fitdla: This is broken!!  JXP -- 28 Aug 2013'
         ;                       ; Return state to Base
         ;WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
         ;return
         x_fitdla, state.fx, wave=state.wave, FWHM=state.fwhm, $
                   ZIN=state.lines[0].zabs, /BETA, /BLOCK, INFLG=4
      end
      'CRUDEHI': begin
          state.crude_val[1] = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDELO': begin
          state.crude_val[0] = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDE': begin
          state.crude = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'EXACT': begin
          state.exact = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'SPECTRUM': begin
          ;; Save
          fd_all_conti[*,state.cur_spec] = state.conti
          fd_all_cstr[state.cur_spec] = state.cstr
          fd_svxymnx[*,state.cur_spec] = state.xymnx

          ;; Update
          state.cur_spec = ev.value
          x_fitdla_SetState, state, state.cur_spec

          ;; Proceed
          flg_plt = 1
      end
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].zabs = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CONTI' : begin
          widget_control, state.conti_id, get_value=tmp
          state.conti = tmp
      end
      'BVAL' : begin
          widget_control, state.bval_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'NCOLM' : begin
          widget_control, state.Ncolm_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = Zoom
                          x_speczoomreg, state
                          if state.flg_zoom NE 2 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state,  /no_copy
                              return
                          endif else state.flg_zoom = 0
                      end 
                      4 : begin
                          state.lines[state.curlin].zabs = $
                            (state.xpos / state.wrest) - 1.
                          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
                          widget_control, state.zabs_id, $
                            set_value=state.lines[state.curlin].zabs
                      end
                      else: 
                  endcase
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
                  widget_control, state.xpos_id, set_value=state.xpos
                  widget_control, state.ypos_id, set_value=state.ypos
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          case eventch of
             'q': begin
                OUT_CONTI=state.conti[0:state.npix-1]
                widget_control, ev.top, /destroy
                return
             end
             '_': begin
                state.zro_lvl = state.ypos
                flg_plt = 1
             end
             else: x_fitline, state, eventch, FLG_PLT=flg_plt ;; FIT LINE
          endcase
      end
      ;; BUTTONS
      'IDLOUT' : begin
         if fd_nspec GT 1 then begin
             fd_all_conti[*,state.cur_spec] = state.conti
             fd_all_cstr[state.cur_spec] = state.cstr
             fd_svxymnx[*,state.cur_spec] = state.xymnx
             zro_lvl = state.zro_lvl
             save, fd_nspec, fd_all_cstr, fd_all_npix, fd_all_conti, $
                   fd_all_fx, fd_all_sig, fd_all_wave, fd_svxymnx, zro_lvl,  $
                   filename='multi_fitdla.idl', /compres
          endif
         x_fitline_idlout, state, /ERR
      end
      'PSFILE' : x_fitdla_psfile, state
      'SPECPLT' : begin
          if state.nlin NE 0 then zin = state.lines[state.curlin].zabs $
            else zin = 0
          x_specplot, state.fx, state.sig, wave=state.wave, inflg=4, $
            zin=zin, /lls, /block, zout=zout
          ;; Push new zabs and show new fit -- Not working...
          ;state.lines[state.curlin].zabs = zout
          ;widget_control, state.zabs_id, set_value=zout
          ;x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
          ;stop
      end
      'DONE' : begin
         OUT_CONTI=state.conti[0:state.npix-1]
         widget_control, ev.top, /destroy
         return
      end
  endcase

  ;; PXMNX
  if state.xymnx[0] NE state.old_xymnx[0] OR $
    state.xymnx[2] NE state.old_xymnx[2] then begin
      state.old_xymnx = state.xymnx
      state.xpmnx = x_getxpmnx(state)
  endif

; Update Plot
  if (flg_plt MOD 2) EQ 1 then x_fitdla_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_fitdla_UpdatePlot, state
  
; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

;  if state.flg_dum EQ 1 then stop
;  state.flg_dum = 0
;  stop
  if state.psfile NE 1 then begin
     titl = state.title
     thick=2
     ppos = state.pos 
     csz = 1.7
  endif else begin
     thick = 2
     xmrg = [7.5,0.5]
     ymrg = [4.5,0.5]
     csz = 1.4
  endelse


  plot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.fx[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, $
        position=ppos, xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
        xtitle='!17Wavelength', ytitle='Flux', $
        title=titl, $
        background=clr.white, $
        thick=thick, $
        color=clr.black, $
        xmarg=xmrg, ymarg=ymrg, $
        xcharsize=csz, $
        ycharsize=csz

  if file_test('template.txt') EQ 1 then begin
     ; scale qso template to roughly y value of data
     good = where(state.wa_template GT 0)
     wa_template = (1 + state.zqso) * state.wa_template[good]
     ind = where(wa_template GT state.xymnx[0] AND wa_template LT state.xymnx[2])
     if (n_elements(ind) gt 1) then begin
        temp = 0.7 * (state.xymnx[1] + state.xymnx[3])
        temp2 = state.template[good]
        template = state.template[good] * temp / median(temp2[ind])
        oplot, wa_template, template,  color=clr.cyan, thick=0.5
     endif
     ;stop
  endif

  ; Plot Error array
  if state.flg_sig EQ 1 then begin
     oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
            state.sig[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, color=clr.red
  endif
  oplot, [-10000., 10000.], replicate(state.zro_lvl,2), color=clr.black, linestyle=5

  ;; FIT
  if state.nlin NE 0 then begin
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        (state.conti[state.xpmnx[0]:state.xpmnx[1]]-state.zro_lvl)* $
        state.fit[state.xpmnx[0]:state.xpmnx[1]] + state.zro_lvl, color=clr.darkgreen, thick=3
      ;; Mark all lines
      for i=0L,state.nlin-1 do begin
          oplot, replicate( (state.lines[i].zabs+1.)*$
                            state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.blue, linestyle=1
      endfor
      ;; Mark current set
      gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin, ngd)
      for i=0L,ngd-1 do begin
          oplot, replicate( (state.lines[gdlin[i]].zabs+1.)*$
                            state.lines[gdlin[i]].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.red, linestyle=1
      endfor
  endif

; CONTINUUM
  ;; Line
  oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.conti[state.xpmnx[0]:state.xpmnx[1]], color=clr.blue, $
    linestyle=2, thick=3
  ;; Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.orange, symsize=5
  endif

; Crude ERROR
  if state.crude NE 0 and state.nlin NE 0 then begin
     ;; HIGH
     tmp = state.lines[0:state.nlin-1]
      mx = max(tmp.N, imx)
      tmp[imx].N = tmp[imx].N + state.crude_val[1]
      crude_hi = x_allvoigt(state.wave[state.xpmnx[0]:state.xpmnx[1]], tmp, SIGMA=state.FWHM/2.35)
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
             state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
             crude_hi, color=clr.red, linestyle=2
      ;; LOW
      tmp = state.lines[0:state.nlin-1]
      tmp[imx].N = tmp[imx].N + state.crude_val[0]
      crude_lo = x_allvoigt(state.wave[state.xpmnx[0]:state.xpmnx[1]], tmp, SIGMA=state.FWHM/2.35)
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
             crude_lo, color=clr.red, linestyle=2
  endif

end

;;;;;;;;;;;;;;;;;;;;
;  Set
;;;;;;;;;;;;;;;;;;;;

pro x_fitdla_SetState, state, qq

  common x_fitdla_common
  state.npix = fd_all_npix[qq]
  state.wave = fd_all_wave[*,qq]
  state.fx = fd_all_fx[*,qq]
  state.sig = fd_all_sig[*,qq]
  
  state.conti = fd_all_conti[*,qq]
  state.cstr = fd_all_cstr[qq]
  if fd_svxymnx[0,qq] LE 1. then begin 
     ydat = state.fx[0:state.npix-1]
     state.svxymnx = [0., $                                              ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $         ; ymin
                      float(n_elements(ydat)-1), $                       ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))]           ; ymax
     state.svxymnx[0] = min(state.wave[0:state.npix-1])
     state.svxymnx[2] = max(state.wave[0:state.npix-1])
  endif else state.svxymnx = fd_svxymnx[*,qq] 
  
  state.xymnx = state.svxymnx
  state.xpmnx = x_getxpmnx(state)

  state.zqso = fd_zqso[qq]
  if state.nlin GT 0 then $
     x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact

  return

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_fitdla_Reset, state

  ;; Plotting
  state.xymnx = state.svxymnx
  state.old_xymnx = state.svxymnx

end

pro x_fitdla_psfile, state

  common x_fitdla_common

  device, get_decomposed=svdecomp

  x_psopen, state.psfilename, /maxs
  !p.multi=[0,2,2]
  state.psfile = 1

  sv_curspec = state.cur_spec
  for qq=0L,fd_nspec-1 do begin

     ;; State
     x_fitdla_SetState, state, qq
     x_fitdla_UpdatePlot, state
  endfor


  x_psclose
  !p.multi=[0,1,1]

  device, decomposed=svdecomp
  state.psfile = 0
  x_fitdla_SetState, state, sv_curspec

  return
end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro x_fitdla_inifit, state, inifit


  ;; File?
  a = findfile(inifit, count=na)
  if na EQ 0 then begin
      print, 'x_fitdla: FILE ', inifit, ' does not exist!'
      stop
  endif

  ;; Restore
  restore, inifit, /relaxed_structure_assignment

  ;; Error
  if keyword_set( SIG ) then begin
      if sig[0] NE 0. then begin
          state.crude_val = sig
          state.crude = 1
          widget_control, state.crude_id, set_value=1
      endif
  endif

  ;; Continuum
  if keyword_set(conti) then begin
     nct = n_elements(conti)
     state.conti[0:nct-1] = conti 
  endif else state.conti[*] = 1.

  ;; conti_str
  if keyword_set( CSTR ) then state.cstr = cstr

  ;; Lines
  a = where(lines.zabs GT 0., nlin)
  state.nlin = nlin
  tmp = state.lines[0:state.nlin-1]
  copy_struct, lines[0:nlin-1], tmp
  state.lines[0:state.nlin-1] = tmp
  state.curlin = 0L

  ;; Zoom in 
  cen = (1.+state.lines[state.curlin].zabs)*state.wrest
  gd = where(abs(state.wave-cen) LT 200.)
  mx = max(state.conti[gd])
  state.xymnx = [cen-200., -0.1*mx, cen+200., mx*2.]

  ;; Zabs, N
  widget_control, state.Ncolm_id, set_value=state.lines[state.curlin].N
  widget_control, state.bval_id, set_value=state.lines[state.curlin].b
  widget_control, state.zabs_id, set_value=state.lines[state.curlin].zabs

  ;; nset
  state.nset = state.lines[state.nlin-1].set

  ;; 
  x_fitline_updfit, state, EXACT=state.exact

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fitdla, yin, XSIZE=xsize, YSIZE=ysize, TITLE=title, $
              WAVE=wave, BLOCK=block, FWHM=fwhm, INIFIT=inifit, $
              MULTI_IN=multi_in,fil_sig=fil_sig, $
              INFLG=inflg, ZIN=zin, BETA=beta, IN_CONTI=in_conti, $
              SMOOTH=smooth, LIMIT=limit, OUTFIT=outfit, OUTPS=outps, $
              FIN_CONTI=fin_conti, BIN=bin

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_fitdla, fx, XSIZE=,YSIZE=, TITLE=, WAVE=, '
    print, '            INIFIT=, INFLG=, FWHM=, /BLOCK, ZIN=, /BETA) [v1.1]'
    return
 endif 

  ;;  Optional Keywords
  common x_fitdla_common

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( FWHM ) then    fwhm = 2. ;; Pixels

  
  ; Set output names
  if not keyword_set(OUTFIT) then  outfit='fort.idl'
  if not keyword_set(OUTPS) then   outps='idl.ps'


  ;; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0

  if inflg EQ 4 then nspec = 1 else nspec = n_elements(yin)
  x_fitdla_initcommon, nspec, INIFIL=MULTI_IN


  if not keyword_set(MULTI_IN) then begin
     if file_test('./redshifts.txt') then $
        readcol, 'redshifts.txt', names_zqso, zqso, format='A,F'
     for kk=0L,fd_nspec-1 do begin
        if INFLG NE 4 then begin
           print, 'Reading: ', yin[kk]
           ydat = x_readspec(yin[kk], INFLG=inflg, head=head, $
                             NPIX=npix, WAV=xdat, SIG=ysig, fil_sig=fil_sig)
        endif else begin
           ydat = x_readspec(yin, INFLG=inflg, head=head, $
                             NPIX=npix, SIG=ysig, fil_sig=fil_sig)
           xdat = wave
        endelse
        
        if file_test('./redshifts.txt') then begin
                                ; read redshift
           i1 = strpos(yin[kk], '/', /reverse_search)
           name = strmid(yin[kk], i1 + 1, strlen(yin[kk]) - 1)
           ind = where(names_zqso EQ name)
           print, kk, name, ind, zqso[ind]
           fd_zqso[kk] = zqso[ind]
        endif
       
        if keyword_set(bin) then begin
           delta_x=median(xdat-shift(xdat,1))
           delta_x_new=bin*delta_x
           xnew_size=ceil((max(xdat)-min(xdat))/delta_x_new)
           xnew=findgen(xnew_size)*delta_x_new+min(xdat)
           ydat=rebin_spectrum(ydat,xdat,xnew)
           npix = n_elements(ydat)
           if keyword_set(ysig) then ysig=rebin_spectrum(ysig,xdat,xnew)/sqrt(bin)
           xdat=xnew
        endif

 
        ;;if no error
        if n_elements(ysig) eq 0 then ysig=ydat-ydat

        ;; Sort
        srt = sort(xdat)
        xdat = xdat[srt]
        ydat = ydat[srt]
        ysig = ysig[srt]
		;print, ysig
        
        if keyword_set(smooth) then begin
           ydat = smooth(ydat, smooth, /NAN)
           ysig = ysig/sqrt(smooth)
        endif

        ;; Save
        fd_all_npix[kk] = npix
        fd_all_fx[0:npix-1,kk] = ydat
        if mean(ysig) ne 0 then fd_all_sig[0:npix-1,kk] = ysig $
        else fd_all_sig[0:npix-1,kk] = replicate(1., npix) 
        fd_all_wave[0:npix-1,kk] = xdat
     endfor
  endif

  tmp1 = { newabslinstrct }

  tmp2 = { contistrct }

;    STATE
  state = { fx: fltarr(300000L), $
            wave: dblarr(300000L), $
            sig: fltarr(300000L), $
            npix: 0L, $
            cur_spec: 0, $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            wrest: 1215.6701d, $
            beta: 0, $
            nlin: 0, $ ; DLA flag
            EXACT: 1,$
            lines: replicate(tmp1,20), $ ; DLA flag
            fit: fltarr(300000L) + 1., $
            conti: replicate(1.,300000L), $
            cstr: tmp2, $
            zro_lvl: 0., $ ;; Zero level for the spectrum
            curlin: 0L, $
            nset: -1, $
            crude: 0, $  ; Crude Error
            crude_val: [-0.1, 0.1], $
            xpos: 0.d, $
            ypos: 0.d, $
            flg_dum: 1, $
            psfile: 0, $ ; Postscript
            psfilename: outps, $  ; mendez
            outfilename:outfit,$ ; mendez
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            old_xymnx:fltarr(4), $
            tmpxy: fltarr(4), $
            xpmnx: lonarr(2), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            bval_id: 0L, $
            Ncolm_id: 0L, $
            conti_id: 0L, $
            lines_id: 0L, $
            crude_id: 0L, $
            spec_id: 0L, $
            exact_id: 0L, $
            crudehi_id: 0L, $
            crudelo_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L, $
            zqso: 0.0, $
            template: fltarr(20000L),$
            wa_template: fltarr(20000L) $
          }

  ;; WAVE
  if keyword_set(BETA) then begin
      state.wrest = 1025.7223d
      state.beta = 1
  endif

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE

	if mean(ysig) ne 0 then state.flg_sig=1L

  ;; Initialize
  state.fx = fd_all_fx[*, 0]
  state.wave = fd_all_wave[*, 0]
  state.sig = fd_all_sig[*, 0]
  state.npix = fd_all_npix[0]
  state.zqso = fd_zqso[0]
  if keyword_set(IN_CONTI) then state.conti[0:state.npix-1] = IN_CONTI

  if file_test('./template.txt') then begin
     readcol, 'template.txt', wa, template
     for j=0L, n_elements(wa)-1 do begin
        state.wa_template[j] = wa[j]
        state.template[j] = template[j]
     endfor
  endif
  if not keyword_set(MULTI_IN) then begin
     ; Setting initial x & y limits
     ydat = state.fx[0:state.npix-1]
     ysorted = ydat[sort(ydat)]
     y95 = ysorted[long(0.95 * n_elements(ydat))]
     state.svxymnx = [0., $                                ; xmin
                      -0.1 * abs(y95), $                   ; ymin
                      float(n_elements(ydat) - 1), $       ; xmax
                      1.5 * abs(y95)]                      ; ymax

; Set svxymnx[0,2]
     state.svxymnx[0] = min(state.wave[0:state.npix-1])
     state.svxymnx[2] = max(state.wave[0:state.npix-1])
     fd_svxymnx[*,state.cur_spec] = state.svxymnx
  endif else begin
     state.svxymnx = fd_svxymnx[*,0]
     state.xymnx = state.svxymnx
     state.conti = fd_all_conti[*,state.cur_spec]
     state.cstr = fd_all_cstr[state.cur_spec]
  endelse

;    Title
  if size(yin[0], /type) EQ 7 then state.title = yin[0]
  if keyword_set( TITLE ) then state.title = title


;  HEY ALEX MAKE SURE TO CORRECT THIS BEFORE CVS UPDATE
;  CHANGE FOR CVS UPDATE 
; this is same as a character, but I am leaving it in since I like to call it from the cmd line.
;  if keyword_set( LIMIT ) then begin   ; it should be this
;  if (keyword_set( LIMIT )) OR (not keyword_set( LIMIT )) then begin
;	state.svxymnx[1] = -0.2  ; reset the limits on the graph 1/0 error
;	state.svxymnx[3] = 1.5   ; upper limit
;  endif

;    WIDGET
  base = WIDGET_BASE( title = 'x_fitdla', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['x_fitdla', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'T -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = 'A -- Al III doublet'
  state.help[19] = 'C -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'N -- AODM'
  state.help[22] = 'q -- Quit '


;;;;;;;;;
;  Toolbar


  ;; Spectra
  state.spec_id = CW_BGROUP(toolbar, strtrim(lindgen(fd_nspec),2), label_left='Spectrum', $
                              row=4, /exclusive, /no_release, $
                              set_value=0,  uvalue='SPECTRUM')

  ;; zabs
  linbar = WIDGET_BASE( toolbar, /row, /frame, /base_align_center,$
                           /align_center)
  state.zabs_id = cw_field(linbar, title='zabs', value=0., /floating, $
                           /column, /return_events, xsize=10, uvalue='ZABS')
  state.Ncolm_id = cw_field(linbar, title='Ncolm', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='NCOLM')
  state.bval_id = cw_field(linbar, title='bval', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='BVAL')


; continuum
  state.conti_id = cw_field(toolbar, title='conti', value=state.conti, /floating, $
                           /column, xsize=12, /return_events, uvalue='CONTI')


; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Crude error
  crude_err = widget_base(toolbar, /row, /align_center, frame=2)
  state.crude_id = CW_BGROUP(crude_err, ['No', 'Yes'], $
                              row=2, /exclusive, /no_release, $
                              set_value=0,  uvalue='CRUDE')
  crude_val = widget_base(crude_err, /column, /align_center, frame=2)
  state.crudehi_id = CW_FIELD(crude_val, value=state.crude_val[1], xsize=5, $
                              title='Crude High', /floating, $
                              /return_events, uvalue='CRUDEHI')
  state.crudelo_id = CW_FIELD(crude_val, value=state.crude_val[0], xsize=5, $
                              title='Crude Low', /floating, $
                              /return_events, uvalue='CRUDELO')
; Exact
  exbase_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.exact_id = CW_BGROUP(exbase_id, ['No', 'Yes'], label_top='Exact?', $
                              row=2, /exclusive, /no_release, $
                              set_value=state.exact,  uvalue='EXACT')

  
  
;;;;;;;;;;;;
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')

;      Done
  buttbase = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  psfile = WIDGET_BUTTON(buttbase, value='PSFILE',uvalue='PSFILE', $
                         /align_right)
  specplt = WIDGET_BUTTON(buttbase, value='Specplt',uvalue='SPECPLT', $
                          /align_right)
  lyb = WIDGET_BUTTON(buttbase, value='Lyb fit',uvalue='LYB', /align_right)
  idlout = WIDGET_BUTTON(buttbase, value='IDL Out',uvalue='IDLOUT', /align_right)
  buttbase2 = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  done = WIDGET_BUTTON(buttbase2, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
;  if not keyword_set( NOCTB ) then loadct, 2, /silent

; HELP
  i = 0
  state.help[0] = '  :::Help Menu for x_fitdla::: '
  state.help[i++] = 'RMB -- Recenter line'
  state.help[i++] = 's/s -- Set region'
  state.help[i++] = 'lrbt -- Set Left Right Bottom Top'
  state.help[i++] = 'T -- Set ymax to 1.1'
  state.help[i++] = 'zz -- Zoom corners'
  state.help[i++] = 'io -- Zoom in/out'
  state.help[i++] = '{}[] -- Pan'
  state.help[i++] = 'w -- Reset the screen'
  state.help[i++] = 'H -- This widget'
  state.help[i++] = 'c -- Add new Lya line'
  state.help[i++] = 'L -- Add new LLS'
  state.help[i++] = 'B -- Add new Beta'
  state.help[i++] = 'd -- Delete current line'
  state.help[i++] = 'C -- Set continnuum' 
  state.help[i++] = 'nN -- Adjust colm'
  state.help[i++] = 'vV -- Adjust bvalue'
  state.help[i++] = '=- -- Loop through lines'
  state.help[i++] = 'P -- Print to postrscipt'
  state.help[i++] = 'I -- IDL output'
  state.help[i++] = '3,4 -- Add/move continuum point'

; Update
  x_fitdla_Reset, state

  ;; Init Fit and ZIN
  if keyword_set(INIFIT) then x_fitdla_inifit, state, inifit
  if keyword_set(INIFIT) and keyword_set(ZIN) then begin
     print, 'Do not set both INIFIT and ZIN at once!'
     stop
  endif
  if keyword_set( ZIN ) then begin
      state.xpos = (1.+zin)*state.wrest
      x_fitline_newlin, state, BETA=state.beta
  endif

  ;; Set pmnx
  state.xpmnx = x_getxpmnx(state)
  
  ;; Plot
  x_fitdla_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'x_fitdla', base, /no_block $
  else xmanager, 'x_fitdla', base
  
  if not keyword_set(OUT_CONTI) then out_conti = 1
  fin_conti = OUT_CONTI

return
end
