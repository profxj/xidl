;+ 
; NAME:
; sdss_chksiew
;    Version 1.0
;
; PURPOSE:
;   Visually check SDSS CIV absorption detected with sdss_fndciv
;
; CALLING SEQUENCE:
;   
;   sdss_chksiew, 
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
;   sdss_chksiew, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Jun-2007 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_chksiew_icmmn, dla, si2fil, ROOT=root

  common sdss_chksiew_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    fit, $
    dlastr, $
    si2str

  strct = { $
          plate: 0L, $
          fiber: 0L, $
          RA: 0., $
          DEC: 0., $
          zabs: 0., $
          NHI: 0., $
          EW: 0., $
          sigEW: 0., $
          flgEW: 0 }

  ;; DLA
  dlastr = dla
  ndla = n_elements(dlastr)
  
  ;;
  a = findfile(si2fil, count=nfil)
  if nfil NE 0 then begin
      print, 'sdss_chksiew: Using ', si2fil, ' as input.'
      si2str = xmrdfits(si2fil,1,/sile)
  endif else begin
      si2str = replicate(strct, ndla)
      si2str.plate = dlastr.sdss_plate
      si2str.fiber = dlastr.sdss_fibid
      x_radec, dlastr.qso_ra, dlastr.qso_dec, rad, decd
      si2str.ra = rad
      si2str.dec = decd
      si2str.zabs = dlastr.zabs
      si2str.NHI = dlastr.NHI
  endelse

  return
end
  

;;;;
; Events
;;;;

pro sdss_chksiew_event, ev

  common sdss_chksiew_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SPLT': x_specplot, fx, sig, wave=wv, zin=state.zabs, /lls, /block, $
        inflg=4
      'GOOD': begin
          si2str[state.cursi2].flgEW = 1
          sdss_chksiew_next, state
      end
      'LOWER': begin
          si2str[state.cursi2].flgEW = 2
          sdss_chksiew_next, state
      end
      'UPPER': begin
          si2str[state.cursi2].flgEW = 3
          sdss_chksiew_next, state
      end
      'FIX': begin
          si2str[state.cursi2].flgEW = 99
          sdss_chksiew_next, state
      end
      'NG': begin
          si2str[state.cursi2].flgEW = -1
          sdss_chksiew_next, state
      end
      'DONE' : begin
          sdss_chksiew_svsi2, state
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          sdss_chksiew_svsi2, state
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_chksiew_Metals, state
  
  common sdss_chksiew_cmm

  ;; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.mdraw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)
  
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)

;  ny = ngd / 2 + (ngd MOD 2 EQ 1)
  ny = ngd 

  ;; Good lines
;  !p.multi = [0,2,ny,0,1]
  !p.multi = [0,1,ny,0,1]
  for j=0L,ngd-1 do begin
      i = gdlin[j]
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]

      ;; Plot
;      if (j NE ny-1) AND (j NE ngd-1 ) then begin
      if (j NE ny-1) then begin
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

  ;; SiII Plot
  !p.multi = [0,1,1]
  widget_control, state.si2_id, get_value=wind
  wset, wind
  wcen = 1526.7066 * (1+state.zabs)
  mn = min(abs(wv-wcen),cpix)
  ppx = lindgen(51) + cpix - 25
  plot, wv[ppx], fx[ppx], color=clr.black, background=clr.white, psym=10
  oplot, wv[ppx], conti[ppx], color=clr.red, linesty=1

  oplot, wv[ppx], (1-[replicate(0.,18),fit,replicate(0.,18)])* $
         conti[ppx], color=clr.blue

  oplot, replicate(wcen,2), [-9e9,9e9], color=clr.black, linestyle=2, thick=3

  ;; CIV
  contam = [1548.195, 1550.770, 2796.352, 2803.531]
  ncon = n_elements(contam)
  for ss=0L,(ncon/2)-1,2 do begin
      ;; Lower
      zc = wcen / contam[ss] - 1
      w2 = (1+zc)*contam[ss+1]
      oplot, replicate(w2,2), [-9e9,9e9], color=clr.orange, linestyle=2, thick=3
      ;; Upper
      zc = wcen / contam[ss+1] - 1
      w2 = (1+zc)*contam[ss]
      oplot, replicate(w2,2), [-9e9,9e9], color=clr.orange, linestyle=2, thick=3
  endfor

  ;; EW
  xyouts, wcen+3., (min(fx[ppx])-0.5)>0., 'EW: '+$
          string(si2str[state.cursi2].ew,format='(f7.3)')+' '+$
          string(si2str[state.cursi2].sigew,format='(f5.3)'), color=clr.black, $
          charsiz=2.

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro sdss_chksiew_update, state, FLG=flg
  flg = sdss_chksiew_setup( state )
  if flg EQ 1 then begin
      sdss_chksiew_fitsi2, state
      sdss_chksiew_updinfo, state
      sdss_chksiew_Metals, state
  endif
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
function sdss_chksiew_setup, state
  common sdss_chksiew_cmm
    
  state.zabs = dlastr[state.curdla].zabs

  ;; Read in data
  plate = dlastr[state.curdla].sdss_plate
  fibid = dlastr[state.curdla].sdss_fibid
  sdss_objinf, [plate, fibid], filnm=datfil, FLGSDSS=flgsdss
  if FLGSDSS EQ 0 then begin
      print, 'You will need to get this one yourself (not DR5)!'
      si2str[state.cursi2].flgEW = 99
      return, 0
  endif

  ;; Read
  parse_sdss, datfil, fx, wv, conti, sig=sig, ZQSO=zem, CDIR=state.con_dir, $
              NPIX=npix
  state.zqso = zem

  ;; xymnx
;  state.xymnx[0] = 1215.6701*(1.+state.zabs) - 200.
;  state.xymnx[2] = 1215.6701*(1.+state.zabs) + 200.
;  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
;  if ngd GT 1 then begin
;      srt = sort(fx[gd])
;      ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
;  endif else ymd = 0.
;  state.xymnx[1] = -1.
;  state.xymnx[3] = ymd*1.5
;  state.svxymnx = state.xymnx

  ;; Set flg
  gd = where(state.velplt.wrest*(state.zabs+1) GT min(wv) AND $
             state.velplt.wrest*(state.zabs+1) GT (state.zqso+1.)*1215.6701 AND $
             state.velplt.wrest*(state.zabs+1) LT max(wv), na)
  if na EQ 0 then return,0
  state.velplt[*].flg = 0
  state.nplt = na
  state.velplt[gd].flg = 1

  ;; Just do the plots
  state.all_velo[*,gd] = x_allvelo(wv, state.zabs, $
                                   state.velplt[gd].wrest,$
                                   state.vmnx, $
                                   all_pmnx=all_pmnx, NPIX=5000L)
  state.all_pmnx[*,gd] = all_pmnx

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set Lines
pro sdss_chksiew_llist, state

  llist = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civ.lst'
  lines = x_setllst(llist, 0)
  ;; Taking only the first 10
  state.ntrans = 10
  state.velplt[0:state.ntrans-1].wrest = lines[0:state.ntrans-1].wave
  state.velplt[0:state.ntrans-1].name = lines[0:state.ntrans-1].name
  delvarx, lines
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.39]

  weak = where(abs(state.velplt.wrest - 1808.0130d) LT 0.01 OR $
               abs(state.velplt.wrest - 1611.2005d) LT 0.01 OR $
               abs(state.velplt.wrest - 2026.136d) LT 0.01 OR $
               abs(state.velplt.wrest - 2260.7805d) LT 0.01, nwk)
  
  if nwk NE 0 then state.velplt[weak].ymnx = [0.7, 1.1]
               

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chksiew_next, state
  common sdss_chksiew_cmm

  jump1: state.curdla = state.curdla + 1 

  if state.curdla EQ n_elements(dlastr) then begin
    sdss_chksiew_svsi2, state
    widget_control, ev.top, /destroy
    return
  endif

  ;; Check for repeat
  if state.flg_redo NE 1 then begin
      ;; Match?
      mt = where(dlastr[state.curdla].sdss_plate EQ si2str.plate AND $
                 dlastr[state.curdla].sdss_fibid EQ si2str.fiber AND $
                 abs(dlastr[state.curdla].zabs-si2str.zabs) LT 0.001 AND $
                 si2str.flgew NE 0, nmt)
      case nmt of 
          0: begin
              ;; If si2str is archived, this DLA may be new.  If so we
              ;; need to add it in
              mt = where(dlastr[state.curdla].sdss_plate EQ si2str.plate AND $
                         dlastr[state.curdla].sdss_fibid EQ si2str.fiber AND $
                         abs(dlastr[state.curdla].zabs-si2str.zabs) LT 0.001, $
                         nmt2)
              if nmt2 EQ 0 then begin
                  tmp = si2str[0]
                  tmp.plate = dlastr[state.curdla].sdss_plate
                  tmp.fiber = dlastr[state.curdla].sdss_fibid
                  x_radec, dlastr[state.curdla].qso_ra, $
                           dlastr[state.curdla].qso_dec, rad, decd
                  tmp.ra = rad
                  tmp.dec = decd
                  tmp.zabs = dlastr[state.curdla].zabs
                  tmp.NHI = dlastr[state.curdla].NHI
                  si2str = [si2str,tmp]
                  ;;
                  state.cursi2 = n_elements(si2str)-1
              endif else state.cursi2 = mt[0]
              flgup = 1
          end
          1: begin
              print, 'Skipping ', dlastr[state.curdla].qso, $
                     dlastr[state.curdla].zabs
              goto, jump1
          end
          else: stop
      endcase
  endif else flgup = 1

  sdss_chksiew_update, state, FLG=flg
  if flg EQ 0 then sdss_chksiew_next, state
  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chksiew_updinfo, state
  common sdss_chksiew_cmm

  ;; Name
  widget_control, state.name_id, $
    set_value=strtrim(dlastr[state.curdla].qso,2)

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs
  widget_control, state.ew_id, $
                  set_value=string(si2str[state.cursi2].ew,format='(f7.3)')+' '+$
                  string(si2str[state.cursi2].sigew,format='(f5.3)')

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chksiew_svsi2, state

  common sdss_chksiew_cmm

  mwrfits, si2str, state.si2fil, /create
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chksiew_fitsi2, state

  common sdss_chksiew_cmm


  rwv = 1526.7066d
  mn = min(abs(wv-(state.zabs+1)*rwv),cpix)
  px = lindgen(15) + cpix - 7
  subpx = lindgen(9) + cpix - 4

  if cpix EQ (npix-1) then begin
      print, 'Redshift too high! ', state.zabs
      return
  endif

  ;; Profile
  gprof = -1.*(fx[px] / conti[px] - 1.)

  wcen = wv[cpix]
  dwv = abs(wv[cpix]-wv[cpix+1])
  gsssig = 2.

  ;; FIT
  fit = x_gaussfit(wv[px], gprof, acoeff, $
                    estimates=[max(gprof), wcen, gsssig], $
                    sigma=sigma, nterms=3, COVAR=covar, $
                    measure_errors=sig[px]/conti[px])

  ;; Save EW
  ewval = acoeff[0] * acoeff[2] * sqrt(!pi*2.) ; A
  sigew1 = sqrt(total( ((sig[subpx]/conti[subpx])*dwv)^2)) ; A
  sigew2 = sqrt(!pi*2.) * sqrt( (acoeff[0]*sigma[2])^2 + $
                                (acoeff[2]*sigma[0])^2 )
;      sigew = sigew1 
;      sigew = sigew2 
;      sigew = sigew1 > sigew2
  sigew = sigew1 < sigew2

  ;; Error checking (for a bad fit)
  ew_chk = total((1. - fx[px]/conti[px])*abs(dwv[px])) ;; Obs Ang
  if abs(ew_chk-ewval) GT 5.*sigew then begin
      print, 'Bad EW value!', ew_chk, ewval
      print, 'Probably a bad Gaussian fit'
      print, 'Setting EW to the lower value', ew_chk, ewval
      ewval = (ewval < ew_chk) > 0.
  endif

  si2str[state.cursi2].ew = ewval / (1.+state.zabs)
  si2str[state.cursi2].sigew = sigew / (1.+state.zabs)

;  if ewval LT -10 then stop
  print, si2str[state.cursi2].ew,  si2str[state.cursi2].sigew 
      
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
; sdss_chksiew, 'sdss_DR3_CIV.fits', 'civ_good.fits', 'civ_bad.fits',
; getenv('SDSSPATH')+'DR3_QSO/ABSLIN/'
pro sdss_chksiew, dla, si2fil, IQSO=iqso, FNEW=fnew, $
              XSIZE=xsize, YSIZE=ysize, ROOT=root, REDO=redo

  common sdss_chksiew_cmm
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_chksiew, dla, si2fil, /FNEW I_YSIZE=, S_YSIZE= [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(CDIR) then cdir = getenv('SDSSPATH')+'DR5_QSO/ABSLIN/'
  if not keyword_set( PATH ) then path = getenv('SDSSPATH')+'/DR5_QSO/' 
  if not keyword_set( SUMMF ) then summf = path+'dr5_qso.fits'
  if not keyword_set( DATDIR ) then datdir = path+'spectro/1d_23/'

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
;  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = -1L

; Initialize the common blcok
  sdss_chksiew_icmmn, dla, si2fil, ROOT=root

  tmp = { velpltstrct }
  tmp2 = { abslinstrct }

; STATE

  state = {             $
          nqal: 0L, $
          curdla: iqso, $
          cursi2: 0L, $
          zabs: 0., $
          zqso: 0., $
          si2fil: si2fil, $
          con_dir: cdir, $
          datdir: datdir, $
          ntrans: 0L, $         ; PLOTTING LINES
          vmnx: [-800., 800.], $
          nplt: 0, $
          flg_redo: keyword_set(redo), $
          civ_lin: tmp2, $
          civ_conti: 0., $
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
          base_id: 0L, $        ; Widgets
          ldraw_id: 0L, $       ; Lya
          ltext_id: 0L, $       ; 
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
          lines_id: 0L, $
          si2_id: 0L, $
          tdraw_id: 0L, $
          left_id: 0L, $
          right_id: 0L, $
          rhs_id: 0L, $
          info_id: 0L, $
          quality_id: 0L, $
          ew_id: 0L, $
          help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES
  sdss_chksiew_llist, state

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_chksiew: Check spectra', /row, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.left_id = WIDGET_BASE( state.base_id, /column, $
                              /base_align_center,/align_center, $
                              xsize=round(3*xsize/4), $
                              ysize=ysize, $
                              uvalue='TOP_BASE', frame=2)
  state.si2_id = widget_draw(state.left_id, xsize=round(3*xsize/4.), $
                              ysize=round(4*ysize/5), /frame, retain=2, $
                              uvalue='SI2DRAW')
  state.right_id = WIDGET_BASE( state.base_id, /row, $
                              /base_align_center,/align_center, $
                              xsize=round(xsize/4.), ysize=ysize, $
                              uvalue='TOP_BASE', frame=2)
  state.mdraw_id = widget_draw(state.right_id, xsize=round(xsize/4.), $
                               ysize=ysize, /frame, retain=2, $
                               uvalue='MDRAW')

;;;;;; Info window ;;;;;;;;;;;
  state.info_id = $
    WIDGET_BASE( state.left_id, /column, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, xsize=xsize/3.)
;               ysize=round(ysize/3.))
  ;; Info
  civinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(civinf, title='Obj ', value=' ', xsize=18)
  state.zabs_id = cw_field(civinf, title='zabs: ', value=state.zabs, xsize=7)
  state.ew_id = cw_field(civinf, title='EW: ', value=' ', xsize=13)

;      BUTTONS
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  good = WIDGET_BUTTON(butbase, value='GOOD',uvalue='GOOD')
  maybe = WIDGET_BUTTON(butbase, value='UPPER',uvalue='UPPER')
  bad  = WIDGET_BUTTON(butbase, value='LOWER', uvalue='LOWER')
  bad  = WIDGET_BUTTON(butbase, value='FIX', uvalue='FIX')
  bad  = WIDGET_BUTTON(butbase, value='NG', uvalue='NG')
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  splt = WIDGET_BUTTON(butbase2, value='SPLT',uvalue='SPLT')


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
;  if qalstr[0].zabs[0] EQ 0. then sdss_chksiew_next, state

  ;; PLOT
  sdss_chksiew_next, state
  sdss_chksiew_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chksiew', base
  delvarx, fx, wv, npix, sig, qalstr

  return
end
	
