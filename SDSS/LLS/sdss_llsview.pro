;+ 
; NAME:
; sdss_llsview
;    Version 1.1
;
; PURPOSE:
;   GUI used to check for (and designate) BAL quasars from SDSS
;
; CALLING SEQUENCE:
;   sdss_llsview, qalfil, IQSO=, XSIZE=, YSIZE=, /RECHK
;
; INPUTS:
;  qalfil  -- Name of FITS file containing the QAL structure
;
; RETURNS:
;
; OUTPUTS:
;  qalfil --  With the BAL info modified accordingly
;
; OPTIONAL KEYWORDS:
;   IQSO  -- Number of the first QSO to start with [default: 0L]
;  /RECHK -- Only plot thos QSOs identified as BAL previously
;   MATCHFIL -- to skip QSOs that have already been checked 
;               (see sdss_balupdate.pro)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_llsview, 'DR2_QSO/sdss_DR2_QAL.fits', MATCHFIL='DR2_QSO/matched.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Feb-2004 Written by JXP
;   07-Jun-2004 Revised by SHF
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_llsview_icmmn, gzfil, qsofil, MATCHFIL=matchfil, RECHK=rechk, $
  CHKFIL=chkfil, LLS=lls, REBUT=rebut

  common sdss_llsview_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    gzstr, $
    chkstr, $
    indx, $
    alldr, $
    allra, $
    alldec, $
    qsostr, $
    nindx

  ;; gzstr
  gzstr = xmrdfits(gzfil,1)
  qsostr = xmrdfits(qsofil, 1)

  ;; All fitted
  alldr = lls
  ndr = n_elements(alldr)
  x_radec, alldr.qso_ra, alldr.qso_dec, rad, decd
  allra = rad
  alldec = decd

  ;; QALSTR
  if not keyword_set( ZMIN ) then zmin = 2.1
  if not keyword_set( RMAX ) then rmax = 22.5
  if keyword_set( RECHK ) then flg_rechk = 1 else flg_rechk = 0

  a = findfile(chkfil+'*',count=nfil)
  if nfil NE 0 then chkstr = xmrdfits(chkfil, /silent) 

  if nfil EQ 0 or keyword_set(REBUT) then begin
      chkstr = replicate(0L, n_elements(gzstr.z1))
      for jj=0L,ndr-1 do begin
          isame=where(abs(gzstr.ra - allra[jj]) LT 0.001 AND $
                      abs(gzstr.dec - alldec[jj]) LT 0.001 AND $
                      abs(gzstr.zem - alldr[jj].qso_zem) LT .05 $
                      , nisame)
          if nisame GT 0 then chkstr[isame] = 1
      endfor
  endif

  return
end
  
  
;;;;
; Events
;;;;

pro sdss_llsview_event, ev

  common sdss_llsview_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SAVE': sdss_llsview_svqal, state
      'NEXT': begin
          sdss_llsview_next, state
          sdss_llsview_plot, state
      end
      'PREV': begin
          sdss_llsview_prev, state
          sdss_llsview_plot, state
      end
;;;;;;;;; Metals ;;;;;;;;;;
      'BALLST' : begin
          widget_control, state.ballst_id, get_value=tmp
          for ii=0L,state.nplt-1 do begin
              if tmp[ii] EQ 1 then  chkstr[state.curqso + ii] = 2 
          endfor
          sdss_llsview_plot, state
      end
      'PLOTLST' : begin
          if ev.select EQ 1 then sdss_llsview_splot, state
      end
      'DONE' : begin
          sdss_llsview_svqal, state
          widget_control, ev.top, /destroy
          return
      end
      'DNNSV' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase


  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
; Plot
;;;;;;;;;;;;;;;;;;;;


pro sdss_llsview_plot, state
  
  common sdss_llsview_cmm

  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  csize = 2.5
  clr = getcolor(/load)
  !p.multi = [0,1,state.nplt,0,1]
  jend = (state.curqso+state.nplt < state.nqso)-1
  
  for kk=state.curqso, jend do begin
      jj = kk
      sdss_objinf, [gzstr.plate[jj], gzstr.fib[jj]], SDSS=QSOSTR, FILNM=fil
      ;; Read data
      parse_sdss, fil, fx, wv, SIG=sig, NPIX=npix
      fx = smooth(fx,3)

      ;; Set zqso
      state.zqso = gzstr.zem[jj]

      ;; xymnx
      wv = wv / (1.+state.zqso)
      if wv[0] LT 1100 then begin
          px = where(abs(wv-1215.) LT 50)
          state.xymnx[3] = median(fx[px])*1.5
      endif else state.xymnx[3] = 5.*median(fx)
      
      ;; Plot
      if (kk NE jend) then begin
          spaces = replicate('!17 ',30)
          plot, wv, fx, $
            xrange=[state.xymnx[0],state.xymnx[2]], $            
            yrange=[state.xymnx[1],state.xymnx[3]], $
            charsize=csize, psym=10, background=clr.white, color=clr.black, $
            xmargin=[7,1], ymargin=[0,0], xstyle=1, ystyle=1, $
            xtickn=spaces
      endif else begin
          plot, wv, fx, $
            xrange=[state.xymnx[0],state.xymnx[2]], $
            yrange=[state.xymnx[1],state.xymnx[3]], $
            charsize=csize, psym=10, background=clr.white, color=clr.black, $
            xtitle='Wavelength', xmargin=[7,1], ymargin=[4,0], xstyle=1, ystyle=1
      endelse
      
      ;; Lines
      swv = ((gzstr.z1[jj]+1)*911.7) / (gzstr.zem[jj]+1.)  
      oplot, [swv,swv], [-9e9,9e9], color=clr.red

      ;; QSO limit
      swv = 911.7
      oplot, [swv,swv], [-9e9,9e9], color=clr.green

      oplot, [-9e9,9e9], [0.,0.], color=clr.gray, linesty=2

      ;; LLS
      isame=where(abs(gzstr.ra[jj] - allra) LT 0.001 AND $
                  abs(gzstr.dec[jj] - alldec) LT 0.001 AND $
                  abs(gzstr.zem[jj] - alldr.qso_zem) LT .05 $ ;; Allows for modified zem
                  , nisame)
;      if gzstr.plate[jj] EQ 1958 and gzstr.fib[jj] EQ 449 then stop
      if nisame NE 0 then begin
          for ss=0L,nisame-1 do begin
              wvd = (alldr[isame[ss]].zabs + 1.)*911.4 / $
                    (gzstr.zem[jj]+1.)
              cd = clr.blue
              oplot, [wvd,wvd], [-9e9,9e9], color=cd
              wvd = (alldr[isame[ss]].zabs + 1.)*1215.67 / $
                    (gzstr.zem[jj]+1.)
              cd = clr.blue
              oplot, [wvd,wvd], [-9e9,9e9], color=cd, linesty=2
          endfor
      endif

      ;; Error
      oplot, wv, sig, psym=10, color=clr.orange

      ;; Labels
      xyouts, state.xymnx[0]+20., state.xymnx[3]*0.85, $
              string(gzstr.plate[jj],format='(i4)')+' '+$
              string(gzstr.fib[jj],format='(i4)'), color=clr.red, charsize=1.5
  endfor 


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_llsview_next, state

  common sdss_llsview_cmm

  state.curqso = state.curqso + state.nplt
  widget_control, state.indx_id, set_value=state.curqso
  if (state.curqso + state.nplt) GE state.nqso then begin
      state.curqso = state.nqso - state.nplt 
      print, 'You have reached the end!!! '
      spc = strarr(state.nplt)
      spc[*] = '  '
;      printcol, state.curqso+lindgen(state.nplt), spc, $
  endif
  sdss_llsview_resetbut, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Specplot
pro sdss_llsview_splot, state

  common sdss_llsview_cmm

  !p.multi = [0,1,1,0,0]
  widget_control, state.plotlst_id, get_value=tmp
  sdss_objinf, [gzstr.plate[indx[state.curqso+tmp]],$
                gzstr.fib[indx[state.curqso+tmp]]], /plot
  !p.multi = [0,1,state.nplt,0,1]
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro sdss_llsview_prev, state

  state.curqso = state.curqso - state.nplt
  widget_control, state.indx_id, set_value=state.curqso
  if state.curqso LT 0 then state.curqso = 0
  sdss_llsview_resetbut, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_llsview_resetbut, state
  common sdss_llsview_cmm

  tmpv = lindgen(state.nplt)
  for ii=0L,state.nplt-1 do begin
      tmpv[ii] = (chkstr[state.curqso+ii] GE 1)
  endfor
  widget_control, state.ballst_id, set_value=tmpv
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_llsview_svqal, state

  common sdss_llsview_cmm

  widget_control, /hourglass
  
  mwrfits, chkstr, state.chkfil, /create
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
;  lls_struct, lls, '~/LLS/Lists/dr5_lls.lst', ROOT=getenv('SDSSPATH'), /NOSYS
; sdss_llsview, 'DR5_QSO/LLS/dr5_lls_goz_sn2.fits', lls,
; '~/SDSS/DR5_QSO/dr5_qso.fits', CHKFIL='chk_lls.fits', IQSO=100
pro sdss_llsview, gzfil, lls, qsofil,  IQSO=iqso, REBUT=rebut, $
                 XSIZE=xsize, YSIZE=ysize, CHKFIL=chkfil

  common sdss_llsview_cmm
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'sdss_llsview, qalfil, IQSO=, XSIZE=, YSIZE=, /RECHK, CHKFIL=, /REBUT [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-100
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
  if not keyword_set( IQSO ) then iqso = 0L
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'

; Initialize the common blcok
  sdss_llsview_icmmn, getenv('SDSSPATH')+gzfil, qsofil, RECHK=rechk, $
    CHKFIL=chkfil, LLS=lls, REBUT=rebut

  if not keyword_set( CHKFIL ) then chkfil = 'chk_lls.fits'

; STATE
;sdss_goz, 'DR5_QSO/tst_goz.fits', 'DR5_QSO/out_tstgz.fits', snrlmt=4.
  state = {             $
            nqso: 0L, $
            curqso: iqso, $
            chkfil: chkfil, $
            zqso: 0., $
            ntrans: 0L, $       ; PLOTTING LINES
            nplt: 5, $
            psfile: 0, $
            flg_rechk: 0, $
            sdsspath: sdsspath, $            
            help: strarr(50), $
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            size: lonarr(2), $
            base_id: 0L, $      ; Widgets
            draw_id: 0L, $    
            indx_id: 0L, $    
            ballst_id: 0L, $       ; BAL List
            plotlst_id: 0L, $       ; Plot List
            ngballst_id: 0L, $       ; BAL List
            help_text_id: 0L $
          }

  if keyword_set( RECHK ) then state.flg_rechk = 1

; Other setup
  state.nqso = n_elements(gzstr.z1)
  state.xymnx[0] = 700.
  state.xymnx[1] = -1.
  state.xymnx[2] = 1300.
  state.svxymnx = state.xymnx

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_llsview: Check spectra', /row, $
                    UNAME='BASE', /tlb_size_events)
  state.base_id = base
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  LLS DRAW
  drawbase = $
    WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center, $
               /tracking_events, uvalue='LDRAW_BASE', frame=2)

  state.draw_id = widget_draw(drawbase, xsize=xsize-100, $
                              ysize=ysize, /frame, uvalue='DRAW')

;      BUTTONS
  butbase = widget_base(state.base_id, /column, /align_center, frame=2)
  state.indx_id = cw_field(butbase, title='I: ', value=0., xsize=7)
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  dnsv = WIDGET_BUTTON(butbase, value='DNNSV',uvalue='DNNSV')
  state.ballst_id = cw_bgroup(butbase, sindgen(state.nplt), $
                              label_top='LLS?', $
                              row=state.nplt, UVALUE='BALLST', /frame, $
                             /return_index, /nonexclusive)
  state.plotlst_id = cw_bgroup(butbase, sindgen(state.nplt), $
                                label_top='PLOT', $
                                row=state.nplt, UVALUE='PLOTLST', /frame, $
                                /return_index, /exclusive)


  widget_control, state.indx_id, set_value=state.curqso


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ; PLOT
  sdss_llsview_plot, state
  sdss_llsview_resetbut, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_llsview', base
  delvarx, fx, wv, npix, sig
  !p.multi = [0,1,1]

  return
end
	
