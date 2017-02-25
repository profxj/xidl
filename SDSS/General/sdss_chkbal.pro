;+ 
; NAME:
; sdss_chkbal
;    Version 1.1
;
; PURPOSE:
;   GUI used to check for (and designate) BAL quasars from SDSS
;
; CALLING SEQUENCE:
;   sdss_chkbal, qalfil, IQSO=, XSIZE=, YSIZE=, /RECHK
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
;   sdss_chkbal, 'DR2_QSO/sdss_DR2_QAL.fits', MATCHFIL='DR2_QSO/matched.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Feb-2004 Written by JXP
;   07-Jun-2004 Revised by SHF
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_chkbal_icmmn, qalfil, MATCHFIL=matchfil, RECHK=rechk, CHKFIL=chkfil, GZFIL=gzfil

  common sdss_chkbal_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    gzstr, $
    qalstr, $
    chkstr, $
    indx, $
    alldr, $
    allra, $
    alldec, $
    nindx

  ;; gzstr
  if keyword_set(GZFIL) then gzstr = xmrdfits(gzfil,1)

  ;; All fitted
  sdss_dlastrct, alldr, /all
  ndr = n_elements(alldr)
  allra = dblarr(ndr)
  alldec = dblarr(ndr)
  for jj=0L,ndr-1 do begin
      if strlen(strtrim(alldr[jj].qso_ra)) EQ 0 then stop
      x_radec, alldr[jj].qso_ra, alldr[jj].qso_dec, rad, decd
      allra[jj] = rad
      alldec[jj] = decd
  endfor

  ;; QALSTR
  if not keyword_set( ZMIN ) then zmin = 2.1
  if not keyword_set( RMAX ) then rmax = 22.5
  if keyword_set( RECHK ) then flg_rechk = 1 else flg_rechk = 0

  if keyword_set( CHKFIL ) then chkstr = xmrdfits(chkfil, 1, /silent)

  if keyword_set( MATCHFIL ) then begin
      match = xmrdfits(getenv('SDSSPATH')+matchfil, 0, /silent)
      qalstr = xmrdfits(qalfil, 1, /silent)
      indx = where(qalstr.z_qso GT zmin AND qalstr.qso_mag LT Rmax AND $
                   qalstr.flg_bal GE flg_rechk, nindx)
  
      print, 'old nindx', nindx
      
      for i=0, nindx-1 do begin
          imatch = where(indx[i] EQ match, nimatch)
          if nimatch GT 0 then indx[i] = -1
      endfor
      
      indx=indx[where(indx NE -1)]
      nindx=n_elements(indx)
      print, 'new nindx', nindx
  
  endif else begin 
  
      qalstr = xmrdfits(qalfil, 1, /silent)
      indx = where(qalstr.z_qso GT zmin AND qalstr.qso_mag LT Rmax AND $
                   qalstr.flg_bal GE flg_rechk, nindx)
      
      print, 'nindx', nindx
  
  endelse
  return
end
  
  
;;;;
; Events
;;;;

pro sdss_chkbal_event, ev

  common sdss_chkbal_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SAVE': sdss_chkbal_svqal, state
      'NEXT': begin
          sdss_chkbal_next, state
          sdss_chkbal_plot, state
      end
      'PREV': begin
          sdss_chkbal_prev, state
          sdss_chkbal_plot, state
      end
;;;;;;;;; Metals ;;;;;;;;;;
      'BALLST' : begin
          widget_control, state.ballst_id, get_value=tmp
          widget_control, state.ngballst_id, get_value=tmp2
          for ii=0L,state.nplt-1 do begin
              if tmp2[ii] EQ 1 then $
                qalstr[indx[state.curqso + ii]].flg_bal = 2 $
              else begin
                  if tmp[ii] EQ 1 then $
                    qalstr[indx[state.curqso + ii]].flg_bal = 1 else $
                    qalstr[indx[state.curqso + ii]].flg_bal = 0
              endelse
          endfor
          sdss_chkbal_plot, state
      end
      'PLOTLST' : begin
          if ev.select EQ 1 then sdss_chkbal_splot, state
      end
      'NGBALLST' : begin
          widget_control, state.ballst_id, get_value=tmp
          widget_control, state.ngballst_id, get_value=tmp2
          for ii=0L,state.nplt-1 do begin
              if tmp2[ii] EQ 1 then $
                qalstr[indx[state.curqso + ii]].flg_bal = 2 $
              else begin
                  if tmp[ii] EQ 1 then $
                    qalstr[indx[state.curqso + ii]].flg_bal = 1 else $
                    qalstr[indx[state.curqso + ii]].flg_bal = 0
              endelse
          endfor
          sdss_chkbal_plot, state
      end
      'DONE' : begin
          sdss_chkbal_svqal, state
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


pro sdss_chkbal_plot, state
  
  common sdss_chkbal_cmm

  ; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  csize = 2.5
  clr = getcolor(/load)
  !p.multi = [0,1,state.nplt,0,1]
  jend = (state.curqso+state.nplt < nindx)-1
  
  for kk=state.curqso, jend do begin
      jj = indx[kk]
      ;; Read data
      parse_sdss, getenv(state.sdsspath) $
        +strtrim(qalstr[jj].file_name,2), fx, wv,$
        SIG=sig, NPIX=npix

      ;; Set zqso
      state.zqso = qalstr[jj].z_qso

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
      if not keyword_set(GZSTR) then $
        swv = qalstr[jj].start_wave / (qalstr[jj].z_qso+1.)  $
      else begin
          pp = where(qalstr[jj].plate EQ gzstr.plate AND $
                     qalstr[jj].fiberid EQ gzstr.fib, npp)
          if npp NE 1 then swv = 1300. else $
            swv = ((gzstr.z1[pp]+1)*1215.6701) / (qalstr[jj].z_qso+1.)  
;          if qalstr[jj].plate EQ 272 and qalstr[jj].fiberid EQ 88 then stop
      endelse
      oplot, [swv,swv], [-9e9,9e9], color=clr.red

      ;; DLA
      if keyword_set( CHKSTR ) then begin
          isame=where(abs(qalstr[jj].ra - chkstr.ra) LT 0.00001 AND $
                      abs(qalstr[jj].dec - chkstr.dec) LT 0.0001 AND $
                      abs(qalstr[jj].z_qso - chkstr.z_qso) LT .01 AND $
                      chkstr.zabs GT 2.2, nisame)
          if nisame NE 0 then begin
              for ss=0L,nisame-1 do begin
                  wvd = (chkstr[isame[ss]].zabs + 1.)*1215.67 / $
                    (qalstr[jj].z_qso+1.)
                  case chkstr[isame[ss]].flg_NHI of
                      0: cd = clr.brown
                      1: cd = clr.orange
                      2: cd = clr.green
                      3: cd = clr.blue
                      else: stop
                  endcase
                  oplot, [wvd,wvd], [-9e9,9e9], color=cd
              endfor
          endif
      endif else begin
          if qalstr[jj].ndla1 NE 0 then begin
              for ss=0L,qalstr[jj].ndla1 - 1 do begin
                  wvd = (qalstr[jj].dla_zabs1[ss] + 1.)*1215.67 / $
                    (qalstr[jj].z_qso+1.)
                  oplot, [wvd,wvd], [-9e9,9e9], color=clr.red
              endfor
          endif
      endelse

      ;; Fitted?
;      if qalstr[jj].plate EQ 569 and qalstr[jj].fiberid EQ 23 then stop
      isame=where(abs(qalstr[jj].ra - allra) LT 0.0001 AND $
                  abs(qalstr[jj].dec - alldec) LT 0.0001 AND $
                  abs(qalstr[jj].z_qso - alldr.qso_zem) LT .01 AND $
                  alldr.zabs GT 2.2, nisame)
      if nisame NE 0 then begin
          for ss=0L,nisame-1 do begin
              wvd = (alldr[isame[ss]].zabs + 1.)*1215.67 / $
                    (alldr[isame[ss]].qso_zem+1.)
              if alldr[isame[ss]].NHI GE 20.3 then cd = clr.green $
              else cd = clr.orange
              oplot, [wvd,wvd], [-9e9,9e9], color=cd, linest=2, thick=5
          endfor
      endif
          
      ;; Error
      oplot, wv, sig, psym=10, color=clr.orange

      ;; Labels
      xyouts, state.xymnx[0]+20., state.xymnx[3]*0.85, qalstr[jj].qso_name, $
        color=clr.red, charsize=1.5
      case qalstr[jj].flg_bal of
          1: begin
              xyouts, state.xymnx[0]+20., state.xymnx[3]*0.7, 'BAL', $
                color=clr.cyan, charsize=1.5
              xyouts, state.xymnx[2]-100., state.xymnx[3]*0.7, 'BAL', $
                color=clr.cyan, charsize=1.5
          end
          2: begin
              xyouts, state.xymnx[0]+20., state.xymnx[3]*0.7, 'NGBAL', $
                color=clr.brown, charsize=1.5
              xyouts, state.xymnx[2]-100., state.xymnx[3]*0.7, 'NGBAL', $
                color=clr.brown, charsize=1.5
          end
          else:
      endcase
  endfor 


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chkbal_next, state

  common sdss_chkbal_cmm

  state.curqso = state.curqso + state.nplt
  widget_control, state.indx_id, set_value=state.curqso
  if (state.curqso + state.nplt) GE state.nqso then begin
      state.curqso = state.nqso - state.nplt 
      print, 'You have reached the end!!! '
      spc = strarr(state.nplt)
      spc[*] = '  '
      printcol, state.curqso+lindgen(state.nplt), spc, $
        qalstr[state.curqso+lindgen(state.nplt)].qso_name
  endif
  sdss_chkbal_resetbut, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Specplot
pro sdss_chkbal_splot, state

  common sdss_chkbal_cmm

  !p.multi = [0,1,1,0,0]
  widget_control, state.plotlst_id, get_value=tmp
  sdss_objinf, [qalstr[indx[state.curqso+tmp]].plate,$
                qalstr[indx[state.curqso+tmp]].fiberid], /plot
  !p.multi = [0,1,state.nplt,0,1]
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro sdss_chkbal_prev, state

  state.curqso = state.curqso - state.nplt
  widget_control, state.indx_id, set_value=state.curqso
  if state.curqso LT 0 then state.curqso = 0
  sdss_chkbal_resetbut, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkbal_resetbut, state
  common sdss_chkbal_cmm

  tmpv = lindgen(state.nplt)
  tmpv2 = lindgen(state.nplt)
  for ii=0L,state.nplt-1 do begin
      tmpv[ii] = (qalstr[indx[state.curqso+ii]].flg_bal GE 1)
      tmpv2[ii] = (qalstr[indx[state.curqso+ii]].flg_bal GE 2)
  endfor
  widget_control, state.ballst_id, set_value=tmpv
  widget_control, state.ngballst_id, set_value=tmpv2
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkbal_svqal, state

  common sdss_chkbal_cmm

  widget_control, /hourglass
  
  mwrfits, qalstr, state.qalfil, /create
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

pro sdss_chkbal, qalfil, MATCHFIL=matchfil, IQSO=iqso, GZFIL=gzfil, $
                 XSIZE=xsize, YSIZE=ysize, RECHK=rechk, CHKFIL=chkfil

  common sdss_chkbal_cmm
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'sdss_chkbal, qalfil, IQSO=, XSIZE=, YSIZE=, /RECHK, CHKFIL= [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-100
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
  if not keyword_set( IQSO ) then iqso = 0L
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'

; Initialize the common blcok
  sdss_chkbal_icmmn, getenv('SDSSPATH')+qalfil, RECHK=rechk, $
    MATCHFIL=matchfil, CHKFIL=chkfil, GZFIL=gzfil

; STATE
;sdss_goz, 'DR5_QSO/tst_goz.fits', 'DR5_QSO/out_tstgz.fits', snrlmt=4.
  state = {             $
            nqso: 0L, $
            curqso: iqso, $
            qalfil: getenv('SDSSPATH')+qalfil, $
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
  state.nqso = nindx
  state.xymnx[0] = 850.
  state.xymnx[2] = 1600.
  state.svxymnx = state.xymnx

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_chkbal: Check spectra', /row, $
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
                              label_top='BALLST', $
                              row=state.nplt, UVALUE='BALLST', /frame, $
                             /return_index, /nonexclusive)
  state.ngballst_id = cw_bgroup(butbase, sindgen(state.nplt), $
                                label_top='NGBAL', $
                                row=state.nplt, UVALUE='NGBALLST', /frame, $
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
  sdss_chkbal_plot, state
  sdss_chkbal_resetbut, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chkbal', base
  delvarx, fx, wv, npix, sig, qalstr
  !p.multi = [0,1,1]

  return
end
	
