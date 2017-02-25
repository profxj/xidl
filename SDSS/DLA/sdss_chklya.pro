;+ 
; NAME:
; sdss_chklya
;    Version 1.1
;
; PURPOSE:
;   GUI used to check for missing DLA
;
; CALLING SEQUENCE:
;   sdss_chklya, qalfil, IQSO=, XSIZE=, YSIZE=, /RECHK, ZMIN=
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
;   ZMIN= -- Minimum redshift [default: 2.]
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
;   sdss_chklya, '/DR2_QSO/sdss_DR2_QAL.fits', MATCHFIL='DR2_QSO/matched.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Feb-2004 Written by JXP
;   07-Jun-2004 Revised by SHF
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_chklya_icmmn, qalfil, chkfil, RECHK=rechk, ZMIN=zmin, GZZ=gzz, $
                       RMIN=rmin

  common sdss_chklya_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    qalstr, $
    chkstr, $
    indx, $
    nindx

  if not keyword_set( ZMIN ) then zmin = 2.1
  if not keyword_set( RMIN ) then rmin = 0.
  if not keyword_set( RMAX ) then rmax = 21.5
  if keyword_set( RECHK ) then flg_rechk = 1 else flg_rechk = 0

  ;; CHKD
  chkstr = xmrdfits(chkfil, 1, /silent)

  ;; QALSTR
  qalstr = xmrdfits(qalfil, 1, /silent)
  indx = where(qalstr.z_qso GT zmin AND qalstr.qso_mag LT Rmax AND $
               qalstr.flg_bal NE 2 and qalstr.qso_mag GT Rmin, nindx)

  if keyword_set( GZZ ) then begin
      subi = where(((qalstr[indx].start_wave / 1215.67) - 1) LT gzz[1] AND $
                   qalstr[indx].z_qso GT gzz[0])
      indx = indx[subi]
  endif

  print, 'nqso = ', n_elements(indx)
  
  return
end
  
  
;;;;
; Events
;;;;

pro sdss_chklya_event, ev

  common sdss_chklya_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'INDX': begin
          widget_control, state.indx_id, get_value=tmp
          while(tmp[0] GT state.curqso) do begin
              sdss_chklya_next, state
          endwhile
          sdss_chklya_plot, state
      end
          
      'SAVE': sdss_chklya_svqal, state
      'NEXT': begin
          sdss_chklya_next, state
          sdss_chklya_plot, state
      end
      'PREV': begin
          sdss_chklya_prev, state
          sdss_chklya_plot, state
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
          sdss_chklya_plot, state
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
          sdss_chklya_plot, state
      end
      'DONE' : begin
          sdss_chklya_svqal, state
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


pro sdss_chklya_plot, state
  
  common sdss_chklya_cmm

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
      wv = (wv / 1215.67) - 1.

      state.xymnx[1] = -0.5
      state.xymnx[0] = state.gzz[0]
      state.xymnx[2] = state.gzz[1]
;      state.xymnx[0] = ( (qalstr[jj].start_wave / 1215.67)  - 1. ) > state.gzz[0]
;      state.xymnx[2] = state.zqso < state.gzz[1]

      pix = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], npix)
      if npix NE 0 then state.xymnx[3] = 2.5*median(fx[pix]) else state.xymnx[3] = 1.
      
      ;; Plot
      plot, wv, fx, $
        xrange=[state.xymnx[0],state.xymnx[2]], $
        yrange=[state.xymnx[1],state.xymnx[3]], $
        charsize=csize, psym=10, background=clr.white, color=clr.black, $
        xmargin=[7,1], ymargin=[4,1], xstyle=1, ystyle=1
      axis, xrange=[(state.xymnx[0]+1)*1215.67, $
                    (state.xymnx[2]+1)*1215.67], xaxis=1, charsize=csize, $
        color=clr.black, xstyle=1
      
      ;; DLA
      if qalstr[jj].ndla1 NE 0 then begin
          for ss=0L,qalstr[jj].ndla1 - 1 do begin
              wvd = qalstr[jj].dla_zabs1[ss] 
              if wvd LE 0. then continue
              mtch = where(abs(chkstr.ra - qalstr[jj].ra) LT 0.001 AND $
                           abs(chkstr.dec - qalstr[jj].dec) LT 0.001 AND $
                           abs(chkstr.zabs-wvd) LT 0.02, nmtch)
              if nmtch NE 1 then stop
              case chkstr[mtch].flg_NHI of
                  0: dclr = clr.red
                  1: dclr = clr.green
                  else: dclr = clr.blue
              endcase
              oplot, [wvd,wvd], [-9e9,9e9], color=dclr
          endfor
      endif

      ;; Error
      oplot, wv, sig, psym=10, color=clr.orange
      oplot, [-99,99.], [0.,0.], color=clr.cyan, linestyle=2
      oplot, [state.zqso, state.zqso], [-99., 1e9], color=clr.red, linestyle=1
      strtz = ( (qalstr[jj].start_wave / 1215.67)  - 1. ) > state.gzz[0]
      oplot, [strtz,strtz], [-99., 1e9], color=clr.black, linestyle=2

      ;; Labels
      xyouts, state.xymnx[0]+0.03, state.xymnx[3]*0.85, $
        string(kk,format='(i4)')+' '+qalstr[jj].qso_name, $
        color=clr.red, charsize=1.5
      case qalstr[jj].flg_bal of
          1: begin
              xyouts, state.xymnx[0]+.03, state.xymnx[3]*0.7, 'BAL', $
                color=clr.cyan, charsize=1.5
          end
          2: begin
              xyouts, state.xymnx[0]+.03, state.xymnx[3]*0.7, 'NGBAL', $
                color=clr.brown, charsize=1.5
          end
          else:
      endcase
  endfor 


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chklya_next, state

  common sdss_chklya_cmm

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
  sdss_chklya_resetbut, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Prev
pro sdss_chklya_prev, state

  state.curqso = state.curqso - state.nplt
  widget_control, state.indx_id, set_value=state.curqso
  if state.curqso LT 0 then state.curqso = 0
  sdss_chklya_resetbut, state

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklya_resetbut, state
  common sdss_chklya_cmm

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
pro sdss_chklya_svqal, state

  common sdss_chklya_cmm

;  mwrfits, qalstr, state.qalfil, /create
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

pro sdss_chklya, qalfil, chkfil, MATCHFIL=matchfil, IQSO=iqso, XSIZE=xsize, $
                 YSIZE=ysize, ZMIN=zmin, GZZ=gzz, RMIN=rmin

  common sdss_chklya_cmm
;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_chklya, qalfil, chkfil, IQSO=, XSIZE=, YSIZE=, ZMIN=, RMIN= ' + $
      'GZZ= [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-100
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
  if not keyword_set( IQSO ) then iqso = 0L
  if not keyword_set( ZMIN ) then zmin = 2.
  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'

; Initialize the common blcok
  sdss_chklya_icmmn, getenv('SDSSPATH')+qalfil, getenv('SDSSPATH')+chkfil, $
    RECHK=rechk, ZMIN=zmin, GZZ=gzz, RMIN=rmin

  if not keyword_set( GZZ ) then gzz = [0., 10.]


; STATE

  state = {             $
            nqso: 0L, $
            curqso: iqso, $
            qalfil: getenv('SDSSPATH')+qalfil, $
            zqso: 0., $
            gzz: gzz, $
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
  base = WIDGET_BASE( title = 'sdss_chklya: Check spectra', /row, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
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
  state.indx_id = cw_field(butbase, title='I: ', value=0., xsize=7, $
                          /return_events, UVALUE='INDX')
  save = WIDGET_BUTTON(butbase, value='SAVE',uvalue='SAVE')
  prev = WIDGET_BUTTON(butbase, value='PREV',uvalue='PREV')
  next = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')
  dnsv = WIDGET_BUTTON(butbase, value='DNNSV',uvalue='DNNSV')
;  state.ballst_id = cw_bgroup(butbase, sindgen(state.nplt), $
  state.ballst_id = cw_bgroup(butbase, sindgen(state.nplt), $
                              row=state.nplt, UVALUE='BALLST', /frame, $
                             /return_index, /nonexclusive)
  state.ngballst_id = cw_bgroup(butbase, sindgen(state.nplt), $
                              row=state.nplt, UVALUE='NGBALLST', /frame, $
                             /return_index, /nonexclusive)

  widget_control, state.indx_id, set_value=state.curqso


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Realize
  WIDGET_CONTROL, base, /realize

  ; PLOT
  sdss_chklya_plot, state
  sdss_chklya_resetbut, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chklya', base
  delvarx, fx, wv, npix, sig, qalstr

  return
end
	
