;+ 
; NAME:
; hst_chkciv
;    Version 1.0
;
; PURPOSE:
;   Visually check HST CIV absorption detected with civ_find
;
; CALLING SEQUENCE:
;   
;   hst_chkciv,civfil,outfil
;
; INPUTS:
;   civfil - civcandstct FITS file of candidates
;
; RETURNS:
;
; OUTPUTS:
;   outfil - civcandstrct FITS file of candidates w/ comments
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;   SIGAOD     - plot error bars on AOD profile
;   inspec     - inputing simple structure inspec={wave:wv,flux:fx,error:er}
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hst_chkciv, 'civcand.fits','civcand_chk.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   civ_aodm_mtch()
;
; REVISION HISTORY:
;   Apr-2008     Written by JXP, adapted from SDSS MgII codes by GEP/SHF
;   25-Apr-2008  clean up, debug, add stuff, KLC
;   30-Oct-2009  Remove CIV-specific magic numbers so works for SiIV+, KLC
;-
;------------------------------------------------------------------------------
@civ_aodm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro hst_chkciv_icmmn, civfil, inspec=inspec

  common hst_chkciv_cmm, $
    npix, $
    fx, $
    wv, $
    sig, $
    conti, $
    civstr

  ;; QALSTR
  civstr = xmrdfits(strtrim(civfil,2), 1, /silent)

  if keyword_set(inspec) then begin
     ;; Handle input of single spectrum
     ;; (snglspec=1/state.flg_snglspec=1)
     ;; note: conti not actually used (lines commented out below)
     wv = inspec.wave
     npix = n_elements(wv)
     if npix ne n_elements(inspec.flux) or $
        npix ne n_elements(inspec.error) then $
           stop,'hst_chkciv_icmmn stop: number of elements in inspec.wave/flux/error must match'
     fx = inspec.flux
     sig = inspec.error
  endif

  return
end
  

;;;;
; Events
;;;;

pro hst_chkciv_event, ev

  common hst_chkciv_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'RATE_GROUP': civstr[state.curciv].rating_eye = ev.value
      'COMM_GROUP': begin 
          widget_control, state.comm_id, get_value=tmp
          civstr[state.curciv].comment_eye[0:n_elements(tmp)-1]=tmp
      end
      'NOTES_GROUP': begin 
          widget_control, state.notes_id, get_value=tmp
          civstr[state.curciv].comment_eye[10+lindgen(n_elements(tmp))]=tmp
      end
      'OTHER': begin 
          widget_control, state.other_id, get_value=tmp
          civstr[state.curciv].other_comments = tmp
      end
      'SPLT': x_specplot, fx, sig, wave=wv, zin=state.zabs, /lls, /block, $
;        inflg=4,xrange=1549.5*(1+state.zabs)+[-10,10],$
        inflg=4,xrange=civstr[state.curciv].wrest[0]*(1+state.zabs)+[-10,10],$
                          yrange=state.velplt[0].ymnx
      'SKIP': begin
          hst_chkciv_next, state, /SKIP, KILL=kill
          if not keyword_set(KILL) then hst_chkciv_update, state $
          else begin 
              widget_control, ev.top, /destroy
              return
          endelse
      end
      'NEXT': begin
          hst_chkciv_next, state, KILL=kill
          if not keyword_set(KILL) then hst_chkciv_update, state $
          else begin
              widget_control, ev.top, /destroy
              return
          endelse
       end
      'BACK': begin
         hst_chkciv_next, state, /BACK, KILL=kill
          if not keyword_set(KILL) then hst_chkciv_update, state $
          else begin
              widget_control, ev.top, /destroy
              return
          endelse
      end
      'DONE' : begin
          hst_chkciv_svciv, state
          widget_control, ev.top, /destroy
          return
      end
      'SAVE' : begin
          hst_chkciv_svciv, state
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Metals Plot
;;;;;;;;;;;;;;;;;;;;


pro hst_chkciv_Metals, state
  
  common hst_chkciv_cmm

  ;; Set plot window
  if state.psfile NE 1 then begin
     widget_control, state.draw_id, get_value=wind
     wset, wind
  endif

  clr = getcolor(/load)
  
  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)

  ;; Adjust velocity bounds for optimal zoom
  sfvmnx = state.vmnx           ;must restore at end for next iteration
  mx = max(2.998e5*(civstr[state.curciv].wv_lim[gdlin,0]/$
                    (civstr[state.curciv].wrest[gdlin]*$
                     (1+civstr[state.curciv].zabs[gdlin]))-1.),/absolute)
  state.vmnx = [-abs(mx)-25,abs(mx)+25]
  mx = max(2.998e5*(civstr[state.curciv].wv_lim[gdlin,1]/$
                    (civstr[state.curciv].wrest[gdlin]*$
                     (1+civstr[state.curciv].zabs[gdlin]))-1.),/absolute)
  if abs(mx) gt state.vmnx[1] then state.vmnx = [-abs(mx)-25,abs(mx)+25]
  mx = 600. < max(abs(state.all_velo))
  if state.vmnx[1] gt mx then state.vmnx = [-mx,mx] ;cap
  if state.vmnx[1] lt 150. then state.vmnx = [-150.,150.] ;floor

  ny = ngd + 1                  ; incl AOD profile

  ;; PLOT Good lines
  !p.multi = [0,1,ny,0,1]
  for j=0L,ngd-1 do begin
     i = gdlin[j]
     pixmin = state.all_pmnx[0,i]
     pixmax = state.all_pmnx[1,i]

     ;; Plot
     if (j NE ny-1) AND (j NE ngd-1 ) then begin
        spaces = replicate('!17 ',30)
        plot, state.all_velo[0:state.all_pmnx[2,i],i], $
              fx[pixmin:pixmax], xrange=state.vmnx, $
              yrange=state.velplt[i].ymnx, xtickn=spaces, xmargin=[9,3], $
              ymargin=[0,0], $
              charsize=1.8, psym=10, background=clr.white, color=clr.black, $
              xstyle=1, ystyle=1
     endif else begin
        plot, state.all_velo[0:state.all_pmnx[2,i],i], $
              fx[pixmin:pixmax], xrange=state.vmnx, $
              yrange=state.velplt[i].ymnx, xmargin=[9,3], ymargin=[3,0], $
              charsize=1.8, psym=10, background=clr.white, color=clr.black, $
              xstyle=1, ystyle=1
     endelse

     ;; Display error
     oplot,state.all_velo[0:state.all_pmnx[2,i],i], sig[pixmin:pixmax], $
          psym=10,color=clr.red

     ;; Make darker the region used for measurements
     mt = min(abs(state.velplt[i].wrest - civstr[state.curciv].wrest), imn)
     v1 = 3e5*(civstr[state.curciv].wv_lim[imn,0]/$
               (civstr[state.curciv].wrest[imn]*$
                (1+state.zabs))-1.) ;wrt 1548
     v2 = 3e5*(civstr[state.curciv].wv_lim[imn,1]/$
               (civstr[state.curciv].wrest[imn]*$
                (1+state.zabs))-1.) ;wrt 1548
     gd = where(state.all_velo[0:state.all_pmnx[2,i],i] GE v1 AND $
                state.all_velo[0:state.all_pmnx[2,i],i] LE v2, ngd)

     if ngd NE 0 then begin
        oplot, state.all_velo[gd,i], fx[pixmin+gd], $
               color=clr.black, thick=3, psym=10
        oplot, state.all_velo[gd,i], sig[pixmin+gd], $
               color=clr.red, thick=3, psym=10
     endif 

     ;; Labels
     xyouts, 0.07*(state.vmnx[1]-state.vmnx[0])+state.vmnx[0], $
             state.velplt[i].ymnx[0]+ $
             (state.velplt[i].ymnx[1]-state.velplt[i].ymnx[0])*0.05, $
             strtrim(state.velplt[i].name,2), $
             color=clr.black, charsize=1.5
     
     ;; Lines
     oplot, [0., 0.], state.velplt[i].ymnx, color=clr.blue, linestyle=2
     oplot, [-5000., 5000.], [0.,0.], color=clr.green, linestyle=3
     oplot, [-5000., 5000.], [1.,1.], color=clr.green, linestyle=3
  endfor

  ;; AOD Profile plot
  gd = where(civstr[state.curciv].aodm_civ[*,1] gt 0.,ngd) ;CIV 1550
  mx = max(civstr[state.curciv].aodm_civ[gd,1])
  mn = min(civstr[state.curciv].aodm_civ[gd,1])
  gd = where(civstr[state.curciv].aodm_civ[*,0] gt 0.,ngd) ;CIV 1548
  mn = min(civstr[state.curciv].aodm_civ[gd,0]) < mn
  mx = max(civstr[state.curciv].aodm_civ[gd,0]) > mx
  aodcolmnx = [mn > 10,mx+0.3]

  ;; Tick interval
  if aodcolmnx[1]-aodcolmnx[0] GT 2. then ytint = 0.8 else ytint = 0
  
  ;; Plot (assumes last plot in column)
  plot, state.vmnx, aodcolmnx, xrange=state.vmnx, $
        yrange=aodcolmnx, xmargin=[9,3], ymargin=[3,0], /NODATA, $
        charsize=1.8, psym=10, background=clr.white, $
        color=clr.black, ytickinterval=ytint, $
        xstyle=1, ystyle=1,$
        ytitle='log !8N!X!DAOD!N('+state.ion[0]+'!E+'+state.ion[3]+'!N)'
  oplot,civstr[state.curciv].aodm_vel[*,0],civstr[state.curciv].aodm_civ[*,0],$
        psym=10,color=clr.black
  ;; Offset of vcent(1550) from vcent(1548)
  dv = 2.998e5*(civstr[state.curciv].zabs[1]-$
                civstr[state.curciv].zabs[0])/$
       (1.+civstr[state.curciv].zabs[0])
  oplot,civstr[state.curciv].aodm_vel[*,1]+dv,$
        civstr[state.curciv].aodm_civ[*,1],psym=10,color=clr.red

  ;; CIV 1548
  vmnxciv = 2.998e5*(civstr[state.curciv].wv_lim[0,*]/$
                     (civstr[state.curciv].wrest[0]*$
                      (1+civstr[state.curciv].zabs[0]))-1.)
  mn = min(civstr[state.curciv].aodm_vel[*,0]-vmnxciv[0],imn,/absolute)
  mx = min(civstr[state.curciv].aodm_vel[*,0]-vmnxciv[1],imx,/absolute)
  if keyword_set(state.sigaod) then $
     oploterror,civstr[state.curciv].aodm_vel[imn:imx,0],$
                civstr[state.curciv].aodm_civ[imn:imx,0],$
                civstr[state.curciv].sigaodm_civ[imn:imx,0],$
                psym=10,color=clr.black,nskip=4,errcolor=clr.black
  oplot,civstr[state.curciv].aodm_vel[imn:imx,0],$
        civstr[state.curciv].aodm_civ[imn:imx,0],psym=10,color=clr.black,$
        thick=3

  ;; CIV 1550
  vlim = 2.998e5*(civstr[state.curciv].wv_lim[1,*]/$
                  (civstr[state.curciv].wrest[1]*$
                   (1+civstr[state.curciv].zabs[1]))-1.)
  mn = min(civstr[state.curciv].aodm_vel[*,1]-vlim[0],imn,/absolute)
  mx = min(civstr[state.curciv].aodm_vel[*,1]-vlim[1],imx,/absolute)
  if keyword_set(state.sigaod) then $
     oploterror,civstr[state.curciv].aodm_vel[imn:imx,1]+dv,$
                civstr[state.curciv].aodm_civ[imn:imx,1],$
                civstr[state.curciv].sigaodm_civ[imn:imx,1],$
                psym=10,color=clr.red,nskip=5,errcolor=clr.red
  oplot,civstr[state.curciv].aodm_vel[imn:imx,1]+dv,$
        civstr[state.curciv].aodm_civ[imn:imx,1],psym=10,color=clr.red,$
        thick=3
  vmnxciv = [vmnxciv[0] < vlim[0], vmnxciv[1] > vlim[1]]
  
;  ;; Labels
;  flg_mtch = civ_aodm_mtch(civstr[state.curciv])
;  xyouts,0.05*(state.vmnx[1]-state.vmnx[0])+state.vmnx[0],$
;         0.85*(aodcolmnx[1]-aodcolmnx[0])+aodcolmnx[0],$
;         strtrim(round(flg_mtch*100),2)+'% AOD agree',$
;         color=clr.black,charsize=1.5

  ;; Lines
  oplot, [0., 0.], aodcolmnx, color=clr.blue, linestyle=2 ;vertical

  ;; Label (Lya forest)
  if (civstr[state.curciv].flg_sys[0] and 8) eq 0 then $
     xyouts,0.5,0.98,'Lya Forest',color=clr.red,charsize=1.8,$
            alignment=0.5,/normal

  ;; Reset for next pass
  !p.multi = [0,1,1]
  state.vmnx = sfvmnx
;  widget_control, state.tdraw_id, get_value=wind
;  wset, wind
;  gd = where(wv GT 1400*(1.+state.zabs) AND wv LT (1.+state.zabs)*1700)
;  plot, wv[gd]/(1.+state.zabs), fx[gd]/conti[gd], xrange=[1480., 1620.], $
;    yrange=[-0.2,1.4], xmargin=[9,3], ymargin=[3,0], $
;    charsize=1.8, psym=10, background=clr.white, color=clr.black, $
;    xstyle=1, ystyle=1
;  oplot, wv[gd]/(1+state.zabs), sig[gd]/conti[gd], color=clr.red, psym=10
;  oplot, 1548.19*[1.,1], [-9,9.], color=clr.blue, linestyle=1
;  oplot, 1550.77*[1.,1], [-9,9.], color=clr.blue, linestyle=1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro hst_chkciv_update, state
  hst_chkciv_updinfo, state
  hst_chkciv_Metals, state
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro hst_chkciv_setup, state
  common hst_chkciv_cmm
    
  if not keyword_set(state.flg_snglspec) then begin
     ;; Read data (avoid re-reading if possible)
     file_list = getenv('MLSS_DIR')+'/'+strtrim(civstr[state.curciv].instr_fil,2)
     if not strmatch(file_list, state.sv_instfil) then civ_readspec, file_list, wv, fx, sig
     state.sv_instfil = file_list
  endif

  ;; Set zabs
  state.zabs = civstr[state.curciv].zabs[0] ;; 1548
  state.zqso = civstr[state.curciv].zqso

  ;; Set EW
  state.ew = civstr[state.curciv].ew[0] ;; 1548

  ;; xymnx
;  state.xymnx[0] = 1548.195*(1.+state.zabs) - 200.
;  state.xymnx[2] = 1548.195*(1.+state.zabs) + 200.
  state.xymnx[0] = civstr[state.curciv].wrest[0]*(1.+state.zabs) - 200.
  state.xymnx[2] = civstr[state.curciv].wrest[0]*(1.+state.zabs) + 200.
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  if ngd GT 1 then begin
      srt = sort(fx[gd])
      ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
  endif else ymd = 0.
  state.xymnx[1] = -1.
  state.xymnx[3] = ymd*1.5
  state.svxymnx = state.xymnx

  gd = where(state.velplt.wrest*(state.zabs+1) GT min(wv) AND $
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
pro hst_chkciv_llist, state

  llist = getenv('XIDL_DIR')+'/HST/CIV/hst_civ.lst'
  lines = x_setllst(llist, 0)
  state.ntrans = n_elements(lines)
  state.velplt[0:state.ntrans-1].wrest = lines.wave
  state.velplt[0:state.ntrans-1].name = lines.name
  delvarx, lines
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.39]

  ;; Check that doublet first in CIV structure
  ;; is CIV (leave HI)
  if abs(state.velplt[0].wrest-state.dblt.wvI) gt 1.e-4 then begin
     state.velplt[0:1].wrest = [state.dblt.wvI,state.dblt.wvII]
     state.velplt[0:1].name = state.ion[0] + ' ' + state.ion[1:2]
  endif 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro hst_chkciv_next, state, BACK=back, SKIP=skip, KILL=kill
  common hst_chkciv_cmm

  if keyword_set(back) then begin
     if state.curciv ne 0 then state.curciv = state.curciv - 1 
  endif else state.curciv = state.curciv + 1 
  kill = 0

  ;; Check for the last one
  if state.curciv EQ n_elements(civstr) then begin
      print, 'hst_chkciv:  All done'
      hst_chkciv_svciv, state
      kill = 1
      return
  endif

  if keyword_set(SKIP) then begin
      while(civstr[state.curciv].rating_eye GT 0 $
            AND state.curciv LT n_elements(civstr)) do begin
          state.curciv=state.curciv+1
;          print, state.curciv, civstr[state.curciv].rating_eye
      endwhile
  endif

  ;; Setup
  hst_chkciv_setup, state

  hst_chkciv_setbutt, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hst_chkciv_setbutt, state
  common hst_chkciv_cmm

  ;; Reset buttons
  widget_control, state.rate_id, set_value=civstr[state.curciv].rating_eye
  widget_control, state.comm_id, get_value=tmp
  widget_control, state.comm_id, set_value=$
                  civstr[state.curciv].comment_eye[lindgen(n_elements(tmp))]
  widget_control, state.notes_id, get_value=tmp
  widget_control, state.notes_id, set_value=$
                  civstr[state.curciv].comment_eye[10+lindgen(n_elements(tmp))]

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro hst_chkciv_updinfo, state
  common hst_chkciv_cmm

  ;; Namej
  widget_control, state.name_id, $
    set_value=strtrim(civstr[state.curciv].qso,2)

  ;; RA, DEC
;  widget_control, state.mag_id, set_value=civstr[state.curciv].rmag
  widget_control, state.ra_id, set_value=civstr[state.curciv].ra
  widget_control, state.dec_id, set_value=civstr[state.curciv].dec
  widget_control, state.curciv_id, set_value=state.curciv
  val = strtrim(round(civstr[state.curciv].ew[0]),2)+"+/-"+$
        strtrim(round(civstr[state.curciv].sigew[0]),2)
  widget_control, state.ew1_id, set_value=val
  val = strtrim(round(civstr[state.curciv].ew[1]),2)+"+/-"+$
        strtrim(round(civstr[state.curciv].sigew[1]),2)
  widget_control, state.ew2_id, set_value=val
  widget_control, state.rto_id, set_value=civstr[state.curciv].ew[0]/$
                  civstr[state.curciv].ew[1]
  widget_control, state.sigrto_id, set_value=abs(civstr[state.curciv].ew[0]/$
                  civstr[state.curciv].ew[1])*$
                  sqrt((civstr[state.curciv].sigew[0]/$
                        civstr[state.curciv].ew[0])^2 + $
                       (civstr[state.curciv].sigew[1]/$
                        civstr[state.curciv].ew[1])^2)
  widget_control, state.dv_id, $
                  set_value=2.998e5*(civstr[state.curciv].zabs[1]-$
                                     civstr[state.curciv].zabs[0])/$
                  (1.+civstr[state.curciv].zabs[0])
  widget_control, state.aod_id,$
                  set_value=round(civ_aodm_mtch(civstr[state.curciv])*100.)
  widget_control, state.flg_id,set_value=civstr[state.curciv].flg_sys[0]
  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs[0]

  ;; Comments
  widget_control, state.other_id, set_value=civstr[state.curciv].other_comments

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set MGII structure
pro hst_chkciv_setciv, state

  common hst_chkciv_cmm

  ;; What is this for? (KLC)
;  if civstr[0].zabs NE 0 then begin
;    tmp = { hstmgiistrct }
;    civstr = [civstr, tmp]
; endif

  nciv = n_elements(civstr)
  state.curciv = nciv-1
  state.flg_new = 1

  ;; Name, etc
  civstr[state.curciv].qso_name = qalstr[state.curqso].qso_name
  civstr[state.curciv].ra = qalstr[state.curqso].ra
  civstr[state.curciv].dec = qalstr[state.curqso].dec
  civstr[state.curciv].rmag = qalstr[state.curqso].rmag
  civstr[state.curciv].hst_obs = qalstr[state.curqso].hst_obs
  civstr[state.curciv].z_qso = qalstr[state.curqso].z_qso
  civstr[state.curciv].zabs = qalstr[state.curqso].zabs
  civstr[state.curciv].ew = qalstr[state.curqso].ew
  civstr[state.curciv].sigew = qalstr[state.curqso].sigew
  civstr[state.curciv].wrest = qalstr[state.curqso].wrest

  ;; Flag
  civstr[state.curciv].ew[49] = flg_civ

  return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hst_chkciv_svciv, state

  common hst_chkciv_cmm

  mwrfits, civstr, state.outfil, /create, /silent

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
pro hst_chkciv, civfil, outfil, ICIV=iciv, SIGAOD=sigaod, $
                XSIZE=xsize, YSIZE=ysize, inspec=inspec, _extra=extra

  common hst_chkciv_cmm
;
  if  N_params() LT 2  then begin 
     print,'Syntax - ' + $
           'hst_chkciv, civfil, outfil, [ICIV=, /sigaod'
     print, '        XSIZE=, YSIZE=, inspec=, _extra=] [v1.0]'
     return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then ysize = ssz[1]-100
  if not keyword_set( ICIV ) then iciv = 0L
  if not keyword_set( SIGAOD ) then sigaod = 0 

; Initialize the common blcok
  hst_chkciv_icmmn, civfil, inspec=inspec

  tmp = { velpltstrct }
  tmp2 = { newabslinstrct }
  tmp3 = { doubletstrct }

; STATE

  state = {             $
          ion:strarr(4),$ ; e.g. ['CIV','1548','1550','3']
          dblt:tmp3, $    ; contains some info redundant with state.ion
          outfil: outfil, $
          nciv: n_elements(civstr), $
          curciv: ICIV, $
          sv_instfil: '', $
          flg_snglspec:keyword_set(inspec), $ ; lock spec after hst_chkciv_icmmn
          flg_new: keyword_set(FNEW), $
          zabs: 0., $
          zqso: 0., $
          con_dir: '', $
          ntrans: 0L, $         ; PLOTTING LINES
          vmnx: [-800., 800.], $
          nplt: 0, $
          civ_lin: tmp2, $
          civ_conti: 0., $
          all_velo: dblarr(5000, 300), $  
          all_pmnx: lonarr(3, 300), $  
          velplt: replicate(tmp, 300), $
          xpos: 0.0, $
          ypos: 0.0, $
          ipress: 0L, $
          pos: [0.1,0.1,0.95,0.95], $ ; Plotting
          sigaod: sigaod, $
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
          fxval_id: 0L, $
          iwvval_id: 0L, $
          swvval_id: 0L, $
          zabs_id: 0L, $
          xmax_id: 0L, $
          name_id: 0L, $
          nspec_id: 0L, $
          pmin_id: 0L, $
          pmax_id: 0L, $
          lines_id: 0L, $
          left_id: 0L, $
          draw_id: 0L, $
          right_id: 0L, $
          info_id: 0L, $
          quality_id: 0L, $
          hits_id: 0L, $
          NHI_id: 0L, $
          NHIb_id: 0L, $
          mtl_id: 0L, $
          stat_id: 0L, $
          ra_id: 0L, $
          dec_id: 0L, $
          mag_id: 0L, $
          rate_id: 0L, $
          comm_id: 0L, $
          other_id: 0L, $
          notes_id: 0L, $
          curciv_id: 0L, $
          ew1_id: 0L, $
          ew2_id: 0L, $
          rto_id: 0L, $ 
          sigrto_id: 0L, $
          dv_id: 0L, $
          aod_id:0L, $
          flg_id:0L, $
          help_text_id: 0L, $
          ew: 0L$
  }

;;;;;;;;;;;;;;
; SETUP LINES
  ;; Load Doublet information
  prs = strsplit(civstr[state.curciv].ion[0],/extract)
  state.ion[0:1] = prs
  prs = strsplit(civstr[state.curciv].ion[1],/extract)
  state.ion[2] = prs[1]
  getion, civstr[state.curciv].wrest[0], ion, elm
  state.ion[3] = ion
  state.dblt = dblt_retrieve(state.ion[0])

  hst_chkciv_llist, state
  state.civ_lin = x_setline(1215.6701d)

;    WIDGET
  base = WIDGET_BASE( title = 'hst_chkciv: Check spectra', /row, $
                      UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.left_id = WIDGET_BASE( state.base_id, /row, $
                               /base_align_center,/align_center, $
                               xsize=xsize/2, ysize=ysize, $
                               uvalue='TOP_BASE', frame=2)
  state.draw_id = widget_draw(state.left_id, xsize=xsize/2, $
                              ysize=ysize, /frame, retain=2, $
                              uvalue='DRAW')
  state.info_id = WIDGET_BASE( state.base_id, /column, $
                               /base_align_center,/align_center, $
                               xsize=xsize/2, ysize=ysize, $
                               uvalue='INFO', frame=2)

;;;;;; Info window ;;;;;;;;;;;
  ;; Info
  ;; QSO, RA, DEC, iciv
  radeci = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(radeci, title='Obj: ', value=' ', xsize=18)
  state.ra_id = cw_field(radeci, title='RA: ', value=0., xsize=10)
  state.dec_id = cw_field(radeci, title='DEC: ', value=0., xsize=10)
  state.curciv_id = cw_field(radeci, title='IDX: ', value=state.curciv, $
                             xsize=10)

  civinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.zabs_id = cw_field(civinf, title='zabs: ', value=state.zabs, xsize=7)
  state.flg_id = cw_field(civinf,title='Flag: ', $
                          value=civstr[state.curciv].flg_sys[0],$
                          xsize=3)  
  val = strtrim(round(civstr[state.curciv].ew[0]),2)+"+/-"+$
        strtrim(round(civstr[state.curciv].sigew[0]),2)
  state.ew1_id = cw_field(civinf, title='EW '+state.ion[1]+': ', $
                          value=val, xsize=9)
  val = strtrim(round(civstr[state.curciv].ew[1]),2)+"+/-"+$
        strtrim(round(civstr[state.curciv].sigew[1]),2)
  state.ew2_id = cw_field(civinf, title='EW '+state.ion[2]+': ', $
                          value=val, xsize=9)

  ;; Comparison info
  civinf2 = widget_base(state.info_id, /row, /align_center, frame=2)
  state.rto_id = cw_field(civinf2, title='EW Ratio: ', $
                          value=civstr[state.curciv].ew[0]/$
                          civstr[state.curciv].ew[1], xsize=5)
  state.sigrto_id = cw_field(civinf2, title='Sig(Ratio): ', $
                             value=abs(civstr[state.curciv].ew[0]/$
                             civstr[state.curciv].ew[1])*$
                             sqrt((civstr[state.curciv].sigew[0]/$
                                         civstr[state.curciv].ew[0])^2 + $
                                        (civstr[state.curciv].sigew[1]/$
                                         civstr[state.curciv].ew[1])^2), $
                             xsize=5)
  state.dv_id = cw_field(civinf2,title=state.ion[0]+' dv (km/s): ', $
                         value=2.998e5*(civstr[state.curciv].zabs[1]-$
                                        civstr[state.curciv].zabs[0])/$
                         (1.+civstr[state.curciv].zabs[0]),xsize=5)
  state.aod_id = cw_field(civinf2,title='AOD (%): ', $
                          value=round(civ_aodm_mtch(civstr[state.curciv])*100),$
                          xsize=3)


;  ;; RA, DEC
;  radeci = widget_base(state.info_id, /row, /align_center, frame=2)
;;  state.mag_id = cw_field(radeci, title='MAG: ', value=0., xsize=10)
;  state.ra_id = cw_field(radeci, title='RA: ', value=0., xsize=10)
;  state.dec_id = cw_field(radeci, title='DEC: ', value=0., xsize=10)
;  state.curciv_id = cw_field(radeci, title='IDX: ', value=state.curciv, $
;                             xsize=10)
  

  ;; Lya
  lyainf = widget_base(state.info_id, /row, /align_center, frame=2)

  ;; BUTTONS
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  good = WIDGET_BUTTON(butbase, value='NEXT',uvalue='NEXT')
  back = WIDGET_BUTTON(butbase, value='BACK',uvalue='BACK')
  skip = WIDGET_BUTTON(butbase, value='SKIP', uvalue='SKIP')
  save = WIDGET_BUTTON(butbase, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')

  ;; Rating
  mtlinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.rate_id = cw_bgroup(mtlinf, ['0=??', '1=No', '2','3', '4', '5', '6=Yes!'], $
                            /exclusive, row=7, UVALUE='RATE_GROUP', frame=1, $
                            LABEL_TOP='Rating')

  ;; Comments (primary reason)
  ;; Lya blend = "one or both CIV lines is at least blended with Lya'
  ;; Poor S/N = "Hard to judge because bad data"
  ;; Poor continuum fit = "One or both CIV lines may be due to
  ;;      misplaced continuum (perhaps check continuum"
  ;; Mismatched profiles = "CIV line profiles are inconsistent"
  ;; Wrong EW ratio = "For the strength of the lines, the EW ratio 
  ;;      is wrong (e.g. weak lines should have ratio = 2)"
  state.comm_id = cw_bgroup(mtlinf, ['Lya blend', $
                                     'Poor S/N', $
                                     'Poor continuum fit', $
                                     'Mismatched profiles', $
                                     'Wrong EW ratio', $
                                     'Weak '+state.ion[0], $
                                     'No Lya absorption'], $
                            /nonexclusive, row=6, UVALUE='COMM_GROUP',frame=1, $
                            LABEL_TOP='Comments')

  ;; Description
  ;; Blended = "one or both CIV lines blended"
  ;; Multi-component = "Doublet is multi-component"
  ;; Adjust limits = "There is more or less to these lines and 
  ;;      the limits for inclusion should be changed"
  ;; Repeat = "lines are practically the same as the lines
  ;;      I just saw (due to goofiness of automatic detection"
  state.notes_id = cw_bgroup(mtlinf, ['Blended', $
                                      'Multi-component', $
                                     'Adjust limits',$
                                     'Repeat'], $
                             row=6, UVALUE='NOTES_GROUP', /nonexclusive,$
                            LABEL_TOP=state.ion[0]+' Descript',fram=1)

  mtlinf2 = widget_base(state.info_id, /row, /align_center, frame=2)
  state.other_id = cw_field(mtlinf2, $
                            title='Other Comment: (Remember to Press Return!)', value='', $
                            /column, xsize=85, /return_events, uvalue='OTHER')

  ;; More buttons
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
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
;  if qalstr[0].zabs[0] EQ 0. then hst_chkciv_next, state

                                ; Load data
  hst_chkciv_setup, state
  hst_chkciv_setbutt, state

                                ; PLOT
  hst_chkciv_update, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'hst_chkciv', base
  delvarx, fx, wv, npix, sig, qalstr

  return
end
	
