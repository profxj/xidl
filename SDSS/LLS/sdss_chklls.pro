;+ 
; NAME:
; sdss_chklls
;    Version 1.0
;
; PURPOSE:
;   Visually check SDSS LLS candidates
;
; CALLING SEQUENCE:
;   
;   sdss_chklls, x, maskid, expsr, XSIZE=, YSIZE=
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
;   sdss_chklls, x, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Dec-2005 Written by JXP/JMO
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;

pro sdss_chklls_icmmn, drpath, qso_fil, szstart, szend, templ_fil, llsfil, $
  MOCK_FIL=mock_fil, AMOD_FIL=amod_fil, DATFILS=datfils, GZFIL=gzfil

  common sdss_chklls_cmm, $
    mock_data, $
    npix, $
    fx, $
    wv, $
    sig, $
    nqsos, $
    llssearch, $
    allmod, $
    allnh, $
    zllsarr, $
    template, $
    lls_fstrct, $
    llsstate, $
    gzstr

  ;; GZFIL
  if keyword_set(GZFIL) then gzstr = xmrdfits(gzfil,1)


  ;; QALSTR
  if strlen(strtrim(qso_fil,2)) GT 0 then begin
      dr=xmrdfits(drpath+qso_fil,1,/silent)
      range = where(dr.z GE szstart AND dr.z LE szend,nqsos)
  endif 

  ;; Mocks
  if keyword_set(MOCK_FIL) then begin
      mock_data = xmrdfits(mock_fil,1)
      nqsos = n_elements(mock_data)
  endif

  ;; Existing file?
  if x_chkfil(llsfil) EQ 1 then begin ;; YES
      llssearch = xmrdfits(llsfil, 1, /silent) 
      if nqsos NE n_elements(llssearch) then begin 
          print, 'sdss_chklls: NQSOs in that z range does not ' + $
                 'match the input file.'
          print, 'sdss_chklls: continue only if you are ' + $
                 'doing a subset or MOCKs'
          stop 
      endif
  endif else begin  ;; NO
      tmp = {sdssllsstrct}
      llssearch = replicate(tmp, nqsos)
      llssearch.plate=dr[range].plate
      llssearch.mjd=dr[range].mjd
      llssearch.fiber=dr[range].fiberid
      llssearch.zem=dr[range].z
      llssearch.gmag=dr[range].psf_g
      llssearch.umag=dr[range].psf_u
      llssearch.uerr=dr[range].psf_su

      if not keyword_set(DATFILS) then begin
          ;;create the file name 
          for qq=0L,nqsos-1 do begin
              IF(llssearch[qq].plate LT 1000) THEN BEGIN
                  splate=strcompress('0'+string(llssearch[qq].plate),/remove_all)
              ENDIF ELSE BEGIN
                  splate=strcompress(string(llssearch[qq].plate))
              ENDELSE
              
              smjd=string(llssearch[qq].mjd, format='(i5.5)')
              
              sfiber=string(llssearch[qq].fiber, format='(i3.3)')
              
              filename = strcompress('spSpec-'+smjd+'-'+splate+'-'+sfiber+'.fit.gz', $
                                     /remove_all)
              filedir = strcompress(drpath+'spectro/1d_26/'+splate+'/1d/', $
                                    /remove_all)
              file=filedir+filename
              llssearch[qq].datfil = file
          endfor
      endif else begin
          readcol, datfils, files, format='A'
          llssearch.datfil = strtrim(files,2)
      endelse
      
  endelse
  
  nnhi = 20
  nh1 = xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/nhi16_19b30.fits',0)
  npnh = n_elements(nh1)
  allnh = fltarr(npnh,nnhi)
  for i=0L,nnhi-1 do allnh[*,i] = xmrdfits(getenv('XIDL_DIR')+ $
                                           '/SDSS/LLS/nhi16_19b30.fits',i)
  
  
  ;; Create model image
  nztest = long((alog10(911.7633*6) - alog10(911.7833*4.2))/0.0001) + 11
  zllsarr = (10^(findgen(nztest)*0.0001 + alog10(4.2*911.7633)))/911.7633 - 1.0
  wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
  mn = min(abs(wv_mod-911.7633),mmn)
  
  wv_sdss = 10^(3.58 + dindgen(3900)*1e-4)
  mn = min(abs(4.2*911.7633-wv_sdss),imn)
  shft = imn-mmn

  if keyword_set(AMOD_FIL) then begin
      if x_chkfil(amod_fil) EQ 1 then begin
          if keyword_set(ALLMOD) then delvarx, allmod
          restore, amod_fil 
      endif else begin
          allmod = fltarr(npnh, nztest, nnhi)
          for rr=0l,nnhi-1 do begin
              for qq=0L,nztest-1 do begin
                  allmod[*,qq,rr] = shift(allnh[*,rr], shft+qq)
                  allmod[(npnh+qq+shft)<(npnh-1):*,qq,rr] = 1.
              endfor
          endfor
          save, allmod, filename=amod_fil
      endelse
  endif
  ;; Template
  if strlen(strtrim(templ_fil,2)) GT 0 then $
    template=xmrdfits(templ_fil,/silent)

  return
end
  

;;;;
; Events
;;;;

pro sdss_chklls_event, ev

  common sdss_chklls_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
;      'ZABS' : begin
;          widget_control, state.zabs_id, get_value=tmp
;          state.zabs = tmp
;          sdss_chklls_updplt, state
;      end
      'SPLT': sdss_chklls_specplot, state
      'FIDDLE': sdss_chklls_fiddle, state
      'SEARCH': begin
          sdss_chklls_search, state
          sdss_chklls_updinfo, state
          sdss_chklls_updplt, state
      end
      'DONE' : begin
          sdss_chklls_updstr, state
          sdss_chklls_svlls, state
          widget_control, ev.top, /destroy
          return
      end
      'NOSV' : begin
          widget_control, ev.top, /destroy
          return
      end
      'SMTH' : begin
          state.smooth = state.smooth +1
          sdss_chklls_updplt, state
      end
      'NOSMTH' : begin
          state.smooth = 0
          sdss_chklls_updplt, state
      end
      'NEXT' : begin
          widget_control, /hourglass
          sdss_chklls_updstr, state
          state.curqso = state.curqso + 1
          if state.curqso GE n_elements(llssearch) then begin
              state.curqso = state.curqso - 1
              sdss_chklls_updstr, state
              sdss_chklls_svlls, state
              widget_control, ev.top, /destroy
              return
          endif
          sdss_chklls_setup, state
          sdss_chklls_updinfo, state
          sdss_chklls_updplt, state
       end
      'PREV' : begin
         sdss_chklls_updstr, state
         if state.curqso GT 0 then begin
            state.curqso = state.curqso -1
            sdss_chklls_setup, state
            sdss_chklls_updinfo, state
            sdss_chklls_updplt, state
         endif else begin
            state.curqso = state.curqso
            sdss_chklls_setup, state
            sdss_chklls_updinfo, state
            sdss_chklls_updplt, state
         endelse
      end
      'SAVE' : begin
          sdss_chklls_updstr, state
          sdss_chklls_svlls, state
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklls_updstr, state
  common sdss_chklls_cmm

;  widget_control, state.nhi_id, get_value=v_nhi
;  llssearch[state.curqso].flg_nhi = v_nhi

  widget_control, state.clls_id, get_value=v_lls
  llssearch[state.curqso].flg_lls = v_lls

  widget_control, state.tlls_id, get_value=v_lls
  llssearch[state.curqso].flg_tweak = v_lls

  widget_control, state.plls_id, get_value=v_lls
  llssearch[state.curqso].flg_extra = v_lls


;  widget_control, state.zabs_id, get_value=z_lls
;  llssearch[state.curqso].zlls = z_lls

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklls_updplt, state
  common sdss_chklls_cmm

  clr = getcolor(/load)
  ;; Lya
  widget_control, state.lya_id, get_value=wind
  wset, wind

  plot,lls_fstrct.wv, lls_fstrct.spec, $
       xrange=(1.0+state.zabs)*1215.67+[-50.,50],xstyle=1, $
       charsize=0.75,ymin=0,ystyle=1, color=clr.black, background=clr.white,$
       psym=10
  oplot,[(1.0+state.zabs)*1215.67,(1.0+state.zabs)*1215.67],$
        [0,1000.],color=clr.green

  ;; LLS
  widget_control, state.lls_id, get_value=wind
  wset, wind
  mx = max(lls_fstrct.savmod)

  plotwavlls=(1.0 + state.zabs)*912.0
  plot,lls_fstrct.wv,lls_fstrct.speci, $ 
    xrange=[plotwavlls-100.,plotwavlls+100],xstyle=1, $
    charsize=0.75,ystyle=1, color=clr.black, background=clr.white, $
    yrange=3.*mx*[-0.3,1], psym=10
  oplot,lls_fstrct.wv,lls_fstrct.savmod,color=clr.red
  oplot,[-1e9,1e9], [0.,0.],color=clr.green,linestyle=2

  smthf = median(lls_fstrct.speci,30)
  sub = where(lls_fstrct.wv LT plotwavlls, nsub)
  if nsub LE 1 then stop
  oplot,lls_fstrct.wv[sub],smthf[sub], color=clr.orange, linesty=1, thick=3

  ;; CHI
  widget_control, state.chi_id, get_value=wind
  wset, wind

  plot,lls_fstrct.zlls,lls_fstrct.chi2, $
    yrange=[0.95*min(lls_fstrct.chi2),1.05*max(lls_fstrct.chi2)],ystyle=1
  oplot,[lls_fstrct.zem,lls_fstrct.zem], $
    [0.95*min(lls_fstrct.chi2),1.05*max(lls_fstrct.chi2)],linestyle=2
  oplot,[state.zabs,state.zabs],[-1e9,1e9],linestyle=3,color=clr.red
  xyouts,state.zabs+0.001,0.9*max(lls_fstrct.chi2),'LLS',color=clr.red
  xyouts,lls_fstrct.zem+0.001,0.9*max(lls_fstrct.chi2),'zem'

  ;; Full spec
  widget_control, state.fspec_id, get_value=wind
  wset, wind

  plot,lls_fstrct.wv,smooth(lls_fstrct.spec,2*state.smooth+1), $ 
    color=clr.black, background=clr.white, $
    xmin=lls_fstrct.xmin,xstyle=1,yrange=[-0.5,max(lls_fstrct.model)],ystyle=1
  oplot,lls_fstrct.wv, lls_fstrct.model, color=clr.red
  oplot,lls_fstrct.wv, lls_fstrct.sig, color=clr.green, linesty=1

  ;; GZ
  if state.flg_gz EQ 1 then begin
      mt = where(llssearch[state.curqso].plate EQ gzstr.plate AND $
                 llssearch[state.curqso].fiber EQ gzstr.fib, nmt)
      case nmt of 
          0: print, 'QSO not in statistical sample'
          1: oplot, replicate((gzstr.z1[mt]+1)*912.,2), $
                 [-1e9,1e9], color=clr.purple, linesty=2
          else: stop
      endcase
  endif

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; Setup data
pro sdss_chklls_setup, state
  common sdss_chklls_cmm
    
  qq = state.curqso
  ;;
  widget_control, state.plls_id, set_value=(llssearch[qq].flg_extra >0)
  widget_control, state.clls_id, set_value=(llssearch[qq].flg_lls > 0)
  widget_control, state.tlls_id, set_value=(llssearch[qq].flg_tweak)

  print,'Now doing qso '+llssearch[qq].datfil+' (#'+ $
        strcompress(string(qq),/remove_all)+' of '+ $
        strcompress(string(nqsos-1),/remove_all)+')'

  ;; Search if need be
  if llssearch[qq].flg_lls EQ 0 OR state.flg_mock then sdss_chklls_search, state $
  else begin  ;; Create the model
      ;; Data
      sdss_objinf, [llssearch[state.curqso].plate, $
                    llssearch[state.curqso].fiber], filnm=fil
      parse_sdss, fil, fx, wave, sig=sig, NPIX=npix
      specwav0 = alog10(wave[0])
      logwav = findgen(npix)*0.0001 + specwav0

      ;; Continuum
      templ_spec = sdss_qsotempl(template, 2.8, logwav, fx, $
                                 sig, llssearch[qq].zem, EPX=epx, FPX=fpx)
      if llssearch[qq].conti_scale LT 1e-4 then llssearch[qq].conti_scale = 1.
      templ_spec = templ_spec * llssearch[qq].conti_scale
      ;; Model
      model = allnh[*,llssearch[qq].taulls]
      
      npnh = n_elements(model)
      wv_mod = 10^(2.7d + dindgen(npnh)*1e-4)
      mn = min(abs(wv_mod-911.7633),mmn)
      
      mn = min(abs((1+llssearch[qq].zlls)*911.7633-wave),imn)  
      
      i1 = mmn-imn
      i2 = (npnh-1) < (i1+npix-1) 
      fit = model[i1:i2]

      ;; Pad as need be
      if n_elements(fit) LT n_elements(templ_spec) then begin
          fit = [fit, replicate(1., n_elements(templ_spec)-n_elements(fit))]
      endif

      lls_fstrct = { $
                 zlls: findgen(100), $
                 chi2: fltarr(100), $
                 zem: llssearch[qq].zem, $
                 model: fit*templ_spec, $
                 wv: wave[0:epx-fpx], $
                 spec: fx[0:epx-fpx], $
                 sig: sig[0:epx-fpx], $
                 xmin: specwav0, $
                 speci: fx[0:epx-fpx]/templ_spec, $
                 savmod: fit[0:epx-fpx] $
               }
  endelse

  ;; Set zabs
  state.zabs = llssearch[qq].zlls

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chklls_search, state
  common sdss_chklls_cmm

  qq = state.curqso

  ;; Grab template if necessary
  if keyword_set(MOCK_DATA) then begin
      zem = fix(llssearch[qq].zem*10)
      ;; Grab the template file
      temp_files = findfile(state.mock_path+'mock_templ*fits', $
                            count=nfil)
;      temp_files = findfile(getenv('LLSPAP')+ $
;                            '/SDSS/DR7/Analysis/Mocks/mock_templ*fits', $
      if nfil EQ 0 then begin
          print, 'sdss_chklls: Should not get here! Continue at your own risk'
          stop
          return
      endif
      flg_tmpl = 0
      tfil = ''
      for tt=0L,nfil-1 do begin
          ;; Parse
          prs = strsplit(temp_files[tt],'_',/extrac)
          np = n_elements(prs)
          z2 = long(strmid(prs[np-1],0,2))
          z1 = long(prs[np-2])
          if zem LT z2 and zem GE z1 then tfil = temp_files[tt]
      endfor
      if strlen(tfil) EQ 0 then stop
      template = xmrdfits(tfil,/silen)
      mock_tmp = mock_data[qq]
  endif

  ;; Read data
  gdmod = where(zllsarr LT llssearch[qq].zem+0.1)
  llsstate=sdss_findlls(template,2.8,llssearch[qq].datfil,$
                        allmod=allmod[*,gdmod,*], FSTRCT=lls_fstrct, $
                        ZEM=llssearch[qq].zem, $
                        MOCK_DATA=mock_tmp)

  llssearch[qq].llsflg=llsstate.llsflg
  llssearch[qq].zlls=llsstate.zlls
  llssearch[qq].taulls=llsstate.taulls
  llssearch[qq].blls=llsstate.blls
  llssearch[qq].zstart=llsstate.zstart
  llssearch[qq].zend=llsstate.zend
  llssearch[qq].uminusg=llsstate.uminusg
  llssearch[qq].svchi=llsstate.svchi
  llssearch[qq].svz=llsstate.svz

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chklls_updinfo, state
  common sdss_chklls_cmm

  ;; Namej
  widget_control, state.tau_id, $
    set_value=(16. + llssearch[state.curqso].taulls*0.2)

  ;; RA, DEC
  widget_control, state.idx_id, set_value=state.curqso
  widget_control, state.plate_id, set_value=llssearch[state.curqso].plate
  widget_control, state.fib_id, set_value=llssearch[state.curqso].fiber

  ;; Magnitudes
  widget_control, state.gmag_id, $
                  set_value=(llssearch[state.curqso].umag- $
                             llssearch[state.curqso].gmag)
  widget_control, state.umag_id, set_value=llssearch[state.curqso].umag
  widget_control, state.uerr_id, set_value=llssearch[state.curqso].uerr

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs

  return

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklls_svlls, state

  common sdss_chklls_cmm

  mwrfits, llssearch, state.llsfil, /create
;  mwrfits, badlls, state.badllsfil, /create
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklls_fiddle, state

  common sdss_chklls_cmm
  tmp = llssearch[state.curqso]
;          stop
  ;; Mock?
  if state.flg_mock then mock_tmp = mock_data[state.curqso]
  sdss_llsfit, [tmp.plate, tmp.fiber], $
               tmp, template, OUTSTR=tmp2, MOCK_DATA=mock_tmp

  ;; Need to save continuum (renorm) here too
  llssearch[state.curqso] = tmp2
  state.zabs = tmp2.zlls
  widget_control, state.zabs_id, set_value=state.zabs
;  stop
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chklls_specplot, state

  common sdss_chklls_cmm
  sdss_objinf, [llssearch[state.curqso].plate, $
                llssearch[state.curqso].fiber], filnm=fil
  parse_sdss, fil, flux, wave, sig=sig
  x_specplot, flux, sig, wave=wave, inflg=4, zin=state.zabs, /lls, $
              /block, ZOUT=zout
  llssearch[state.curqso].zlls = zout
  state.zabs = zout
  widget_control, state.zabs_id, set_value=state.zabs
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
; sdss_chklls, 'dr5_qso.fits', 'LLS/hiztemplate_32_34.fits', 'lls_search.fits',
;  SZSTART=3.2, SZEND=3.4, AMOD_FIL='LLS/all_mod.idl',
;  gzfil='LLS/dr5_lls_s2n2.fits'
pro sdss_chklls, qso_fil, templ_fil, outfil, IQSO=iqso, $
                 XSIZE=xsize, L_YSIZE=i_ysize, M_YSIZE=s_ysize, $
                 YSIZE=ysize, OBJNM=objnm, XMAX=xmax, LLIST=llist, $
                 szstart=SZSTART, szend=SZEND, $
                 drpath=DRPATH, AMOD_FIL=amod_fil, DATFILS=datfils, $
                 GZFIL=gzfil, MOCK_FIL=mock_fil, MOCK_PATH=mock_path

  common sdss_chklls_cmm
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'sdss_chklls, qso_fil, templ_fil, outfil, IQSO=, [v1.0]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-100
  if not keyword_set( MINVAL ) then minval = 6.
  if not keyword_set( IQSO ) then iqso = 0L
  if not keyword_set( DRPATH ) then drpath = getenv('SDSSPATH')+'/DR7_QSO/'

; Initialize the common blcok
  sdss_chklls_icmmn, drpath, qso_fil, szstart, szend, templ_fil, outfil, $
                     AMOD_FIL=amod_fil, DATFILS=datfils, GZFIL=gzfil, $
                     MOCK_FIL=mock_fil
  
;  tmp = { velpltstrct }
;  tmp2 = { abslinstrct }

; STATE

  state = {             $
          nqal: 0L, $
          curqso: iqso, $
          curqal: 0, $
          curlls: 0L, $
          llsfil: outfil, $
          drpath: drpath, $
          flg_new: 0, $
          flg_gz: 0, $
          zabs: 0., $
          zqso: 0., $
          minval: minval, $
          ntrans: 0L, $         ; PLOTTING LINES
          smooth: 0, $
          flg_mock: 0, $
          mock_path: '', $
          nplt: 0, $
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
          tau_id: 0L, $
          nspec_id: 0L, $
          pmin_id: 0L, $
          pmax_id: 0L, $
          lines_id: 0L, $
          top_id: 0L, $
          lls_id: 0L, $
          lya_id: 0L, $
          bottom_id: 0L, $
          rhs_id: 0L, $
          info_id: 0L, $
          quality_id: 0L, $
          scr1_id: 0L, $
          scr2_id: 0L, $
          hits_id: 0L, $
          NHI_id: 0L, $
          chi_id: 0L, $
          fspec_id: 0L, $
          stat_id: 0L, $
          plate_id: 0L, $
          gmag_id: 0L, $
          umag_id: 0L, $
          uerr_id: 0L, $
          fib_id: 0L, $
          idx_id: 0L, $
          cLLS_id: 0L, $
          tLLS_id: 0L, $
          pLLS_id: 0L, $
          help_text_id: 0L $
          }

;;;;;;;;;;;;;;
; SETUP LINES
;  sdss_chklls_llist, state
;  state.lls_lin = x_setline(1215.6701d)

  if keyword_set(GZFIL) then state.flg_gz = 1
  if keyword_set(MOCK_FIL) then begin
      state.flg_mock = 1
      if not keyword_set(MOCK_PATH) then $
        state.mock_path = getenv('LLSPAP')+'/SDSS/DR7/Analysis/Mocks/'
  endif

; Other setup
;  state.nqal = n_elements(qalstr)
;  state.con_dir = con_dir

;    WIDGET
  base = WIDGET_BASE( title = 'sdss_chklls: Check spectra', /column, $
                    UNAME='BASE', /tlb_size_events, xoffset=200L)
  state.base_id = base
  

  state.top_id = WIDGET_BASE( state.base_id, /row, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(1*ysize/3.), $
                              uvalue='TOP_BASE', frame=2)
  state.lya_id = widget_draw(state.top_id, xsize=round(xsize/4), $
                              ysize=round(ysize/3), /frame, retain=2, $
                              uvalue='LYADRAW')
  state.lls_id = widget_draw(state.top_id, xsize=round(3*xsize/4), $
                              ysize=round(ysize/3), /frame, retain=2, $
                              uvalue='LLSDRAW')
  state.bottom_id = WIDGET_BASE( state.base_id, /row, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(2*ysize/3.), $
                              uvalue='TOP_BASE', frame=2)

;;;;;; Info window ;;;;;;;;;;;
  state.info_id = $
    WIDGET_BASE( state.bottom_id, /column, /base_align_center,/align_center, $
                 uvalue='INFO_BASE', frame=2, xsize=xsize/3.)
;               ysize=round(ysize/3.))
  ;; Info
  llsinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.tau_id = cw_field(llsinf, title='N(HI) ', value=0., xsize=5)
  state.zabs_id = cw_field(llsinf, title='zabs', value=0., $
                           /floating, uvalue='ZABS', xsize=10)
;  state.ew_id = cw_field(llsinf, title='EW: ', value=state.ew, xsize=7)

  ;; RA, DEC
  radeci = widget_base(state.info_id, /row, /align_center, frame=2)
  state.idx_id = cw_field(radeci, title='IDX: ', value=0., xsize=5)
  state.plate_id = cw_field(radeci, title='PLATE: ', value=0L, xsize=4)
  state.fib_id = cw_field(radeci, title='FIB: ', value=0L, xsize=4)
  magi = widget_base(state.info_id, /row, /align_center, frame=2)
  state.gmag_id = cw_field(magi, title='u-g: ', value=0., xsize=5)
  state.umag_id = cw_field(magi, title='u mag: ', value=0., xsize=5)
  state.uerr_id = cw_field(magi, title='u err: ', value=0., xsize=5)

  ;; Lya
  lyainf = widget_base(state.info_id, /row, /align_center, frame=2)

;      BUTTONS
  NHIBase = widget_base(state.info_id, /row, /align_center, frame=2)
  state.tLLS_id = cw_bgroup(NHIBase, ['z', 'N', 'C'], LABEL_LEFT='Tweak?: ', $
                           column=3, UVALUE='TWKLLS', /NONEXCLUS)
  LLSBase = widget_base(state.info_id, /column, /align_center, frame=2)
  state.cLLS_id = cw_bgroup(LLSBase, ['?', 'N', 'Y', 'YP', 'M', 'X','S','B'], $
                            LABEL_LEFT='LLS?: ', $
                           /exclusive, column=8, UVALUE='LLSBGROUP')
  state.pLLS_id = cw_bgroup(LLSBase, ['N', 'Y', 'M'], $
                            LABEL_LEFT='Additional (P)LLS?: ', $
                           /exclusive, column=3, UVALUE='PLLS')
;  defi = WIDGET_BUTTON(butbase, value='DEFINITE',uvalue='DEFINITE')
;  good = WIDGET_BUTTON(butbase, value='GOOD',uvalue='GOOD')
;  maybe = WIDGET_BUTTON(butbase, value='MAYBE',uvalue='MAYBE')
;  bad  = WIDGET_BUTTON(butbase, value='BAD', uvalue='BAD')

  butbase3 = widget_base(state.info_id, /row, /align_center, frame=2)
  smth = WIDGET_BUTTON(butbase3, value='SMTH', uvalue='SMTH')
  smth = WIDGET_BUTTON(butbase3, value='NOSMTH', uvalue='NOSMTH')
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
  next = WIDGET_BUTTON(butbase2, value='NEXT', uvalue='NEXT')
  prev = WIDGET_BUTTON(butbase2, value='PREV', uvalue='PREV')
  save = WIDGET_BUTTON(butbase2, value='SAVE', uvalue='SAVE')
  done = WIDGET_BUTTON(butbase2, value='DONE',uvalue='DONE')
  done = WIDGET_BUTTON(butbase2, value='NOSV',uvalue='NOSV')
  splt = WIDGET_BUTTON(butbase2, value='SPLT',uvalue='SPLT')
  splt = WIDGET_BUTTON(butbase2, value='FIDDLE',uvalue='FIDDLE')
  srch = WIDGET_BUTTON(butbase2, value='SEARCH',uvalue='SEARCH')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Metals DRAW
  state.mdrawbase_id = $
    WIDGET_BASE( state.bottom_id, /column, /base_align_center,/align_center, $
               uvalue='SDRAW_BASE', frame=2)

  state.chi_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./3), $
                              ysize=round(1.*ysize/3), /frame, retain=2, $
                              uvalue='CHIDRAW')
  state.fspec_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./3), $
                              ysize=round(1.*ysize/3), /frame, retain=2, $
                              uvalue='FSPDRAW')


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
;  if qalstr[0].zabs[0] EQ 0. then sdss_chklls_next, state

  ; Load data
  sdss_chklls_setup, state

  sdss_chklls_updinfo, state
  sdss_chklls_updplt, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'sdss_chklls', base, /no_block
;  delvarx, fx, wv, npix, sig, qalstr
;  delvarx, allmod

  return
end
	
