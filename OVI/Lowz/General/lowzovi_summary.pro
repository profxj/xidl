;+ 
; NAME:
; lowzovi_summary   
;   Version 1.1
;
; PURPOSE:
;    Launches a GUI which can be used to plot and print some 
;  simple info from the galaxy survey.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   ovi  -- Galaxy survey
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /REDSHIFT -- List only those galaxies with redshift measurements
;  MXY=      -- Maximum number of galaxies to show without scrolling
;               [Default: 30L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lowzovi_summary, ovi
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lowzovi_summary_initlzovin

common lowzovi_summary_common, lzovi_gal, lzovi_fspec, lzovi_flg, lzovi_indx, $
  lzovi_prsstate

return
end

;;;;
; Events
;;;;

pro lowzovi_summary_event, ev

common lowzovi_summary_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'BUTTON': begin
          widget_control, state.button_id, get_value=buttons
          state.flg_psfil = buttons[0]
          state.flg_survey = buttons[1]
      end
      'ZIMG': lowzovi_summary_zimg, state
      'RADEC': lowzovi_summary_radec, state
      'PLT': lowzovi_summary_pltobj, state
      'PARSE': lowzovi_summary_updlist, state
      'ZHIST': lowzovi_summary_zhist, state
      'LIST': begin
          state.lindx = widget_info(state.list_id, /list_select)
          lowzovi_summary_setlist, state
      end
      'OBJNMLIST': $
        state.objnm_sel = widget_info(state.listobjnm_id, /list_select)
      'FSPECLIST': begin
          state.fspec_sel = widget_info(state.listfspec_id, /list_select)
          lowzovi_summary_fspeclst, state, state.flg_fspec
      end
      'FANLY': lowzovi_summary_chnganly, state
      'CHGTYPE': lowzovi_summary_chgtype, state
      'DONE' : begin
          delvarx, lzovi_indx, lzovi_gal, lzovi_prsstate
          widget_control, ev.top, /destroy
          return
      end
      else:
  endcase

;
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
pro lowzovi_summary_mklist, state, list
common lowzovi_summary_common

  if lzovi_indx[0] EQ -1 then list = ' No Obj satisfy those criterion!! ' $
  else begin
      ;; Make array
      list = strarr(n_elements(lzovi_indx))

      case state.nfilt of 
          2: fmt = '(a6,1x,i1,1x,a23,1x,a3,1x,f8.5,1x,f6.2,1x,2f6.0)'
          else: stop
      endcase

      ;; Impact parameter (arcmin)
      rho = sqrt(lzovi_gal.ra^2 + lzovi_gal.dec^2)/60.
      
      ;; Loop
      for q=0L,n_elements(list)-1 do begin
          jj = lzovi_indx[q]
          ;; Mag
          if lzovi_gal[jj].flg_anly MOD 2 EQ 0 then magstr = ' ' else begin
              magstr = ' '
              for i=0L,state.nfilt-1 do $
                magstr = magstr+string(lzovi_gal[jj].mag[i],FORMAT='(f5.2)')$
                +string(lzovi_gal[jj].magerr[i],FORMAT='(f6.2)')+' '
              magstr = strtrim(magstr,2)
          endelse
          ;; List
          list[q] = string(strtrim(lzovi_gal[jj].id,2)+ $
                           strtrim(lzovi_gal[jj].obj_id,2),$
                           lzovi_gal[jj].flg_anly, $
                           magstr, $
                           lzovi_gal[jj].gal_type, $
                           lzovi_gal[jj].z, rho[jj], $
                           lzovi_gal[jj].xypix, $
                           FORMAT=fmt)
      endfor
  endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro lowzovi_summary_updlist, state, NOINTER=nointer
common lowzovi_summary_common

  ;; Update index
  lowzovi_summary_setcindx, state, NOINTER=nointer
  ;; Make List
  lowzovi_summary_mklist, state, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=state.lindx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Flip Anly
pro lowzovi_summary_chnganly, state
common lowzovi_summary_common

  ;; Check indx
  if lzovi_indx[0] LT 0 then begin
      print, 'lowzovi_summary_flipanly: Select entry first!'
      return
  endif

  stop

  ;; Update
  lowzovi_summary_updlist, state, /NOINTER

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro lowzovi_summary_chgtype, state
common lowzovi_summary_common

  ;; Check indx
  if lzovi_indx[0] LT 0 then begin
      print, 'lowzovi_summary_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['E', 'Sa', 'Sb','Sc'])

  ;; Change
  lzovi_gal[lzovi_indx].type = type

  ;; Update
  lowzovi_summary_updlist, state, /NOINTER

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro lowzovi_summary_setcindx, state, NOINTER=nointer
common lowzovi_summary_common

  lowzovi_parse, lzovi_gal, lzovi_prsstate, INDX=indx, NOINTER=nointer
  lzovi_indx = temporary(indx)

  ;; Nobj
  nobj = n_elements([lzovi_indx])
  widget_control, state.nobj_id, set_value=nobj

  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Create fspec list
pro lowzovi_summary_fspeclst, state, flg
common lowzovi_summary_common

  gdfspec = where( strlen(strtrim(lzovi_gal[state.indx].fspec_fil,2)) NE 0,ngd)
  if ngd NE 0 then begin
      list = strtrim(lzovi_gal[state.indx].fspec_fil[gdfspec],2)
      widget_control, state.listfspec_id, set_value=list
      state.fspec_sel = state.fspec_sel < (ngd-1)
      widget_control, state.listfspec_id, set_list_select=state.fspec_sel 
      ;; Check to see if already selected
      if strtrim(state.fspec_fil[state.fspec_sel],2) $
        NE list[state.fspec_sel] then begin
          ;; Open
          lzovi_fspec = xmrdfits(list[state.fspec_sel], 1, $
                               structyp='lwdfspecstrct', /silent)
      endif
      ;; Setup Obj Names
      lowzovi_summary_objfillst, state
      state.fspec_fil[0:ngd-1] = list
  endif

  return
end
          

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Create Obj file list
pro lowzovi_summary_objfillst, state, flg
common lowzovi_summary_common

  ;; Grab the names
  tmp = x_getobjnm(lzovi_fspec, /lst)
  ;; Allow for 10000L
  indx = where(tmp EQ state.tobj)
  ;; Obj files
  gd = where(strlen(strtrim(lzovi_fspec[indx].obj_fil,2)) GT 0,ngd)
  if ngd EQ 0 then stop
  tmp = lzovi_fspec[indx].obj_fil[gd]
  nlst = n_elements(tmp)
  state.objnm_list[0:nlst-1] = tmp
  widget_control, state.listobjnm_id, set_value=tmp
  widget_control, state.listobjnm_id, set_list_select=state.objnm_sel

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; pltobj -- Plot
pro lowzovi_summary_pltobj, state
common lowzovi_summary_common


  widget_control, /hourglass   
  ;; Run plt obj
  wfccd_pltobj, state.fspec_fil[state.fspec_sel], state.tobj, $
    state.objnm_sel, /fspec

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; List
pro lowzovi_summary_setlist, state
common lowzovi_summary_common

  ;; INDEX
  state.indx = lzovi_indx[state.lindx]
  state.obj = strtrim(lzovi_gal[state.indx].id,2)+ $
    strtrim(lzovi_gal[state.indx].obj_id,2)
  state.tobj = strtrim(lzovi_gal[state.indx].id MOD 10000L,2)+ $
    strtrim(lzovi_gal[state.indx].obj_id,2)
  lowzovi_summary_fspeclst, state, state.flg_fspec

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; zhist -- zhist
pro lowzovi_summary_zhist, state
common lowzovi_summary_common

; Optional Keywords
  if not keyword_set( BIN ) then bin = 0.005

; PSFIL

  if keyword_set( PSFIL ) then begin
      device, decompose=0
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 4
      ps_open, file=psfil, font=1, /color, bpp=8
;      !y.margin = [5,2]
  endif

  ; PLOT
  clr = getcolor(/load)
  wset, 0
  plothist, lzovi_gal[lzovi_indx].z, bin=bin, xrange=zmnx, $
    xtitle='!17z', ymargin=[5,0], background=clr.white, color=clr.black, $
    /fill, fcolor=clr.green, ystyle=1

  if keyword_set( ZQSO ) then oplot, [zqso, zqso], $
    [0., 1e5], color=clr.red, linestyle=2

; CLOSE PSFIL

  if keyword_set(PSFIL) then begin
      ps_close, /noprint, /noid
      device, decomposed=1
      !x.thick = 1
      !y.thick = 1
      !p.charthick = 1
;      !y.margin = [4,2]
  endif



  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; zimg -- Image + Redshifts
pro lowzovi_summary_zimg, state
common lowzovi_summary_common

  ;; Hourglass
  widget_control, /hourglass   

  ;; Check for xatv
;  ans = x_guinum(2, TITLE='Have you set up xatv? (1/0)')

  xatverase
  xatvplot, lzovi_gal[lzovi_indx].xypix[0]-1., $
    lzovi_gal[lzovi_indx].xypix[1]-1., psym=1
 
;  Label

  for q=0L,n_elements(lzovi_indx)-1 do $
    xatvxyouts, lzovi_gal[lzovi_indx[q]].xypix[0]-40., $
    lzovi_gal[lzovi_indx[q]].xypix[1]+10., $
    string(lzovi_gal[lzovi_indx[q]].z,FORMAT='(f7.4)'),$
    charsize=2.0, color='red'

;  if not keyword_set(NONM) then begin
;      for q=0L,n_elements(lzovi_indx)-1 do begin
;          xatvxyouts, lzovi_gal[lzovi_indx[q]].xyimg[0]+2., $
;            lzovi_gal[lzovi_indx[q]].xyimg[1]+12., string(list[q],FORMAT='(a5)'),$
;        charsize=2.0, color='blue'
;      endfor 
;  endif


  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; radec -- Spatial dist of the galaxies
pro lowzovi_summary_radec, state
common lowzovi_summary_common


 
  ;; Max min
  mxx = max(lzovi_gal[lzovi_indx].ra, min=mnx)
  mxy = max(lzovi_gal[lzovi_indx].dec, min=mny)

 if (mxx-mnx) GT (mxy-mny) then begin
     diff = (mxx-mnx) - (mxy-mny)
     mxy = mxy + diff/2.
     mny = mny - diff/2.
 endif else begin
     diff = -(mxx-mnx) + (mxy-mny)
     mxx = mxx + diff/2.
     mnx = mnx - diff/2.
 endelse

 if keyword_set( state.flg_psfil ) then begin
     device, decompose=0
     ps_open, filename='summ.ps', /color, bpp=8, /maxs
     !p.thick = 5
     !p.charthick = 3
     !x.thick = 3
     !y.thick = 3
 endif
  clr = getcolor(/load)
  plot, [mnx,mxx]/60., [mny,mxy]/60., /nodata, color=clr.black, $
    background=clr.white, charsize=2.5, xmargin=[6.5,1], ymargin=[3,1],$
    xtitle='!17delta RA (arcmin)', ytitle='delta Dec(arcmin)'

  oplot, lzovi_gal[lzovi_indx].ra/60., lzovi_gal[lzovi_indx].dec/60., $
    color=clr.blue, psym=1, symsize=2.

  oplot, [0.], [0.], psym=2, color=clr.red, symsize=2.

 if keyword_set( state.flg_psfil ) then begin
      ps_close, /noprint, /noid
      device, decompose=1
      !p.thick = 1
      !p.charthick = 1
     !x.thick = 1
     !y.thick = 1
 endif
  return
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

pro lowzovi_summary, ovi, MXY=mxy, REDSHIFT=redshift, INDX=indx

common lowzovi_summary_common

;
;  if  N_params() LT 1  then begin 
;    print,'Syntax - ' + $
;             'lowzovi_summary, ovi [v1.0]'
;    return
;  endif 

  if not keyword_set(OVI) then begin
      files= findfile('*_gal.fits', count=nfil)
      if nfil EQ 0 then begin
          print, 'lowzovi_summary: Input a structure please'
          return
      endif
      ovi = xmrdfits(files[0], 1)
  endif

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( MXY ) then mxy = 30L
  if not keyword_set( LST_FONT ) then lst_font = '6x10'

; Init common
  lowzovi_summary_initlzovin
  if size(ovi, /type) EQ 7 then $
    lzovi_gal = xmrdfits(ovi, 1, structyp='galsurveystrct',/silent) $
  else lzovi_gal = ovi

;    STATE
  state = { $
            indx: -1L, $    ; True index
            lindx: -1L, $   
            flg_sci: 0L, $  ; 0L means include serendip
            flg_spec: 0L, $ ; 1L means spec only
            flg_z: 0L, $    ; 1L means good z only
            flg_psfil: 0L, $
            flg_survey: 0L, $
            nfilt: 0L, $
            obj: ' ', $
            tobj: ' ', $
            objnm: 0L, $
            objnm_list: strarr(100), $
            objnm_sel: 0L, $
            flg_fspec: 0L, $
            fspec_sel: 0L, $
            fspec_fil: strarr(100), $
            list_id: 0L, $
            listfspec_id: 0L, $
            listobjnm_id: 0L, $
            button_id: 0L, $
            nobj_id: 0L, $
            base_id: 0L $
          }

  if keyword_set( REDSHIFT ) then state.flg_z = 1L

;    WIDGET
  base = WIDGET_BASE( title = 'lowzovi_summary', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='ovi_summary')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parse Button List
;  spec_but = CW_BGROUP(toolbar, ['Spectra', 'Redshift', 'Science'], $
;                       FONT='courier', $
;                       row=3, /nonexclusive, $
;                       uvalue='PARSE', $
;                       LABEL_TOP="Parse List", /frame)
  spec_but = widget_button(toolbar, value='PARSE', uvalue='PARSE')
  state.nobj_id = cw_field(toolbar, /column, title='Nobj', value=0L, xsize=5L)
;  widget_control, spec_but, set_value=[state.flg_spec,state.flg_z,state.flg_sci]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Buttons
  butbase = widget_base(toolbar, /column, /align_center)
  plt = WIDGET_BUTTON(butbase, value='PLT',uvalue='PLT')
  chg_anly = WIDGET_BUTTON(butbase, value='Chng Anly',uvalue='FANLY', /align_right)
  chg_type = WIDGET_BUTTON(butbase, value='Chng Type',uvalue='CHGTYPE', $
                           /align_right)
  zhist = WIDGET_BUTTON(butbase, value='ZHIST',uvalue='ZHIST', FONT='courier')
  zimg = WIDGET_BUTTON(butbase, value='ZIMG',uvalue='ZIMG', FONT='courier')
  radec = WIDGET_BUTTON(butbase, value='RADEC',uvalue='RADEC', FONT='courier')

  state.button_id = CW_BGROUP(toolbar, ['PSFIL', $
                                         'flg_survey'], $
                              row=2, set_value=[state.flg_psfil, $
                                                state.flg_survey], $
                              /nonexclusive, font=stdfont, uvalue='BUTTON')

; Find number of Filters
  a = where(lzovi_gal.flg_anly MOD 2 EQ 1, na)
  if na NE 0 then begin
      b = where( strlen(strtrim(lzovi_gal[a[0]].filter,2)) EQ 0, nb)
      if nb EQ 0 then state.nfilt = 10 else state.nfilt = b[0]
  endif else state.nfilt = 0
  
; Get the index list
  if keyword_set( INDX ) then begin
      lzovi_indx = indx 
      widget_control, state.nobj_id, set_value=nobj
  endif else lowzovi_summary_setcindx, state
  lowzovi_summary_mklist, state, list

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lists
  ysz = mxy < n_elements(list)
  listbase = WIDGET_BASE(state.base_id, /row, /frame, $
                         /align_top, uvalue='LISTS')
  state.list_id = widget_list(listbase, FONT=lst_font, $
                              value=list, xsize=70L, ysize=ysz, $
                              uvalue = 'LIST')
  state.listfspec_id = WIDGET_LIST(listbase, FONT=lst_font, $
                               VALUE=[' '], xsize=25, $
                               uvalue='FSPECLIST', ysize = 3)
  state.listobjnm_id = WIDGET_LIST(listbase, FONT=lst_font, $
                               VALUE=[' '], xsize=25, $
                               uvalue='OBJNMLIST', ysize = 5)
  widget_control, state.listobjnm_id, set_list_select=-1

;   DONE
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

; Realize
  WIDGET_CONTROL, base, /realize
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  xmanager, 'lowzovi_summary', base, /no_block

; Finish
 
  return

end
