;+ 
; NAME:
; ovi_summary   
;   Version 1.0
;
; PURPOSE:
;    Launches a GUI which can be used to plot and print some 
;  simple info from the galaxy survey.
;
; CALLING SEQUENCE:
; ovi_summary, ovi, MXY=mxy, REDSHIFT=redshift
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
;   string = ovi_summary( 'Enter filename: ')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ovi_summary_initcmmn

common ovi_summary_common, cmm_ovi, cmm_fspec, cmm_flg, cmm_indx

return
end

;;;;
; Events
;;;;

pro ovi_summary_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'ZIMG': ovi_summary_zimg, state
      'PLT': ovi_summary_pltobj, state
      'PARSE': begin
          WIDGET_CONTROL, ev.id, get_value = val
          state.flg_spec = val[0]
          state.flg_z = val[1]
          state.flg_sci = val[2]
          ovi_summary_updlist, state
      end
      'ZHIST': ovi_summary_zhist, state
      'LIST': begin
          state.lindx = widget_info(state.list_id, /list_select)
          ovi_summary_setlist, state
      end
      'OBJNMLIST': $
        state.objnm_sel = widget_info(state.listobjnm_id, /list_select)
      'FSPECLIST': begin
          state.fspec_sel = widget_info(state.listfspec_id, /list_select)
          ovi_summary_fspeclst, state, state.flg_fspec
      end
      'FANLY': ovi_summary_chnganly, state
      'CHGTYPE': ovi_summary_chgtype, state
      'DONE' : begin
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
pro ovi_summary_mklist, state, list
common ovi_summary_common

  if cmm_indx[0] EQ -1 then list = ' No Obj satisfy those criterion!! ' $
  else begin
      ;; Make array
      list = strarr(n_elements(cmm_indx))

      case state.nfilt of 
          2: fmt = '(a6,1x,i1,1x,a23,1x,a3,1x,f8.5)'
          else: stop
      endcase

      ;; Loop
      for q=0L,n_elements(list)-1 do begin
          jj = cmm_indx[q]
          ;; Mag
          if cmm_ovi[jj].flg_anly MOD 2 EQ 0 then magstr = ' ' else begin
              magstr = ' '
              for i=0L,state.nfilt-1 do $
                magstr = magstr+string(cmm_ovi[jj].mag[i],FORMAT='(f5.2)')$
                +string(cmm_ovi[jj].magerr[i],FORMAT='(f6.2)')+' '
              magstr = strtrim(magstr,2)
          endelse
          ;; List
          list[q] = string(strtrim(cmm_ovi[jj].id,2)+ $
                           strtrim(cmm_ovi[jj].obj_id,2),$
            cmm_ovi[jj].flg_anly, $
            magstr, $
            cmm_ovi[jj].gal_type, $
            cmm_ovi[jj].z, $
            FORMAT=fmt)
      endfor
  endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro ovi_summary_updlist, state
common ovi_summary_common

  ;; Update index
  ovi_summary_setcindx, state
  ;; Make List
  ovi_summary_mklist, state, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=state.lindx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Flip Anly
pro ovi_summary_chnganly, state
common ovi_summary_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'ovi_summary_flipanly: Select entry first!'
      return
  endif

  stop

  ;; Update
  ovi_summary_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro ovi_summary_chgtype, state
common ovi_summary_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'ovi_summary_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['E', 'Sa', 'Sb','Sc'])

  ;; Change
  cmm_ovi[cmm_indx].type = type

  ;; Update
  ovi_summary_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro ovi_summary_setcindx, state
common ovi_summary_common

  ;; Science
  if state.flg_sci EQ 0 then cmm_indx = lindgen(n_elements(cmm_ovi)) $
  else cmm_indx = where(cmm_ovi.obj_id EQ 'a')

  ;; Spec
  if state.flg_spec EQ 1L then begin
      a = where( cmm_ovi[cmm_indx].flg_anly MOD 4 GT 1, na)
      if na EQ 0 then cmm_indx = -1 else cmm_indx = cmm_indx[a]
  endif

  ;; Redshift
  if state.flg_z EQ 1L then begin
      a = where( cmm_ovi[cmm_indx].flg_anly MOD 8 GT 3, na)
      if na EQ 0 then cmm_indx = -1 else cmm_indx = cmm_indx[a]
  endif

  ;; Nobj
  nobj = n_elements([cmm_indx])
  widget_control, state.nobj_id, set_value=nobj

  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Create fspec list
pro ovi_summary_fspeclst, state, flg
common ovi_summary_common

  gdfspec = where( strlen(strtrim(cmm_ovi[state.indx].fspec_fil,2)) NE 0,ngd)
  if ngd NE 0 then begin
      list = strtrim(cmm_ovi[state.indx].fspec_fil[gdfspec],2)
      widget_control, state.listfspec_id, set_value=list
      state.fspec_sel = state.fspec_sel < (ngd-1)
      widget_control, state.listfspec_id, set_list_select=state.fspec_sel 
      ;; Check to see if already selected
      if strtrim(state.fspec_fil[state.fspec_sel],2) $
        NE list[state.fspec_sel] then begin
          ;; Open
          cmm_fspec = xmrdfits(list[state.fspec_sel], 1, $
                               structyp='lwdfspecstrct', /silent)
      endif
      ;; Setup Obj Names
      ovi_summary_objfillst, state
      state.fspec_fil[0:ngd-1] = list
  endif

  return
end
          

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Create Obj file list
pro ovi_summary_objfillst, state, flg
common ovi_summary_common

  ;; Grab the names
  tmp = x_getobjnm(cmm_fspec, /lst)
  indx = where(tmp EQ state.obj)
  ;; Obj files
  gd = where(strlen(strtrim(cmm_fspec[indx].obj_fil,2)) GT 0,ngd)
  if ngd EQ 0 then stop
  tmp = cmm_fspec[indx].obj_fil[gd]
  nlst = n_elements(tmp)
  state.objnm_list[0:nlst-1] = tmp
  widget_control, state.listobjnm_id, set_value=tmp
  widget_control, state.listobjnm_id, set_list_select=state.objnm_sel

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; pltobj -- Plot
pro ovi_summary_pltobj, state
common ovi_summary_common


  widget_control, /hourglass   
  ;; Run plt obj
  stop
  wfccd_pltobj, state.fspec_fil[state.fspec_sel], state.obj, $
    state.objnm_sel, /fspec

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; List
pro ovi_summary_setlist, state
common ovi_summary_common

  ;; INDEX
  state.indx = cmm_indx[state.lindx]
  state.obj = strtrim(cmm_ovi[state.indx].id,2)+ $
    strtrim(cmm_ovi[state.indx].obj_id,2)
  ovi_summary_fspeclst, state, state.flg_fspec

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; zhist -- zhist
pro ovi_summary_zhist, state
common ovi_summary_common

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
  plothist, cmm_ovi[cmm_indx].z, bin=bin, xrange=zmnx, $
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
pro ovi_summary_zhist, state
common ovi_summary_common

  ;; Check for xatv
  ans = x_guinum(2, TITLE='Have you set up xatv? (1/0)')

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

pro ovi_summary, ovi, MXY=mxy, REDSHIFT=redshift

common ovi_summary_common

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ovi_summary, ovi, MXY=, REDSHIFT= [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( MXY ) then mxy = 30L
  if not keyword_set( LST_FONT ) then lst_font = '6x10'

; Init common
  ovi_summary_initcmmn
  cmm_ovi = xmrdfits(ovi, 1, structyp='galsurveystrct',/silent)

;    STATE
  state = { $
            indx: -1L, $    ; True index
            lindx: -1L, $   
            flg_sci: 0L, $  ; 0L means include serendip
            flg_spec: 0L, $ ; 1L means spec only
            flg_z: 0L, $    ; 1L means good z only
            nfilt: 0L, $
            obj: ' ', $
            objnm: 0L, $
            objnm_list: strarr(100), $
            objnm_sel: 0L, $
            flg_fspec: 0L, $
            fspec_sel: 0L, $
            fspec_fil: strarr(100), $
            list_id: 0L, $
            listfspec_id: 0L, $
            listobjnm_id: 0L, $
            nobj_id: 0L, $
            base_id: 0L $
          }

  if keyword_set( REDSHIFT ) then state.flg_z = 1L

;    WIDGET
  base = WIDGET_BASE( title = 'ovi_summary', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='ovi_summary')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parse Button List
  spec_but = CW_BGROUP(toolbar, ['Spectra', 'Redshift', 'Science'], $
                       FONT='courier', $
                       row=3, /nonexclusive, $
                       uvalue='PARSE', $
                       LABEL_TOP="Parse List", /frame)
  state.nobj_id = cw_field(toolbar, /column, title='Nobj', value=0L, xsize=5L)
  widget_control, spec_but, set_value=[state.flg_spec,state.flg_z,state.flg_sci]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Buttons
  butbase = widget_base(toolbar, /column, /align_center)
  plt = WIDGET_BUTTON(butbase, value='PLT',uvalue='PLT')
  chg_anly = WIDGET_BUTTON(butbase, value='Chng Anly',uvalue='FANLY', /align_right)
  chg_type = WIDGET_BUTTON(butbase, value='Chng Type',uvalue='CHGTYPE', $
                           /align_right)
  zhist = WIDGET_BUTTON(butbase, value='ZHIST',uvalue='ZHIST', FONT='courier')
  zimg = WIDGET_BUTTON(butbase, value='ZIMG',uvalue='ZIMG', FONT='courier')


; Find number of Filters
  a = where(cmm_ovi.flg_anly MOD 2 EQ 1, na)
  if na NE 0 then begin
      b = where( strlen(strtrim(cmm_ovi[a[0]].filter,2)) EQ 0, nb)
      if nb EQ 0 then state.nfilt = 10 else state.nfilt = b[0]
  endif else state.nfilt = 0
  
; Create the list
  ovi_summary_setcindx, state
  ovi_summary_mklist, state, list

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lists
  ysz = mxy < n_elements(list)
  listbase = WIDGET_BASE(state.base_id, /row, /frame, $
                         /align_top, uvalue='LISTS')
  state.list_id = widget_list(listbase, FONT=lst_font, $
                              value=list, xsize=60L, ysize=ysz, $
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

  xmanager, 'ovi_summary', base, /no_block

; Finish
;  delvarx, cmm_indx, cmm_ovi
  
  return

end
