;+ 
; NAME:
; dla_analyse
;  V1.1
;
; PURPOSE:
;    Launches a GUI which enables simple plotting and inspection of
;     individual DLA
; CALLING SEQUENCE:
;   dla_analyse, [list]
;
; INPUTS:
;   [list] -- List of DLA to inspect  
;       [default: /u/xavier/DLA/Lists/tot_dla.lst]
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ROOT=  Path to the DLA tree
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_analyse, '/u/xavier/DLA/Lists/tot_dla.lst'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Events
;;;;

pro dla_analyse_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DLALIST' : begin
          state.curdla=ev.index
          dla_analyse_update, state
          if state.flg_dla MOD 2 EQ 0 then state.flg_dla = state.flg_dla+1
      end
      'SPECPLT': dla_analyse_specplt, state
      'VELPLT': dla_analyse_velplt, state
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;
;  UPDATE

pro dla_analyse_update, state

  ;; ZABS
  widget_control, state.zabs_id, set_value=state.dla[state.curdla].zabs
  ;; Name
  widget_control, state.name_id, set_value=strtrim(state.dla[state.curdla].qso,2)
  ;; NHI
  widget_control, state.nhi_id, set_value=state.dla[state.curdla].NHI
  ;; Dv
  widget_control, state.delv_id, set_value=state.dla[state.curdla].lwfvel
  ;; FeH
  widget_control, state.fehflg_id, set_value=state.dla[state.curdla].flgfe
  widget_control, state.feh_id, set_value=state.dla[state.curdla].feh
  ;; AlphH
  widget_control, state.alphflg_id, set_value=state.dla[state.curdla].flgalpha
  widget_control, state.alph_id, set_value=state.dla[state.curdla].alpha
  ;; CII
  widget_control, state.ciiflg_id, set_value=state.dla[state.curdla].flgCII
  widget_control, state.cii_id, set_value=state.dla[state.curdla].CII
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;
;  SPEC PLOT

pro dla_analyse_specplt, state

  ;; Check that one is selected
  if state.flg_dla MOD 2 EQ 0 then begin
      print, 'dla_analyse_plt: Need to select a dla first!'
      return
  end

  ;; Check for Low
  if state.dla[state.curdla].flglw NE 0 then begin
      datfil = strtrim(state.dla[state.curdla].lwfil,2)
      ipos = strpos(datfil, 'f.fits')
      errfil = datfil
      strput, errfil, 'e', ipos
      x_specplot, datfil, errfil, /QAL, /NRM, zin=state.dla[state.curdla].zabs
  endif else print, 'dla_analyse_plt: No data file!'
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;
;  VEL PLOT

pro dla_analyse_velplt, state

  ;; Check that one is selected
  if state.flg_dla MOD 2 EQ 0 then begin
      print, 'dla_analyse_velplt: Need to select a dla first!'
      return
  end

  ;; Check for Low
  if state.dla[state.curdla].flglw NE 0 then begin
      datfil = strtrim(state.dla[state.curdla].lwfil,2)
      x_velplt, datfil, state.dla[state.curdla].zabs
  endif else print, 'dla_analyse_plt: No data file!'
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_analyse, list, ROOT=root

; Optional Keywords
  if not keyword_set( list ) then begin
      print, 'dla_analyse: Using tot_dla.lst as the list'
      list = '~/DLA/Lists/all_mtl.lst'
  endif

  if not keyword_set(ROOT) then root = getenv('DLA')

; Parse DLA
  parse_dlalst, dla, list, /noelm, ROOT=root

; Files
  readcol, list, files, format='A'

; STATE

  state = {             $
            dla: dla, $
            curdla: 0L, $
            list: list, $
            flg_dla: 0, $
            files: files, $ 
            ndla: n_elements(files), $
            base_id: 0L, $      ; Widgets
            dla_id: 0L, $
            name_id: 0L, $
            nhi_id: 0L, $
            delv_id: 0L, $
            feh_id: 0L, $
            fehflg_id: 0L, $
            alph_id: 0L, $
            alphflg_id: 0L, $
            cii_id: 0L, $
            ciiflg_id: 0L, $
            zabs_id: 0L $
          }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    WIDGET
  base = WIDGET_BASE( title = 'dla_analyse', /column, $
                    UNAME='BASE', xoffset=300L, yoffset=300L)
	  state.base_id = base
  
;;;;;;;;;;;;;;
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='dla_analyse', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.1)', /align_center)

;;;;;;;;;;;;;;;;;;
;  SUMMARY

  ;; Name + zabs
  firstpair = WIDGET_BASE( toolbar, /column, /base_align_left,$
                         /align_center)
  state.name_id = cw_field(firstpair, value='', $
                             xsize=15, title='Name:', UVALUE='Name')
  state.zabs_id = cw_field(firstpair, value=0., $
                             xsize=10, title='zabs:', UVALUE='Zabs')
  ;; NHI+delv
  nhipair = WIDGET_BASE( toolbar, /column, /base_align_left,$
                         /align_center)
  state.nhi_id = cw_field(nhipair, value=0., $
                             xsize=7, title='N(HI): ', UVALUE='NHI')
  state.delv_id = cw_field(nhipair, value=0., $
                             xsize=7, title='Dv: ', UVALUE='Dv')
  ;; Fe/H
  fehpair = WIDGET_BASE( toolbar, /column, /base_align_left,$
                         /align_center)
  state.fehflg_id = cw_field(fehpair, value=0, $
                             xsize=3, title='[Fe/H]_flg:', UVALUE='FeHflg')
  state.feh_id = cw_field(fehpair, value=0., $
                             xsize=8, title='[Fe/H]: ', UVALUE='FeH')
  ;; Alpha/H
  alphpair = WIDGET_BASE( toolbar, /column, /base_align_left,$
                         /align_center)
  state.alphflg_id = cw_field(alphpair, value=0, $
                             xsize=3, title='[A/H]_flg:', UVALUE='AHFLG')
  state.alph_id = cw_field(alphpair, value=0., $
                             xsize=8, title='[Aa/H]:', UVALUE='AH')
  ;; CII
  ciipair = WIDGET_BASE( toolbar, /column, /base_align_left,$
                         /align_center)
  state.ciiflg_id = cw_field(ciipair, value=0, $
                             xsize=2, title='CII*_flg:', UVALUE='CIIFLG')
  state.cii_id = cw_field(ciipair, value=0., $
                             xsize=8, title='CII*:    ', UVALUE='CII')
;  state.msk_id = cw_field(toolbar, value='', $
;                             /return_events, xsize=5,$
;                             title='Mask', UVALUE='MASK')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  listbase = WIDGET_BASE(state.base_id, /row, /frame, /base_align_center,$
                              uvalue='LISTS')

; OBJ List
  state.dla_id = WIDGET_LIST(listbase, $
                             VALUE=state.files, $
                             uvalue='DLALIST', ysize = 20)
  widget_control, state.dla_id, set_list_select=-1

;      BUTTONS
  butbase = widget_base(listbase, /column, /align_center)
  specplt = WIDGET_BUTTON(butbase, value='SPECPLT',uvalue='SPECPLT')
  velplt = WIDGET_BUTTON(butbase, value='VELPLT',uvalue='VELPLT')
;  zhist = WIDGET_BUTTON(butbase, value='ZHIST',uvalue='ZHIST')
;  chk = WIDGET_BUTTON(butbase, value='CHK',uvalue='CHK')
;  edit = WIDGET_BUTTON(butbase, value='EDIT',uvalue='EDIT')
;  img = WIDGET_BUTTON(butbase, value='IMG',uvalue='IMG')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')

; Realize
  WIDGET_CONTROL, base, /realize

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'dla_analyse', base, /no_block

  return
end
