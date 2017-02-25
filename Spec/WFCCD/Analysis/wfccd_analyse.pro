;+ 
; NAME:
; wfccd_analyse
;    Version 1.0
;
; PURPOSE:
;   Calls various routines to analyse the data
;
; CALLING SEQUENCE:
;   
;   wfccd_analyse, wfccd, maskid, expsr, XSIZE=, YSIZE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_analyse, wfccd
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_analyse_initcmm, i_wfccd

common wfccd_analyse_cmm, $
  a_wfccd, $
  a_wfspec

  a_wfccd = i_wfccd

 return
end


;;;;
; Events
;;;;

pro wfccd_analyse_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'OBJLIST' : begin
          case state.flg_spec of 
              0: begin
                  ;; Parse obj #
                  fil = strtrim(state.obj_list[ev.index],2)
                  ipos = strpos(fil, 'c_')
                  state.obj=long(strmid(fil,ipos+2,2))
              end
              1: state.obj=long(state.obj_list[ev.index])
              else:
          endcase
          state.obj_fil = strtrim(state.obj_list[ev.index],2)
          wfccd_analyse_setobj, state
      end
      'OBJNMLIST' : begin
          state.objnm=ev.index
          state.obj_nm = strtrim(state.objnm_list[ev.index],2)
          wfccd_analyse_setobjnm, state
      end
      'EXPLIST' : state.expsr=ev.index
      'PLT': wfccd_analyse_pltobj, state
      'ZHIST': wfccd_analyse_zhist, state
      'EDIT': wfccd_analyse_edit, state
      'CHK': wfccd_analyse_chk, state
      'IMG': wfccd_analyse_img, state
      'REDO': wfccd_analyse_redo, state
      ;; ZFIND
      'ZFIND': wfccd_analyse_zfind, state
      'WVMNVAL': state.wvmnval = ev.value
      'WVMXVAL': state.wvmxval = ev.value
      'ZMINVAL': state.zmin = ev.value
      'ZMAXVAL': state.zmax = ev.value
      'ZFINDBUT': begin
          case ev.value of 
              0: begin  ; WVMNX
                  if state.flg_zfind MOD 2 EQ 1 then $
                    state.flg_zfind = state.flg_zfind - 1 else $
                    state.flg_zfind = state.flg_zfind + 1 
              end
              1: begin  ; Plot
                  if state.flg_zfind MOD 4 GE 2 then $
                    state.flg_zfind = state.flg_zfind - 2 else $
                    state.flg_zfind = state.flg_zfind + 2 
              end
              2: begin  ; Zmin
                  if state.flg_zfind MOD 8 GE 4 then $
                    state.flg_zfind = state.flg_zfind - 4 else $
                    state.flg_zfind = state.flg_zfind + 4 
              end
              3: begin  ; Zmax
                  if state.flg_zfind MOD 16 GE 8 then $
                    state.flg_zfind = state.flg_zfind - 8 else $
                    state.flg_zfind = state.flg_zfind + 8 
              end
          endcase
      end
      ;; EDIT OBJSTR
      'OBJBUT': begin
          case ev.value of 
              0: begin  ; WVMNX
                  if state.flg_objstr MOD 2 EQ 1 then $
                    state.flg_objstr = state.flg_objstr - 1 else $
                    state.flg_objstr = state.flg_objstr + 1 
              end
              1: begin  ; Plot
                  if state.flg_objstr MOD 4 GE 2 then $
                    state.flg_objstr = state.flg_objstr - 2 else $
                    state.flg_objstr = state.flg_objstr + 2 
              end
              else:
          endcase
      end
      'EDITOBJ': wfccd_analyse_editobjstr, state
      ;
      'FSPECOBJ': begin
          if ev.select EQ 1 then begin
              state.flg_spec = ev.value
              wfccd_analyse_fspecobj, state
          endif
      end
      'NSPEC': begin
          widget_control, state.nspec_id, get_value=nplt
          state.nplt = nplt
          state.npg =  (state.nobj/state.nplt) + (state.nobj MOD state.nplt NE 0)
          state.curpg = state.curpg < (state.npg-1)
          wfccd_analyse_Plot, state
          widget_control, state.pg_id, set_value=state.curpg
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; getobjlst -- Grabs the list of 'Obj' or 'Fspec' values
pro wfccd_analyse_getobjlst, state

common wfccd_analyse_cmm

  case state.flg_spec of 
      0: begin
          a = findfile('Extract/Fspec*fits*', count=cnt)
          if cnt EQ 0 then begin
              print, 'wfccd_analyse_getobjlst: No Fspec files found!'
              return
          endif
          state.nobj_list = cnt
          state.obj_list[0:cnt-1] = strtrim(a,2)
          ;; Truncate the gz
          for qq=0L,cnt-1 do begin
              slen = strlen(state.obj_list[qq])
              if strmid(state.obj_list[qq],slen-2,2) EQ 'gz' then $
                state.obj_list[qq] = strmid(state.obj_list[qq],0,slen-3)
          endfor
      end
      1: begin
          a = a_wfccd[uniq(a_wfccd.mask_id,sort(a_wfccd.mask_id))].mask_id
          b = a[where(a GE 0L)]
          state.nobj_list = n_elements(b)
          state.obj_list[0:state.nobj_list-1] = strtrim(b,2)
      end
      else:
  endcase

  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; getobjnmlst -- Grabs the list of Object names
pro wfccd_analyse_getobjnmlst, state

common wfccd_analyse_cmm

  ;; Get list
  tmp = x_getobjnm(a_wfspec, /LST)
  state.nobjnm_list = n_elements(tmp)
  state.objnm_list[0:state.nobjnm_list-1] = tmp

  ;; Set GUI
  widget_control, state.objnm_id, $
    set_value=state.objnm_list[0:state.nobjnm_list-1]
  widget_control, state.objnm_id, set_list_select=-1

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; setobjvar -- Sets the obj var
pro wfccd_analyse_setobj, state

common wfccd_analyse_cmm

  case state.flg_spec of
      0: wfccd_wrfspec, a_wfspec, state.obj_fil, /read
      1: a_wfspec = mrdfits(state.obj_fil,1,STRUCTYP='specobjstrct',/silent)
      else:
  endcase

  ;; Field
  widget_control, state.field_id, set_value=strtrim(a_wfspec[0].field,2)
  ;; Grab the obj names
  wfccd_analyse_getobjnmlst, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; setobjnm -- Sets the obj_nm and updates things
pro wfccd_analyse_setobjnm, state

common wfccd_analyse_cmm

  ;; Set index
  state.objnm = x_getobjnm(a_wfspec, state.obj_nm)

  ;; Set exposure
  case state.flg_spec of
      0: begin
          widget_control, state.exp_id, set_value=a_wfspec[state.objnm].texp
          state.expsr = 0L
          widget_control, state.expsr_id, $
            set_value=a_wfspec[state.objnm].obj_fil[0:a_wfspec[state.objnm].nexp-1]
          widget_control, state.expsr_id, set_list_select=state.expsr
      end
      1: widget_control, state.exp_id, set_value=a_wfspec[state.objnm].exp
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; pltobj -- Plot
pro wfccd_analyse_pltobj, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Set exposure
  case state.flg_spec of
      0: wfccd_pltobj, state.obj_fil, state.obj_nm, state.expsr, /fspec
      1: wfccd_pltobj, a_wfccd, state.obj, state.obj_nm
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; zhist -- zhist
pro wfccd_analyse_zhist, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Set exposure
  case state.flg_spec of
      0: wfccd_fspeczhist, state.obj_fil
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; zfind -- zfind
pro wfccd_analyse_zfind, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Plot
  if state.flg_zfind MOD 4 GE 2 then plot = 1 else plot=0
  ;; WVMNX
  if state.flg_zfind MOD 2 EQ 1 then wvmnx=[state.wvmnval,state.wvmxval] else $
    wvmnx=0
  ;; ZMIN, ZMAX
  if state.flg_zfind MOD 8 GE 4 then zmin = state.zmin else zmin = 0
  if state.flg_zfind MOD 16 GE 8 then zmax = state.zmax else zmax = 0
  ;; Set exposure
  case state.flg_spec of
      0: wfccd_zfind, state.obj_fil, state.obj_nm, PLOT=plot, WVMNX=wvmnx, $
        zmin=zmin, zmax=zmax
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; edit -- edit
pro wfccd_analyse_edit, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Set exposure
  case state.flg_spec of
      0: wfccd_editspec, state.obj_fil, state.obj_nm, /fspec
      1: wfccd_editspec, state.obj_fil 
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; chk -- Plot the spectra
pro wfccd_analyse_chk, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Set exposure
  case state.flg_spec of
      0: wfccd_chkfspec, state.obj_fil, /zfind
      1: wfccd_chkspec, a_wfccd, state.obj, state.expsr
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; img -- Plot the spectra
pro wfccd_analyse_img, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Set exposure
  case state.flg_spec of
      0: begin
          xatv, a_wfspec[0].img_fil, /invert
          wfccd_fspeczimg, state.obj_fil
      end
      else:
  endcase

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; fspecobj -- FSPEC vs. OBJ
pro wfccd_analyse_fspecobj, state

common wfccd_analyse_cmm

  wfccd_analyse_getobjlst, state
  widget_control, state.obj_id, $
    set_value=state.obj_list[0:state.nobj_list-1]

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; pltobj -- Plot
pro wfccd_analyse_editobjstr, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Options
  if state.flg_objstr MOD 2 NE 0 then NOSLIT = 1 else NOSLIT = 0
  if state.flg_objstr MOD 4 GE 2 then NOOBJ = 1 else NOOBJ = 0
  ;; Call program
  wfccd_editobjstr, a_wfccd, state.obj, state.expsr, NOSLIT=noslit, $
    NOOBJ=noobj

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; REDO -- Redo
pro wfccd_analyse_redo, state

common wfccd_analyse_cmm

  widget_control, /hourglass   
  ;; Call program
  extrct_obj, a_wfccd, state.obj, lindgen(a_wfspec[state.objnm].nexp), $
    state.obj_nm, /mkall

  return
end

          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro wfccd_analyse, wfccd, OBJ=obj

common wfccd_analyse_cmm

;
  if  N_params() LT 1 then begin 
    print,'Syntax - ' + $
      'wfccd_pltobj, wfccd, /OBJ  [v1.0]'
    return
  endif 

;  Optional Keywords

  if N_params() EQ 0 then begin
      fils = findfile('Extract/Fspec*fits', count=nfil)
      if nfil EQ 0 then return 
      wfccd = x_guilist(fils)
  endif
      
  ;; Initialize the Common block
  wfccd_analyse_initcmm, wfccd

; STATE

  state = {             $
            wfccd: wfccd, $
            path: '', $
            flg_spec: 0, $  ; 1=obj, 0=Fspec
            obj_list: strarr(100), $
            obj_fil: '',$
            obj: 0L, $   ; OBJ or FSPEC Index
            objnm: 0L, $
            obj_nm: '', $
            expsr: 0L, $
            nobj_list: 0L, $
            objnm_list: strarr(100), $
            nobjnm_list: 0L, $
            wvmnval: 3800., $      ; Zfind
            wvmxval: 9000., $      
            zmin: -0.03, $      
            zmax: 0.5, $      
            flg_zfind: 0, $       ;; Binary flag
            flg_objstr: 0, $       ;; Binary flag
            base_id: 0L, $      ; Widgets
            obj_id: 0L, $
            objnm_id: 0L, $
            field_id: 0L, $
            msk_id: 0L, $
            fspecobj_id: 0L, $
            expsr_id: 0L, $
            exp_id: 0L $
          }

;;; 
; Set Obj
  if keyword_set(OBJ) then state.flg_spec = 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    WIDGET
  base = WIDGET_BASE( title = 'wfccd_analyse', /column, $
                    UNAME='BASE', xoffset=300L, yoffset=300L)
	  state.base_id = base
  
;;;;;;;;;;;;;;
; FSPEC vs OBJ
  fspecbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
  state.fspecobj_id = CW_BGROUP(fspecbar, ['Fspec', 'Obj'], $
                              row=1, /exclusive, /return_index, $
                              set_value=state.flg_spec,  uvalue='FSPECOBJ')
;;;;;;;;;;;;;;
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

;        Version + Name + trace
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='wfccd_analyse', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.0)', /align_center)

;;;;;;;;;;;;;;;;;;
; Field and Mask
  state.field_id = cw_field(toolbar, value='', $
                             /return_events, xsize=15,$
                             title='Field', UVALUE='FIELD')
  state.exp_id = cw_field(toolbar, value=0., /float, $
                             /return_events, xsize=6,$
                             title='Exp', UVALUE='EXP')
;  state.msk_id = cw_field(toolbar, value='', $
;                             /return_events, xsize=5,$
;                             title='Mask', UVALUE='MASK')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      DRAW
  listbase = WIDGET_BASE(state.base_id, /row, /frame, /base_align_center,$
                              uvalue='LISTS')

; OBJ List
  wfccd_analyse_getobjlst, state
  state.obj_id = WIDGET_LIST(listbase, VALUE=state.obj_list[0:state.nobj_list-1], $
                              uvalue='OBJLIST', ysize = 20)
  widget_control, state.obj_id, set_list_select=-1

; OBJNM List
  state.objnm_id = WIDGET_LIST(listbase, $
                               VALUE=[' '], xsize=10, $
                               uvalue='OBJNMLIST', ysize = 20)
  widget_control, state.objnm_id, set_list_select=-1

; EXSPR/OBJ List
  state.expsr_id = WIDGET_LIST(listbase, $
                               VALUE=[' '], xsize=17, $
                               uvalue='EXPLIST', ysize = 20)
  widget_control, state.objnm_id, set_list_select=-1


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      BUTTONS
  butbase = widget_base(listbase, /column, /align_center)
  plt = WIDGET_BUTTON(butbase, value='PLT',uvalue='PLT')
  zhist = WIDGET_BUTTON(butbase, value='ZHIST',uvalue='ZHIST')
  chk = WIDGET_BUTTON(butbase, value='CHK',uvalue='CHK')
  edit = WIDGET_BUTTON(butbase, value='EDIT',uvalue='EDIT')
  img = WIDGET_BUTTON(butbase, value='IMG',uvalue='IMG')
  redo = WIDGET_BUTTON(butbase, value='REDO',uvalue='REDO')
  done = WIDGET_BUTTON(butbase, value='DONE',uvalue='DONE')

; ZHIST
  zfindbase = widget_base(listbase, /column, /align_center, /frame)
  wvmnxbase = widget_base(zfindbase, /row, /align_center)
  wvmnxbut = CW_BGROUP(wvmnxbase, ['WVMNX','PLOT', 'ZMIN', 'ZMAX'], $
                       /nonexclusive, $
                       row=4, /return_index, $
                       uvalue='ZFINDBUT')
  wvmnxvalb = widget_base(wvmnxbase, /column, /align_center)
  wvmnvalid = cw_field(wvmnxvalb, value=3800., /float, $
                             /return_events, xsize=8,$
                             title='WVMN:', UVALUE='WVMNVAL')
  wvmxvalid = cw_field(wvmnxvalb, value=9000., /float, $
                             /return_events, xsize=8,$
                             title='WVMX:', UVALUE='WVMXVAL')
  zminid = cw_field(wvmnxvalb, value=0., /float, $
                             /return_events, xsize=8,$
                             title='ZMIN:', UVALUE='ZMINVAL')
  zmaxid = cw_field(wvmnxvalb, value=0.5, /float, $
                             /return_events, xsize=8,$
                             title='ZMAX:', UVALUE='ZMAXVAL')
  zfind = WIDGET_BUTTON(zfindbase, value='ZFIND',uvalue='ZFIND', xsize=10)

  
;;;;;;;;;
;;  EDIT OBJSTR
  objstrbase = widget_base(listbase, /column, /align_center, /frame)
  editobj = WIDGET_BUTTON(objstrbase, value='EDITOBJ',uvalue='EDITOBJ')
  objbut = CW_BGROUP(objstrbase, ['NOSLIT','NOOBJ'], $
                       /nonexclusive, $
                       row=2, /return_index, $
                       uvalue='OBJBUT')
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
;  strhelp = strarr(50)
;  strhelp = ['   Help Menu   ',$
;             'LMB - Truncate/Extend trace', $ 
;             'RMB - Contrast/Brightness', $
;             'CMB/CMB - Zoom' $ 
;             ]
;  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=4,$
;                             /scroll)

; Realize
  WIDGET_CONTROL, base, /realize

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'wfccd_analyse', base, /no_block

  !P.MULTI= [0,1,1]
  return
end

