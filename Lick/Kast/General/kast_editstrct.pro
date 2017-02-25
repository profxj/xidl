;+ 
; NAME:
; kast_editstrct   
;   Version 1.1
;
; PURPOSE:
;    Launches a gui to edit the Kast IDL structure
;
; CALLING SEQUENCE:
;   kast_editstrct, kast
;
; INPUTS:
;  kast  --  Kast IDL structure
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
;   kast_editstrct, kast
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_editstrct_initcmmn

common kast_editstrct_common, cmm_kast, cmm_flg, cmm_indx
cmm_indx = -1

end

;;;;
; Events
;;;;

pro kast_editstrct_event, ev

common kast_editstrct_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'LIST': cmm_indx = widget_info(state.list_id, /list_select)
      'FANLY': kast_editstrct_flipanly, state
      'CHGNAME': kast_editstrct_chgname, state
      'CHGTYPE': kast_editstrct_chgtype, state
      'CHGMODE': kast_editstrct_chgmode, state
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
pro kast_editstrct_mklist, list
common kast_editstrct_common

  ;; Make array
  list = strarr(n_elements(cmm_kast))

  ;; Loop
  for q=0L,n_elements(list)-1 do begin
      if cmm_kast[q].obj_id GE 0L then begin
          if cmm_kast[q].obj_id LT 10 then $
            id_str = '0'+strtrim(cmm_kast[q].obj_id,2) $
          else id_str = strtrim(cmm_kast[q].obj_id,2)
      endif else id_str = '-1'
      if cmm_kast[q].setup GE 0L then begin
          if cmm_kast[q].setup LT 10 then $
            setup_str = '0'+strtrim(cmm_kast[q].setup,2) $
          else setup_str = strtrim(cmm_kast[q].setup,2)
      endif else id_str = '-1'
      case cmm_kast[q].mode of
          0: mode = 'IMG'
          1: mode = 'SPEC'
;          2: mode = 'ECH'
          else:
      endcase
      list[q] = string(cmm_kast[q].frame, $
        cmm_kast[q].img_root, $
        cmm_kast[q].side, $
        cmm_kast[q].flg_anly, $
        cmm_kast[q].Obj, $
        setup_str, $
        id_str, $
        mode, $
        cmm_kast[q].type, $
        cmm_kast[q].grising, $
        long(cmm_kast[q].exp), $
        cmm_kast[q].slit, $
        FORMAT='(i4,1x,a8,i2,1x,i1,1x,a12,1x,a3,1x,a3,1x,a4,a5,a10,i5,f8.0)')
  endfor

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro kast_editstrct_updlist, state
common kast_editstrct_common

  kast_editstrct_mklist, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=cmm_indx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Flip Anly
pro kast_editstrct_flipanly, state
common kast_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'kast_editstrct_flipanly: Select entry first!'
      return
  endif

  ;; Flip
  for q=0L,n_elements(cmm_indx)-1 do begin
      if  cmm_kast[cmm_indx[q]].flg_anly EQ 1 then cmm_kast[cmm_indx[q]].flg_anly = 0 $
      else cmm_kast[cmm_indx[q]].flg_anly = 1
  endfor

  ;; Update
  kast_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro kast_editstrct_chgtype, state
common kast_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'kast_editstrct_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['OBJ','STD','DRK','TWI','QTZ','ARC','ZRO','ACQ','SLT'])

  ;; Change
  cmm_kast[cmm_indx].type = type

  ;; Update
  kast_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change name
pro kast_editstrct_chgname, state
common kast_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'kast_editstrct_chgname: Select entry first!'
      return
  endif

  ;; Get type
  name = x_guistring('Name')

  ;; Change
  cmm_kast[cmm_indx].Obj = name

  ;; Update
  kast_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change mode
pro kast_editstrct_chgmode, state
common kast_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'kast_editstrct_chgmode: Select entry first!'
      return
  endif

  ;; Get type
  mode = x_guilist(['IMG','LWD','ECH'], indx=indx)

  ;; Change
  cmm_kast[cmm_indx].mode = indx

  ;; Update
  kast_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_editstrct, kast

common kast_editstrct_common

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'kast_editstrct, kast   [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( LSTFONT ) then lstfont = '6x10'

; Init common
  kast_editstrct_initcmmn
  cmm_kast = kast

;    STATE
  state = { $
            indx: -1L, $
            list_id: 0L, $
            base_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'kast_editstrct', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='kast_editstrct')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

; Buttons
  flip_anly = WIDGET_BUTTON(toolbar, value='Flip Anly',uvalue='FANLY', /align_right)
  chg_type = WIDGET_BUTTON(toolbar, value='Chng Type',uvalue='CHGTYPE', /align_right)
  chg_name = WIDGET_BUTTON(toolbar, value='Chng Name',uvalue='CHGNAME', /align_right)
  chg_mode = WIDGET_BUTTON(toolbar, value='Chng Mode',uvalue='CHGMODE', /align_right)

; Create the list
  kast_editstrct_mklist, list

; PD Lists
  ysz = 30 < n_elements(list)
  state.list_id = widget_list(base, value=list, xsize=80L, ysize=ysz, $
                              FONT=lstfont, uvalue = 'LIST', /MULTIPLE)

;   DONE
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

; Realize
  WIDGET_CONTROL, base, /realize
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  xmanager, 'kast_editstrct', base

; Finish
  kast = temporary(cmm_kast)
  delvarx, cmm_indx
  
  return

end
