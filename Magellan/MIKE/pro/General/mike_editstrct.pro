;+ 
; NAME:
; mike_editstrct   
;   Version 1.1
;
; PURPOSE:
;    Launches a gui to edit the MIKE structure
;
; CALLING SEQUENCE:
;   
;   mike_editstrct, mike
;
; INPUTS:
;  mike  --  ESI structure
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
;   mike_editstrct, mike
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_editstrct_initcmmn

common mike_editstrct_common, cmm_mike, cmm_flg, cmm_indx
cmm_indx = -1

end

;;;;
; Events
;;;;

pro mike_editstrct_event, ev

common mike_editstrct_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'LIST': cmm_indx = widget_info(state.list_id, /list_select)
      'FANLY': mike_editstrct_flipanly, state
      'CHGNAME': mike_editstrct_chgname, state
      'CHGTYPE': mike_editstrct_chgtype, state
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
pro mike_editstrct_mklist, list
common mike_editstrct_common

  ;; Make array
  list = strarr(n_elements(cmm_mike))

  ;; Loop
  for q=0L,n_elements(list)-1 do begin
      if cmm_mike[q].obj_id GE 0L then begin
          if cmm_mike[q].obj_id LT 10 then id_str = '0'+strtrim(cmm_mike[q].obj_id,2) $
          else id_str = strtrim(cmm_mike[q].obj_id,2)
      endif else id_str = '-1'
;      case cmm_mike[q].mode of
;          0: mode = 'IMG'
;          1: mode = 'LWD'
;          2: mode = 'ECH'
;          else:
;      endcase
      list[q] = string(q,$
        cmm_mike[q].img_root, $
        cmm_mike[q].flg_anly, $
        cmm_mike[q].Obj, $
        id_str, $
        cmm_mike[q].setup, $
        cmm_mike[q].type, $
        long(cmm_mike[q].exp), $
        cmm_mike[q].colbin, $
        cmm_mike[q].rowbin, $
        cmm_mike[q].slit, $
        cmm_mike[q].arclamp, $
        cmm_mike[q].qtzlamp, $
        FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a4,i3,a5,i5,2i2,f5.2,i2,i2)')
  endfor

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro mike_editstrct_updlist, state
common mike_editstrct_common

  mike_editstrct_mklist, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=cmm_indx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Flip Anly
pro mike_editstrct_flipanly, state
common mike_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mike_editstrct_flipanly: Select entry first!'
      return
  endif

  ;; Flip
  for q=0L,n_elements(cmm_indx)-1 do begin
      if  cmm_mike[cmm_indx[q]].flg_anly EQ 1 then cmm_mike[cmm_indx[q]].flg_anly = 0 $
      else cmm_mike[cmm_indx[q]].flg_anly = 1
  endfor

  ;; Update
  mike_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro mike_editstrct_chgtype, state
common mike_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mike_editstrct_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['OBJ','STD','DRK','TWI','MFLT','TFLT','ARC','ZRO'])

  ;; Change
  cmm_mike[cmm_indx].type = type

  ;; Update
  mike_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change name
pro mike_editstrct_chgname, state
common mike_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mike_editstrct_chgname: Select entry first!'
      return
  endif

  ;; Get type
  name = x_guistring('Name')

  ;; Change
  cmm_mike[cmm_indx].Obj = name

  ;; Update
  mike_editstrct_updlist, state

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

pro mike_editstrct, mike

common mike_editstrct_common

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'mike_editstrct, mike   [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( LSTFONT ) then lstfont = '6x10'

; Init common
  mike_editstrct_initcmmn
  cmm_mike = mike

;    STATE
  state = { $
            indx: -1L, $
            list_id: 0L, $
            base_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'mike_editstrct', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='mike_editstrct')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

; Buttons
  flip_anly = WIDGET_BUTTON(toolbar, value='Flip Anly',uvalue='FANLY', /align_right)
  chg_type = WIDGET_BUTTON(toolbar, value='Chng Type',uvalue='CHGTYPE', /align_right)
  chg_name = WIDGET_BUTTON(toolbar, value='Chng Name',uvalue='CHGNAME', /align_right)

; Create the list
  mike_editstrct_mklist, list

; PD Lists
  ysz = 30 < n_elements(list)
  state.list_id = widget_list(base, value=list, xsize=60L, ysize=ysz, $
                              FONT=lstfont, uvalue = 'LIST', /MULTIPLE)

;   DONE
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

; Realize
  WIDGET_CONTROL, base, /realize
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  xmanager, 'mike_editstrct', base

; Finish
  mike = temporary(cmm_mike)
  delvarx, cmm_indx
  
  return

end
