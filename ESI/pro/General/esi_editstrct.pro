;+ 
; NAME:
; esi_editstrct   
;   Version 1.1
;
; PURPOSE:
;    Launches a gui to edit the ESI structure
;
; CALLING SEQUENCE:
;   
;   esi_editstrct, esi
;
; INPUTS:
;  esi  --  ESI structure
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
;   esi_editstrct, esi
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

pro esi_editstrct_initcmmn

common esi_editstrct_common, cmm_esi, cmm_flg, cmm_indx
cmm_indx = -1

end

;;;;
; Events
;;;;

pro esi_editstrct_event, ev

common esi_editstrct_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'LIST': cmm_indx = widget_info(state.list_id, /list_select)
      'FANLY': esi_editstrct_flipanly, state
      'CHGNAME': esi_editstrct_chgname, state
      'CHGTYPE': esi_editstrct_chgtype, state
      'CHGMODE': esi_editstrct_chgmode, state
      'CHGSLIT': esi_editstrct_chgslit, state
      'CHGOBJID': esi_editstrct_chgobjid, state
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
pro esi_editstrct_mklist, list
common esi_editstrct_common

  ;; Make array
  list = strarr(n_elements(cmm_esi))

  ;; Loop
  for q=0L,n_elements(list)-1 do begin
      if cmm_esi[q].obj_id GE 0L then begin
          if cmm_esi[q].obj_id LT 10 then id_str = '0'+strtrim(cmm_esi[q].obj_id,2) $
          else id_str = strtrim(cmm_esi[q].obj_id,2)
      endif else id_str = '-1'
      case cmm_esi[q].mode of
          0: mode = 'IMG'
          1: mode = 'LWD'
          2: mode = 'ECH'
          else:
      endcase
      list[q] = string(cmm_esi[q].frame, $
                       cmm_esi[q].img_root, $
                       cmm_esi[q].flg_anly, $
                       cmm_esi[q].Obj, $
                       id_str, $
                       mode, $
                       cmm_esi[q].cbin, $
                       cmm_esi[q].rbin, $
                       cmm_esi[q].type, $
                       long(cmm_esi[q].exp), $
                       cmm_esi[q].slit, $
                       cmm_esi[q].arclamp, $
                       cmm_esi[q].qtzlamp, $
                       cmm_esi[q].imfilt, $
                       FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a3,1x,a4,1x,' + $
                       'i2,1x,i2,a5,i5,f5.2,i2,i2,a2,f5.2,2i2)')
  endfor

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro esi_editstrct_updlist, state
common esi_editstrct_common

  esi_editstrct_mklist, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=cmm_indx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Flip Anly
pro esi_editstrct_flipanly, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_flipanly: Select entry first!'
      return
  endif

  ;; Flip
  for q=0L,n_elements(cmm_indx)-1 do begin
      if  cmm_esi[cmm_indx[q]].flg_anly EQ 1 then cmm_esi[cmm_indx[q]].flg_anly = 0 $
      else cmm_esi[cmm_indx[q]].flg_anly = 1
  endfor

  ;; Update
  esi_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro esi_editstrct_chgtype, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['OBJ','STD','DRK','TWI','IFLT','DFLT','ARC','ZRO'])

  ;; Change
  cmm_esi[cmm_indx].type = type

  ;; Update
  esi_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change name
pro esi_editstrct_chgname, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_chgname: Select entry first!'
      return
  endif

  ;; Get type
  name = x_guistring('Name')

  ;; Change
  cmm_esi[cmm_indx].Obj = name

  ;; Update
  esi_editstrct_updlist, state

  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change name
pro esi_editstrct_chgobjid, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_chgobjid: Select entry first!'
      return
  endif

  ;; Get type
  objid = x_guinum('Obj_id')

  ;; Change
  cmm_esi[cmm_indx].Obj_id = objid

  ;; Update
  esi_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change mode
pro esi_editstrct_chgmode, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_chgmode: Select entry first!'
      return
  endif

  ;; Get type
  mode = x_guilist(['IMG','LWD','ECH'], indx=indx)

  ;; Change
  cmm_esi[cmm_indx].mode = indx

  ;; Update
  esi_editstrct_updlist, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change slit
pro esi_editstrct_chgslit, state
common esi_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'esi_editstrct_chgslit: Select entry first!'
      return
  endif

  ;; Get type
  slit = x_guilist(['0.33','0.5','0.75','1.0','1.25','6.0'], indx=indx)

  ;; Change
  case indx of 
      0: cmm_esi[cmm_indx].slit = 0.33
      1: cmm_esi[cmm_indx].slit = 0.50
      2: cmm_esi[cmm_indx].slit = 0.75
      3: cmm_esi[cmm_indx].slit = 1.00
      4: cmm_esi[cmm_indx].slit = 1.25
      5: cmm_esi[cmm_indx].slit = 6.00
      else: stop
  endcase

  ;; Update
  esi_editstrct_updlist, state

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

pro esi_editstrct, esi

common esi_editstrct_common

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'esi_editstrct, esi   [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( LSTFONT ) then lstfont = '6x10'

; Init common
  esi_editstrct_initcmmn
  cmm_esi = esi

;    STATE
  state = { $
            indx: -1L, $
            list_id: 0L, $
            base_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'esi_editstrct', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='esi_editstrct')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)

; Buttons
  flip_anly = WIDGET_BUTTON(toolbar, value='Flip Anly',uvalue='FANLY', /align_right)
  chg_type = WIDGET_BUTTON(toolbar, value='Chng Type',uvalue='CHGTYPE', /align_right)
  chg_name = WIDGET_BUTTON(toolbar, value='Chng Name',uvalue='CHGNAME', /align_right)
  chg_mode = WIDGET_BUTTON(toolbar, value='Chng Mode',uvalue='CHGMODE', /align_right)
  chg_slit = WIDGET_BUTTON(toolbar, value='Chng Slit',uvalue='CHGSLIT', /align_right)
  chg_obj = WIDGET_BUTTON(toolbar, value='Chng Objid',uvalue='CHGOBJID', /align_right)

; Create the list
  esi_editstrct_mklist, list

; PD Lists
  ysz = 30 < n_elements(list)
  state.list_id = widget_list(base, value=list, xsize=60L, ysize=ysz, $
                              FONT=lstfont, uvalue = 'LIST', /MULTIPLE)

;   DONE
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

; Realize
  WIDGET_CONTROL, base, /realize
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  xmanager, 'esi_editstrct', base

; Finish
  esi = temporary(cmm_esi)
  delvarx, cmm_indx
  
  return

end
