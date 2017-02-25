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

pro mage_editstrct_initcmmn

common mage_editstrct_common, cmm_mage, cmm_flg, cmm_indx
cmm_indx = -1

end

;;;;
; Events
;;;;

pro mage_editstrct_event, ev

common mage_editstrct_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'LIST': cmm_indx = widget_info(state.list_id, /list_select)
      'CHGOBJID': mage_editstrct_chgobjid, state
      'CHGNAME': mage_editstrct_chgname, state
      'CHGTYPE': mage_editstrct_chgtype, state
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
pro mage_editstrct_mklist, list
common mage_editstrct_common

  ;; Make array
  list = strarr(n_elements(cmm_mage))

  ;; Loop
  for q=0L,n_elements(list)-1 do begin

      list[q] = string(q,$
        cmm_mage[q].fitsfile, $
        cmm_mage[q].Object, $
        cmm_mage[q].exptype, $
        long(cmm_mage[q].exptime), $
        cmm_mage[q].slit, $
        cmm_mage[q].obj_id, $
        FORMAT='(i4,1x,a15,1x,a20,1x,a9,i4,2x,a5,i5)')
   endfor

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro mage_editstrct_updlist, state
common mage_editstrct_common

  mage_editstrct_mklist, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=cmm_indx
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Flip Anly
;; pro mage_editstrct_flipanly, state
;; common mage_editstrct_common

;;   Check indx
;;   if cmm_indx[0] LT 0 then begin
;;       print, 'mike_editstrct_flipanly: Select entry first!'
;;       return
;;   endif

;;   Flip
;;   for q=0L,n_elements(cmm_indx)-1 do begin
;;       if  cmm_mage[cmm_indx[q]].flg_anly EQ 1 then cmm_mage[cmm_indx[q]].flg_anly = 0 $
;;       else cmm_mage[cmm_indx[q]].flg_anly = 1
;;   endfor

;;   Update
;;   mage_editstrct_updlist, state

;;   return
;; end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro mage_editstrct_chgtype, state
common mage_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mike_editstrct_chgtype: Select entry first!'
      return
  endif

  ;; Get type
  type = x_guilist(['XE-FLASH','DOMEFLT','SCIENCE','ARC','STD','BRIGHT','ZRO','TRASH'])

  ;; Change
  cmm_mage[cmm_indx].exptype = type

  ;; Update
  mage_editstrct_updlist, state

  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change objid
pro mage_editstrct_chgobjid, state
common mage_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mage_editstrct_objid: Select entry first!'
      return
  endif

  ;; Get type
  obj_id = x_guinum(0,Title='Object ID')

  ;; Change
  cmm_mage[cmm_indx].obj_id = obj_id

  ;; Update
  mage_editstrct_updlist, state

  return
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change name
pro mage_editstrct_chgname, state
common mage_editstrct_common

  ;; Check indx
  if cmm_indx[0] LT 0 then begin
      print, 'mage_editstrct_chgname: Select entry first!'
      return
  endif

  ;; Get type
  name = x_guistring('Name')

  ;; Change
  cmm_mage[cmm_indx].Object = name[0]

  match = where(strtrim(cmm_mage.object) EQ strtrim(name[0]) $
                AND cmm_mage.obj_id NE cmm_mage[cmm_indx].obj_id, nmatch)

  if (nmatch GT 0) then begin
     newid = cmm_mage[match].obj_id
     print, newid[0]
     cmm_mage[cmm_indx].obj_id = newid[0]
     print, "Found another exposure of this object, amending obj_id"
  endif

  ;; Update
  mage_editstrct_updlist, state

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

pro mage_editstrct, filename

common mage_editstrct_common   ;common block for editing mage structure

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'mage_editstrct, fitsfile   [v0.1]'
    return
  endif 

;  Optional Keywords
;  if not keyword_set( XOFFSET ) then xoffset = 200
;  if not keyword_set( YOFFSET ) then yoffset = 200
;  if not keyword_set( LSTFONT ) then lstfont = '6x10'

; Init common
  mage = xmrdfits(filename, 1)
  mage_editstrct_initcmmn
  cmm_mage = mage

;    STATE
  state = { $
            indx: -1L, $
            list_id: 0L, $
            base_id: 0L $
          }

;    WIDGET
  base = WIDGET_BASE( title = 'mage_editstrct', /column, $
                      xoffset=xoffset,yoffset=yoffset)
  state.base_id = base

; TITLE
;  titl_id = widget_label(base, value='mike_editstrct')

; Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center, /align_center)

; Buttons
  chg_objid = WIDGET_BUTTON(toolbar, value='Change ObjID',uvalue='CHGOBJID', /align_right)
  chg_type = WIDGET_BUTTON(toolbar, value='Change Type',uvalue='CHGTYPE', /align_right)
  chg_name = WIDGET_BUTTON(toolbar, value='Change Name',uvalue='CHGNAME', /align_right)

  wlabel = WIDGET_LABEL(state.base_id, /align_left, value='Index   Filename          Object           Type       Exp. Slit  Obj.ID')
; Create the list
  mage_editstrct_mklist, list

; PD Lists
  ysz = 30 < n_elements(list)    ;This takes whatever is smaller, the list or 30
  state.list_id = widget_list(base, value=list, xsize=70L, ysize=ysz, uvalue = 'LIST', /MULTIPLE)

;   DONE
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)

; Realize
  WIDGET_CONTROL, base, /realize
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  xmanager, 'mage_editstrct', base
; Finish
  mage = temporary(cmm_mage)
  delvarx, cmm_indx
  mwrfits, mage, filename, /create
  return

end
