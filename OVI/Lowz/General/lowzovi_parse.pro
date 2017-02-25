;+ 
; NAME:
; lowzovi_parse   
;   Version 1.1
;
; PURPOSE:
;  Launches a GUI to allow the user to interactively parse all of the
;  galaxies in a specific structure file.
; 
; CALLING SEQUENCE:
; lowzovi_parse, gal, prs_state, RMAX=, ZMNX=, SPEC=, FLGSCI=, 
;    SURVEY=, /NOINTER, INDX=, MAX_RHO=
;
; INPUTS:
;   gal       -- Galaxy structure
; [prs_state] --
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  RMAX=  -- Maximum magnitude to consider [default: 99.9]
;  ZMNX=  -- 2-element array specifying max min of redshift interval
;  /SPEC  -- Parse all objects with a spectrum
;  FLGSCI= -- Sci flag (1=science target)
;  /SURVEY -- In the main survey?
;  /NOINTER -- Do not launch the GUI
;  MAX_RHO= -- Maximum impact parameter
;
; OPTIONAL OUTPUTS:
;  INDX=  -- Indices of galaxies matching the criteria
;
; COMMENTS:
;
; EXAMPLES:
;   lowzovi_parse, ovi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lowzovi_parse_initcmmn

common lowzovi_parse_common, lzoviparse_gal, lzoviparse_indx, lzoviparse_state

end

;;;;
; Events
;;;;

pro lowzovi_parse_event, ev

common lowzovi_parse_common

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'LIST': 
      'BUTTON': begin
          widget_control, state.button_id, get_value=buttons
          state.flg_rmag = buttons[0]
          state.flg_impact = buttons[1]
          state.flg_zmnx = buttons[2]
          state.flg_spec = buttons[3]
          state.flg_z = buttons[4]
          state.flg_sci = buttons[5]
          state.flg_survey = buttons[n_elements(buttons)-1]
          lowzovi_parse_updlist, state
      end
      'RHO': begin
          widget_control, state.Rho_id, get_value=str_rho
          state.max_rho = float(str_rho)
          if state.flg_impact EQ 1 then lowzovi_parse_updlist, state
      end
      ; zmin, zmax
      'ZMIN': begin
          widget_control, state.zmin_id, get_value=str_zmin
          state.zmin = float(str_zmin)
          if state.flg_zmnx EQ 1 then lowzovi_parse_updlist, state
      end
      'ZMAX': begin
          widget_control, state.zmax_id, get_value=str_zmax
          state.zmax = float(str_zmax)
          if state.flg_zmnx EQ 1 then lowzovi_parse_updlist, state
      end
      'RMAG': begin
          widget_control, state.rmag_id, get_value=str_rmag
          state.rmax = float(str_rmag)
          if state.flg_rmag EQ 1 then lowzovi_parse_updlist, state
      end
      'DONE' : begin
          lzoviparse_state = state
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
pro lowzovi_parse_mklist, state, list
common lowzovi_parse_common

  if lzoviparse_indx[0] EQ -1 then list = ' No Obj satisfy those criterion!! ' $
  else begin
      ;; Make array
      list = strarr(n_elements(lzoviparse_indx))

      case state.nfilt of 
          2: fmt = '(a6,1x,i1,1x,a23,1x,f6.1,a3,1x,f8.5)'
          else: stop
      endcase

      ;; Loop
      for q=0L,n_elements(list)-1 do begin
          jj = lzoviparse_indx[q]
          ;; Mag
          if lzoviparse_gal[jj].flg_anly MOD 2 EQ 0 $
            then magstr = ' ' else begin
              magstr = ' '
              for i=0L,state.nfilt-1 do $
                magstr = magstr+string(lzoviparse_gal[jj].mag[i],FORMAT='(f5.2)')$
                +string(lzoviparse_gal[jj].magerr[i],FORMAT='(f6.2)')+' '
              magstr = strtrim(magstr,2)
          endelse
          ;; Calculate rho
          rho = sqrt(lzoviparse_gal[jj].ra^2 + lzoviparse_gal[jj].dec^2)
          ;; List
          list[q] = string(strtrim(lzoviparse_gal[jj].id,2)+ $
                           strtrim(lzoviparse_gal[jj].obj_id,2),$
            lzoviparse_gal[jj].flg_anly, $
            magstr, $
            rho, $
            lzoviparse_gal[jj].gal_type, $
            lzoviparse_gal[jj].z, $
            FORMAT=fmt)
      endfor
  endelse

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Update List
pro lowzovi_parse_updlist, state
common lowzovi_parse_common

  lowzovi_parse_setcindx, state
  lowzovi_parse_mklist, state, list
  widget_control, state.list_id, set_value=list
  widget_control, state.list_id, set_list_select=state.lindx
  nobj = n_elements([lzoviparse_indx])
  widget_control, state.nobj_id, set_value=nobj
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Change type
pro lowzovi_parse_setcindx, state
common lowzovi_parse_common

  ;; Science
  if state.flg_sci EQ 0 then lzoviparse_indx = lindgen(n_elements(lzoviparse_gal)) $
  else lzoviparse_indx = where(lzoviparse_gal.obj_id EQ 'a')

  ;; Spec
  if state.flg_spec EQ 1L then begin
      a = where( lzoviparse_gal[lzoviparse_indx].flg_anly MOD 4 GT 1, na)
      if na EQ 0 then lzoviparse_indx = -1 else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

  ;; Rmag
  if state.flg_rmag EQ 1 then begin
      case state.nfilt of 
          2: a = where( lzoviparse_gal[lzoviparse_indx].mag[1] LT state.rmax, na)
          else: stop
      endcase
      if na EQ 0 then lzoviparse_indx = -1 $
      else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

  ;; Redshift
  if state.flg_z EQ 1L then begin
      a = where( lzoviparse_gal[lzoviparse_indx].flg_anly MOD 8 GT 3, na)
      if na EQ 0 then lzoviparse_indx = -1 else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

  ;; Survey
  if state.flg_survey EQ 1 then begin
      a = where( lzoviparse_gal[lzoviparse_indx].flg_survey MOD 2 EQ 1, na)
      if na EQ 0 then lzoviparse_indx = -1 $
      else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

  ;; Impact
  if state.flg_impact EQ 1 then begin
      rho = sqrt(lzoviparse_gal[lzoviparse_indx].ra^2 + $
             lzoviparse_gal[lzoviparse_indx].dec^2)
      a = where( rho LT state.max_rho, na)
      if na EQ 0 then lzoviparse_indx = -1 $
      else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

  ;; Zmin, zmax
  if state.flg_zmnx EQ 1 then begin
      a = where( lzoviparse_gal[lzoviparse_indx].z LT state.zmax $
                 AND lzoviparse_gal[lzoviparse_indx].z GT state.zmin, na)
      if na EQ 0 then lzoviparse_indx = -1 $
      else lzoviparse_indx = lzoviparse_indx[a]
  endif
  if lzoviparse_indx[0] EQ -1 then return

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

pro lowzovi_parse, gal, prs_state, RMAX=rmax, ZMNX=zmnx, $
                   SPEC=spec, FLGSCI=flgsci, SURVEY=survey, $
                   NOINTER=nointer, INDX=indx, MAX_RHO=max_rho

common lowzovi_parse_common

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'lowzovi_parse, ovi, prs_state, RMAX=, ZMNX=, SPEC=, FLGSCI=, SURVEY=, /NOINTER, INDX=, MAX_RHO= [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 200
  if not keyword_set( YOFFSET ) then yoffset = 200
  if not keyword_set( MXY ) then mxy = 30L
  if not keyword_set( STDFONT ) then stdfont = '8x13'
  if not keyword_set( LST_FONT ) then lst_font = '6x10'

; Init common
  lowzovi_parse_initcmmn
  lzoviparse_gal = gal

;    STATE
  state = { $
            indx: -1L, $
            lindx: -1L, $
            flg_sci: 1L, $  ; 0L means include serendip
            flg_spec: 0L, $ ; 1L means spec only
            flg_survey: 0L, $ ; 1L means survey
            flg_z: 0L, $    ; 1L means good z only
            flg_zmnx: 0L, $    ; 1L means set within zmin,zmax
            flg_impact: 0L, $
            max_rho: 0., $
            zmax: 99.9, $
            zmin: 0.0, $
            nfilt: 2L, $
            flg_rmag: 0, $
            rmax: 0., $
            list_id: 0L, $
            button_id: 0L, $
            zmax_id: 0L, $
            zmin_id: 0L, $
            Rho_id: 0L, $
            Rmag_id: 0L, $
            nobj_id: 0L, $
            base_id: 0L $
          }

;; INIT VALUES
  if keyword_set( PRS_STATE ) then state = prs_state
  if keyword_set( RMAX ) then begin
      state.flg_rmag = 1
      state.rmax = rmax
  endif else state.rmax = 99.9

  if keyword_set( MAX_RHO ) then begin
      state.flg_impact = 1
      state.max_rho = max_rho
  endif else state.max_rho = 9999.9

  if keyword_set( SPEC ) then state.flg_spec = spec
  if keyword_set( FLGSCI ) then state.flg_sci = flgsci
  ;; Zmin max
  if keyword_set( ZMNX ) then begin
      state.zmin = zmnx[0]
      state.zmax = zmnx[1]
      state.flg_zmnx = 1L
  endif
  


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lists
;;
  ;; Setup initial list
  lowzovi_parse_setcindx, state
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    WIDGET
  if not keyword_set( NOINTER ) then begin
      base = WIDGET_BASE( title = 'lowzovi_parse', /row, $
                          xoffset=xoffset,yoffset=yoffset)
      state.base_id = base
; TITLE
;  titl_id = widget_label(base, value='lowzovi_parse')
      
      lowzovi_parse_mklist, state, list
;
      ysz = mxy < n_elements(list)
      state.list_id = widget_list(base, value=list, xsize=60L, ysize=ysz, $
                                  uvalue = 'LIST', FONT=lst_font)
;  widget_control, state.objnm_id, set_list_select=-1
      
; Toolbar
      toolbar = WIDGET_BASE( state.base_id, /frame, /row, title='Parsing') 
      parsebar = WIDGET_BASE( toolbar, /frame, /row, title='Parsing') 
      

; Buttons
      state.button_id = CW_BGROUP(parsebar, ['Rmag', 'Impact', 'zmin,zmax', $
                                            'flg_spec', $
                                            'flg_z', $
                                            'flg_sci', $
                                            'flg_survey'], $
                                  row=7, set_value=[state.flg_rmag, $
                                                    state.flg_impact, $
                                                    state.flg_zmnx, $
                                                    state.flg_spec, $
                                                    state.flg_z, $
                                                    state.flg_survey], $
                                  /nonexclusive, font=stdfont, uvalue='BUTTON')
; Fields
      field_id = WIDGET_BASE(parsebar, /column) 
      state.Rmag_id = widget_text(field_id, font=stdfont, value=strtrim(state.rmax,2),$
                                  /editable, XSIZE=6L, UVALUE='RMAG')
      state.Rho_id = widget_text(field_id, font=stdfont, $
                                 value=strtrim(state.max_rho,2),$
                                 /editable, XSIZE=6L, UVALUE='RHO')
      base_zmnx_id = WIDGET_BASE(field_id, /column) 
      state.zmin_id = widget_text(base_zmnx_id, font=stdfont, $
                                 value=strtrim(state.zmin,2),$
                                 /editable, XSIZE=6L, UVALUE='ZMIN')
      state.zmax_id = widget_text(base_zmnx_id, font=stdfont, $
                                 value=strtrim(state.zmax,2),$
                                 /editable, XSIZE=6L, UVALUE='ZMAX')
      
;   DONE
      state.nobj_id = cw_field(toolbar, /column, title='Nobj', value=0L, xsize=5L)
      done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
      
; Realize
      WIDGET_CONTROL, base, /realize
      WIDGET_CONTROL, base, set_uvalue = state, /no_copy
      
      xmanager, 'lowzovi_parse', base
      prs_state = temporary(lzoviparse_state)
  endif else prs_state = temporary(state)

; Finish
  
  indx = temporary(lzoviparse_indx)
  delvarx, lzoviparse_gal
  
  return
end
