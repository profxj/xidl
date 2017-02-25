;+ 
; NAME:
; x_helpwidg   
;   Version 1.1
;
; PURPOSE:
;    Creates a widget to display a help list which is simply an 
;   array of strings that describe a GUI [or anything].
;
; CALLING SEQUENCE:
;   x_helpwidg, help
;
; INPUTS:
;   help       - String array of help statements
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_helpwidg, ['Help', 'h: Show this']
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
	; REVISION HISTORY:
;   28-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Events
;;;;

pro x_helpwidg_event, ev


  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
      else:
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_helpwidg, help, XSIZE=xsize, YSIZE=ysize

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_helpwidg, help, XSIZE=, YSIZE= [V1.1]'
    return
  endif 


;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 50
  if not keyword_set( YSIZE ) then    ysize = 10

;    STATE
  state = { $
            base_id: 0L, $      ; Widgets
            help_id: 0L $
          }
      
;    WIDGET
  base = WIDGET_BASE( title = 'x_help: Help', /column, /align_center)
  state.base_id = base
  
;        Help
  help_text_id = widget_text(base, value=help, xsize=xsize, ysize=ysize,$
                             /scroll)
;      Done
  done = WIDGET_BUTTON(base, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  

; Update
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  
; Send to the xmanager
  xmanager, 'x_helpwidg', base, /no_block

  return

end
