;+ 
; NAME:
; x_guistring   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs a string from the user
;
; CALLING SEQUENCE:
;   value = x_guistring(title)
;
; INPUTS:
;   title - Title
;
; RETURNS:
;   value - String
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
;   string = x_guistring( 'Enter filename: ')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_guistring, title

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'value = x_guistring(title) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  device, get_screen_size=ssz
  if not keyword_set( XOFFSET ) then    xOFFSET = ssz[0]*0.3
  if not keyword_set( YOFFSET ) then    yOFFSET = ssz[1]*0.2

;    WIDGET
  base = WIDGET_BASE( title = 'x_guistring', /column, $
                      xoffset=xoffset,yoffset=yoffset)

  
; PD Lists
  cw_id = cw_field(base, value='', title=title, uvalue = 'cwstring', $
                  /column, /return_events)

; Realize
  WIDGET_CONTROL, base, /realize
  
; LOOP
  string = ''
  repeat begin
      ev = widget_event(base)
      string = ev.value
  end until string NE ''
  
  WIDGET_CONTROL, /destroy, base

  return, string

end
