;+ 
; NAME:
; x_guistring   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs input from the user
;
; CALLING SEQUENCE:
;   
;   string = x_guistring(title)
;
; INPUTS:
;   title - Title
;
; RETURNS:
;   string - String
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_guistring, title

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'value = x_guistring(title) [v1.0]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 500
  if not keyword_set( YOFFSET ) then yoffset = 500

;    

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
