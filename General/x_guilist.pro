;+ 
; NAME:
; x_guilist   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs input from the user
;
; CALLING SEQUENCE:
;   
;   string = x_guilist(list, title)
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
;   string = x_guilist( 'Enter filename: ')
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

function x_guilist, list, title, INDX=indx

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'value = x_guistring(list, [title], INDX=) [v1.0]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 500
  if not keyword_set( YOFFSET ) then yoffset = 500

;    

;    WIDGET
  base = WIDGET_BASE( title = 'x_guilist', /column, $
                      xoffset=xoffset,yoffset=yoffset)

; TITLE
  if keyword_set( TITLE ) then $
    titl_id = widget_label(base, value=title)

  
; PD Lists
  ysz = 20 < n_elements(list)
  wdgt_id = widget_list(base, value=list, ysize=ysz, uvalue = 'widglist')

; Realize
  WIDGET_CONTROL, base, /realize
  
; LOOP
  click = 0L
  repeat begin
      ev = widget_event(base)
      indx = ev.index
      click = ev.clicks
  end until click NE 0L
  
  WIDGET_CONTROL, /destroy, base

  return, list[indx]

end
