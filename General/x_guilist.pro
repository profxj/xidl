;+ 
; NAME:
; x_guilist   
;   Version 1.1
;
; PURPOSE:
;    Launches a cw_field, waits for the user to make a choice and
;    then returns the index and value of the choice.  Routine will
;    block any other GUI.
;
; CALLING SEQUENCE:
;   
;   value = x_guilist(list, title)
;
; INPUTS:
;   list  - String array
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
;   value = x_guilist( ['a', 'b', 'c'], INDX=indx )
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

function x_guilist, list, title, INDX=indx, MAXY=maxy

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'value = x_guistring(list, [title], INDX=, MAXY=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( MAXY ) then    maxy = 20
  device, get_screen_size=ssz
  if not keyword_set( XOFFSET ) then    xOFFSET = ssz[0]*0.3
  if not keyword_set( YOFFSET ) then    yOFFSET = ssz[1]*0.2

;    WIDGET
  base = WIDGET_BASE( title = 'x_guilist', /column, $
                      xoffset=xoffset,yoffset=yoffset)

; TITLE
  if keyword_set( TITLE ) then $
    titl_id = widget_label(base, value=title)

  
; PD Lists
  ysz = MAXY < n_elements(list)
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
