;+ 
; NAME:
; x_guinum   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs a number from the user
;
; CALLING SEQUENCE:
;   
;   num = x_guinum(flg)
;
; INPUTS:
;   flg = Return (0: float, 1: double, 2: Long)
;
; RETURNS:
;   num - Number
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  TITLE= -- Title for the GUI [default: 'Num']
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   num = x_guinum( 0 )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_guinum, flg, TITLE=title

;common x_slctline_ans

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'line = x_guinum(flg, TITLE=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  device, get_screen_size=ssz
  if not keyword_set( XOFFSET ) then    xOFFSET = round(ssz[0]*0.5)
  if not keyword_set( YOFFSET ) then    yOFFSET = round(ssz[1]*0.5)
  if not keyword_set( TITLE ) then title = 'Num'

;    

;    WIDGET
  base = WIDGET_BASE( title = 'x_guinum', /column, $
                      xoffset=xoffset,yoffset=yoffset)

  
; PD Lists
  cw_id = cw_field(base, value=-999., title=title, uvalue = 'cwnum', $
                  /return_events, /floating)

; Realize
  WIDGET_CONTROL, base, /realize
  
; LOOP
  num = 0.d
  repeat begin
      ev = widget_event(base)
      num = double(ev.value)
  end until ev.update EQ 1
  
  WIDGET_CONTROL, /destroy, base

; FIND THE RIGHT WAV

  case flg of 
      0: return, float(num)
      1: return, num
      2: return, long(num)
      else: return, num
  endcase
end
