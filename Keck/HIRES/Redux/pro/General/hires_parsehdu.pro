;+ 
; NAME:
; hires_parsehdu
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

function hires_parsehdu, hires, head, keywd

  ;;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'guess = hires_parsehdu(hires, head, keywd) [v1.1]'
      return, -1
  endif 

  ;; Structure
  pos = n_elements(hires)


  ;; Pane
  pane = strtrim(sxpar(head,keywd.pane),2)
  if pane NE 'PANE' then stop


  ;; CCDs
  mosmod = sxpar(head,keywd.mosmod)
  if strpos(mosmod, 'B') NE -1 then begin
      tmp = { hiresstrct }
      tmp.chip = 1
      case strpos(mosmod, 'B') of 
          0: tmp.exten = 1 
          3: tmp.exten = 2 
          6: tmp.exten = 3
          else: stop
      endcase
      hires = [hires,tmp]
      idx = [pos]
      pos = pos + 1
  endif
  if strpos(mosmod, 'G') NE -1 then begin 
      tmp = { hiresstrct }
      tmp.chip = 2
      case strpos(mosmod, 'G') of 
          0: tmp.exten = 1 
          3: tmp.exten = 2 
          6: tmp.exten = 3
          else: stop
      endcase
      hires = [hires,tmp]
      if keyword_set(IDX) then idx = [idx,pos] else idx = [pos]
      pos = pos + 1
  endif
  if strpos(mosmod, 'R') NE -1 then begin
      tmp = { hiresstrct }
      tmp.chip = 3
      case strpos(mosmod, 'R') of 
          0: tmp.exten = 1 
          3: tmp.exten = 2 
          6: tmp.exten = 3
          else: stop
      endcase
      hires = [hires,tmp]
      if keyword_set(IDX) then idx = [idx,pos] else idx = [pos]
      pos = pos + 1
  endif

  ;; AMP
  amode = strtrim(sxpar(head,keywd.ampmod),2)
  hires[idx].ampmode = amode
  case amode of
      'SINGLE:A': begin
          hires[idx].amp = 1
      end
      'SINGLE:B': begin
          hires[idx].amp = 2
      end
      'DUAL:A+B': stop
      else: stop
  endcase
          
  return, idx

end
    
