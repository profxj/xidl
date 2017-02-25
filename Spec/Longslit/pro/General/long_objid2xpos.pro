;+
; NAME:
;   long_objid2xpos
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;
; BUGS:
;   
; REVISION HISTORY:
;   11-Jan-2011  Written by JXP
;-  
;-----------------------------------------------------------------------------
function long_objid2xpos, scifil, objid, USAGE=usage, REVERSE=reverse

  if  KEYWORD_SET(USAGE) THEN BEGIN
      print,'Syntax - ' + $
        'long_reduce, planfile, /REV [v1.0]' 
      return, -1
  endif 

  strct = xmrdfits(scifil, 5)

  if keyword_set(REVERSE) then begin
     ny = n_elements(strct[0].ypos)
     mn = min(abs(strct.xpos[ny/2] - objid), imn)
     return, strct[imn].objid
  endif else begin
     idx = where(strct.objid EQ objid, nidx)
     if nidx NE 1 then begin
        print, 'long_objid2xpos:  No match!!', objid
        return, -1.
     endif
     ostr = strct[idx]
     ;; Half way up
     ny = n_elements(ostr.ypos)
     ;; Return
     return, ostr.xpos[ny/2]
  endelse
end
;------------------------------------------------------------------------------
