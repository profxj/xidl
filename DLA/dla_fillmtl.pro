;+ 
; NAME:
; dla_fillmtl
;  V1.0
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   dla_fillmtl, dla
;
; INPUTS:
;
; RETURNS:
;   dla      - IDL structure
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   02-Jan-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro dla_fillmtl, dla, OFF=off

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_fillmtl, dla, OFF= (v1.0)' 
    return
  endif 

; Keywords
  if not keyword_set( OFF ) then off = 0.4
  if not keyword_set( SIGOFF ) then sigoff = 0.1

;  LOOP
  for q=0L,n_elements(dla)-1 do begin
      ;;  Take Alpha if first
      if dla[q].flgAlpha EQ 1 or dla[q].flgAlpha GT 3 then begin
          dla[q].mtl.flg = 1
          dla[q].mtl.val = dla[q].Alpha
          dla[q].mtl.sig = dla[q].sigAlpha
          continue
      endif
      ;;  Try Zn second
      if dla[q].flgZn EQ 1 then begin
          dla[q].mtl.flg = 1
          dla[q].mtl.val = dla[q].ZnH
          dla[q].mtl.sig = dla[q].sigZnH
          continue
      endif
      ;;  Try Fe last
      if dla[q].flgFe EQ 1 OR dla[q].flgFe GT 3 then begin
          dla[q].mtl.flg = 2
          dla[q].mtl.val = dla[q].FeH + off
          dla[q].mtl.sig = dla[q].sigFeH + sigoff
          continue
      endif
  endfor

  return
end
