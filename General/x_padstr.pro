;+ 
; NAME:
; x_padstr
;
; PURPOSE:
;    Pads an input string to a input length
;
; CALLING SEQUENCE:
;   
;   padstr = x_padstr(instr, len, /TRIM)
;
; INPUTS:
;   instr -- Input string (or array)
;   len   -- Desired length of string
;  [PAD_CHAR] -- Input charcater for padding
;
; RETURNS:
;   padstr -- Padded string.  Default is the back
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /TRIM --  First trim the string of extraneous blank spaces
;  /REVERSE -- Pad from the front
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   pad = x_padstr('Blah', 15L)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Nov-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function x_padstr, istr, len, pad_char, TRIM=trim, REVERSE=reverse

; Check

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'padstr = x_padstr(istr, len, [PAD_CHAR], /TRIM, /REVERSE) v(1.1)'
    return, -1
  endif 

  if keyword_set(PAD_CHAR) then begin
      if strlen(pad_char) GT 1 then begin
          message, 'x_padstr: You can only use a character.'
          return, -1
      endif
      padc = pad_char
  endif else padc = ' '

  ;; Array
  nstr = n_elements(istr)
  if nstr GT 1 then begin
     outarr = strarr(nstr)
     for ii=0L,nstr-1 do $
        outarr[ii] = x_padstr(istr[ii], len, pad_char, TRIM=trim, REVERSE=reverse)
     return, outarr
  endif

  ;; TRIM
  if keyword_set( TRIM ) then pstr = strtrim(istr,2) else pstr = istr

  ;; Pad
  olen = strlen(pstr)

  ;; 
  if olen GE len then return, pstr else begin
      pads = replicate(padc, len-olen)
      pads = strjoin(pads)
      if keyword_set(REVERSE) then return, pads+pstr else return, pstr+pads
  endelse


end
