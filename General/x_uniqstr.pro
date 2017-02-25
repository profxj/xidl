;+ 
; NAME:
; x_uniqstr   
;    Version 1.1
;
; PURPOSE:
;    Obsolote:  Superseded by UNIQ
;
; CALLING SEQUENCE:
;   
; uniq = x_uniqstr(strings, COUNT=count)
;
; INPUTS:
;   strings - Array of strings
;
; RETURNS:
;   uniq  - Array of unique members
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   sort - Sort the output
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   uniq = x_uniqstr( lbls, count=count)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_uniqstr, strngs, COUNT=count, SORT=sort

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'uniq = x_uniqstr, strngs, COUNT= (v1.0)'
      return, -1
  endif 

  nstr = n_elements(strngs)
  tmp = strarr(nstr)
  tmp[0] = strngs[0]
  count = 0L

  for i=1L,nstr-1 do begin
      flg_match = 0
      for j=0L,count do begin
          if strngs[i] EQ tmp[j] then begin
              flg_match = 1
              break
          endif
      endfor
;  Match?
      if flg_match NE 1 then begin
          count = count + 1
          tmp[count] = strngs[i]
      endif
  endfor
  count = count + 1
  if not keyword_set( SORT ) then return, tmp[0:count-1] else begin
      srt_tmp = tmp[0:count-1]
      ans = sort(srt_tmp)
      return, srt_tmp[ans]
  endelse
    

end
  
      
