;+
; NAME:
;   vprint
;
; PURPOSE:
;   decides whether to print based on verbosity level
;
; CALLING SEQUENCE:
;   vprint, verblevel, string
; 
; INPUTS:
;   verblevel  - level of verbosity
;   string     - string to print
;
; KEYWORDS:
;   same as for print
;
; COMMON BLOCKS:
;   verbosity_level
;
; EXAMPLES:
;   vprint, 1, 'vprint is easy to use!'
;
; COMMENTS:
;   Use vprint for output from scripted programs. 
;   high vlevel means one must set high verbosity (with verbset)
;   in order for the statement to print. 
;
;   We indent by a number of spaces equal to the verbosity. 
;
;   This is analogous to R. Lupton's verbecho in dervish
;
; REVISION HISTORY:
;
;       Sat Feb 23 17:23:22 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
pro vprint, vlevel, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, _extra=extra

  common verbosity_level, verblevel

  if n_params() GT 11 then message, 'too many params!'

  if size(vlevel, /tname) EQ 'STRING' then begin 
     message, 'vprint, vlevel, string', /info
     string = vlevel
     vlevel = 0
  endif 

  if n_elements(verblevel) NE 0 then begin 
     loquacious = verblevel GE vlevel
  endif else begin 
     message, 'Please set verbosity level with verbset, verblevel', /info
     loquacious = 1B
  endelse 

  spc = strmid('                        ', 0, vlevel)

  if loquacious then begin 
     if n_params() eq 2 then print, spc, s0, _extra=extra
     if n_params() eq 3 then print, spc, s0, s1, _extra=extra
     if n_params() ge 4 then begin 
        exstr = "print"
        for i=0, n_params()-2 do begin
           exstr = exstr+', s'+string(i, format='(I1)')
        endfor 

        foo = execute(exstr)
     endif 
  endif 
end

