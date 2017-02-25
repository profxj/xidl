;+ 
; NAME:
; vpparse
;   Version 1.1
;
; PURPOSE:
;   parses a vpfit fort.26 output file and puts data into a structure of type
;   'vpstrct' which holds the relavent data.
; 
;
; CALLING SEQUENCE:
;   vpparse, vpfil, VPSTR=vpstr
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Written by GEP
;   17 Sep 2010  modified by KLC
;-
;------------------------------------------------------------------------------
function find_last_num, str0
  ;; Analyze input string (str0) to see what of the characters are
  ;; numbers. Assumes the letters come at the end.
  ;; Also tests for ill-defined numbers (which print as *****)

  len = strlen(str0)
  l0 = strmid(str0,0,1)         ; first character
  if l0 eq '*' then begin
     ;; The rest of them will be *
     return,-1
  endif 

  done = 0
  
  while not done do begin
     lf = strmid(str0,len-1,1)  ; last character
     if stregex(lf,'[A-Za-z]',/boolean) then begin
        ;; Still a string
        len = len - 1
        if len lt 0 then done = 1
     endif else begin
        ;; String done, return present location as index of first
        ;; letter character 
        done = 1
     endelse 
  endwhile

  ;; Then the following command should yield a number:
  ;; num = double(strmid(str0,0,len)) if len > 0 
  return,len
  
end                             ; find_last_num()


function pop_char, str0, flt=flt, dbl=dbl, first=first, both=both, $
                   multiple=multiple
  ;; Remove last (default), first or both last and first characters
  ;; of a string and return it as a string (default), float (/flt), 
  ;; or double (/dbl).
  ;; Also enable with '/multiple' to find several last characters
  catch, error_flag
  str = strtrim(str0,2)
  if keyword_set(flt) or keyword_set(dbl) then num = 1

  if keyword_set(multiple) and keyword_set(num) then begin
     ;; Only work on numbers
     indx = find_last_num(str0)
     if indx eq -1 then val = -1.d $
     else val = double(strmid(str0,0,indx))
     if keyword_set(flt) then return,float(val) $
     else return,val
  endif else begin
     ;; Do the normal thing
     nwstr2 = strmid(str,0,strlen(str)-1)
     if keyword_set(first) or keyword_set(both) then begin
        if keyword_set(first) then nwstr = strmid(str,1) $
        else nwstr = strmid(nwstr2,1)
     endif else nwstr = nwstr2
     if keyword_set(flt) then val = float(nwstr)
     if keyword_set(dbl) then val = double(nwstr)
     if error_flag ne 0 and keyword_set(num) then $
        stop, 'pop_char(): cannnot convert ',str,' to number'  
     catch, /cancel
     if keyword_set(val) then return, val else return, nwstr
  endelse 
end                             ; pop_char()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro g_vpparse, vpfil, VPSTR=vpstr

    vpstr = {vpstrct}
    nlin = file_lines(vpfil)
    openr, 1, vpfil

    dumc=''
    nreg = 0
    nion = 0

    for i=0,nlin-1 do begin
        readf, 1, dumc
        prs = strsplit(dumc,' ', /extract)
        if prs[0] eq '!' then continue ; comment line
        if prs[0] EQ '%%' then begin   ;  a Region line
            vpstr.fluxfil[nreg] = prs[1]
            vpstr.reg_beg[nreg] = double(prs[3])
            vpstr.reg_end[nreg] = double(prs[4])
            nreg = nreg + 1
        endif else begin
            if strlen(prs[0]) EQ 1 then begin  ; the ion is like 'H I'
                vpstr.ion[nion] = prs[0]+' '+prs[1]
                vpstr.z[nion] = pop_char(prs[2],/dbl,/multiple)
                vpstr.zerr[nion] = pop_char(prs[3],/dbl,/multiple)
                vpstr.z_str[nion] = prs[2]
                vpstr.b[nion] = pop_char(prs[4],/dbl,/multiple)
                vpstr.berr[nion] = pop_char(prs[5],/dbl,/multiple)
                vpstr.b_str[nion] = prs[4]
                vpstr.n[nion] = pop_char(prs[6],/dbl,/multiple)
                vpstr.nerr[nion] = pop_char(prs[7],/dbl,/multiple)
                vpstr.n_str[nion] = prs[6]
            endif else begin                       ; the ion is like 'SiIV'
                vpstr.ion[nion] = prs[0]
                vpstr.z[nion] = pop_char(prs[1],/dbl,/multiple)
                vpstr.zerr[nion] = pop_char(prs[2],/dbl,/multiple)
                vpstr.z_str[nion] = prs[1]
                vpstr.b[nion] = pop_char(prs[3],/dbl,/multiple)
                vpstr.berr[nion] = pop_char(prs[4],/dbl,/multiple)
                vpstr.b_str[nion] = prs[3]
                vpstr.n[nion] = pop_char(prs[5],/dbl,/multiple)
                vpstr.nerr[nion] = pop_char(prs[6],/dbl,/multiple)
                vpstr.n_str[nion] = prs[5]
            endelse
            nion = nion + 1
        endelse                            ;  a transition line
    endfor
    close, 1

    vpstr.nreg = nreg
    vpstr.nion = nion

end
