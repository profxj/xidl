;+ 
; NAME:
; x_setllst
;  Version 1.1
;
; PURPOSE:
;    Returns a line list sturcture gven a file and a flag
;
; CALLING SEQUENCE:
;   
;   llist = x_setllst(file, fmt)
;
; INPUTS:
;   file       - Line list filename
;   fmt        - Format of line file  
;                    0='lin.dat', 1=Gal or QSO line list
;
; RETURNS:
;   llist - Structure of linelists
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
;   llist = x_setllst('/u/xavier/bin/junk/lin.dat', 0)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Sat-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setllst, file, fmt

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'llist = x_setllst(file, fmt) [v1.1]'
    return, -1
  endif 

;

  dumc = ' '
  dumc2 = ' '
  dumf1 = 0.0
  dumd1 = 0.0d
  dumd2 = 0.0d
  dumi1 = 0
  nlin = 0L

  ; Open the file
  close, 9
  openr, 9, file, ERROR = err
  if(err NE 0) then begin
      print, 'File does not exist', err
      return, -1
  endif
  
  tmp = { lliststrct }
  case fmt of 
      0: begin   ; QAL line list
          readf, 9, nlin
          tmplst = replicate(tmp, nlin)
          for i=0L,nlin-1 do begin
              readf, 9, format='(f9.4,a12,1x,f10.5,2x,i2)', $
                dumd1, dumc, dumd2, dumi1
              tmplst[i].wave = dumd1
              tmplst[i].name = strtrim(dumc,2)
              tmplst[i].fval = dumd2
              tmplst[i].ref = dumi1
          endfor
          close, 9
      end
      1: begin ; Galaxy line list + QSO line list
          readf, 9, nlin
          tmplst = replicate(tmp, nlin)
          for i=0L,nlin-1 do begin
              readf, 9, format='(f8.2,2x,i1,3x,a2,1x,a5)', $
                dumf1, dumi1, dumc, dumc2
              tmplst[i].wave = dumf1
              tmplst[i].elm = strtrim(dumc,2)
              tmplst[i].ion = strtrim(dumc2,2)
              tmplst[i].flg = dumi1
              tmplst[i].name = tmplst[i].elm+tmplst[i].ion+' '+$
                strmid(strtrim(tmplst[i].wave,2), 0, 6)
          endfor
          close, 9
       end
      2: begin    ;molecular line
         readcol, file, dumf1, dumi1, dumc, dumc2, format='D,I,A,A'
         nlin=n_elements(dumi1)
         tmplst = replicate(tmp,nlin)
         for i=0L,nlin-1 do begin
            tmplst[i].wave = dumf1[i]
            tmplst[i].elm = dumc[i]
            tmplst[i].ion = dumc2[i]
            tmplst[i].flg = dumi1[i]
            tmplst[i].name = tmplst[i].elm+tmplst[i].ion+' '
            ;+$
            ;                 strmid(strtrim(tmplst[i].wave,2), 0, 6)
         endfor
         close, 9
      end
      else: return, -1
  endcase

  return, tmplst
end

