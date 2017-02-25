;+ 
; NAME:
; getabnd
;
; PURPOSE:
;    Given the name or Z, get abnd and or Z.  Reads the file
;    abnd_dat (default: abnd.dat in $XIDL_DIR/Abund (Meteoritic where
;    possible)) 
;
; CALLING SEQUENCE:
;   
;   getabnd, name, Z, abnd, flag=, abnd_dat=
;
; INPUTS:
;   name (Z)    - Element name
;
; RETURNS:
;   Z       - Atomic number
;   abnd    - Meteoritic Abundance
;   [name]  - Element name [returned if flag = 1]
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   flag    - 0=Name input, 1=Z input
;   abnd_dat - Element  Abundance  Atomic Number
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES: ;   getabnd, 'C', Z, abnd
;   getabnd, elm, 14, abnd, flag=1
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro getabnd, nam, Zval, abnd, FLAG=flag, abnd_dat=abnd_dat

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'getabnd, nam, Zval, abnd, [flag=, abnd_dat=] [v1.1]'
    return
  endif 

;
  if not keyword_set( FLAG ) then    flag    = 0 

;  Parse down the element name allowing for whitespace

  if(flag EQ 0) then begin
      c1 = 0
      c2 = 0
      
      for i=0,strlen(nam)-1 do begin
          c1 = i
          if(strmid(nam,c1,1) NE ' ') then break
      endfor
      
      if(strlen(nam) GT c1+1) then $
        if(strmid(nam,c1+1,1) EQ 'I' OR strmid(nam,c1+1,1) EQ 'V' OR $
           strmid(nam,c1+1,1) EQ ' ') then c2 = c1 else c2 = c1+1
      
      nnam = ' '
      if(c2 EQ c1) then nnam = strmid(nam,c1,1) + ' ' $
      else nnam = strmid(nam,c1,c2-c1+1)
  endif
  
  
; Initialize

  dumc = ' '
  dumr = 0.0
  dumi = 0

  if keyword_set(abnd_dat) then fil = abnd_dat $
  else fil = getenv('XIDL_DIR')+'/Abund/abnd.dat'
  close, 9
  openr, 9, fil, ERROR = err
  if(err NE 0) then stop, 'getabnd: File does not exist', err
  readf, 9, dumc  ;; Reference
;
  while not eof(9) do begin
      readf, 9, format='(a2, 1x, f5.2, 1x, i2)', dumc, dumr, dumi
      case flag of 
          0: begin
              if(strcmp(nnam,dumc,2) EQ 1) then begin
                  Zval = dumi
                  abnd = dumr
                  goto, finish
              endif
          end
          1: begin
              if(Zval EQ dumi) then begin
                  nam = dumc
                  abnd = dumr
                  goto, finish
              endif
          end
          else: stop
      endcase
   endwhile

  
  if not keyword_set(NAM) then nam = ' '
  print, 'getabnd: Warning No match', nam
  print, 'getabnd: Assuming a value of 0.00'
  abnd = 0.00

  finish: close, 9
  nam = strtrim(nam,2)

  return
end

