;+ 
; NAME:
; x_chkfil   
;    Version 1.1
;
; PURPOSE:
;    Check for a file and return the number of entries matching
;
; CALLING SEQUENCE:
;   
; flag = x_chkfil(fil, COUNT=, /SILENT)
;
; INPUTS:
;   fil - Filename 
;
; RETURNS:
;   flag  - 0: No file; 1: One file; 2: More than one
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SILENT -- Suppress messages to screen
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   flg = x_chkfil('junk.dat', count=cnt, /SILENT)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;   08-Jun-2011 Use IDL's file_search() instead of findfile()
;               for speed, KLC
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_chkfil, fil, COUNT=count, SILENT=silent

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'flg = x_chkfil(fil, COUNT=, /SILENT) (v1.1)'
      return, -1
  endif 

  a = file_search(fil, COUNT=count)
  case count of 
      0: begin
          if not keyword_set( SILENT ) then $
            print, 'x_chkfil: '+strtrim(fil,2)+' not found!'
          return, 0
      end
      1: return, 1
      else: return, 2
  endcase

end
