;+ 
; NAME:
; mike_taghead   
;     Version 1.1
;
; PURPOSE:
;   Examine the mike structure to determine if a given setup has the
;   requisite calibration files.  It then sets a number of tag names
;   accordingly.   Finally, it gives the object frames a running index
;   based on the object name.
;
;   It is recommended to rerun this program if: (1) you have bombed out of
;   IDL without saving the MIKE structure, (2) you have changed the
;   setup values or calibration files.  I doubt it will ever hurt to
;   rerun.
;
;
; CALLING SEQUENCE:
;   
;  mike_setup, mike
;
; INPUTS:
;   mike --  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   OUTFIL -- Name for the ASCII file output (default:
;             'mike_summ.txt')
;
; COMMENTS:
;
; EXAMPLES:
;   mike_taghead, head
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21-Jun-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_taghead, head

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_taghead, head [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( READFIL ) then readfil = getenv('MIKE_DIR')+'/doc/README'
  
  close, 78
  openr, 78, readfil
  cmmnt = ''
  readf, 78, cmmnt
  close, 78

  sxaddpar, head, 'MIKE_TAG', cmmnt, 'Version #'

  return
end
