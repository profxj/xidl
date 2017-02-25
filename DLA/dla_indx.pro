;+ 
; NAME:
; dla_indx
;  V1.1
;
; PURPOSE:
;    Returns the index for a DLA given the name and z (if necessary)
;     Defaults to the first entry if multiple.  This is a simple and
;     rather uninteresting program.
;
; CALLING SEQUENCE:
;   
;   idx = dla_indx(struct, name, [z]) 
;
; INPUTS:
;   struct - dla structure
;   name - Name of the quasar (string; need not be complete)
;   [z] - redshift (optional)
;
; RETURNS:
;   index
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
;   idx = dla_indx(sdla, 'Q0000')
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   24-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------
function dla_indx, sdla, name, z

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'idx = dla_indx(sdla, name, [z]) [v1.1]'
    return, -1
  endif 
;


  if not keyword_set( z ) then $
    all = where(strmid(sdla.qso, 0, strlen(name)) EQ name, nall) $
  else $
    all = where(strmid(sdla.qso, 0, strlen(name)) EQ name AND $
                abs(sdla.zabs-z) LT 0.1, nall) 

  if nall EQ 0 then begin
      print, 'DLA '+name+' not found'
      return, -1
  endif

  return, all[0]

end
