;+ 
; NAME:
; dla_indx
;  V1.0
;
; PURPOSE:
;    Returns the index for a DLA given the name and z (if necessary)
;     Defaults to the first entry if multiple
;
; CALLING SEQUENCE:
;   
;   idx = dla_indx(struct, name, [z]) 
;
; INPUTS:
;   struct - dla structure
;   name - string
;   z - redshift (optional)
;
; RETURNS:
;   index
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   IDX - the same index
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
function dla_indx, sdla, name, z, IDX=idx

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'idx = dla_indx(sdla, name, [z], IDX=) [v1.0]'
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

  if arg_present( IDX ) then idx = all[0]

  return, all[0]

end
