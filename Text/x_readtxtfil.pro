;+ 
; NAME:
; x_readtxtfil
;    Version 1.1
;
; PURPOSE:
;  Returns the x-value corresponding to the x-pixel on a GUI.
;   Requires pos, xymnx or the appropriate Structure.
;
; CALLING SEQUENCE:
;  xval = xgetx_plt(xtmp, pos, xymnx, size, /STRCT)
;   
; INPUTS:
;  xtmp -- x-pixel value
;  pos  -- Fraction of plot window covered by plot (2 element array)
;  xymnx -- x-y limits of plot window (4 element array: x0,y0,x1,y1)
;
; RETURNS:
;   xval -- 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STRCT -- xtmp contains a structure with the relevant tags
;
; OPTIONAL OUTPUTS:
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
;   04-Apr-2008 Written by JXP
;-
;------------------------------------------------------------------------------

function x_readtxtfil, fil

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = x_readtextfil(fil)'
    return, -1
  endif 

  ;; Get lines
  nlin = FILE_LINES(fil)
  lines = strarr(nlin)

  ;; Open fil
  get_lun, fili
  openr, fili, fil
  
  readf, fili, lines
  return, lines
end
  
