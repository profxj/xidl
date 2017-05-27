;+
; NAME:
;   trans_tset
;
; PURPOSE:
;   translate a tset to be appropriate for a different range of x
;
; CALLING SEQUENCE:
;   tset = trans_tset(aset, dx)
;
; INPUTS:
;   aset       - input traceset to be translated
;   dx         - number of pixels to translate independent variable
;
; OUTPUTS:
;   tset       - traceset with "outliers" rejected
;          tset = $
;            { func    :    'legendre'  , $
;              xmin    :    0  , $
;              xmax    :    4095   , $
;              coeff   :    array[ncoeff, ntrace] $
;            }
;
; COMMENTS: 
;   This does the "dumbest" thing, but it works. 
;
; BUGS:
;
; EXAMPLES:
;   If for example, one has a "red" solution for 0,4095, and wants
;   that tset shifted to -4096,-1 (for the "blue" solution) , one
;   could write:
;   tset2 = trans_tset(tset, -4096)
; 
;   This transformation would be a bit more trivial if wavelengths
;   were not always stored as log10(lambda).  Or maybe not. 
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;-
;----------------------------------------------------------------------------
function trans_tset, aset, dx
  
  traceset2xy, aset, xtemp, ytemp
  traceset2xy, aset, xtemp+dx, yval

  xy2traceset, xtemp, yval, tset, /silent

  return, tset
end


