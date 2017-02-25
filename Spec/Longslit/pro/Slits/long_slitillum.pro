;+
; NAME:
;   long_slitillum
;
; PURPOSE:
;  Generate an illumination image for all the slits in a mask
;
; CALLING SEQUENCE:
;  image = long_slitillum( illumflatfile, tset_slits)
;
; INPUTS:
;  illumflatfile  -- Name of illumination flat image
;  tset_slits     -- slit structure
;
; OUTPUTS:
;  image -- Illumination image
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

function long_slitillum, illumflatfile, slitmask, ximg, edgmask
    
    nslit = max(slitmask)
    slitillum = 1.0*(slitmask GT 0)

     for slitid=1L, nslit do begin
       igood = where(slitmask EQ slitid, ngood)
       illum_set = xmrdfits(illumflatfile, slitid, /silent)
       if keyword_set(illum_set) then begin
          eval_illum = bspline_valu(ximg[igood], illum_set)
          if keyword_set(eval_illum) then slitillum[igood] = eval_illum
       endif
     endfor
     unitinds = WHERE(slitmask LE 0.0 OR edgmask, nunit)
     IF nunit GT 0 THEN slitillum[unitinds] = 1.0
     return, slitillum
end

    

