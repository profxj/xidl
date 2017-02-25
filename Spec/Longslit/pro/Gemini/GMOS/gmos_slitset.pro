;+
; NAME:
;   gmos_slitmask
;
; PURPOSE:
;   Generate slitmask structure for GEMINI CCD
;
; CALLING SEQUENCE:
;   gmos_slitmask, nx,ny
;
; INPUTS:
;
; OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This will break if the binning is not 2x4
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

FUNCTION GMOS_SLITSET, nx, ny

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 3), $
    dims    :    long([nx, ny]) $
  }

tset_slits = replicate(tset_proto, 2)

; This hard wires the GEMINI slits. This treats the gaps in the
; gemini Longslit as 3 different slits. 
tset_slits[0].coeff[0, *] = [11.0, 393.0, 774.0]
tset_slits[1].coeff[0, *] = [382.0, 763.0, 1144.0]


RETURN, tset_slits
END
