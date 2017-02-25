;+
; NAME:
;   mmt_slitmask
;
; PURPOSE:
;   Generate slitmask structure for GEMINI CCD
;
; CALLING SEQUENCE:
;   niri_slitset,nx,ny
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
;   06-Jun-2006  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

FUNCTION NIRI_SLITSET, nx, ny

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 1), $
    dims    :    long([nx, ny]) $
  }

;oscan_left = 20
;oscan_right = 975
oscan_left = 11
oscan_right = 980
tset_slits = replicate(tset_proto, 2)
tset_slits[0].coeff[0, 0] = oscan_left
tset_slits[1].coeff[0, 0] = oscan_right

RETURN, tset_slits
END
