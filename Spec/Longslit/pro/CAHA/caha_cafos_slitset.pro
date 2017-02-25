;+
; NAME:
;   mmt_slitmask
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

FUNCTION caha_cafos_slitset, nx, ny

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 1), $
    dims    :    long([nx, ny]) $
  }

tset_slits = replicate(tset_proto, 2)

tset_slits[0].coeff[0, 0] = 0.0
tset_slits[1].coeff[0, 0] = double(nx-1L)

RETURN, tset_slits
END
