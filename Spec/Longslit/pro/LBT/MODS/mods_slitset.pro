;+
; NAME:
;   gmos_slitmask
;
; PURPOSE:
;   Generate slitmask structure for MODS CCD
;
; CALLING SEQUENCE:
;   mods_slitmask, nx,ny
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
;   This will break if the binning is not 1x1
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

FUNCTION MODS_SLITSET, nx, ny, CHANNEL= CHANNEL

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 5), $
    dims    :    long([nx, ny]) $
  }

tset_slits = replicate(tset_proto, 2)

; This hard wires the MODS longslit. This treats the gaps in the
; MODS Longslit as 5 different slits. 

if CHANNEL eq 'MODS1B' then begin
   tset_slits[0].coeff[0, *] = [635.0, 1167.0, 1691.0, 2217.0, 2749.0]+30.0
   tset_slits[1].coeff[0, *] = [1129.0, 1657.0, 2181.0, 2701.0, 3071.0]-30.0
endif else if CHANNEL eq 'MODS1R' then begin
   tset_slits[0].coeff[0, *] = [285.0, 798.0, 1310.0, 1824.0, 2337.0]
   tset_slits[1].coeff[0, *] = [774.0, 1283.0, 1799.0, 2309.0, 2826.0]
endif 

RETURN, tset_slits
END
