;+
; NAME:
;   mmt_slitset
;
; PURPOSE:
;   Generate slitmask structure for the BCS CCD
;
; CALLING SEQUENCE:
;   tset_slit = mmt_slitset( nx,ny )
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
;   10-Apr-2009  Added case for updgraded red Channel
;-
;------------------------------------------------------------------------------

FUNCTION CAHA_SLITSET, nx, ny, hdr

caha_hh = caha_hdr(hdr)
path = strcompress(sxpar(caha_hh, 'PATH'), /rem)

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 1), $
    dims    :    long([nx, ny]) $
  }
;
;tset_slits = replicate(tset_proto, 2)

;print, 'mmt_slitset ', channel 
CASE PATH OF
   'BLUE': BEGIN
      tset_slits = replicate(tset_proto, 2)
      tset_slits[0].coeff[0, 0] = 80.0
      tset_slits[1].coeff[0, 0] = 430.0 
   END
   'RED': BEGIN
      tset_slits = replicate(tset_proto, 2)
      tset_slits[0].coeff[0, 0] = 30.0
      tset_slits[1].coeff[0, 0] = 529.0 
   END
   ELSE: message, 'unrecognized channel'
ENDCASE
RETURN, tset_slits
END
