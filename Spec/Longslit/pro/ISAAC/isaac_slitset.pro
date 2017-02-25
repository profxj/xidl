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

FUNCTION ISAAC_SLITSET, nx, ny, NSLIT = NSLIT1

;; default to use a single side of the detector for now. 
  IF KEYWORD_SET(NSLIT1) THEN NSLIT = NSLIT1 ELSE NSLIT = 1L
  
  CASE NSLIT OF
     1: BEGIN 
        tset_proto = $
           { func    :    'legendre', $
             xmin    :    0.0, $
             xmax    :    float(ny-1), $
             coeff   :    dblarr(3, 1), $
             dims    :    long([nx, ny]) $
           }
        oscan_lef1 = 95
        oscan_rig1 = 511
        tset_slits = replicate(tset_proto, 2)
        tset_slits = replicate(tset_proto, 2)
        tset_slits[0].coeff[0, 0] = oscan_lef1
        tset_slits[1].coeff[0, 0] = oscan_rig1
     END
     2: BEGIN  
        tset_proto = $
           { func    :    'legendre', $
             xmin    :    0.0, $
             xmax    :    float(ny-1), $
             coeff   :    dblarr(3, 2), $
             dims    :    long([nx, ny]) $
           }
        
        oscan_lef1 = 95
        oscan_rig1 = 511
        oscan_lef2 = 512
        oscan_rig2 = 913
        tset_slits = replicate(tset_proto, 2)
        tset_slits[0].coeff[0, 0] = oscan_lef1
        tset_slits[0].coeff[0, 1] = oscan_lef2
        tset_slits[1].coeff[0, 0] = oscan_rig1
        tset_slits[1].coeff[0, 1] = oscan_rig2
     END
     ELSE: message, 'unsupported NSLIT'
  ENDCASE
RETURN, tset_slits
END
