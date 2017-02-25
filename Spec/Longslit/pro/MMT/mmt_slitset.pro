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

FUNCTION MMT_SLITSET, nx, ny, CHANNEL = CHANNEL

; defaults to blue
IF NOT KEYWORD_SET(CHANNEL) THEN CHANNEL = 'BLUE'

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
CASE CHANNEL OF
    'BLUE': BEGIN
        ;; Translate the standard trace
        std_file= getenv('LONGSLIT_DIR') +  $
          '/calib/standards/std_trace/mmt_blue_300_stdtrace.sav'
        restore, std_file
        pos_set = struct_addtags(pos_set, create_struct('DIMS', long([nx, ny])))
        tset_slits = replicate(pos_set, 2)
        tset_slits[0].coeff[0, 0] = 5.0D
        tset_slits[1].coeff[0, 0] = (float(nx-1)-5.0) ;; JXP (handles trim)
    END
    'RED': BEGIN
        tset_slits = replicate(tset_proto, 2)
        tset_slits[0].coeff[0, 0] = 30.0
        tset_slits[1].coeff[0, 0] = 529.0 
    END
    'NEWRED': BEGIN
        print,'mmt_slitset: nx, ny ', nx,' ',ny
        tset_slits = replicate(tset_proto, 2)
        ;; Edited by JXP (6/2010)
;        tset_slits[0].coeff[0, 0] =  (float(nx)/2.0) - 30.0 ;; x position where slit mask starts
;        tset_slits[1].coeff[0, 0] =  (float(nx)/2.0) + 30.0 ;; y position where slit mask ends
        tset_slits[0].coeff[0, 0] =  5.  ;; Trim just the edges
        tset_slits[1].coeff[0, 0] = (float(nx-1)-5.0) ;; JXP (handles trim)
    END
    ELSE: message, 'unrecognized channel'
ENDCASE
;print,'mmt_slitset :',tset_slits[0].coeff[0,0]
;print,'mmt_slitset :',tset_slits[0]
;print,'mmt_slitset :',tset_slits[1].coeff[0,0]
;print,'mmt_slitset :',tset_slits[1]
RETURN, tset_slits
END
