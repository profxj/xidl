;+
; NAME:
;   gmos_slitset
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

FUNCTION XGMOS_SLITSET, nx, ny, XBIN, CHANNEL = CHANNEL

; defaults to blue
IF NOT KEYWORD_SET(XBIN) THEN XBIN=1

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 1), $
    dims    :    long([nx, ny]) $
  }
;
;tset_slits = replicate(tset_proto, 2)
; 820, 1490

print,'xgmos_slitset: nx, ny ', nx,' ',ny
tset_slits = replicate(tset_proto, 2)

xstart = round(2*820./xbin)  
xend = round(2*1490./xbin) 
IF xstart GT nx THEN BEGIN
   ;; This is a rastered standard star image, just using 0,nx as slit
   ;; boundary
   ;; x position where slit mask starts
   tset_slits[0].coeff[0, 0] =  170.0
;; x position where slit mask ends
   tset_slits[1].coeff[0, 0] =  500.0
ENDIF ELSE BEGIN
;; x position where slit mask starts
   tset_slits[0].coeff[0, 0] =  xstart
;; x position where slit mask ends
   tset_slits[1].coeff[0, 0] =  xend
ENDELSE   
RETURN, tset_slits
END
