;+
; NAME:
;   bok_slitset
;
; PURPOSE:
;   Generate slitmask structure for the Bok B&C spectrograph.
;
; CALLING SEQUENCE:
;   tset_slit = mmt_slitset( nx,ny,xtrim=xtrim)
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
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   17-Dec-2009  J. Moustakas, based on Hennawi's GMOS_SLITSET()
;-
;------------------------------------------------------------------------------

FUNCTION bok_SLITSET, nx, ny, xtrim=xtrim

tset_proto = $
  { func    :    'legendre', $
    xmin    :    0.0, $
    xmax    :    float(ny-1), $
    coeff   :    dblarr(3, 1), $
    dims    :    long([nx, ny]) $
  }
;
tset_slits = replicate(tset_proto, 2)

if (n_elements(xtrim) eq 0) then begin
   xtrim1 = 0.0
   xtrim2 = 0.0
endif else begin
   if (strtrim(xtrim,2) ne '') then begin
      xxtrim = float(strsplit(xtrim,',',/extract))
      xtrim1 = xxtrim[0]
      xtrim2 = xxtrim[1]
   endif else begin
      xtrim1 = 5.0
      xtrim2 = 5.0
   endelse
endelse

; one single longslit
tset_slits[0].coeff[0, 0] = xtrim1
tset_slits[1].coeff[0, 0] = (float(nx-1)-xtrim2)

RETURN, tset_slits
END
