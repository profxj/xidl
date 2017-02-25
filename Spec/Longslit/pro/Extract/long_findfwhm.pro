;+
; NAME:
;   long_findfwhm
;
; PURPOSE:
;   Calculate the spatial FWHM of all the objects in the slit.  The
;   routine assumes a simple algorithm.
;
; CALLING SEQUENCE:
;  long_findfwhm, model, x, peak, peak_x, lwhm, rwhm
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   traceset2xy
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;------------------------------------------------------------------------------
pro long_findfwhm, model, x, peak, peak_x, lwhm, rwhm

    peak = max(model * (abs(x) LT 1.), pl)
    peak_x =  x[pl]

    lh = min(where(reverse((x LT peak_x) AND (model LT 0.5*peak))))
    rh = min(where((x GT peak_x) AND (model LT 0.5*peak)))
    if lh[0] NE -1 then lwhm = (reverse(x))[lh] else lwhm = -0.5*2.3548
    if rh[0] NE -1 then rwhm = x[rh] else rwhm =  0.5*2.3548

    return
end

