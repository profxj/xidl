;+
; NAME:
;   clean_tset
;
; PURPOSE:
;   clean traceset by removing traces that are "outliers" in shape. 
;
; CALLING SEQUENCE:
;   clean_tset, tset
;
; INPUTS:
;   tset       - traceset to be cleaned
;
; OUTPUTS:
;   tset       - traceset with "outliers" rejected
;          tset = $
;            { func    :    'legendre'  , $
;              xmin    :    0  , $
;              xmax    :    4095   , $
;              coeff   :    array[ncoeff, ntrace] $
;            }
;
; COMMENTS: 
;   We examine traceset polynomial coeffs of order 1..N for Nth order
;     polynomial (N = ncoeff-1) and fit them with a quadratic function
;     of x.  Outliers greater than 10 sigma are rejected. 
;
;   The (saturated) alignment star boxes are sometimes detected by
;     this, because of the sharp transition from saturated to
;     unsaturated regions as a function of DEIMOS-Y.  This appears to
;     do no harm downstream. 
;
; BUGS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;   14-Oct-2002  Reject based on the tset coeffs themselves - DPF
;-
;----------------------------------------------------------------------------
pro clean_tset, tset

  ntrace = (size(tset.coeff, /dimens))[1]
  ncoeff = (size(tset.coeff, /dimens))[0]

; -------- Loop over coefficients and remove anything that doesn't fit
  good = 1B
  x = tset.coeff[0, *]
  ndeg = 2    ; cubic fits better but may be less stable.  Use 10 sigma below:
  nsig = 3
  for i=1, ncoeff-1 do begin 
     poly_iter, x, tset.coeff[i, *], ndeg, nsig, yfit, coeff=coeff
     thresh = 10*djsig(tset.coeff[i, *]-yfit)
     good = good AND (abs(tset.coeff[i, *]-yfit) LT thresh)
  endfor 
  
  coeff = tset.coeff[*, where(good, ngood)]
  print, 'Rejected by CLEAN: ', ntrace-ngood
; make new tset
  cleantset = $
    { func    :    tset.func   , $
      xmin    :    tset.xmin   , $
      xmax    :    tset.xmax   , $
      coeff   :    coeff $
    }

; overwrite old tset  
  tset = cleantset

  return
end
