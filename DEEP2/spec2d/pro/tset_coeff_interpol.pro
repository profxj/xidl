
; analogous to IDL interpol, but for coeffs of a tset
; this actually does a linear fit, but COULD call interpol?
function tset_coeff_interpol, tset, X, U

  cf   = tset.coeff
  ncf  = (size(cf, /dimens))[0]
  nset = (size(cf, /dimens))[1]

  if nset NE n_elements(X) then message, 'dimension mismatch!'
  
  newcoeff = dblarr(ncf, n_elements(U))
  for i=0, ncf-1 do begin 
     lf = linfit(X, cf[i, *])
     newcoeff[i, *] = lf[0] + lf[1]*U
  endfor 

  return, newcoeff
end
