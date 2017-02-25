;
;  Xmeas are the measured centroids of the good arclines
;  Xerr are the corresponding error estimates from trace_crude/trace_fweight
;  Xbase is the fiducial line position at the center of the order
;  nxcoeff is the number of free parameters to fit the slope (default: 2)
;
pro fit_singleorder, xmeas, xerr, ybase, xbase, xfit=xfit, slopefit=slopefit, $
       nxcoeff=nxcoeff

    if NOT keyword_set(nxcoeff) then nxcoeff=2

    nrows = (size(xmeas))[1]
    nlines = (size(xmeas))[0] EQ 1 ? 1 : (size(xmeas))[2]
    ; setup the basis for each line center

    nparams = nlines + nxcoeff
    work1 = reform(replicate(1,nrows) # ((identity(nlines))[*]), $
                  nrows*nlines, nlines)
    tmp_basis = (fchebyshev(xbase, nxcoeff))[*]
    workslope = reform(tmp_basis[*] ## replicate(1,nrows), $
                 nrows*nlines, nxcoeff) * (ybase[*] # replicate(1,nxcoeff))

    invvar = 1.0/(xerr^2 + (xerr EQ 0)) * (xerr GT 0 AND xerr LT 90)
;    invvar = invvar[*,good_lines]

    a1 = [[work1],[workslope]]  
    a2 = transpose(a1 * (invvar[*] # replicate(1,nparams)))

    alpha = a2 # a1
    beta  = a2 # (xmeas[*])
    choldc, alpha, p, /double
    res = cholsol(alpha,p,beta)
    xfit = xmeas*0.0
    xfit[*] = a1 # res
    slopefit = res[nlines:*] 
  
    return
end 


