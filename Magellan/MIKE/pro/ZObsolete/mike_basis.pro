; NAME:
; mike_basis   
;     Version 1.1
;
; PURPOSE:
;    Created a 'smoothed' (quasi-2D) version of the traces given the
;    polynomial coefficients for a set of 1D traces.  The code
;    performs a PCA analysis the polynomial coefficients.  
;    By default the code ignores the first order term in the PCA analysis.
;    
; CALLING SEQUENCE:
;   
;
; INPUTS:
;  tset - Trace set produced by xy2traceset
;  xcen - Trace values
;  outmask -  Mask of good trace values
;
; RETURNS:
;   Updated fit based on the PCA analysis
;
; OUTPUTS:
;   
;
; OPTIONAL KEYWORDS:
;  NCOEFF - Number of coefficients to perfom PCA (default: tset.coeff)
;  MSKTRC - Mask out specific traces in the analysis 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  The code is particularly useful for interpolating through orders
;  with large blank regions (e.g. DLA) or for extrapolating in poor
;  S/N orders.
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_pca
;
; REVISION HISTORY:
;   ??-2003 Written by SMB
;-
;------------------------------------------------------------------------------

function mike_basis, tset, xcen, outmask, ncoeff=ncoeff, eigenvec=eigenvec, $
                     outstr=outstr, msktrc=msktrc, skipx0=skipx0, $
                     poly_ncoeff=poly_ncoeff, x0in=x0in

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'fit = mike_basis(tset, xcen, outmask, NCOEFF=, EIGENVEC=,' + $
        'OUTSTR=, MSKTRC=) [v1.1]'
      return, -1
  endif 
   
  ;; Make the model
  traceset2xy, tset, yset, xfit
  if NOT keyword_set(ncoeff) then ncoeff=(size(tset.coeff))[2]
  
  nrow = n_elements(yset[*,0]) 
  ntrace = n_elements(yset[0,*]) 

  msk = replicate(1B, ntrace)
  if n_elements( MSKTRC ) GT 0 then $
     if msktrc[0] NE -1 then msk[msktrc] = 0B
  
     
  if NOT keyword_set(x0in) then x0in = findgen(ntrace)
 
  usetrace = where(total(outmask EQ 0,1) LT nrow/2 $
                   AND msk EQ 1B, ngoodtrace)

      
  
  if ngoodtrace LT ncoeff-1 then begin
      print, "Less than 7 good traces,", ngoodtrace, "  PCA won't really work"
      print, "Using Trace_guess for now, standard star guide would help"
      return, -1L
  endif 
  
  ;; PCA analysis
  mike_pca, tset.coeff[1:ncoeff-1,usetrace], ncoeff-1, eigencoeff, hidden
  
  ;; find the median x value.... 
  
  ynorm = 2.0 * yset[*,0]/(nrow-1) - 1.
  legendre_poly = flegendre(ynorm,ncoeff)
  eigenvec = legendre_poly[*,1:*] # eigencoeff
  
  med_hidden = djs_median(hidden,2)
  med_highorder = med_hidden
  med_highorder[0] = 0
  
  
  high_order_matrix = rebin(transpose(med_highorder), ntrace, ncoeff-1)
  
  ;; Linear first to first two PCA coefficients
  if NOT keyword_set(poly_ncoeff) then $
    poly_ncoeff = (long(3.3 * ngoodtrace / ntrace ) > 1 ) < 3
  print, 'mike_basis: fitting first order with polynomial order ', $
    poly_ncoeff

  mask0 = usetrace*0+1
  y = (hidden[0,*])[*]
  for jj=0,ngoodtrace-5 do begin & $
    coeff0 = func_fit(x0in[usetrace], y, poly_ncoeff, $
             function_name='poly', invvar = mask0, yfit=fit0) & $
    qdone = djs_reject(y, fit0, outmask=mask0, inmask=mask0, $
                 maxdev=1.0, maxrej=1) & $
    if qdone EQ 1 then break & $
  endfor


  mask1 = mask0
  y = (hidden[1,*])[*]
  for jj=0,ngoodtrace-5 do begin & $
    u = where(mask1) & $
    coeff1 = ladfit(x0in[usetrace[u]], y[u]) & $
    fit1 = poly(x0in[usetrace], coeff1) & $
    qdone = djs_reject(y, fit1, outmask=mask1, inmask=mask1, $
                 maxdev=0.3, maxrej=1) & $
    if qdone EQ 1 then break & $
  endfor


  rej0  = where(mask0 EQ 0)
  rej1  = where(mask1 EQ 0)
  rejpt = { $
            rej0: rej0, $
            rej1: rej1 $
          }

  high_order_matrix[*,0] = poly(x0in, coeff0)
  high_order_matrix[*,1] = poly(x0in, coeff1)
  high_fit = high_order_matrix

  ;;
;  x_splot, usetrace, hidden[1,*], ytwo=fit1, /block, psym1=1
;  stop
  
  high_order_fit = eigenvec # transpose(high_order_matrix)
  sub = ( xcen - high_order_fit ) * outmask
;      vec_rec = eigenvec[*,0] # replicate(1,ntrace)
;      vec_rec_t = transpose(vec_rec)
  
  denom = total(outmask, 1)
  numer = total(sub,1)
  
  x0 = fltarr(ntrace)
  mask = x0
  x0fit = x0
  reduced_chi2 = 0.0
  svx0 = numer/(denom + (denom EQ 0))
 
  if NOT keyword_set(skipx0) then begin
    mask = (abs(denom) GT 10 )
;    mask = (abs(denom) GT 10 AND msk EQ 1B)
    for ii = 1,2 do begin 
        good = where(mask)
        x0[good] = numer[good]/denom[good]  ; Average position
        x0res = poly_fit(x0in[good], x0[good], 3)
        x0fit = poly(x0in, x0res)
        chi2 = (x0 - x0fit)^2*mask
        mask = mask * (chi2 LT total(chi2)/2.)
        reduced_chi2 =  total(chi2)/total(mask)
        print, 'mike_basis: reduced_chi2 --', reduced_chi2
    endfor
    if reduced_chi2 GT 2.0 then $
      print, 'MIKE_BASIS: Very large residuals, check traces' $
    else if reduced_chi2 GT 0.5 then $
      print, 'MIKE_BASIS: Large residuals, check traces' 

    bad =where(mask EQ 0)
    if bad[0] NE -1 then x0[bad] = x0fit[bad]

  endif

  x3fit = eigenvec # transpose(high_order_matrix) + x0 ## replicate(1,nrow)

  ;; Output structure
  if arg_present(OUTSTR) then begin
      outstr = { $
                 hidden: hidden, $
                 high_fit: high_fit, $
                 usetrc: usetrace, $
                 rejpt: rejpt, $
                 red_chi2: reduced_chi2, $
                 x0: svx0, $
                 x0fit: x0fit, $
                 x0msk: mask $
               }
  endif
  
  return, x3fit
end
