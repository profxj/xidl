;+
; NAME:
;   bspline_workit
;
; PURPOSE:
;   Internal matrix manipulation for bspline extraction
;
;-
;------------------------------------------------------------------------------
function bspline_workind, xdata, ydata, invvar, action, $
            alpha=alpha, beta=beta, lower=lower, upper=upper

   bw = (size(action))[2] 
   nfull = n_elements(lower)

   ;  The next line is REQUIRED to fill a1

   a2 = action * sqrt(invvar # replicate(1.0d,bw))

   alpha = dblarr(bw,nfull)
   beta = dblarr(nfull)

   bi = lindgen(bw)
   bo = lindgen(bw)
   for i=1L, bw-1 do bi = [bi, lindgen(bw-i)+(bw+1)*i]
   for i=1L, bw-1 do bo = [bo, lindgen(bw-i)+bw*i]

   
   for i=0L, nfull-1 do begin

      itop = i 
      ibottom = (itop < (nfull-bw)) +bw - 1
       
      ict = upper[i] - lower[i] + 1
  
      if (ict GT 0) then begin

         work = a2[lower[i]:upper[i],*] ## transpose(a2[lower[i]:upper[i],*])
         wb   =  (ydata[lower[i]:upper[i]]*sqrt(invvar[lower[i]:upper[i]])) $
                            # a2[lower[i]:upper[i],*] 

         alpha[bo+itop*bw] = alpha[bo+itop*bw] + work[bi]
         beta[itop:ibottom] = beta[itop:ibottom] + wb
      endif
   endfor

   ; Drop break points where minimal influence is located

   min_influence = 1.0e-10 * total(invvar) / nfull

   ; This call to cholesky_band operates on alpha and changes contents

   errb = cholesky_band(alpha, mininf=min_influence) 
   if (errb[0] NE -1) then begin 
      return, errb
   endif
 
   ; this changes beta to contain the solution

   errs = cholesky_solve(alpha, beta)   
   if (errs[0] NE -1) then begin
      return, errs
   endif

;   sc = size(sset.coeff)
;   if (sc[0] EQ 2) then begin
;      sset.icoeff[*,goodbk] = reform(alpha[0,lindgen(nfull)],npoly,nn)
;      sset.coeff[*,goodbk] = reform(beta[lindgen(nfull)], npoly, nn)
;   endif else begin
;      sset.icoeff[goodbk] = alpha[0,lindgen(nfull)]
;      sset.coeff[goodbk] = beta[lindgen(nfull)]
;   endelse

;   if (arg_present(yfit)) then $
;    yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)

   return, -1L 
end
