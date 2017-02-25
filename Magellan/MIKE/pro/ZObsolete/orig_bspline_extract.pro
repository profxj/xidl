;+
; NAME:
;   bspline_extract
;
; PURPOSE:
;   Calculate a B-spline in the least squares sense with rejection, 
;     using a model profile.
;
; CALLING SEQUENCE:
;   sset = bspline_extract(xdata, ydata, invvar, slit_profile, obj_profile, $
;     yfit=yfit, numiter=iiter, bkpt=bkpt, maxiter=maxiter,relative=relative, $
;     upper=upper, lower=lower, outmask=outmask, fullbkpt=fullbkpt, $
;     inmask=inmask, _EXTRA=EXTRA, RED_CHI=reduced_chi)
;
; INPUTS:
;   xdata      - Data x values
;   ydata      - Data y values
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance of y; if not set, then set to be
;                consistent with the standard deviation.  This only matters
;                if rejection is being done.
;   slit_profile - spatial profile of the background at each pixel
;   obj_profile  - spatial profile of the object at each pixel
;   nord       - Order for spline fit; default to 4.
;   x2         - 2nd dependent variable for 2-D spline fitting.
;   npoly      - Polynomial order to fit over 2nd variable (X2); default to 2.
;   xmin       - Normalization minimum for X2; default to MIN(XDATA).
;   xmax       - Normalization maximum for X2; default to MAX(XDATA).
;   oldset     - If set, then use values of FULLBKPT, NORD, XMIN, XMAX, NPOLY
;                from this structure.
;   funcname   - If OLDSET is not specified and this is a 2-D B-spline,
;                then the function for the second variable may be passed.
;                The default is 'legendre' in the call to CREATE_BSPLINESET().
;   maxiter    - Maximum number of rejection iterations; default to 10;
;                set to 0 to disable rejection.
;   upper      - Upper rejection threshhold; default to 5 sigma.
;   lower      - Lower rejection threshhold; default to 5 sigma.
;   _EXTRA     - Keywords for BSPLINE_BKPTS() and/or DJS_REJECT().
;
; OUTPUTS:
;   sset       - Structure describing spline fit.
;                Return 0 if too few good (INVVAR NE 0) points are passed.
;
; OPTIONAL OUTPUTS:
;   outmask    - Output mask, set =1 for good points, =0 for bad points.
;   fullbkpt   - If OLDSET is not specified, then the break points are
;                chosen with a call to BSPLINE_BKPTS() which can be returned
;                with this keyword.
;   yfit       - B-spline fit evaluated at each data point.
;   numiter    - Last iteration (0-indexed)
;
; COMMENTS:
;
;   Wavelengths must be sorted!
;
;   Data points can be masked either by setting their weights to zero
;   (INVVAR[]=0), or by using INMASK and setting bad elements to zero.
;   INMASK is passed to DJS_REJECT().
;
;   If OLDSET is used, then the output structure SSET will be a structure
;   with the same name as OLDSET.  This will allow the two structures to
;   be concatented, i.e.
;     > junk = [oldset, sset]
;
;   Although I'm not sure how to treat data points which fall outside
;   minmax(bkpt), now I will set them equal to minmax with invvar = 0
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   bspline_bkpts()
;   bspline_fit()
;   create_bsplineset()
;   djs_reject()
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------
function bspline_workit, xdata, ydata, invvar, action, sset, $
            lower=lower, upper=upper, npoly=npoly, nord=nord, yfit=yfit, $
            covariance=covariance, alpha=alpha

   nord = sset.nord
   npoly = sset.npoly
   goodbk = where(sset.bkmask[nord:*] NE 0, nbkpt)

   if (nbkpt LT nord) then begin
      if (arg_present(yfit)) then yfit = fltarr(n_elements(ydata))
      return, -2L
   endif

   nn = nbkpt 
   nfull = nn * npoly
   bw = npoly * nord   ; this is the bandwidth

   ;  The next line is REQUIRED to fill a1

   a2 = action * sqrt(invvar # replicate(1,bw))

   alpha = dblarr(bw,nfull+bw)
   beta = dblarr(nfull+bw)

   bi = lindgen(bw)
   bo = lindgen(bw)
   for i=1L, bw-1 do bi = [bi, lindgen(bw-i)+(bw+1)*i]
   for i=1L, bw-1 do bo = [bo, lindgen(bw-i)+bw*i]

   for i=0L, nn-nord do begin

      itop = i * npoly
      ibottom = (itop < nfull + bw) - 1
       
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

   covariance = alpha
   errb = cholesky_band(alpha, mininf=min_influence) 
   if (errb[0] NE -1) then begin 
      if (arg_present(yfit)) then $
        yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)
      return, bspline_maskpoints(sset, errb, npoly)
   endif
 
   ; this changes beta to contain the solution

   errs = cholesky_solve(alpha, beta)   
   if (errs[0] NE -1) then begin
      if (arg_present(yfit)) then $
       yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)
      return, bspline_maskpoints(sset, errs, npoly)
   endif

   sc = size(sset.coeff)
   if (sc[0] EQ 2) then begin
      sset.icoeff[*,goodbk] = reform(alpha[0,lindgen(nfull)],npoly,nn)
      sset.coeff[*,goodbk] = reform(beta[lindgen(nfull)], npoly, nn)
   endif else begin
      sset.icoeff[goodbk] = alpha[0,lindgen(nfull)]
      sset.coeff[goodbk] = beta[lindgen(nfull)]
   endelse

   if (arg_present(yfit)) then $
    yfit = bspline_valu(xdata, sset, x2=xdata, action=action, upper=upper, lower=lower)

   return, 0L
end
;------------------------------------------------------------------------------
function bspline_extract, xdata, ydata, invvar, slit_profile, obj_profile, $
     yfit=yfit, numiter=iiter, bkpt=bkpt, maxiter=maxiter, relative=relative, $
     upper=upper, lower=lower, outmask=outmask, fullbkpt=fullbkpt, $
     inmask=inmask, _EXTRA=EXTRA, RED_CHI=reduced_chi, buff=buff, $
     covariance=covariance, alpha=alpha

   if (n_params() LT 2) then begin
      print, 'Syntax -  sset = bspline_extract( )'
      return, 0
   endif

   ;----------
   ; Check dimensions of inputs

   nx = n_elements(xdata)
   if (n_elements(ydata) NE nx) then $
    message, 'Dimensions of XDATA and YDATA do not agree'

   if (NOT keyword_set(nord)) then nord = 4L
   if n_elements(upper) EQ 0 then upper = 5
   if n_elements(lower) EQ 0 then lower = 5

   if (keyword_set(invvar)) then begin
      if (n_elements(invvar) NE nx) then $
       message, 'Dimensions of XDATA and INVVAR do not agree'
   endif 

   if (n_elements(maxiter) EQ 0) then maxiter = 25L

   yfit = 0 * ydata ; Default return values

   if (NOT keyword_set(invvar)) then begin
      var = variance(ydata, /double)
      if (var EQ 0) then var = 1
      invvar = 0.0 * ydata + 1.0/var
   endif

   maskwork = invvar[*] GT 0
   if keyword_set(inmask) then maskwork = maskwork*inmask
   these = where(maskwork, nthese)
 
   ;----------
   ; Determine the break points and create output structure

   if NOT keyword_set(fullbkpt) then $
      fullbkpt = bspline_bkpts(xdata[these], nord=nord, bkpt=bkpt, $
        _EXTRA=EXTRA)

   
   npoly = 2L
   sset = create_bsplineset(fullbkpt, nord, npoly=npoly) 
   sset.funcname = 'Sky subtraction special'

   if (nthese LT nord) then begin
         print, 'Number of good data points fewer the NORD'
         return, sset
   endif


   ;----------
   ; Iterate spline fit

   iiter = 0
   error = 0
   tempin = 0
   reduced_chi = 0
   if NOT keyword_set(buff) then buff = ''

   while (((error[0] NE 0) OR (keyword_set(qdone) EQ 0)) $
    AND iiter LE maxiter) do begin

      ngood = total(maskwork)
      goodbk = where(sset.bkmask NE 0)

      if (ngood LE 1 OR goodbk[0] EQ -1) then begin
         sset.coeff = 0
         iiter = maxiter + 1; End iterations
      endif else begin
        ; Do the fit.  Return values for ERROR are as follows:
        ;    0 if fit is good
        ;   -1 if all break points are masked
        ;   -2 if everything is screwed

        ;  we'll do the fit right here..............


        bf1 = bspline_action(1.0d*xdata, sset, lower=loweraction, upper=upperaction)
        filler = replicate(1,nord)
        action = dblarr(nx, npoly*nord)
        action[*,lindgen(nord)*npoly] = bf1 * (obj_profile[*] # filler)
        action[*,lindgen(nord)*npoly+1] = bf1 * (slit_profile[*] # filler)

        ; help, where(finite(bf1) EQ 0), where(finite(action) EQ 0)
 
        if total(finite(action) EQ 0) GT 0 then begin
             print, '!! Infinities in action matrix, wavelengths may be very messed up!!!'
            return, 0
        endif

        error = bspline_workit(xdata, ydata, invvar*maskwork, action, sset, alpha=alpha, $
            lower=loweraction, upper=upperaction, yfit=yfit, covariance=covariance)

      endelse

      iiter = iiter + 1

      if (error[0] EQ -2L) then begin
         ; All break points have been dropped.
         return, sset
      endif else if (error[0] EQ 0) then begin
         ; Iterate the fit -- next rejection iteration.
          reduced_chi = total((ydata - yfit)^2  * (invvar*maskwork)) / $
                              (ngood - 2*(n_elements(goodbk)+nord) - 1)


          relative_factor = 1.0
          if keyword_set(relative) then $
                relative_factor = sqrt(reduced_chi) > 1.0

         qdone = djs_reject(ydata, yfit, invvar=invvar, inmask=tempin, $
                            outmask=maskwork, upper=upper*relative_factor, $
                            lower=lower*relative_factor, _EXTRA=EXTRA)

          tempin = maskwork
          print, format='(i4, f7.3,i5)', iiter, reduced_chi, $
            total(maskwork EQ 0)
          print, format='(a, $)', buff

;            total(maskwork EQ 0), replicate(string(8b),16)
      endif

  endwhile

;   print, format='($, i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   print, format='(i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   outmask = maskwork
   ;----------
   ; Re-sort the output arrays OUTMASK and YFIT to agree with the input data.

   return, sset
end
;------------------------------------------------------------------------------
