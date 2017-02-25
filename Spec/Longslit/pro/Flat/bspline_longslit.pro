;+
; NAME:
;   bspline_longslit
;
; PURPOSE:
;   Calculate a B-spline in the least squares sense with rejection, 
;     using a model profile.
;
; CALLING SEQUENCE:
;   sset = bspline_longslit(xdata, ydata, invvar, profile_basis, $
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
;   profile_basis - spatial profiles of each independent b-spline fit
;                  Should have dimensions [NPIX, NBSPLINES]
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
;   bspline_workit()
;   create_bsplineset()
;   djs_reject()
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------
function bspline_longslit, xdata, ydata, invvar, profile_basis, $
     yfit=yfit, numiter=iiter, bkpt=bkpt, maxiter=maxiter, relative=relative, $
     upper=upper, lower=lower, outmask=outmask, fullbkpt=fullbkpt, $
     inmask=inmask, _EXTRA=EXTRA, RED_CHI=reduced_chi, buff=buff, $
     covariance=covariance, alpha=alpha, nord=nord, silent=silent, $
     uaction = uaction, laction=laction, action=action, DEBUG=DEBUG

;  stop                          ; KHRR
   if (n_params() LT 2) then begin
      print, 'Syntax -  sset = bspline_longslit( )'
      return, 0
   endif

   t0 = systime(1)
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

   npoly = n_elements(profile_basis)/nx
   if npoly*nx NE n_elements(profile_basis) then begin
      splog, 'Profile basis size is not a multiple of data points'
      splog, nx, n_elements(profile_basis), npoly   
      return, 0
   endif

   if NOT keyword_set(silent) then $
      splog, npoly, ' basis functions and', nx, ' pixels'

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
     fullbkpt = bspline_bkpts(xdata[these], nord=nord, bkpt=bkpt, /silent, $
          _EXTRA=EXTRA) 

   
   sset = create_bsplineset(fullbkpt, nord, npoly=npoly) 
   sset.funcname = 'Bspline longslit special'

   action_multiple = reform(profile_basis[*] # replicate(1,nord),nx,npoly*nord) 


   if (nthese LT nord) then begin
         print, 'Number of good data points fewer the NORD'
         return, sset
   endif



   ;----------
   ; Iterate spline fit

   iiter = 0
   error = -1L
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


        if error[0] NE 0 then begin
              bf1 = bspline_action(1.d*xdata, sset, lower=laction, upper=uaction) 
;           if keyword_set(DBL) then $
;              bf1 = bspline_action(xdata, sset, lower=laction, upper=uaction) $
;           else $
;              bf1 = bspline_action(float(xdata), sset, lower=laction, upper=uaction) 
          if (n_elements(bf1) NE nx*nord) then begin
             splog, 'ERROR: BSPLINE_ACTION failed!'
             return, 0
          endif

;          if npoly EQ 1 then bf_full = bf1 $
;          else bf_full = transpose(rebin(bf1, nx, nord, npoly) ,[0,2,1])
;          action = action_multiple * bf_full

          ; Same thing as above using less memory
          action = action_multiple
          for ipoly=0L, npoly-1L do $
           action[*,lindgen(nord)*npoly+ipoly] *= bf1
          bf1 = 0  ; clear memory
        endif

        ; help, where(finite(bf1) EQ 0), where(finite(action) EQ 0)
        if total(finite(action) EQ 0) GT 0 then begin
            splog, 'ERROR: Infinities in action matrix, wavelengths may be very messed up!!!'
            return, 0
        endif

        error = bspline_workit(xdata, ydata, invvar*maskwork, action, $
            sset, alpha=alpha, lower=laction, upper=uaction, $
            yfit=yfit, covariance=covariance)
      endelse

      iiter = iiter + 1

      if (error[0] EQ -2L) then begin
         ; All break points have been dropped.
          print, 'bspline_longslit: All break points have been dropped!!'
         return, sset
      endif else if (error[0] EQ 0) then begin
         ; Iterate the fit -- next rejection iteration.
          chi_array = (ydata - yfit) * sqrt(invvar*maskwork)
          reduced_chi = total(chi_array^2) / $
                              (ngood - npoly*(n_elements(goodbk)+nord) - 1)


          relative_factor = 1.0
          if keyword_set(relative) then begin
            nrel  = n_elements(relative)
            if nrel EQ 1 then relative_factor = sqrt(reduced_chi) $
            else begin
                  this_chi2 = total(chi_array[relative]^2) / $
                               (nrel - (n_elements(goodbk)+nord) - 1)
                  relative_factor = sqrt(this_chi2) 
            endelse
            relative_factor = relative_factor > 1
          endif
             
         qdone = djs_reject(ydata, yfit, invvar=invvar, inmask=tempin, $
                            outmask=maskwork, upper=upper*relative_factor, $
                            lower=lower*relative_factor, _EXTRA=EXTRA)

          tempin = maskwork
      if NOT keyword_set(silent) then $
          splog, iiter, reduced_chi, total(maskwork EQ 0), relative_factor, $
                       format='(i4, f8.3, i7, f6.2)'

      endif
      ;; JXP DEBUGGING
      if keyword_set(DEBUG) then begin
         ndbg = n_elements(DEBUG)
         model = fltarr(355L, 2048L)
         model[DEBUG[0:ndbg/2-1]] = yfit
         xatv, model, /bloc
         modelfit = fltarr(n_elements(yfit))
         modelfit[debug[ndbg/2:*]] = yfit
         model2 = fltarr(355L, 2048L)
         model2[DEBUG[0:ndbg/2-1]] = modelfit
         xatv, model2, /bloc
         stop
      endif

  endwhile

;   print, format='($, i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   if NOT keyword_set(silent) then $
     print, format='(i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   outmask = maskwork
   ;----------
   ; Re-sort the output arrays OUTMASK and YFIT to agree with the input data.

   if NOT keyword_set(silent) then splog, 'Elapsed time: ', systime(1)-t0

   return, sset
end
;------------------------------------------------------------------------------
