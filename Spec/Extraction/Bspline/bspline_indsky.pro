;+
; NAME:
;   bspline_indsky
;
; PURPOSE:
;   Bspline extraction with independent breakpoints for OBJ/SKY
;    need a clever way to keep track of sky versus obj coefficients.
;
; CALLING SEQUENCE:
;   sset = bspline_iterfit( )
;
; INPUTS:
;   xdata      - Data x values
;   ydata      - Data y values
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance of y; if not set, then set to be
;                consistent with the standard deviation.  This only matters
;                if rejection is being done.
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
;   bspline_workit()
;   bspline_fit()
;   create_bsplineset()
;   djs_reject()
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------
pro action_merge, act_obj, lowerobj, upperobj, act_sky, lowersky, uppersky, $
                  action, loweraction, upperaction, mapping

    ; goal here is to combine 2 action arrays into one.
    ; need to sort upper and lower limits and retain mapping.

      no = n_elements(lowerobj)
      ns = n_elements(lowersky)
      npix = (size(act_obj))[1]
      objbw = (size(act_obj))[2]
      skybw = (size(act_sky))[2]

      objkey = replicate('O',no)
      skykey = replicate('S',ns)

      lo = [lowerobj, replicate(lowerobj[no-1], objbw-1)]
      ls = [lowersky, replicate(lowersky[ns-1], skybw-1)]
      low = [lo,ls]
      uo = [upperobj, replicate(upperobj[no-1], objbw-1)]
      us = [uppersky, replicate(uppersky[ns-1], skybw-1)]
      
      repl = where(uo EQ -1)
      if repl[0] NE -1 then uo[repl] = ([lo[1:*],(npix)])[repl]-1
      if repl[0] NE -1 then lo[repl] = ([-1,uo])[repl]+1

      repl = where(us EQ -1)
      if repl[0] NE -1 then us[repl] = ([ls[1:*],(npix)])[repl]-1
      if repl[0] NE -1 then ls[repl] = ([-1,us])[repl]+1

      upper = [uo,us]
      upper = upper > [low[1:*] - 1, (npix-1)]
      

      sortkey = sort(low*(no+ns) + upper*(no+ns)^2 + lindgen(n_elements(low)))

      findl = total(histogram(upper, min=0, max=npix-1), /cumul) 

      diff = where(findl - [-1,findl] NE 0,t)
      mapping = sortkey[where(upper[sortkey] EQ diff[1], hgt)]
      loweraction = replicate(0, hgt)
      upperaction = replicate(diff[1], hgt)
      up1= diff[1]
      for i=2, n_elements(diff)-1  do begin & $
         mapping = [mapping, sortkey[where(upper[sortkey] EQ $
                                diff[i], hgt)]] & $
         loweraction =  [loweraction, replicate(up1+1,hgt)] & $
         upperaction =  [upperaction, replicate(diff[i],hgt)] & $
         up1 = diff[i] & $
      endfor

      loweraction[1:*] = loweraction[1:*] > (upperaction + 1)

      ; find bandwidth required...

      ll = where(mapping LT no+objbw-1, complement=uu)
      bw = ((max(ll[objbw-1:*] - [-1,ll])) > (max(uu[skybw-1:*] - [-1,uu]))) 

      action = fltarr(npix, bw)
      nn = n_elements(low)
      diff2 = where(upperaction GT loweraction)

      ladderobj = total(histogram(uo,min=-1,max=npix-2), /cumul) 
      laddersky = total(histogram(us,min=-1,max=npix-2),/cumul) 
      ladder = total(histogram(upperaction,min=-1,max=npix-2),/cumul) 

      col = lindgen(npix) 
      for i=0, objbw -1 do begin 
        row = ll[ladderobj+i] - long(ladder) 
        action[col,row] = act_obj[*,i] 
      endfor

      for i=0, skybw -1 do begin 
        row = uu[laddersky+i] - long(ladder) 
        action[col,row] = act_sky[*,i] 
      endfor

      return
end


pro bspline_indsky, xdata, ydata, invvar, slit_profile, obj_profile, $
     yfit=yfit, numiter=iiter, bkpt=bkpt, maxiter=maxiter, relative=relative, $
     upper=upper, lower=lower, outmask=outmask, fullbkpt=fullbkpt, $
     inmask=inmask, _EXTRA=EXTRA, RED_CHI=reduced_chi, skybkpt=skybkpt, $
     objset=objset, skyset=skyset

   if (n_params() LT 2) then begin
      print, 'Syntax -  sset = bspline_extract( )'
      return
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

   objset = create_bsplineset(fullbkpt, nord) 
;   objset.funcname = 'Object extraction'
   skyset = create_bsplineset(skybkpt, nord) 
;   skyset.funcname = 'Sky extraction'

   if (nthese LT nord) then begin
         print, 'Number of good data points fewer the NORD'
         return
   endif


   ;----------
   ; Iterate spline fit

   iiter = 0
   error = 0
   tempin = 0
   reduced_chi = 0

   while (((error[0] NE 0) OR (keyword_set(qdone) EQ 0)) $
    AND iiter LE maxiter) do begin

      ngood = total(maskwork)
      goodobjbk = where(objset.bkmask[objset.nord:*] NE 0)
      goodskybk = where(skyset.bkmask[skyset.nord:*] NE 0)
      if (ngood LE 1 OR goodobjbk[0] EQ -1 OR goodskybk[0] EQ -1) then begin
         objset.coeff = 0
         skyset.coeff = 0
         iiter = maxiter + 1; End iterations
      endif else begin
        ; Do the fit.  Return values for ERROR are as follows:
        ;    0 if fit is good
        ;   -1 if all break points are masked
        ;   -2 if everything is screwed

        ;  we'll do the fit right here..............


        filler = replicate(1,nord)
        act_obj = bspline_action(xdata, objset, lower=lowerobj,upper=upperobj) $
                * (obj_profile[*] # filler)

        act_sky = bspline_action(xdata, skyset, lower=lowersky,upper=uppersky) $
                * (slit_profile[*] # filler)

        action_merge, act_obj, lowerobj, upperobj, act_sky, lowersky,uppersky,$
                  action, loweraction, upperaction, mapping

        error = bspline_workind(xdata, ydata, invvar*maskwork, action, $
            lower=loweraction, upper=upperaction, alpha=alpha, beta=beta)
        if error[0] EQ -1 then begin
           ll = where(mapping LT n_elements(goodobjbk), complement=uu)
           objset.coeff[goodobjbk] = beta[ll]
           objset.icoeff[goodobjbk] = alpha[0,ll]
           if n_elements(uu) NE n_elements(goodskybk) then $
              error = -1L $
           else begin 
             skyset.coeff[goodskybk] = beta[uu]
             skyset.icoeff[goodskybk] = alpha[0,uu]
             objfit = bspline_valu(xdata, objset) * obj_profile[*]  
             skyfit = bspline_valu(xdata, skyset) * slit_profile[*] 
             yfit = objfit + skyfit 
             ngd_bk = n_elements(beta)
             error = 0L
           endelse
        endif else begin
           ll = where(mapping[error] LT n_elements(objset.coeff),complement=uu)
           if ll[0] NE -1 then oe = bspline_maskpoints(objset, $
                    mapping[error[ll]], 1) $
           else oe = -1L

           if uu[0] NE -1 then se = bspline_maskpoints(skyset, $
                    mapping[error[uu]]-n_elements(objset.coeff), 1) $
           else se = -1L

           error = (oe[0] NE -2 AND se[0] NE -2) ? -1 : -2 
        endelse
      endelse

      iiter = iiter + 1

      if (error[0] EQ -2L) then begin
         ; All break points have been dropped.
         return
      endif else if (error[0] EQ 0) then begin
         ; Iterate the fit -- next rejection iteration.
          reduced_chi = total((ydata - yfit)^2  * (invvar*maskwork)) / $
                              (ngood - ngd_bk - 1)


          relative_factor = 1.0
          if keyword_set(relative) then $
                relative_factor = sqrt(reduced_chi) > 1.0

         qdone = djs_reject(ydata, yfit, invvar=invvar, inmask=tempin, $
                            outmask=maskwork, upper=upper*relative_factor, $
                            lower=lower*relative_factor, _EXTRA=EXTRA)

          tempin = maskwork
          print, format='($, i4, f7.3,i5,63x," ")', iiter, reduced_chi, $
            total(maskwork EQ 0)
          if finite(reduced_chi) EQ 0 then begin
             print, 'NaN in answer'
             stop
          endif
;            total(maskwork EQ 0), replicate(string(8b),16)
      endif

  endwhile

;   print, format='($, i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   print, format='(i4, f7.3, i5)', iiter, reduced_chi, total(maskwork EQ 0)
   outmask = maskwork
   ;----------
   ; Re-sort the output arrays OUTMASK and YFIT to agree with the input data.

   return
end
;------------------------------------------------------------------------------
