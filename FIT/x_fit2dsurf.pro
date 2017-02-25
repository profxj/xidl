;+ 
; NAME:
; x_fit2dsurf
;   Version 1.1
;
; PURPOSE:
;    Fits a 2d surface to a set of x,y data with errors (optional).
;    It is recommended to set the FIT parameters with the structure
;    FITSTR. 
;
; CALLING SEQUENCE:
;   
;   fit = x_fit2dsurf(xy, z,[sig], nx=, ny=, MSK=, FITSTR=, REJPT=,
;   NRM=, /NONRM, /SVDFT, LSIG=, HSIG=, NITER=, MAXREJ=, FUNC= )
;
; INPUTS:
;   xy         - 2d array of xy pairs: [N,2]
;   z          - Values at each xy pair
;   [sig]      - Variance in each data point
;
; RETURNS:
;   fit        - Values at each xdat
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   msk        - Mask  (0 = Do NOT include)
;   func       - 2d function ('POLY')
;   nx         - Order in x
;   ny         - Order in y
;   NONRM      - Don't normalize data before fitting
;   MINPT      - Minimum number of points to keep in fit
;   LSIG,HSIG  - Lower and upper sigma rejection
;   /SVDFT     - Use the IDL routine SVDFT for the fitting.   
;                Otherwise use CHOLDC (recommended for speed)
;
; OPTIONAL OUTPUTS:
;   ffit     - Functional form
;   NRM      - [2,2] array of Normalization coefficients
;   REJPT    - Points rejected in the fit
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_fit2dsurf(xy,z, sig)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  f2dpoly
;  CHOLDC
;  SVDFT
;  x_calc2dfit
;
	; REVISION HISTORY:
;   31-Jan-2002 Written by JXP
;   02-Feb-2002 Added rejection
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fit2dsurf, xydat, zval, sig, NX=nx, NY=ny, MSK=msk, FUNC=func, $
                      NONRM=nonrm, NRM=nrm, FITSTR=fitstr, SVDFT=svdft, $
                      LSIG=lsig, HSIG=hsig, NITER=niter, MAXREJ=maxrej, $
                      MINPT=minpt, REJPT=rejpt

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'fit = x_fit2dsurf(xydat, zval, [sig], NX=, NY=, MSK=,' 
    print, '                /NONRM, NRM=, FITSTR=, /SVDFIT, LSIG=,'
    print, '                HSIG=, NITER=, MAXREJ=, MINPT=) [v1.1]'
    return, -1
  endif 


;  Optional Keywords

  ; Fit structure
  if keyword_set(FITSTR) then begin
      func = fitstr.func
      nx = fitstr.nx
      ny = fitstr.ny
      lsig = fitstr.lsig
      hsig = fitstr.hsig
      niter= fitstr.niter
      minpt= fitstr.minpt
      maxrej= fitstr.maxrej
  endif 

  ; func, nord
  if not keyword_set( FUNC ) then func = 'POLY'
  if not keyword_set( nx ) then nx = 3
  if not keyword_set( ny ) then ny = 3
  if not keyword_set( niter ) then niter=1
  if not keyword_set( LSIG ) then lsig=0.
  if not keyword_set( HSIG ) then hsig=0.
  if not keyword_set( MINPT ) then minpt=0

  ; mask
  if keyword_set( MSK ) then igood = where(msk NE 0) $
  else igood = lindgen(n_elements(xydat[*,0]))

  ; sigma array
  if not keyword_set( SIG ) then sig = fltarr(n_elements(zval))+1.

;;;;;;;;;;;;
; Normalize xydat

  if not keyword_set( NONRM ) then begin
      ; x first
      nrm = dblarr(2,2)
      mnx = min(xydat[*,0], MAX=mxx)
      nrm[0,0] = 0.5 * (mnx + mxx)
      nrm[0,1] = mxx - mnx
      xnrm = 2. * (xydat[*,0] - nrm[0,0])/nrm[0,1]
      ; now y
      nrmy = dblarr(2)
      mny = min(xydat[*,1], MAX=mxy)
      nrm[1,0] = 0.5 * (mny + mxy)
      nrm[1,1] = mxy - mny
      ynrm = 2. * (xydat[*,1] - nrm[1,0])/nrm[1,1]
  endif else begin
      xnrm = xydat[*,0]
      ynrm = xydat[*,1]
  endelse

;;;;;;;;;;;;
;  Set Fit arrays

  xfit = double(xnrm[igood]) ; Normalized
  yfit = double(ynrm[igood])
; IVAR
  if not keyword_set(SVDFT) then ivar = 1./(sig^2)
                                
;;;;
; Setup Rejection
  if not keyword_set( maxrej ) then maxrej = long(0.1*n_elements(zval))
  niter = niter > 1L
  if (LSIG GT 0.) OR (HSIG GT 0.) then begin
      tmpniter = niter+1
      flg_rej = 1 
  endif else begin
     flg_rej = 0
     tmpniter = niter
  endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize Function

  case func of 
      'POLY': begin
          ;; Set flag
          if nx EQ ny then flg = 0 else flg = nx
          ;; Initialize
          tmp = f2dpoly(0,-1,XVAL=xfit,YVAL=yfit,FLG=flg)
          func2d = 'f2dpoly'
      end
      else: begin
          print, 'x_fit2dsurf: Invalid function declaration! '
          return, -1
      end
  endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Let's Fit!

  for iiter=0L, tmpniter-1 do begin
      ; SVDFIT
      if keyword_set(SVDFT) then begin
          ffit = svdfit(igood, zval[igood], nx*ny, measure_errors=sig[igood], $
                        function_name=func2d, $
                        /double, singular=singular)
      endif else begin  ; NORMAL EQUATIONS (Best for few points)
          fit = zval[igood]
          A = call_function(func2d, igood, nx*ny)
          alpha = transpose(A) # (A * (replicate(1,nx*ny) # ivar[igood]))
          beta = transpose((zval[igood] * ivar[igood]) # A)
          
          choldc, alpha, p, /double
          ffit = cholsol(alpha, p, beta)
          fit[*] = A # ffit
      endelse
      
      ; REJECTION
      if(flg_rej EQ 1) then begin
          outmask = 0 
          tmpmask = (n_elements(zval[igood]) - minpt) < maxrej
          qdone = djs_reject(zval[igood], fit, outmask=outmask, $ ;invvar=ivar, $
                             upper=hsig, lower=lsig, /sticky, maxrej=tmpmax)
          
          if keyword_set( IVAR ) then ivar[igood] = ivar[igood]*outmask
          igood = igood[where(outmask NE 0)]
          
          rejpt = where(outmask EQ 0, nrej)
          if qdone EQ 1 then break
      endif else break ; No Rej so break
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Uninitialize Function

  case func of 
      'POLY': tmp = f2dpoly(0,-2)
      else: begin
          print, 'x_fit2dsurf: Invalid function declaration! '
          return, -1
      end
  endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the surface values at all xydat

  ; Release memory
  delvarx, xnrm, xfit, yfit, sigfit, zfit

; FIT EVERYWHERE
  fit = x_calc2dfit(xydat, func, ffit, nx, ny, NRM=nrm)

; RMS
  gd = where(outmask NE 0, ngd)
  rms = sqrt(total( (fit[igood[gd]]-zval[igood[gd]])^2 )/(ngd-1.))
  
; Write to fit structure

  if arg_present(fitstr) then begin
      if not keyword_set(fitstr) then fitstr = { fit2dstrct }
      if not keyword_set( NONRM ) then fitstr.nrm = nrm
      fitstr.ffit = ptr_new(ffit)
      fitstr.nx = nx
      fitstr.ny = ny
      fitstr.func = func
      fitstr.lsig = lsig
      fitstr.hsig = hsig
      fitstr.niter = niter
      fitstr.minpt = minpt
      fitstr.maxrej = maxrej
      fitstr.rms = rms
  endif

  ; Return val as double if xval is double or float otherwise
  if size(xydat,/type) EQ 5 OR  size(zval,/type) EQ 5 then $
    return, double(fit) $
  else return, float(fit)

end
