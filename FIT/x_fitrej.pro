;+ 
; NAME:
; x_fitrej
;   Version 1.2
;
; PURPOSE:
;    Fits a function to a set of x,y data with rejection!
;
; CALLING SEQUENCE:
;   
;   fit =
;   x_fitrej(xdat,ydat,func=,nord=,sig=,reg=,
;            lsigma=,hsigma=,maxrej=,niter=,minpt=)
;
; INPUTS:
;   xdat       - Values along one dimension
;   ydat       - Values along the other
;   [func]     - String for Fitting function (POLY, LEGEND, BSPLIN)
;   [nord]     - Order of the fit
;
; RETURNS:
;   fit        - Values at each xdat
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   inter      - Interactive fitting
;   sig        - Errors in the points
;   reg        - Regions of data to fit
;   lsigma     - Lower sigma threshold
;   hsigma     - Upper sigma threshold
;   maxrej     - Maximum points to reject per iteration
;   niter      - Number of iterations for rejection
;   minpt      - Minimum # of points
;   FITSTR     - Fit structure :: OVERIDES ALL OTHER INPUT VALUES!
;
; OPTIONAL OUTPUTS:
;   bsplin     - Bspline structure
;   rms        - RMS of the fit
;   REJPT      - Rejected points
;   GDPT       - Good (unrejected) points
;   NRM      - Normalization numbers for xdat :: dblarr(2)
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_fitrej(x, y, 'POLY', nord=5)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  POLY_FIT
;  POLY
;  SVDFIT
;  SVLEG
;  DJS_REJECT
;
; REVISION HISTORY:
;   20-Nov-2001 Written by JXP
;   31-Jan-2002 Revised significantly.  Added fit structure, CHEBY
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fitrej, xdat, ydat, func, nord, SIG=sig, REG=reg, $
                   LSIGMA=lsigma, HSIGMA=hsigma, MAXREJ=maxrej, $
                   NITER=niter, MINPT=minpt, FFIT=ffit, REJPT=rejpt, $
                   MSK=msk, RMS=rms, GDPT=gdpt, NRM=nrm, FITSTR=fitstr, $
                   NONRM=nonrm, IVAR=ivar

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'fit = x_fitrej(xdat, ydat, [func, nord], '
    print, '        [SIG=, REG=, LSIGMA=, HSIGMA=, MAXREJ=, NITER=, MINPT=,'
    print, '        REJPT=, FFIT=, MSK=, NRM=, FITSTR=, /NONRM, IVAR=)  [v1.2]'
    return, -1
  endif 

  igood = lindgen(n_elements(xdat))

;  Optional Keywords

  ; Fit structure
  if keyword_set(FITSTR) then begin
      func = fitstr.func
      nord = fitstr.nord
      lsigma = fitstr.lsig
      hsigma = fitstr.hsig
      minpt = fitstr.minpt
      maxrej = fitstr.maxrej
  endif else begin
      if (not keyword_set( LSIGMA ) AND $
          not keyword_set( HSIGMA )) then return, -1
      if not keyword_set( NITER ) then niter = 9999L
      if not keyword_set( MINPT ) then minpt = 1L
      if not keyword_set( MAXREJ ) then maxrej = 100L
      if not keyword_set( FUNC ) then func = 'POLY'
      if not keyword_set( NORD ) then nord = 0
      ; Create a fit structure
      fitstr = { fitstrct }
      fitstr.func = func
      fitstr.nord = nord
      fitstr.niter = niter
      fitstr.minpt = minpt
      fitstr.maxrej = maxrej
      fitstr.lsig = lsigma
      fitstr.hsig = hsigma
  endelse

; Mask
  if keyword_set( MSK ) then igood = where(msk NE 0)

; Normalize xdat

  if not keyword_set( NONRM ) then begin
      nrm = dblarr(2)
      mnx = min(xdat, MAX=mxx)
      nrm[0] = 0.5 * (mnx + mxx)
      nrm[1] = mxx - mnx
      xnrm = 2. * (xdat - nrm[0])/nrm[1]
      fitstr.nrm = nrm
  endif else xnrm = xdat

;   Pre-set Regions
  if keyword_set( REG ) then begin
      regpt = x_fndreg(xdat[igood], reg, NPNT=npnt)  ; Find the regions 
      if npnt EQ 0 then return, -1
      igood = temporary(igood[regpt])
  endif 

  xfit = xnrm[igood]
  yfit = ydat[igood]
  if keyword_set( SIG ) then sigfit = sig[igood]
  if keyword_set( IVAR ) then ivarfit = ivar[igood]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; MAIN LOOP
  
;  First Fit
  fit = x_fit(xfit, yfit, SIG=sigfit, IVAR=ivarfit, FITSTR=fitstr, /NONRM)
  if fit[0] EQ -1 then return, -1

  iter = 0L
  while iter LT fitstr.niter do begin
      
;   Reject  

      tmpmax = (n_elements(xfit) - fitstr.minpt) < fitstr.maxrej
      if tmpmax LE 0 then begin
          outmask = intarr(n_elements(xfit)) + 1
          break
      endif

      qmod = djs_reject(yfit, fit, outmask=outmask, upper=fitstr.hsig, $
                        lower=fitstr.lsig, maxrej=tmpmax, /sticky) ;, $
;                       SIGMA=sigfit)
      ;; Consider the above line more carefully!  JXP  (9/23/03)

;   Check number of good points
      tmp = where(outmask EQ 1, ntmp)
      if (qmod EQ 1 OR ntmp LE fitstr.minpt) then break $
        else iter = iter + 1

;   FIT
      fit = x_fit(xfit, yfit, SIG=sigfit, FITSTR=fitstr, MSK=outmask, /NONRM, $
                 IVAR=ivarfit)
      if fit[0] EQ -1 or n_elements(fit) EQ 1 then return, -1

  endwhile

;  Calculate fit everywhere

  allfit = x_calcfit(xdat, FITSTR=fitstr)

; Rejected points
  if arg_present( REJPT ) then begin
      rej = where(outmask EQ 0, nrej)
      if nrej EQ 0 then rejpt = -1 else rejpt = igood[rej]
      delvarx, rej
  endif

; Good points
  if arg_present( GDPT ) then begin
      gd = where(outmask NE 0, ngd)
      gdpt = igood[gd]
      delvarx, gd
  endif


; RMS

  gd = where(outmask NE 0, ngd)
  rms = sqrt(total( (fit[gd]-yfit[gd])^2 )/(ngd-1.))
  fitstr.rms = rms

  delvarx, xfit, yfit, fit, sigfit

  ; Return fit as double if xdat or ydat is double
  if size(ydat,/type) EQ 5 OR  size(xdat,/type) EQ 5 then $
    return, double(allfit) $
  else return, float(allfit)

end
