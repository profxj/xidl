;+ 
; NAME:
; x_maxschecht   
;    Version 1.1
;
; PURPOSE:
;  Run a Maximum Likelihood fit with the Shecter Luminosity functional
;  form (i.e. a Gamma function)
;  The default is two find the maximum in two iterations of 100 steps
;
; CALLING SEQUENCE:
;  x_maxschect, arr, nmnx, amnx, val, NA=, NN=, STA=, STN=
;
; INPUTS:
;   Array
;
; RETURNS:
;
; OUTPUTS:
;  val == Best fit values  0=Ngamma; 1=alpha
;
; OPTIONAL KEYWORDS:
;  LOG -- Input N range is in logarithmic value
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_maxshect
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function x_maxschect_kscumf, x

common x_maxschect_cmmn, schect_k, schect_alph, schect_nmin, schect_nstr

return, schect_k * schect_nstr * gamma(schect_alph+1)* $
                    (igamma(schect_alph+1.,x/schect_nstr) - $
                    igamma(schect_alph+1.,schect_Nmin/schect_nstr)) 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_maxschecht, arr, inmnx, iamnx, val, NA=na, NN=nn, Nmin=nmin, $
                 LIK=lik, KSPROB=ksprob, NOISE=noise, ERR=err, CL=cl, $
                  LOG=log

common x_maxschect_cmmn
  ; 
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_maxschecht, arr, nmnx, amnx, val, NA=, NN=, LIK= [v1.1]'
      return
  endif 
  if not keyword_set(NMIN) then stop
  if not keyword_set(NN) then nn = 100L
  if not keyword_set(NA) then na = 100L
  if not keyword_set(CL) then cl = 0.683
  if not keyword_set(NITER) then niter = 2L

  nmnx = inmnx
  amnx = iamnx
  for qq=0L,niter-1 do begin
      ;; Set the step size
      stn = (nmnx[1]-nmnx[0])/double(NN)
      sta = (amnx[1]-amnx[0])/double(NA)

      ;; Calculate image of parameters
      nvec = nmnx[0] + dindgen(nn)*stn
      if keyword_set(LOG) then nvec = 10^nvec
      avec = amnx[0] + dindgen(na)*sta
      nimg = nvec # replicate(1.,nn)
      aimg = replicate(1.,na) # avec

      ;; Likelihood
      narr = n_elements(arr)
      lik = aimg*total(alog(arr)) - aimg*alog(Nimg)*narr - $
        total(arr) / Nimg - narr*alog(Nimg) - $
        narr*alog( gamma(aimg+1) * (1. - igamma(aimg+1., Nmin/Nimg)))
      maxL = max(lik,imx)

      val = [ nimg[imx], aimg[imx] ]

      ;; Val
;      if not keyword_set(SILENT) then print, 'x_maxschect: val = ', val

      ;; Reset range
      if keyword_set(LOG) then nmnx = alog10(val[0]) + stn*[-1.,1] $
      else nmnx = val[0] + stn*[-1.,1]
      amnx = val[1] + sta*[-1.,1]
  endfor

  ;; KS Test
  if arg_present(KSPROB) then begin
      schect_k = 1. / val[0] / (gamma(val[1]+1.)* $
                              (1.-igamma(val[1]+1.,Nmin/val[0])))
      schect_alph = val[1]
      schect_nstr = val[0]
      schect_nmin = nmin
      if keyword_set( NOISE ) then darr = arr * 10^(noise*randomu(-2442,narr)) $
      else darr = arr
      ;; 
      ksone, darr, 'x_maxschect_kscumf', d, ksprob, PLOT=plot 
      if keyword_set( DEBUG ) then stop
  endif

  ;; Error -- One parameter uncertainties
  if arg_present(ERR) then begin
      nmnx = inmnx
      amnx = iamnx
      stn = (nmnx[1]-nmnx[0])/double(NN)
      sta = (amnx[1]-amnx[0])/double(NA)
      ;; Nval
      nvec = nmnx[0] + dindgen(nn)*stn
      if keyword_set(LOG) then nvec = 10^nvec
      ;; Likelihood
      lgl = val[1]*total(alog(arr)) - val[1]*alog(nvec)*narr - $
        total(arr) / nvec - narr*alog(nvec) - $
        narr*alog( gamma(val[1]+1) * (1. - igamma(val[1]+1., Nmin/nvec)))
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd EQ 0 then stop
      ;; Error
      errn = x_maxcl(nvec[gd], lgl[gd], val[0], CL=cl)
      ;; aval
      avec = amnx[0] + dindgen(na)*sta
      lgl = avec*total(alog(arr)) - avec*alog(val[0])*narr - $
        total(arr) / val[0] - narr*alog(val[0]) - $
        narr*alog( gamma(avec+1) * (1. - igamma(avec+1., Nmin/val[0])))
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd EQ 0 then stop
      ;; Error
      erra = x_maxcl(avec[gd], lgl[gd], val[1], CL=cl)
      err = [errn, erra]
      err = reform(err,2,2)
  endif
  return

end
