;+ 
; NAME:
; x_maxdblpow   
;    Version 1.1
;
; PURPOSE:
;  Run a Maximum Likelihood fit on a double power law of the form:
;    (N/Nd)^beta where (beta = b1 for N<Nd and beta=b2 for N>Nd)
;
; CALLING SEQUENCE:
;  x_maxdblpow, arr, nmnx, b1mnx, b2mnx, val, NA=, NN=, STA=, STN=
;
; INPUTS:
;   Array
;   nmnx = Logarithmic range of break point
;
; RETURNS:
;
; OUTPUTS:
;  val == Best fit values
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_maxdblpow
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function x_maxdblpow_kscumf, x

common x_maxdblpow_cmmn, dblpow_k, dblpow_nd, dblpow_nmin, dblpow_b1, dblpow_b2

  high = where(x GT dblpow_nd, nhigh, complement=low, ncomplement=nlow)
  val = dblarr(n_elements(x))
  if nhigh NE 0 then begin
      val[high] =  dblpow_k * dblpow_nd * $
        (( 1. - (dblpow_Nmin/dblpow_nd)^(dblpow_b1+1))/(1.+dblpow_b1) + $
        ((x[high]/dblpow_nd)^(dblpow_b2+1) -1)/(1.+dblpow_b2))
  endif
  if nlow NE 0 then begin
      val[low]= dblpow_k * dblpow_nd * $
        ( (x[low]/dblpow_nd)^(dblpow_b1+1) - $
          (dblpow_Nmin/dblpow_nd)^(dblpow_b1+1))/ (1.+dblpow_b1) 
  endif
  if n_elements(x) EQ 1 then val = val[0]
  
  return, val
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_maxdblpow, arr, nmnx, b1mnx, b2mnx, val, NB1=nb1, NN=nn, Nmin=nmin, $
                 LIK=lik, B1STP=b1stp, B2STP=b2stp, NB2=nb2, STN=stn, $
                 KSPROB=ksprob, NOISE=noise, NMAX=nmax, ERR=err

common x_maxdblpow_cmmn

  ; 
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'x_maxdblpow, arr, nmnx, b1mnx, b2mnx, val, NB1=, NB2=, NN=, LIK= [v1.1]'
      return
  endif 
  if not keyword_set(NMIN) then stop
  if not keyword_set(NMAX) then nmax = 1d99

  if not keyword_set(STN) then begin
      if not keyword_set( NN ) then begin
          print, 'x_maxshect: Must set STN or NN'
          return
      endif
      stn = (nmnx[1]-nmnx[0])/float(NN)
  endif else NN = round( (nmnx[1]-nmnx[0]) / stn)

  ;; Beta1
  if not keyword_set(B1STP) then begin
      if not keyword_set( NB1 ) then begin
          print, 'x_maxshect: Must set B1SPT or NB1'
          return
      endif
      b1stp = (b1mnx[1]-b1mnx[0])/float(NB1)
  endif else NB1 = round( (b1mnx[1]-b1mnx[0]) / b1stp)

  ;; Beta2
  if not keyword_set(B2STP) then begin
      if not keyword_set( NB2 ) then begin
          print, 'x_maxshect: Must set B2SPT or NB2'
          return
      endif
      b2stp = (b2mnx[1]-b2mnx[0])/float(NB2)
  endif else NB2 = round( (b2mnx[1]-b2mnx[0]) / b2stp)

  ;; Calculate image of parameters
  beta1 = b1mnx[0] + dindgen(nb1)*b1stp
  beta2 = b2mnx[0] + dindgen(nb2)*b2stp
  b1img = beta1 # replicate(1., nb1)
  b2img = replicate(1., nb2) # beta2

  ;; Liklihood
  narr = n_elements(arr)

  svb1 = dblarr(nn)
  svb2 = dblarr(nn)
  svmx = dblarr(nn)
  Nsval = Nmnx[0] + findgen(nn)*stn
  for ii=0L,nn-1 do begin
      Ns = nsval[ii]
      ns10 = 10^Ns
      ;; Calculate denominator for all beta
      denom = ns10*( (1. - (Nmin/Ns10)^(b1img+1.))/(1.+b1img) + $
                     ((Nmax/Ns10)^(b2img+1) -1.)/(b2img+1.))

      ;; Numerator
      high = where(arr GT Ns10,nhigh,complement=low,ncomplement=nlow)
      if nhigh NE 0 then num1 = b2img * total(alog(arr[high]/ns10)) $
      else num1 = 0.
;      low = where(arr LE Ns,nlow)
      if nlow NE 0 then num2 = b1img * total(alog(arr[low]/ns10)) $
      else num2 = 0.
      num = num1 + num2

      ;; Liklihood  (Signs!)
      lik = (num1+num2) - narr*alog(denom)
      svmx[ii] = max(lik,imx)
      svb1[ii] = b1img[imx]
      svb2[ii] = b2img[imx]
  endfor
  maxL = max(svmx,idx)
  val = [nsval[idx], svb1[idx], svb2[idx]]

  ;; KS Test
  if arg_present(KSPROB) then begin
      ns10 = 10^val[0]
      dblpow_k = 1. / ( ns10*(1. - (Nmin/Ns10)^(val[1]+1))/(1.+val[1]) + $
                        ((Nmax/Ns10)^(val[2]+1) - 1.)/(val[2]+1)) 
      dblpow_b1 = val[1]
      dblpow_b2 = val[2]
      dblpow_nd = ns10
      dblpow_nmin = nmin
      if keyword_set( NOISE ) then darr = arr * 10^(noise*randomu(-2442,narr)) $
      else darr = arr
      ;; 
      ksone, darr, 'x_maxdblpow_kscumf', d, ksprob, PLOT=plot 
      if keyword_set( DEBUG ) then begin
          x_psclose
          ksone, darr, 'x_maxdblpow_kscumf', d, ksprob, /PLOT
          stop
      endif
  endif

  ;; Error -- One parameter uncertainties
  if arg_present(ERR) then begin
      stn = (nmnx[1]-nmnx[0])/double(NN)
      stb1 = (b1mnx[1]-b1mnx[0])/double(Nb1)
      stb2 = (b2mnx[1]-b2mnx[0])/double(Nb2)
      ;;;;;;;;;;;;;;;;;;;;;
      ;; Nsval
      Nsval = Nmnx[0] + findgen(nn)*stn
      ;; Likelihood
      ns10 = 10^Nsval
      denom = ns10*( (1. - (Nmin/Ns10)^(val[1]+1.))/(1.+val[1]) + $
                     ((Nmax/Ns10)^(val[2]+1) -1.)/(val[2]+1.))
      lgl = dblarr(nn)
      for ii=0L,nn-1 do begin
          ;; Numerator
          high = where(arr GT Ns10[ii],nhigh,complement=low,ncomplement=nlow)
          if nhigh NE 0 then num1 = val[2] * total(alog(arr[high]/ns10[ii])) $
          else num1 = 0.
          if nlow NE 0 then num2 = val[1] * total(alog(arr[low]/ns10[ii])) $
          else num2 = 0.
          num = num1 + num2
          lgl[ii] = num - narr*alog(denom[ii])
      endfor
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd LT 20L then stop
      ;; Error
      errn = x_maxcl(nsval[gd], lgl[gd], val[0], CL=cl)
      ;;;;;;;;;;;;;;;;;;;;;
      ;; b1val
      b1val = b1mnx[0] + findgen(nb1)*stb1
      ;; Likelihood
      ns10 = 10^val[0]
      denom = ns10*( (1. - (Nmin/Ns10)^(b1val+1.))/(1.+b1val) + $
                     ((Nmax/Ns10)^(val[2]+1) -1.)/(val[2]+1.))
      ;; Numerator
      high = where(arr GT Ns10,nhigh,complement=low,ncomplement=nlow)
      if nhigh NE 0 then num1 = val[2] * total(alog(arr[high]/ns10)) $
      else num1 = 0.
      if nlow NE 0 then num2 = b1val * total(alog(arr[low]/ns10)) $
      else num2 = 0.
      num = num1 + num2
      lgl = num - narr*alog(denom)
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd LT 20L then stop
      ;; Error
      errb1 = x_maxcl(b1val[gd], lgl[gd], val[1], CL=cl)
      ;;;;;;;;;;;;;;;;;;;;;
      ;; b2val
      b2val = b2mnx[0] + findgen(nb2)*stb2
      ;; Likelihood
      ns10 = 10^val[0]
      denom = ns10*( (1. - (Nmin/Ns10)^(val[1]+1.))/(1.+val[1]) + $
                     ((Nmax/Ns10)^(b2val+1) -1.)/(b2val+1.))
      ;; Numerator
      high = where(arr GT Ns10,nhigh,complement=low,ncomplement=nlow)
      if nhigh NE 0 then num1 = b2val * total(alog(arr[high]/ns10)) $
      else num1 = 0.
      if nlow NE 0 then num2 = val[1] * total(alog(arr[low]/ns10)) $
      else num2 = 0.
      num = num1 + num2
      lgl = num - narr*alog(denom)
      ;; Normalize
      lgl = lgl - maxL
      ;; Truncate
      gd = where(lgl GT -15.,ngd)
      if ngd LT 20L then stop
      ;; Error
      errb2 = x_maxcl(b2val[gd], lgl[gd], val[2], CL=cl)

      ;; Save
      err = [errn, errb1, errb2]
      err = reform(err,2,3)
  endif

  return

end
