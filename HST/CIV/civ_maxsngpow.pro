;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_maxsngpow.pro
; Author: Kathy Cooksey                      Date: 31 Dec 2008
; Project: CIV survey with Xavier Prochaska 
; Description: Fit the column density or EW frequency 
;              distribution with Maximum-Likelihood Method,
;              as outlined in Storrie-Lombarde et al (1996)
; Input: 
;    strct_fil -- civcandstrct of sample to fit
;    sens_fil -- cumulative pathlength X versus N and EW
; Optional:
;    guess -- initial guess [coefficient, power]
;    /siiv -- SiIV survey
;    intlim= -- limits of integration [min, saturated, max]
;    /ew -- fit f(EW), where EW is for CIV 1548
;    /ncolm -- fit f(N), where N is error-weighted column 
;              density
;    data_norm= -- normalization for power-law function (default:
;                  N_0 = 10.^14, W_0 = 400 mA for CIV; 
;                  N_0 = 10.^13.5, W_0 = 150 mA for SiIV)
;    conv_factor_lim= -- limit for convergence
;    infit_fil= -- name of previous fit file for starting point
;                  (also increase requirement for convergence)
;    /erronly -- re-calculate error ellipses only
;    /ksonly -- re-compute the K-S test statistics only
;    /plot -- show K-S and error-ellipses plots
;    /resmplerr -- increase grid resolution for computing error
;                  ellipses
;    nbin_alpha= -- number of exponent (alpha) bins for
;                   maximum likelihood grid (of size dalpha)
;    nbin_coeff= -- number of coefficient (normalization k) bins
;                   for maximum likelihood grid (of size dcoeff)
;    dcoeff= -- size of coeff bins
;    dalpha= -- size of alpha bins
;    int_param= -- 2-element array defining IDL's qromb's integration
;                  parameters ([eps, jmax])
;    /print_time -- print start and end times
;    /redshift -- fit for f(N) = Dnumber / (Dcolumn * Dz)
; Output: 
;    ostrct_fil -- name of output fit structure
; Example:
; History:
;    31 Dec 2008 -- created by KLC
;    27 Jan 2009 -- add effect of saturated features
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_sensitivity                ; need civ_sensitivity_x()
@civ_calcewn                    ; need civ_calcewn_ndblt()
@civ_omciv                      ; need civ_omciv_fit()
function civ_maxsngpow_int_freqdistr, var
  ;; Integral of (var/var_min)^px[1] * X(var) over d(var)
  common civ_maxsngpow_cmmn, sngpow_intlim, sngpow_data, sngpow_sigdata, $
     sngpow_xdata, sngpow_cumx, sngpow_coeff, sngpow_alpha, $
     sngpow_datatype, sngpow_flg, sngpow_datanorm,  qromb_par, sngpow_z

  if sngpow_datatype eq 'NCOLM' then $
     x_var = civ_sensitivity_x(sngpow_cumx,ncolm=alog10(var),z=z_var) $
  else x_var = civ_sensitivity_x(sngpow_cumx,ew=var,z=z_var) 
  if keyword_set(sngpow_z) then x_var = z_var ; f(n) defined over redshift

  return, (var/sngpow_datanorm)^sngpow_alpha * x_var
  
end                             ; civ_maxsngpow_intfn()


function civ_maxsngpow_int_deriv_freqdistr, var
  ;; Integral of ln(var/var_min) * (var/var_min)^alpha * X(var)
  common civ_maxsngpow_cmmn

  if sngpow_datatype eq 'NCOLM' then $
     x_var = civ_sensitivity_x(sngpow_cumx,ncolm=alog10(var),z=z_var) $
  else x_var = civ_sensitivity_x(sngpow_cumx,ew=var,z=z_var)
  if keyword_set(sngpow_z) then x_var = z_var ; f(n) defined over redshift
  
  return,  alog(var/sngpow_datanorm) * $
           (var/sngpow_datanorm)^sngpow_alpha * x_var

end                             ; civ_maxsngpow_intdfn()


function civ_maxsngpow_logL, px, gradient
  ;; px = [coeff, alpha]
  ;; ln L = -exp(px[0])*INT((var/var_0)^px[1] * 
  ;;                        X(var) * dvar, var_min, var_sat)
  ;;        -exp(px[0])*INT((var/var_0)^px[1] * 
  ;;                        X(var) * dvar, var_sat, var_max) 
  ;;        + px[0] * (n_elements(data[unsat]) + n_elements(data[sat]))
  ;;        + px[1] * SUM(ln(data/var_0))
  ;;        + SUM(ln(X(data)) 
  ;;        + n_elements(data[sat]) * 
  ;;          ln( INT((var/var_0)^px[1] * 
  ;;                  X(var) * dvar, var_sat, var_max) )
  ;;        ;- ln( n_elements(data[sat]) ! )
  common civ_maxsngpow_cmmn

  ;; Be smarter for grid searches
  sz = size(px)
  if sz[0] eq 1 then nelem = sz[0] $
  else nelem = sz[1]

  logL = dblarr(nelem)
  int_freqdistr = dblarr(nelem,2)
  if arg_present(gradient) then begin
     int_deriv_freqdistr = dblarr(nelem,2)
     gradient = dblarr(nelem,2)
  endif 

  if nelem eq 1 then begin
     coeff_arr = px[0]
     alpha_arr = px[1]
  endif else begin
     coeff_arr = px[*,0]
     alpha_arr = px[*,1] 
  endelse 

  ;; Determine saturation limit
  unsat = where((sngpow_flg and 2) eq 0 or (sngpow_flg and 6) ge 6,$
                nunsat,complement=sat,ncomplement=nsat)

  for ii=0L,nelem-1 do begin
     ;; Set for frequency distribution integrals
     sngpow_alpha = alpha_arr[ii]
     
     if ii gt 0 then begin
        if sngpow_alpha eq alpha_arr[ii-1] then begin
           ;; Save computational time and dont repeat integration
           int_freqdistr[ii,*] = int_freqdistr[ii-1,*]
           if keyword_set(gradient) then $
              int_deriv_freqdistr[ii,*] = int_deriv_freqdistr[ii-1,*]
           continue
        endif 
     endif 

     int_freqdistr[ii,0] = $
        qromb('civ_maxsngpow_int_freqdistr',sngpow_intlim[0],sngpow_intlim[1],$
              eps=0.5*qromb_par[0], /double, jmax=qromb_par[1]) 

     ;; Save time
     if nsat gt 0 then $
        int_freqdistr[ii,1] = $
        qromb('civ_maxsngpow_int_freqdistr',sngpow_intlim[1],sngpow_intlim[2],$
              eps=qromb_par[0], /double, jmax=qromb_par[1]) $
     else int_freqdistr[ii,1] = 0
     
     if arg_present(gradient) then begin
        ;; d(ln L)/d(px[1]) 
        int_deriv_freqdistr[ii,0] = $
           qromb('civ_maxsngpow_int_deriv_freqdistr',$
                  sngpow_intlim[0],sngpow_intlim[1], eps=0.5*qromb_par[0], $
                 /double, jmax=qromb_par[1])
        ;; Save time
        if nsat gt 0 then $
           int_deriv_freqdistr[ii,1] = $
           qromb('civ_maxsngpow_int_deriv_freqdistr',$
                  sngpow_intlim[1],sngpow_intlim[2], eps=qromb_par[0], $
                 /double, jmax=qromb_par[1]) $
        else int_deriv_freqdistr[ii,1] = 0
     endif                      ; gradient set 
  endfor                        ; loop ii = nelem

  sum_freqdistr = total(alog(sngpow_data[unsat]/sngpow_datanorm))

  logL = -exp(coeff_arr) * int_freqdistr[*,0] $
         -exp(coeff_arr) * int_freqdistr[*,1] + $
         coeff_arr * double(n_elements(sngpow_data)) + $
         alpha_arr * sum_freqdistr + total(alog(sngpow_xdata[unsat]))
  if nsat gt 0 then logL = logL + $
    nsat * alog( int_freqdistr[*,1] ) 

  if keyword_set(gradient) then begin
     ;; gradient = [d(ln L)/d(px[0]), d(ln L)/d(px[1])]
     ;; d(ln L)/d(px[0]) = -exp(px[0]) * INT((var/var_min)^px[1] * 
     ;;                                      X(var) * dvar, var_min, var_sat) 
     ;;                    -exp(px[0]) * INT((var/var_min)^px[1] * 
     ;;                                      X(var) * dvar, var_sat, var_max) 
     ;;                    + n_elements(data)
     ;; d(ln L)/d(px[1]) = -exp(px[1]) * $
     ;;                     INT(ln(var/var_min)*(var/var_min)^px[1] * 
     ;;                          X(var) * dvar, var_min, var_max)
     ;;                    -exp(px[1]) * $
     ;;                     INT(ln(var/var_min)*(var/var_min)^px[1] * 
     ;;                          X(var) * dvar, var_sat, var_max)
     ;;                    + SUM(ln(data/var_min)) 
     ;;                    + n_elements(data[sat]) * 
     ;;                      INT(ln(var/var_min)*(var/var_min)^px[1] * 
     ;;                          X(var) * dvar, var_sat, var_max) / 
     ;;                      INT((var/var_min)^px[1] * X(var) * dvar, 
     ;;                          var_sat, var_max)

     ;; d(ln L)/d(px[0]) 
     gradient[*,0] = -exp(coeff_arr) * int_freqdistr[*,0] $
                     -exp(coeff_arr) * int_freqdistr[*,1] + $
                     double(n_elements(sngpow_data))
     
     
     gradient[*,1] = -exp(coeff_arr)*int_deriv_freqdistr[*,0] $
                     -exp(coeff_arr)*int_deriv_freqdistr[*,1] $
                     + sum_freqdistr 
     if nsat gt 0 then gradient[*,1] = gradient[*,1] $
       + nsat * int_deriv_freqdistr[*,1] / int_freqdistr[*,1]

     if nelem eq 1 then $
        gradient = -gradient[0,*] $
     else gradient = -gradient  ; minimize (so negative)     
;     print,-logL[0],px[0],gradient[0],px[1],gradient[1] 
  endif                         ; gradient set

  if nelem eq 1 then begin
     return, -logL[0]          ; minimize
  endif else begin
     return, -logL
  endelse 
end                             ; civ_maxsngpow_logL()


function civ_maxsngpow_kscumf, var
  ;; CDF = INT( f(varp) * dvarp , var_min, var ) 
  ;;     = exp(coeff_norm)/(1 + alpha) * var_min^(-alpha) *
  ;;       (var^(1+alpha)*x(var) - var_min^(1+alpha)*x(var_min))
  ;; coeff_norm = alog((1+alpha)*var_min^alpha / 
  ;;                   (var_max^(1+alpha)*x(var_max) - 
  ;;                    var_min^(1+alpha)*x(var_min)))
  common civ_maxsngpow_cmmn
  k1 = exp(sngpow_coeff)/(1.+sngpow_alpha) * $
          sngpow_datanorm^(-sngpow_alpha)
  k2 = sngpow_intlim[0]^(1.+sngpow_alpha) 
  rslt = k1 * (var^(1.+sngpow_alpha) - k2[0])
  return, rslt
end


function civ_maxsngpow_omciv, px, intlim=intlim, siiv=siiv, _extra=extra
  ;; CIV mass density 
  ;; Analytic formula
  common civ_maxsngpow_cmmn

  if not keyword_set(intlim) then intlim = [10.^13,10.^15]

  ;; Be smarter for grid searches
  sz = size(px)
  if sz[0] eq 1 then nelem = sz[0] $
  else nelem = sz[1]

  logL = dblarr(nelem)
  int_freqdistr = dblarr(nelem,2)
  if arg_present(gradient) then begin
     int_deriv_freqdistr = dblarr(nelem,2)
     gradient = dblarr(nelem,2)
  endif 

  if nelem eq 1 then begin
     coeff_arr = px[0]
     alpha_arr = px[1]
  endif else begin
     coeff_arr = px[*,0]
     alpha_arr = px[*,1] 
  endelse 

  omciv = dblarr(nelem)
  for ii=0L,nelem-1 do begin
     omciv[ii] = civ_omciv_fit(intlim, exp(coeff_arr[ii]), $
                               alpha_arr[ii], sngpow_datanorm, siiv=siiv,$
                              _extra=extra)
  endfor 
  
  return, omciv
end                             ; civ_maxsngpow_omciv()


function civ_maxsngpow_mcerr,civcand
  ;; Monte-Carlo error analysis
  common civ_maxsngpow_cmmn 

  nciv = n_elements(civcand)

  if sngpow_datatype eq 'NCOLM' then begin
     ncolm = ncolm + randomn(seed,nciv)*signcolm
  endif else begin
     civcand.ew[0] = civcand.ew[0] + $
                     randomn(seed,nciv)*civcand.sigew[0]
  endelse 

end                             ; civ_maxsngpow_mcerr()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_maxsngpow, strct_fil, sens_fil, guess, siiv=siiv, data_norm=data_norm,$
                   intlim=intlim, ew=ew, ncolm=ncolm, ostrct_fil=ostrct_fil, $
                   conv_factor_lim=conv_factor_lim, infit_fil=infit_fil,$
                   erronly=erronly,ksonly=ksonly,plot=plot,resmplerr=resmplerr,$
                   nbin_alpha=nbin_alpha,nbin_coeff=nbin_coeff,$
                   dcoeff=dcoeff,dalpha=dalpha,int_param=int_param,$
                   print_time=print_time,redshift=redshift,_extra=extra

  ;; Instantiate common block
  ;; (sngpow_intlim, sngpow_data, sngpow_sigdata, sngpow_xdata, 
  ;;  sngpow_cumx, sngpow_coeff, sngpow_alpha, sngpow_datatype,
  ;;  sngpow_flg, sngpow_datanorm,  qromb_par, sngpow_z)
  common civ_maxsngpow_cmmn 
  if size(sens_fil,/type) eq 7 then $
     sngpow_cumx = xmrdfits(sens_fil,1,/silent) $
  else sngpow_cumx = sens_fil
  if keyword_set(int_param) then qromb_par = int_param $
  else qromb_par = [1.e-7, 20.] ; eps=, jmax=
  if keyword_set(redshift) then sngpow_z = 1 ; f(N) is over path length

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  if keyword_set(erronly) or keyword_set(ksonly) then begin
     ;; Read in previous results
     fitstrct = xmrdfits(infit_fil,1,/silent)
     civcand = xmrdfits(infit_fil,2,/silent)
     
     ;; Set up common values
     sngpow_datatype = strtrim(fitstrct.datatype,2)
     sngpow_datanorm = fitstrct.datanorm
     sngpow_data = fitstrct.data
     sngpow_sigdata = fitstrct.sigdata
     sngpow_flg = fitstrct.flag
     sngpow_xdata = fitstrct.xdata
     sngpow_intlim = fitstrct.intlim
     sngpow_coeff = alog(fitstrct.coeff)
     sngpow_alpha = fitstrct.alpha

     if keyword_set(ksonly) and not keyword_set(erronly) then $
        errstrct = xmrdfits(infit_fil,3,/silet)
  endif else begin
     ;;;;;;;;;;;;;;
     ;; Solve for coeff and alpha
     ;;;;;;;;;;;;;;
     ;; Data
     if size(strct_fil,/type) eq 7 then $
        strct = xmrdfits(strct_fil,1,/silent) $
     else strct = strct_fil
     ;; Trim strct based on both features being 3 sigma
     civ_group,civcand,strct,$
               flg_sys=384,nciv=nciv,_extra=extra ; rating=,zlim=
     if nciv eq 0 then $
        stop,'civ_maxsngpow: no '+dblt_name

     ;; Data must be in linear space
     ;; Sort for kstest
     ion_indx = where(stregex(civcand.ion,dblt_name,/boolean),dblt_flag)
     ion_indx = ion_indx[0]
     if keyword_set(ew) then begin
        sngpow_datatype = 'EW'
        if keyword_set(data_norm) then sngpow_datanorm = data_norm $
        else begin 
           if keyword_set(siiv) then sngpow_datanorm = 150.d $
           else sngpow_datanorm = 400.d ;130.d
        endelse 
        if keyword_set(resmplerr) then $
           civcand.ew[ion_indx] = civcand.ew[ion_indx] + $
                                  randomn(seed,nciv)*civcand.sigew[ion_indx]
        srt = sort(civcand.ew[ion_indx])
        civcand = civcand[srt]
        sngpow_data = double(civcand.ew[ion_indx])
        sngpow_sigdata = double(civcand.sigew[ion_indx])
        sngpow_flg = replicate(1l,nciv) ; dummy
        sngpow_xdata = civ_sensitivity_x(sngpow_cumx,ew=sngpow_data,z=z_var)

        if not keyword_set(ostrct_fil) then $
           ostrct_fil = 'few_fit.fits'
        scl = 1.5
     endif 
     if keyword_set(ncolm) then begin
        sngpow_datatype = 'NCOLM'
        if keyword_set(data_norm) then sngpow_datanorm = data_norm $
        else begin
           if keyword_set(siiv) then sngpow_datanorm = 10.d^13.5 $
           else sngpow_datanorm = 10.d^14. ;13.5
        endelse 
        ncolm  = civ_calcewn_ndblt(civcand,dblt_name,signcolm=signcolm,$
                                   flg_colm=flg_colm,/silent)
        if keyword_set(resmplerr) then $
           ncolm = ncolm + randomn(seed,nciv)*signcolm
        srt = sort(ncolm)
        civcand = civcand[srt]
        sngpow_data = ncolm[srt] ; error-weighted
        sngpow_sigdata = signcolm[srt]
        sngpow_flg = flg_colm[srt]

        sngpow_xdata = civ_sensitivity_x(sngpow_cumx,$
                                         ncolm=alog10(sngpow_data),z=z_var)
        if not keyword_set(ostrct_fil) then $
           ostrct_fil = 'fn_fit.fits'
        scl = 1.e2              ; 1.e3
     endif 
     if keyword_set(sngpow_z) then sngpow_xdata = z_var ; f(N) defined over z

     ;; Sub-sample, can be used to determine intlim[1] (saturated limit)
     unsat = where((sngpow_flg and 2) eq 0 or (sngpow_flg and 6) ge 6,$
                   nunsat,complement=sat,ncomplement=nsat)

     ;; Set integration limits
     sngpow_intlim = dblarr(3)
     if not keyword_set(intlim) then begin
        ;; Accept all data
        sngpow_intlim[0] = min(sngpow_data,max=mx)
        if nsat eq 0 then sngpow_intlim[1] = scl * mx $ ; something large
        else sngpow_intlim[1] = min(sngpow_data[sat])
        sngpow_intlim[2] = scl * mx ; something large 
     endif else begin
        if n_elements(intlim) ne 3 then $
           stop,'civ_maxsngpow: integration limit = 3-element array'
        sngpow_intlim = intlim
     endelse 

     ;; Trim saturated absorbers from unsaturated limits
     ;; and absorbers below lowest integration limit
     bd = where((sngpow_flg and 2) ge 2 and (sngpow_flg and 6) lt 6 $
                and sngpow_data lt intlim[1],nbd,complement=gd)
     if nbd ne 0 then begin
        print,'civ_maxsngpow: Number of saturated absorbers below intlim[1]',nbd
        sngpow_data = sngpow_data[gd]
        sngpow_sigdata = sngpow_sigdata[gd]
        sngpow_flg = sngpow_flg[gd]
        sngpow_xdata = sngpow_xdata[gd]
        civcand = civcand[gd]
     endif 
     ;; Trim absorbers outside the limits of interest
     bd = where(sngpow_data lt intlim[0] or $
                sngpow_data gt intlim[2],nbd,complement=gd)
     if nbd ne 0 then begin
        print,'civ_maxsngpow: Number of absorbers outside intlim',nbd
        sngpow_data = sngpow_data[gd]
        sngpow_sigdata = sngpow_sigdata[gd]
        sngpow_flg = sngpow_flg[gd]
        sngpow_xdata = sngpow_xdata[gd]
        civcand = civcand[gd]
     endif 
     ;; Set unsaturated absorber flags above saturation limits to
     ;; saturated
     bd = where(sngpow_data ge intlim[1] and (sngpow_flg and 2) lt 2,nbd)
     if nbd ne 0 and sngpow_datatype eq 'NCOLM' then begin
        print,'civ_maxsngpow: Number of absorbers that should be saturated',nbd
        sngpow_flg[bd] = sngpow_flg[bd] or 2 ; saturated
     endif 

     ;; Check that there's still enough data!
     nciv = n_elements(sngpow_data)
     if nciv le 1 then begin
        ;; Terrible!
        print,'civ_maxsngpow: insufficient data to fit',nciv
        print,'civ_maxsngpow: exiting without creating ',ostrct_fil
        return
     endif 

     ndim = 2                   ; dimensions of fit
     if not keyword_set(guess) then begin
        if keyword_set(infit_fil) then begin
           fitstrct = xmrdfits(infit_fil,1,/silent)
           guess = [alog(fitstrct.coeff/fitstrct.datanorm^fitstrct.alpha*$
                         sngpow_datanorm^fitstrct.alpha),$
                    fitstrct.alpha]
           if not keyword_set(conv_factor_lim) then $
              conv_factor_lim = 0.5*fitstrct.conv_factor
        endif else begin
           ;; Use prior knowledge to set these 
           if keyword_set(ncolm) then $
              guess = [alog(4.48e-14/(10.^13.5)^(-1.55)*$
                            sngpow_datanorm^(-1.5)), -1.5d] $
           else guess = [alog(1.08e-2/130.^(-1.63)*sngpow_datanorm^(-1.6)), $
                         -1.6d]
        endelse 
     endif 

     ;; Initialize
     ;; Frequency distribution:
     ;; f(data) = exp(px[0]) * (data/data_min)^px[1]
     ;; where data is either column density or equivalent width (linear).
     ;; Cast with exp(px[0]) to make gradient of maximum likelihood
     ;; function (ln L) with respect to px[0] on order of gradient with
     ;; respect to px[1].
     ;; Keep data within reasonable range of values by dividing by
     ;; "typical" (minimum, in this case) value. Especially useful for
     ;; f(N).
     if keyword_set(print_time) then $
        print,'civ_maxsngpow: starting minf_conj_grad ',systime()
     rslt = dblarr(10L,3)

     minf_conj_grad, guess, f_min, conv_factor, func_name='civ_maxsngpow_logL',$
                     /initialize,_extra=extra

     count = 0
     rslt[count,*] = [exp(guess[0]),guess[1],f_min]
     if not keyword_set(conv_factor_lim) then conv_factor_lim = 1.e-5
     while (conv_factor gt conv_factor_lim) do begin
        if (count mod 10) eq 0 and keyword_set(print_time) then $
           print,'civ_maxsngpow: minf_conj_grad loop = '+$
                 string(count,format='(i4)'),systime()
        minf_conj_grad, guess, f_min, conv_factor, $
                        func_name='civ_maxsngpow_logL',$
                        _extra=extra ; /use_deriv, tolerance=[sqrt(1.e-y)]
        
        count = count + 1
        rslt[count,*] = [exp(guess[0]),guess[1],f_min]
        ihi = n_elements(rslt[*,0])
        if ihi-1 eq count+1 then begin
           ;; Expand
           tmp = dblarr(2*ihi,3)
           tmp[0:ihi-1,0] = rslt[*,0]
           tmp[0:ihi-1,1] = rslt[*,1]
           tmp[0:ihi-1,2] = rslt[*,2]
           rslt = tmp
        endif 

     endwhile 

     ;; Set to store and use later
     sngpow_coeff = guess[0]
     sngpow_alpha = guess[1]

     ;; Save the normal coeff * N^alpha formalism
     fitstrct = {datatype:sngpow_datatype, $ ; store following to use as infit_fil
                 data:sngpow_data, $
                 sigdata:sngpow_sigdata, $
                 flag:sngpow_flg, $
                 xdata:sngpow_xdata, $
                 datanorm:sngpow_datanorm, $
                 intlim:sngpow_intlim, $ ; [lo, "saturated", hi] 
                 zrng:[min(civcand.zabs[ion_indx]),$
                       median(civcand.zabs[ion_indx]),$
                       max(civcand.zabs[ion_indx])], $ ; after trim
                 datalim:[min(sngpow_data,max=mx),mx], $
                 coeff: exp(sngpow_coeff),$ ; best-fit results and error
                 sigcoeff: dblarr(2), $
                 alpha: sngpow_alpha, $
                 sigalpha: dblarr(2), $
                 omciv: civ_maxsngpow_omciv(guess,siiv=siiv), $ 
                 sigomciv: dblarr(2), $ ; defined by following locations
                 sigcoeff_omciv:dblarr(2), $
                 sigalpha_omciv:dblarr(2), $
                 logL:-f_min, $ ; maximum likelihood
                 d_ks:0.d, $      ; K-S test results
                 prob_ks:0.d, $
                 conv_factor:conv_factor, $ ; convergence factor
                 rslt_iter:rslt $           ; info about iterations
                }

     ;; Write-out (safety)
     mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil
     mwrfits,civcand,ostrct_fil,/silent          ; sorted, trimmed

     ;; Set next steps
     ksonly = 1
     erronly = 1
  endelse                                        ; max-L


  if keyword_set(ksonly) then begin
     ;; KS Test
     ;; Normalize function:
     ;; px[0] = ln(1 + px[1]) - ln( (var_max^(1+px[1])/var_min^px[1]) -
     ;;                         var_min)
     ;; Only for unsaturated features
     sub = where(sngpow_data lt sngpow_intlim[1],nsub)
     sngpow_intlim[0] = min(sngpow_data)
     numer = (1. + sngpow_alpha) * sngpow_datanorm^sngpow_alpha
     denom = (sngpow_intlim[1]^(1.+sngpow_alpha)  - $
              sngpow_intlim[0]^(1.+sngpow_alpha) )
     sngpow_coeff = alog(numer / denom)
     ksone_wght, sngpow_data[sub], 'civ_maxsngpow_kscumf', d, ksprob, $
                 wght = 1./sngpow_xdata[sub], plot=plot, _extra=extra

     ;; Save results
     fitstrct.d_ks = d
     fitstrct.prob_ks = ksprob 

     ;; Restore common values 
     sngpow_intlim = fitstrct.intlim
     sngpow_coeff = alog(fitstrct.coeff)

     ;; Write-out again (safety)
     mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil
     mwrfits,civcand,ostrct_fil,/silent          ; sorted, trimmed

     if keyword_set(print_time) then $
        print,'civ_maxsngpow: finished KS Test ',systime()
  endif                         ; /ksonly


  if keyword_set(erronly) then begin
     ;; Errors
     ;; Grid of coeff and alpha
     if not keyword_set(nbin_alpha) then nbin_alpha = 200L
     if not keyword_set(nbin_coeff) then nbin_coeff = 200L
     if not keyword_set(dcoeff) then dcoeff = 0.01
     if not keyword_set(dalpha) then dalpha = 0.005

     ;; Ensure grids have max values "exactly"
     ;; Also know that likelihood turns over faster at high ceff 
     ;; and high alpha
     coeff_grid = sngpow_coeff + dcoeff * (lindgen(nbin_coeff) - nbin_coeff/2)
     alpha_grid = sngpow_alpha + dalpha * (lindgen(nbin_alpha) - nbin_alpha/2)

     ;; Set array of [coeff,alpha] to find likelihood function value
     px = dblarr(nbin_coeff*nbin_alpha,2)
     for ii=0L,nbin_alpha-1 do begin
        ;; Keep grouped by alpha (quicker integration)
        px[lindgen(nbin_coeff)+ii*nbin_coeff,0] = coeff_grid
        px[lindgen(nbin_coeff)+ii*nbin_coeff,1] = alpha_grid[ii]
     endfor  
     ;; Evaluate (could be expensive)
     if keyword_set(print_time) then $
        print,'civ_maxsngpow: calculating grid of log L ',systime()
     logL = -civ_maxsngpow_logL(px) ; want max

     if keyword_set(print_time) then $
        print,'civ_maxsngpow: starting calculating Om(C+3) ',systime()
     omciv = civ_maxsngpow_omciv(px,siiv=siiv,_extra=extra)
     ;; coeff by row; alpha by column 
     ;; coeff fixed when logL[ii,*] == x (in SURFACE)
     ;; alpha fixed when logL[*,jj] == y (in SURFACE)
     logL_grid = reform(logL,nbin_coeff,nbin_alpha)
     omciv_grid = reform(omciv,nbin_coeff,nbin_alpha)

     ;; Orient
     logL_mx = max(logL_grid,imx)
     nn = array_indices(logL_grid,imx)
     if logL_mx ne fitstrct.logL then $
        stop,'civ_maxsngpow: maximum likelihood values do not match',$
             fitstrct.logL,logL_mx
     icoeff_best = nn[0]
     ialpha_best = nn[1]

     ;; Find ellipse bounds (a la /asymerr in fuse_cog)
     dlogl_onesigma = -1.15      ; emprical limit 
     errsurf = where((logL_grid-logL_mx) ge dlogl_onesigma,nerrsurf)
     indx = lonarr(nerrsurf,2)
     for ii=0L,nerrsurf-1 do indx[ii,*] = array_indices(logL_grid,errsurf[ii])
     srt = sort(indx[*,1])      ; lump alpha
     indx[*,0] = indx[srt,0]
     indx[*,1] = indx[srt,1]

     imn = min(indx[*,0],max=imx) ; max extent of coeff
     jmn = min(indx[*,1],max=jmx) ; max extent of alpha

     ;; Check error ellipses closed
     if imn eq 0 or imx eq nbin_coeff-1 then $
        stop,'civ_maxsngpow: error ellipse not closed in k: ',imn,imx
     if imx-imn le 5 then $
        stop,'civ_maxsngpow: error ellipse too small in k'
     if jmn eq 0 or jmx eq nbin_alpha-1 then $
        stop,'civ_maxsngpow: error ellipse not closed in alpha: ',jmn,jmx
     if jmx-jmn le 5 then $
        stop,'civ_maxsngpow: error ellipse too small in alpha'

     ;; Store error results
     fitstrct.sigcoeff = [exp(coeff_grid[icoeff_best])-exp(coeff_grid[imn]),$
                          exp(coeff_grid[imx])-exp(coeff_grid[icoeff_best])]
     fitstrct.sigalpha = [alpha_grid[ialpha_best]-alpha_grid[jmn],$
                          alpha_grid[jmx]-alpha_grid[ialpha_best]]

     ;; Om(CIV) and error
     dum = min(omciv[errsurf],omn,subscript_max=omx)
     nn = array_indices(omciv_grid,errsurf[omn])
     fitstrct.sigcoeff_omciv[0] = fitstrct.coeff - exp(coeff_grid[nn[0]])
     fitstrct.sigalpha_omciv[0] = fitstrct.alpha - alpha_grid[nn[1]]
     fitstrct.sigomciv[0] = fitstrct.omciv - omciv_grid[nn[0],nn[1]]

     nn = array_indices(omciv_grid,errsurf[omx])
     fitstrct.sigcoeff_omciv[1] = exp(coeff_grid[nn[0]]) - fitstrct.coeff
     fitstrct.sigalpha_omciv[1] = alpha_grid[nn[1]] - fitstrct.alpha
     fitstrct.sigomciv[1] = omciv_grid[nn[0],nn[1]] - fitstrct.omciv 


     if keyword_set(plot) then begin
        ;; To plot:
        ;; surface,logL_grid,coeff_grid,alpha_grid,ax=50.,charsize=2
        contour,exp(logL_grid-logL_mx),coeff_grid,alpha_grid,levels=[exp(dlogl_onesigma)]
        oplot,[coeff_grid[icoeff_best]],[alpha_grid[ialpha_best]],psym=4
        oplot,[coeff_grid[imn],coeff_grid[imn]],$
              [alpha_grid[0],alpha_grid[nbin_alpha-1]]
        oplot,[coeff_grid[imx],coeff_grid[imx]],$
              [alpha_grid[0],alpha_grid[nbin_alpha-1]]
        oplot,[coeff_grid[0],coeff_grid[nbin_coeff-1]],$
              [alpha_grid[jmn],alpha_grid[jmn]]
        oplot,[coeff_grid[0],coeff_grid[nbin_coeff-1]],$
              [alpha_grid[jmx],alpha_grid[jmx]]
     endif 

     ;; Store info about error ellipse
     errstrct = {coeff_grid:exp(coeff_grid),alpha_grid:alpha_grid,$
                 logL:logL,omciv:omciv,surf:errsurf}
  endif                         ; /erronly


  ;; Write-out (for real)
  mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil
  mwrfits,civcand,ostrct_fil,/silent          ; sorted, trimmed
  mwrfits,errstrct,ostrct_fil,/silent 
  print,'civ_maxsngpow: created ',ostrct_fil
  if keyword_set(print_time) then $
     print,'civ_maxsngpow: finished ',systime()
  
end
