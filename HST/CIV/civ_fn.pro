;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_fn.pro               
; Author: Kathy Cooksey                      Date: 17 Sep 2008
; Project: HST CIV survey with Xavier Prochaska
; Description: 
; Input: 
; Optional:
; Output: 
; Example:
; History:
;   17 Sep 2008  created by KLC
;   11 Dec 2008  Use dX and the EW cut from g(z)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_calcewn                    ; need civ_calcewn_ndblt()
@civ_sensitivity                ; need civ_sensitivity_x()
pro civ_fn_plot,psfil, strct_fil, sens_fil, siiv=siiv, $
                fnstrct_fil=fnstrct_fil, psym=psym,pclr=pclr,$
                csize=csize,psize=psize,lthick=lthick,$
                xrng=xrng,yrng=yrng, log=log, few=few, $
                fit_fil=fit_fil, label=label, plt_intlim=plt_intlim,$
                plttitle=plttitle, fit_err=fit_err, _extra=extra
  if (N_params() lt 1) then begin
     print,'Syntax - '+$
           'civ_fn_plot,psfil, strct_fil, sens_fil, [/siiv, few=, fnstrct_fil=,'
     print,'            csize=, psize=,lthick=, label='
     print,'            xrng=, yrng=, /log, fit_fil=, /fit_err, _extra=]'
     return
  endif

  ;; Default
  if not keyword_set(plt_intlim) then plt_intlim = 1 ; saturation bound
  if plt_intlim lt 0 then plt_intlim = 0

  if keyword_set(fnstrct_fil) then begin
     ;; Read in f(N)
     if size(fnstrct_fil,/type) eq 8 then fnstrct = fnstrct_fil $
     else fnstrct = xmrdfits(fnstrct_fil,1,/silent)
  endif else begin
     ;; Create f(N)
     civ_fn, strct_fil, sens_fil, siiv=siiv, few=few,$
             ostrct_fil=fnstrct,_extra=extra ; for civ_group
  endelse 

  nbin = n_elements(fnstrct.locn)
  bdfn = where(fnstrct.fn le 0.,nbdfn,complement=gdfn)
  if keyword_set(log) then begin
     fnstrct.sigfn[gdfn,0] = fnstrct.sigfn[gdfn,0]/$
                             abs(alog(10)*fnstrct.fn[gdfn])
     fnstrct.sigfn[gdfn,1] = fnstrct.sigfn[gdfn,1]/$
                             abs(alog(10)*fnstrct.fn[gdfn])
     fnstrct.fn[gdfn] = alog10(fnstrct.fn[gdfn])
     if nbdfn ne 0 then begin
        ;; Set 2-sigma upper limit
        fnstrct.fn[bdfn] = alog10(2.*fnstrct.sigfn[bdfn,1])
        fnstrct.sigfn[bdfn,*] = 0.
     endif 
  endif else begin
     if nbdfn ne 0 then begin
        fnstrct.fn[bdfn] = 1.*fnstrct.sigfn[bdfn,1]
        fnstrct.sigfn[bdfn,*] = 0.
     endif 
  endelse 

  ;; Params
  if not keyword_set(psym) then psym = 4 ; diamond
  if not keyword_set(csize) then csize = 2.0
  if not keyword_set(psize) then psize = 2.
  if not keyword_set(lthick) then lthick = 5.0
  if not keyword_set(xrng) then $
     xrng = [fnstrct.locn[0]-fnstrct.siglocn[0,0]-0.1,$
             fnstrct.locn[nbin-1]+fnstrct.siglocn[nbin-1,1]+0.1]

  if not keyword_set(yrng) then begin
     yrng = [min(fnstrct.fn[gdfn])-max(fnstrct.sigfn[gdfn,0]),$
             max(fnstrct.fn[gdfn])+$
             max(fnstrct.sigfn[gdfn,1])] 
     if keyword_set(log) then yrng = yrng + [-0.02,0.02] $
     else yrng = yrng + [-1e-20,1e-15]
  endif

  if keyword_set(few) then begin
     if keyword_set(siiv) then begin
        if keyword_set(log) then ytitle = 'log !8f!X(!8W!X!D!8r!X,1393!N)'  $
        else ytitle = '!8f!X(!8W!X!D!8r!X,1393!N)'
        xtitle = 'log !8W!X!D!8r!X,1393!N (m'+STRING("305B) +')'
     endif else begin
        if keyword_set(log) then ytitle = 'log !8f!X(!8W!X!D!8r!X,1548!N)'  $
        else ytitle = '!8f!X(!8W!X!D!8r!X,1548!N)'
        xtitle = 'log !8W!X!D!8r!X,1548!N (m'+STRING("305B) +')'
     endelse 
  endif else begin
     if keyword_set(siiv) then begin
        if keyword_set(log) then ytitle = 'log !8f!X(!8N!X(Si!E+3!N))' $
        else ytitle = '!8f!X(!8N!X(Si!E+3!N))'
        xtitle = 'log !8N!X(Si!E+3!N)'
     endif else begin
        if keyword_set(log) then ytitle = 'log !8f!X(!8N!X(C!E+3!N))' $
        else ytitle = '!8f!X(!8N!X(C!E+3!N))'
        xtitle = 'log !8N!X(C!E+3!N)'
     endelse 
  endelse 

  x_psopen,psfil,/maxs,_extra=extra
  !p.multi = [1,1,1]
  !x.margin = [8,3]
  !y.margin = [4,2]

  clr = getcolor(/load)

  plot,xrng,yrng,/nodata,/ystyle,/xstyle,background=clr.white,color=clr.black,$
       ytitle=ytitle,xtitle=xtitle,charsize=csize


  ;; f(N)
  oploterror,fnstrct.locn[gdfn],fnstrct.fn[gdfn],fnstrct.siglocn[gdfn,0],$
             fnstrct.sigfn[gdfn,0],$
             errcolor=clr.black,psym=3,/lobar,color=clr.black,symsize=psize,$
             thick=lthick
  oploterror,fnstrct.locn[gdfn],fnstrct.fn[gdfn],fnstrct.siglocn[gdfn,1],$
             fnstrct.sigfn[gdfn,1],$
             errcolor=clr.black,psym=3,/hibar,color=clr.black,symsize=psize,$
             thick=lthick
  oplot,fnstrct.locn[gdfn],fnstrct.fn[gdfn],color=pclr,psym=psym,$
        symsize=psize
  if nbdfn ne 0 then begin
     plotsym,1,4,thick=2.5,color=clr.black ; down arrow, now psym=8
     oplot,fnstrct.locn[bdfn],fnstrct.fn[bdfn],$
           color=clr.black,psym=8 ; arrow
     oploterror,fnstrct.locn[bdfn],fnstrct.fn[bdfn],fnstrct.siglocn[bdfn,0],$
                replicate(0.,nbdfn),/lobar,errcolor=clr.black,$
                color=clr.black,psym=psym,symsize=psize,thick=lthick ; point
     oploterror,fnstrct.locn[bdfn],fnstrct.fn[bdfn],fnstrct.siglocn[bdfn,1],$
                replicate(0.,nbdfn),/hibar,errcolor=clr.black,$
                color=clr.black,psym=psym,symsize=psize,thick=lthick ; point
  endif 

  ;; Fit
  if keyword_set(fit_fil) then begin
     if size(fit_fil,/type) eq 7 then begin
        fitstrct = xmrdfits(fit_fil,1,/silent) 
        if keyword_set(fit_err) then errstrct = xmrdfits(fit_fil,3,/silent)
     endif else begin
        fitstrct = fit_fil
        if keyword_set(fit_err) then begin
           if size(fit_err,/type) ne 8 then $
              stop,'civ_fn_plot: must set fit_err to fit error structure' $
           else errstrct = fit_err
        endif 
     endelse 

     ;; Linear space first
     xfit = 10.^xrng[0] + dindgen(1000L) * (10.^xrng[1]-10.^xrng[0])/1000L
     yfit = fitstrct.coeff * (xfit/fitstrct.datanorm)^fitstrct.alpha
     if keyword_set(fit_err) then begin
        ;; Potentially time-consuming if do all 1000L
        xerr = [xfit[lindgen(10)*100L],xfit[999L]]
        civ_fn_fit, rslt, fitstrct, in_array=xerr, errstrct=errstrct, log=log
;        ;; Test strict error propagation (too small! b/c now covar term?)
;        yerrlo = sqrt(fitstrct.sigcoeff[0]^2 * $
;                      (xfit/fitstrct.datanorm)^(2.*fitstrct.alpha) + $
;                      (fitstrct.coeff*fitstrct.alpha*$
;                       (xfit/fitstrct.datanorm)^(fitstrct.alpha-1))^2 * $
;                      fitstrct.sigalpha[0]^2)
;        yerrhi = sqrt(fitstrct.sigcoeff[1]^2 * $
;                      (xfit/fitstrct.datanorm)^(2.*fitstrct.alpha) + $
;                      (fitstrct.coeff*fitstrct.alpha*$
;                       (xfit/fitstrct.datanorm)^(fitstrct.alpha-1))^2 * $
;                      fitstrct.sigalpha[1]^2)
     endif 
     
     if keyword_set(log) then begin
        oplot,alog10(xfit),alog10(yfit),$
              linestyle=0,color=clr.red,thick=lthick 
        oplot,alog10([fitstrct.intlim[plt_intlim],$
                      fitstrct.intlim[plt_intlim]]), yrng,$
              linestyle=1,color=clr.black,thick=lthick
        if keyword_set(fit_err) then begin
           ;; log(N+sigN) = log(N + sigLogN * ln(10) * N)
           ;;             = logN + log(1+sigLogN*ln(10))
           ;; which means "1 sigma" in log space is +/-log(1+/-sigLogN*ln(10))
           oplot,alog10(xerr),rslt[*,0]-abs(alog10(1.-rslt[*,1]*alog(10.))),$
                 linestyle=2,color=clr.red ; dashed
           oplot,alog10(xerr),rslt[*,0]+alog10(1.+rslt[*,1]*alog(10.)),$
                 linestyle=2,color=clr.red ; dashed
;           ;; Test
;           oplot,alog10(xfit),alog10(yfit)-abs(alog10(1.-yerrlo*alog(10.))),$
;                 linestyle=2,color=clr.slategray
;           oplot,alog10(xfit),alog10(yfit)+alog10(1.+yerrhi*alog(10.)),$
;                 linestyle=2,color=clr.slategray
;           stop
        endif 

     endif else begin
        oplot,alog10(xfit),yfit,linestyle=0,color=clr.red,thick=lthick
        oplot,[fitstrct.intlim[plt_intlim],$
               fitstrct.intlim[plt_intlim]], yrng,$
              linestyle=1,color=clr.black,thick=lthick
        if keyword_set(fit_err) then begin
           oplot,xerr,rslt[*,0]-rslt[*,1],linestyle=2,color=clr.red ; dashed
           oplot,xerr,rslt[*,0]+rslt[*,2],linestyle=2,color=clr.red ; dashed
        endif 
     endelse                    ; log or not

     ;; Legend (upper right)
     if keyword_set(label) then begin
        dx = xrng[1] - xrng[0]
        dy = yrng[1] - yrng[0]
        if keyword_set(few) then begin
           xyouts, 0.97*dx + xrng[0], 0.95*dy + yrng[0], $
                   '!9a!X!D!8W!X!N = '+$
                   strtrim(string(fitstrct.alpha,format='(f5.2)'),2),$
                   charsize=csize, color=clr.black,alignment=1
           xyouts, 0.97*dx + xrng[0], 0.9*dy + yrng[0], $
                   '!8k!X!D3!N = '+strtrim(string(fitstrct.coeff*1.e3,$
                                                  format='(f5.2)'),2) + $
                   ' m'+STRING("305B)+'!E-1!N',$
                   charsize=csize, color=clr.black,alignment=1
        endif else begin
           xyouts, 0.97*dx + xrng[0], 0.95*dy + yrng[0], $
                   '!9a!X!D!8N!X!N = '+$
                   strtrim(string(fitstrct.alpha,format='(f5.2)'),2),$
                   charsize=csize, color=clr.black,alignment=1
           xyouts, 0.97*dx + xrng[0], 0.9*dy + yrng[0], $
                   '!8k!X!D14!N = '+strtrim(string(fitstrct.coeff*1.e14,$
                                                   format='(f5.2)'),2) + $
                   ' cm!E2!N',$
                   charsize=csize, color=clr.black,alignment=1
        endelse 
        if fitstrct.prob_ks le 1.e-3 then $
           xyouts, 0.97*dx + xrng[0], 0.85*dy + yrng[0], $
                   'P!DKS!N = '+strtrim(string(fitstrct.prob_ks*1.e3,$
                                               format='(f5.3)'),2)+$
                   "x10!E-3!N",charsize=csize, color=clr.black,alignment=1 $
        else $
           xyouts, 0.97*dx + xrng[0], 0.85*dy + yrng[0], $
                   'P!DKS!N = '+strtrim(string(fitstrct.prob_ks,$
                                               format='(f5.3)'),2),$
                   charsize=csize, color=clr.black,alignment=1 
     endif                      ; /label
  endif                         ; fit_fil=

  if keyword_set(plttitle) then begin
     ;; Upper right or lower left
     if keyword_set(label) then $
        xyouts,0.2,0.2,plttitle,color=clr.black,/normal,charsize=csize,$
               alignment=0 $    ; align left side
     else $
        xyouts,0.9,0.88,plttitle,color=clr.black,/normal,charsize=csize,$
               alignment=1.     ; align right side
  endif

  x_psclose
  print,'civ_fn_plot: created ',psfil

end                             ; civ_fn_plot


pro civ_fn_plotcdf,psfil,strct_fil,sens_fil, siiv=siiv, few=few, $
                   csize=csize,psize=psize,lthick=lthick,$
                   xrng=xrng,yrng=yrng, fit_fil=fit_fil, label=label, $
                   redshift=redshift,_extra=extra
  if (N_params() lt 3) then begin
     print,'Syntax - '+$
           'civ_fn_plot,psfil,strct_fil,sens_fil,few=,csize=,lthick=,'+$
           '          xrng=,yrng=, fit_fil=, /label'
     return
  endif

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  if size(sens_fil,/type) eq 7 then sens = xmrdfits(sens_fil,1,/silent) $
  else sens = sens_fil  
  civ_group,strct,strct_fil,nciv=nciv,_extra=extra ;flg_sys=,rating=,zlim=

  if keyword_set(few) then begin
     ncolm = alog10(strct.ew[0])
     signcolm = strct.sigew[0]/abs(alog(10.)*strct.ew[0])
     xpath = civ_sensitivity_x(sens,ew=10.^ncolm,z=zpath) ; 95%, interpolate
     if keyword_set(siiv) then begin
        ytitle = '!8f!X(!8W!X!D!8r!X,1393!N) CDF'
        xtitle = 'log !8W!X!D!8r!X,1393!N (m'+STRING("305B) +')'
     endif else begin
        ytitle = '!8f!X(!8W!X!D!8r!X,1548!N) CDF'
        xtitle = 'log !8W!X!D!8r!X,1548!N (m'+STRING("305B) +')'
     endelse 
  endif else begin
     ;; Error-weighted column density
     ncolm = civ_calcewn_ndblt(strct,dblt_name,signcolm=signcolm,/log,/silent)
     xpath = civ_sensitivity_x(sens,ncolm=ncolm,z=zpath) ; 95%, interpolate
     if keyword_set(siiv) then begin
        ytitle = '!8f!X(!8N!X(Si!E+3!N)) CDF'
        xtitle = 'log !8N!X(Si!E+3!N)'     
     endif else begin
        ytitle = '!8f!X(!8N!X(C!E+3!N)) CDF'
        xtitle = 'log !8N!X(C!E+3!N)'     
     endelse 
  endelse 
  if keyword_set(redshift) then xpath = zpath

  srt = sort(ncolm)             ; order
  ncolm = ncolm[srt]
  strct = strct[srt]
  xpath = xpath[srt]
  nmin = min(ncolm,max=nmax)
  
  ;; CDF
  cdf = total(1./xpath,/cumulative)/total(1./xpath)

  ;; Params
  if not keyword_set(csize) then csize = 1.5
  if not keyword_set(psize) then psize = 2.
  if not keyword_set(lthick) then lthick = 2
  if not keyword_set(xrng) then $
     xrng = [ncolm[0]-signcolm[0]-0.1,$
             ncolm[nciv-1]+signcolm[nciv-1]+0.1]
  if not keyword_set(yrng) then yrng = [0.,1.01] 

  ;; Plot
  x_psopen,psfil,/maxs,_extra=extra
  !p.multi = [1,1,1]
  !x.margin = [8,3]
  !y.margin = [4,2]

  clr = getcolor(/load)

  plot,xrng,yrng,/nodata,/ystyle,/xstyle,background=clr.white,color=clr.black,$
       ytitle=ytitle,xtitle=xtitle,charsize=csize

  ;; CDF
  oplot,[ncolm[0],ncolm],[0,cdf],color=clr.black,psym=10.,thick=lthick

  ;; Fit (normalize to intlim[0] to infinity)
  if keyword_set(fit_fil) then begin
     if size(fit_fil,/type) eq 7 then $
        fitstrct = xmrdfits(fit_fil,1,/silent) $
     else fitstrct = fit_fil

     ;; Linear space first
     xfit = 10^xrng[0] + dindgen(1000L) * (10^xrng[1]-10^xrng[0])/1000L
     yfit = (xfit^(1.+fitstrct.alpha) - $
             fitstrct.datalim[0]^(1.+fitstrct.alpha))/$
            (fitstrct.datalim[1]^(1.+fitstrct.alpha) - $
             fitstrct.datalim[0]^(1.+fitstrct.alpha) )
     oplot,alog10(xfit),yfit, linestyle=0,color=clr.red,thick=lthick 
     ;; Itegration limit
     oplot,alog10([fitstrct.intlim[1],fitstrct.intlim[1]]), yrng,$
           linestyle=1,color=clr.black,thick=lthick

     ;; Legend (upper left)
     if keyword_set(label) then begin
        dx = xrng[1] - xrng[0]
        dy = yrng[1] - yrng[0]
        xyouts, 0.05*dx + xrng[0], 0.95*dy + yrng[0], $
                '!9a!X = '+strtrim(string(fitstrct.alpha,format='(f5.2)'),2),$
                charsize=csize, color=clr.black
        xyouts, 0.05*dx + xrng[0], 0.9*dy + yrng[0], $
                'P!DKS!N = '+strtrim(string(fitstrct.prob_ks,$
                                            format='(f5.3)'),2),$
                charsize=csize, color=clr.black
     endif                      ; /label
  endif                         ; fit_fil=

  x_psclose
  print,'civ_fn_plotcdf: created ',psfil

end                             ; civ_fn_plotcdf


pro civ_fn_fit, rslt, fit_fil, in_array=in_array, $
                errstrct=errstrct, sigcoeff_fn=sigcoeff_fn,$
                sigalpha_fn=sigalpha_fn,log=log
  ;; Calculate the error on the f(N) or f(EW) fits given the column
  ;; density or EW values (as through in_array=). Default is to
  ;; calculate for N_0 or EW_0 as stored in the fit structure.
  ;; Modeled after civ_dndx

  if N_params() lt 2 then begin
     print,'Syntax - civ_fn_fit, rslt, fit_fil, [in_array=,'
     print,'                     errstrct=, sigcoeff_fn=,sigalpha_fn=,/log]'
     print,'         '
  endif 

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'
  
  ;; Read in and check fit and error structures
  if size(fit_fil,/type) eq 7 then begin
     fitstrct = xmrdfits(fit_fil,1,/silent)
     if not keyword_set(errstrct) then $
        errstrct = xmrdfits(fit_fil,3,/silent)
  endif else begin
     fitstrct = fit_fil
     if not keyword_set(errstrct) then $
        stop,'civ_Fn_fit: must input error structure for frequency distribution fit'
  endelse 

  ;; Set value(s) for which to estimate error
  if not keyword_set(in_array) then in_array = fitstrct.datanorm

  ;; Prepare output
  nbin = n_elements(in_array)
  rslt = dblarr(nbin,3)
  rslt[*,0] = fitstrct.coeff * (in_array / fitstrct.datanorm)^ fitstrct.alpha
  sigcoeff_fn = dblarr(nbin,2)
  sigalpha_fn = dblarr(nbin,2)
  
  freqdistr = dblarr(n_elements(errstrct.surf)) ; minimize computation
  alpha_grid = freqdistr
  coeff_grid = freqdistr

  for ii=0,nbin-1 do begin
     ;; Reset
     freqdistr[*] = 0.
     alpha_grid[*] = 0.
     coeff_grid[*] = 0.

     logl_grid = reform(errstrct.logL,n_elements(errstrct.coeff_grid),$
                        n_elements(errstrct.alpha_grid)) ; for reference

     for nn=0L,n_elements(errstrct.surf)-1 do begin
        mm = array_indices(logl_grid,errstrct.surf[nn])
        alpha_grid[nn] = errstrct.alpha_grid[mm[1]]
        coeff_grid[nn] = errstrct.coeff_grid[mm[0]]
     endfor                     ; loop nn=n_elements(errstrct.surf)
     
     ;; f(N) = coeff * (N/N_0) ^ alpha
     freqdistr = coeff_grid * (in_array[ii] / fitstrct.datanorm)^alpha_grid
     
     ;; Find extrema of error ellipse values
     mn = min(freqdistr,imn,max=mx,subscript_max=imx)

     rslt[ii,1] = rslt[ii,0] - mn
     rslt[ii,2] = mx - rslt[ii,0]

     ;; Return locations of these limits
     sigcoeff_fn[ii,*] = [fitstrct.coeff - coeff_grid[imn],$
                         coeff_grid[imx] - fitstrct.coeff]
     sigalpha_fn[ii,*] = [fitstrct.alpha - alpha_grid[imn],$
                         alpha_grid[imx] - fitstrct.alpha]
  endfor                        ; loop ii=nbin

  ;; Convert
  if keyword_set(log) then begin
     rslt[*,1] = rslt[*,1]/(alog(10.) * rslt[*,0])
     rslt[*,2] = rslt[*,2]/(alog(10.) * rslt[*,0])
     rslt[*,0] = alog10(rslt[*,0])
  endif 
end                             ; civ_fn_fit


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_fn, strct_fil, sens_fil, siiv=siiv, nmax=nmax, $
            few=few, psfil=psfil, binsize=binsize, hdf5_fil=hdf5_fil, $
            nlim=nlim, ostrct_fil=ostrct_fil, redshift=redshift,_extra=extra
  if not keyword_set(nmax) then nmax = 16.1 ; defines saturated bin

  if not keyword_set(binsize) then binsize = 0.2
  if not keyword_set(nlim) then begin
     if keyword_set(few) then nlim = [alog10(20.),alog10(3.e4)] $ ; log(mA)
     else nlim = [13.,nmax]                                        ; 14.3
  endif 

  if keyword_set(siiv) then dblt_name = 'SiIV' $
  else dblt_name = 'CIV'

  ;; Trim strct based on both features being 3 sigma
  civ_group,strct,strct_fil,flg_sys=384,nciv=nstrct,$
            _extra=extra        ; zlim=,rating=,/unsat
  if nstrct eq 0 then stop,'civ_fn: no '+dblt_name+$
                           ' with both lines detected at 3 sigma'

  ;; g(z) (unblocked path = sens.zpath; convert to dX)
  if keyword_set(sens_fil) then begin
     if size(sens_fil,/type) eq 7 then sens = xmrdfits(sens_fil,1,/silent) $
     else sens = sens_fil
  endif else begin
     if keyword_set(hdf5_fil) then begin
        if size(hdf5_fil,/type) eq 7 then $ 
           hdf5_sensitivity,hdf5_fil,dblt_name,sens,_extra=extra $ ; cosmology=, zlim=
        else hdf5_sensitivity,strct[0].instr_fil,dblt_name,sens,_extra=extra 
     endif else stop,'civ_fn: sens_fil not set'
;     print,max(sens.cumx_ncolm[*,1])
;     stop
  endelse 

  ;; Determine index
  civ = where(stregex(strct[0].ion,dblt_name,/boolean),nciv)
  if nciv ne 2 then stop,'civ_fn: not a doublet'
  civ = civ[sort(strct[0].wrest[civ])]

  if keyword_set(few) then begin
     ncolm = alog10(strct.ew[civ[0]])
  endif else begin
     ;; Error-weighted average of AODM column density
     ncolm = civ_calcewn_ndblt(strct,dblt_name,flg_colm=flg_colm,/log,/silent)
  endelse    
  if n_elements(ncolm) eq 1 then begin
     ;; Make all bins
     num = ceil((nlim[1]-nlim[0])/binsize)
     stop
  endif else begin
     histn = histogram(ncolm,loc=locn,binsize=binsize,min=nlim[0])
     locn = locn + 0.5*binsize
  endelse 
  nbin = n_elements(histn)
  siglocn = fltarr(nbin,2)
  siglocn[*,0] = 0.5*binsize
  siglocn[*,1] = 0.5*binsize

  ;; Saturated measurement bin
  gd = where(locn lt nlim[1],ngd,complement=bd,ncomplement=nbd)
  if nbd ne 0 then begin
     ;; Re-sample histogram; but have to have defined saturated bin
     tmp = histn[0:ngd]
     test = where(ncolm ge nlim[1] and ncolm le nmax,ntest)
     tmp[ngd] = ntest           ;total(histn[bd])
     histn = tmp
     
     ;; Reset bin locations 
     cent = 0.5*(nmax+nlim[1])  ;mean(locn[bd])
     tmp = fltarr(ngd+1,2)
     tmp[*,0] = siglocn[0:ngd,0]
     tmp[*,1] = siglocn[0:ngd,1]
     tmp[ngd,1] = nmax - cent ;locn[nbin-1]+siglocn[nbin-1,1] - cent
     tmp[ngd,0] = cent - nlim[1] ;(locn[bd[0]]-siglocn[bd[0],0])
     siglocn = tmp
     tmp = locn[0:ngd]
     tmp[ngd] = cent
     locn = tmp

     ;; New bin number
     nbin = ngd + 1
  endif
  locnhi = locn + siglocn[*,1]
  locnlo = locn - siglocn[*,0]

  ;; Instantiate structures
  fn = fltarr(nbin)
  sigfn = fltarr(nbin,2)
  sigdX = dblarr(nbin,2)
  sigdN = fltarr(nbin,2)        ; dN = histn
  if keyword_set(few) then begin
     dX = civ_sensitivity_x(sens,sigx=sigdX,ew=10.^locn,$
                            sigew=siglocn[*,0]*alog(10.)*10.^locn,$
                           z=dz,sigz=sigdz)    ; 95% limit, centered
  endif else begin
     dX = civ_sensitivity_x(sens,sigx=sigdX,ncolm=locn,$
                            sign=siglocn[*,0],$
                           z=dz,sigz=sigdz)    ; 95% limit, centered
  endelse 
  if keyword_set(redshift) then begin
     ;; Define f(N) as number / (dz * dN(ion))
     dX = dz
     sigdX = sigdz
  endif 

  for ii=0,nbin-1 do begin
     ;; f(N)
     denom = (10^locnhi[ii]-10^locnlo[ii])*dX[ii]
     if denom ne 0. then begin
        fn[ii] = histn[ii]/denom
        p = x_poisscl(histn[ii],0.683,/silent) 
        ;; Should round number?
        sigdN[ii,0] = histn[ii]-p[1] ; low
        sigdN[ii,1] = p[0]-histn[ii] ; hi
        ;; Errors add in quadrature .r 
        sigfn[ii,0] = sqrt((sigdN[ii,0]/denom)^2 + $
                           (sigdX[ii,1]*fn[ii]/dX[ii])^2)
        sigfn[ii,1] = sqrt((sigdN[ii,1]/denom)^2 + $
                           (sigdX[ii,0]*fn[ii]/dX[ii])^2)
     endif 
  endfor                        ; loop bins

  fnstrct = {fn:fn, sigfn:sigfn, locn:locn, siglocn:siglocn, $
             dX:dX, sigdX:sigdX, dN:histn, sigdN:sigdN}

  if keyword_set(ostrct_fil) then begin
     if size(ostrct_fil,/type) eq 7 then mwrfits,fnstrct,ostrct_fil,/create $
     else ostrct_fil = fnstrct 
  endif else ostrct_fil = fnstrct

  if keyword_set(psfil) then begin
     civ_fn_plot,psfil,siiv=siiv,fnstrct_fil=fnstrct,few=few,_extra=extra ; fit_fil=, /label
  endif                                                           ; psfil=

end
