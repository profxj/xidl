;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_fitfxw.pro
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
;    /dz -- fit for f(N) = Dnumber / (Dcolumn * Dz)
; Output: 
;    ostrct_fil -- name of output fit structure
; Example:
; History:
;    31 Dec 2008 -- created by KLC
;    27 Jan 2009 -- add effect of saturated features
;    15 Dec 2011 -- adopted from civ_maxsngpow, KLC
;    25 Jan 2012 -- Include Schechter function, KLC
;    05 Oct 2012 -- Working to fitting f(N), KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sdss_fitfxw_calcfiterr, fit_fil, cl, sigma=sigma, errstrct=errstrct,$
                                 delta=delta, mindlogl=mindlogl, logl=logl,$
                                 maxdlogl=maxdlogl, alpha_grid=alpha_grid,$
                                 coeff_grid=coeff_grid, norm_grid=norm_grid,$
                                 nalpha_bin=nalpha_bin,$
                                 ncoeff_bin=ncoeff_bin, nnorm_bin=nnorm_bin,$
                                 sigalpha=sigalpha,sigcoeff=sigcoeff, $
                                 signorm=signorm,ofitstr=ofitstr,oerrstr=oerrstr,$
                                 clobber=clobber,debug=debug

  if size(fit_fil,/type) eq 7 then begin
     fitstrct = xmrdfits(fit_fil,1,/silent)
     errstrct = xmrdfits(fit_fil,2,/silent)
  endif else $
     fitstrct = fit_fil
  
  if keyword_set(errstrct) then begin
     if size(errstrct,/type) eq 8 then begin
        alpha_grid = errstrct.alpha_grid
        coeff_grid = errstrct.coeff_grid
        norm_grid = errstrct.norm_grid
        ncoeff_bin = (size(coeff_grid,/dim))[0]
        nalpha_bin = (size(alpha_grid,/dim))[0]
        nnorm_bin = (size(norm_grid,/dim))[0]
        logl = errstrct.logl
     endif                      ; is structure
  endif                         ; errstrct=
  if not keyword_set(nnorm_bin) then begin
     norm_grid = 0
     nnorm_bin = 0
  endif

  if keyword_set(nnorm_bin) then $
     logl_grid = reform(logl,ncoeff_bin,nalpha_bin,nnorm_bin) $
  else logl_grid = reform(logl,ncoeff_bin,nalpha_bin)
  
  if not keyword_set(logl) then $
     stop,'sdss_fitfxw_calcfiterr(): must have logl= set somehow'

  logl_mx = max(logl_grid,imx)
  nn = array_indices(logL_grid,imx)
  if keyword_set(fitstrct) then begin
     if logl_mx ne fitstrct.logl then $
        print,'sdss_fitfxw_calcfiterr(): maximum likelihood values do not match',$
              fitstrct.logL,logL_mx $
     else begin
        ;; Keep the max on the grid
        gd = where(logl_grid eq fitstrct.logl)
        nn = array_indices(logl_grid,gd[0])
     endelse 
  endif 
  icoeff_best = nn[0]
  ialpha_best = nn[1]
  if keyword_set(nnorm_bin) then inorm_best = nn[2]

  dlogl = logL_grid-logL_mx
  dlogl_tot = total(exp(dlogl))
  npix_tot = float(n_elements(dlogl))
  
  ;; Parameters
  if not keyword_set(delta) then delta = 0.05
  if not keyword_set(mindlogl) then mindlogl = -1. ; empirical = -1.15 for 1sig
  if not keyword_set(maxdlogl) then maxdlogl = 0.5*max(dlogl) ; 50% down

  if not keyword_set(cl) and not keyword_set(sigma) then sigma = 1.
  if keyword_set(sigma) then begin
     if keyword_set(cl) then $
        print,'sdss_fitfxw_calcfiterr(); defaulting to sigma= value, not cl'
     cl = 1. - 2*(1-gauss_pdf(sigma)) ; area enclosed
  endif 

  nloop = ceil((maxdlogl-mindlogl)/delta)
  dlogl_limit = mindlogl - delta*findgen(nloop) ; emprical limit 
  
  ;; Actually a volume
  area = fltarr(nloop,2,/nozero)
  sigcoeff = dblarr(nloop,2,/nozero)
  sigalpha = dblarr(nloop,2,/nozero)
  if keyword_set(nnorm_bin) then signorm = dblarr(nloop,2,/nozero)

  for ii=0,nloop-1 do begin
     ;; Find ellipse bounds (a la /asymerr in fuse_cog)
     errsurf = where(dlogl ge dlogl_limit[ii],nerrsurf)
     
     if errsurf[0] eq -1 then begin
        area[ii,*] = 0.         ; must instantiate
        continue 
     endif 
     area[ii,0] = total(exp(dlogl[errsurf]))/dlogl_tot ; volume
     area[ii,1] = nerrsurf / npix_tot

     icoeff = errsurf mod ncoeff_bin
     ialpha = (errsurf - icoeff) / ncoeff_bin

     imn = min(icoeff,max=imx) ; max extent of coeff
     jmn = min(ialpha,max=jmx) ; max extent of alpha

     ;; Check error ellipses closed
     if imn eq 0 or imx eq ncoeff_bin-1 then $
        print,'sdss_fitfxw_calcfiterr(): WARNING! error ellipse not closed in k: ',$
             imn,imx
     if imx-imn le 5 then $
        print,'sdss_fitfxw_calcfiterr(): error ellipse too small in k'
     if jmn eq 0 or jmx eq nalpha_bin-1 then $
        print,'sdss_fitfxw_calcfiterr(): WARNING! error ellipse not closed in alpha: ',$
             jmn,jmx
     if jmx-jmn le 5 then $
        print,'sdss_fitfxw_calcfiterr(): error ellipse too small in alpha'

     ;; Store error results
     sigcoeff[ii,*] = [coeff_grid[icoeff_best]-coeff_grid[imn],$
                       coeff_grid[imx]-coeff_grid[icoeff_best]]
     sigalpha[ii,*] = [alpha_grid[ialpha_best]-alpha_grid[jmn],$
                       alpha_grid[jmx]-alpha_grid[ialpha_best]]

     if keyword_set(nnorm_bin) then begin
        inorm = (errsurf - ialpha) / nalpha_bin
        kmn = min(inorm,max=kmx)
        if kmn eq 0 or kmx eq nalpha_bin-1 then $
           stop,'sdss_fitfxw_calcfiterr(): error ellipse not closed in alpha: ',$
                kmn,kmx
        if kmx-kmn le 5 then $
           stop,'sdss_fitfxw_calcfiterr(): error ellipse too small in alpha'
        signorm[ii,*] = [norm_grid[inorm_best]-norm_grid[jmn],$
                         norm_grid[jmx]-norm_grid[inorm_best]]
     endif 

  endfor                        ; loop ii=nloop

  ;; Interpolate the best final value
  dlogl_cl_lim = interpol(dlogl_limit,area[*,0],cl)
  frac_cl_lim = interpol(area[*,1],dlogl_limit,dlogl_cl_lim)
  if frac_cl_lim gt 0.5 then $
     stop,'sdss_fitfxw_calcfiterr(): error ellipse greater than 50% of grid',$
          frac_cl_lim
  if keyword_set(debug) then begin
     print,cl,dlogl_cl_lim,frac_cl_lim,format='("cl=",f6.4,1x,"dlog=",f7.4,1x,"frac=",f6.4)'
     
     if keyword_set(nnorm_bin) then begin
        print,'dLogL','cl','Frac','klo','khi','alo','ahi','nlo','nhi',$
              format='(a7,2x,a6,7(2x,a8))'
        printcol,dlogl_limit,area[*,0],area[*,1],sigcoeff[*,0],sigcoeff[*,1],$
                 sigalpha[*,0],sigalpha[*,1],signorm[*,0],signorm[*,1],$
                 format='(f7.4,2x,f6.4,7(2x,f8.5))'
     endif else begin
        print,'dLogL','cl','Frac','klo','khi','alo','ahi',$
              format='(a7,2x,a6,5(2x,a8))'
        printcol,dlogl_limit,area[*,0],area[*,1],sigcoeff[*,0],sigcoeff[*,1],$
                 sigalpha[*,0],sigalpha[*,1],$
                 format='(f7.4,2x,f6.4,5(2x,f8.5))'
     endelse 
  endif                         ; /debug

  ;; Same with errors
  siglo = interpol(sigcoeff[*,0],dlogl_limit,dlogl_cl_lim)
  sighi = interpol(sigcoeff[*,1],dlogl_limit,dlogl_cl_lim)
  sigcoeff = [siglo,sighi]
  siglo = interpol(sigalpha[*,0],dlogl_limit,dlogl_cl_lim)
  sighi = interpol(sigalpha[*,1],dlogl_limit,dlogl_cl_lim)
  sigalpha = [siglo,sighi]
  if keyword_set(nnorm_bin) then begin
     siglo = interpol(signorm[*,0],dlogl_limit,dlogl_cl_lim)
     sighi = interpol(signorm[*,1],dlogl_limit,dlogl_cl_lim)
     signorm = [siglo,sighi]
  endif else begin
     signorm = replicate(!values.d_nan,2)
  endelse 

  if keyword_set(fitstrct) then begin
     ;; Save
     if keyword_set(debug) then begin
        print,fitstrct.sigcoeff[0],fitstrct.sigcoeff[1],$
              fitstrct.sigalpha[0],fitstrct.sigalpha[1],$
              fitstrct.sigdatanorm[0],fitstrct.sigdatanorm[1],$
              format='("Old errors: ",6(f8.4,1x))'
        print,sigcoeff[0],sigcoeff[1],sigalpha[0],sigalpha[1],$
              signorm[0],signorm[1],format='("New errors: ",6(f8.4,1x))'
     endif 
     fitstrct.sigalpha = sigalpha
     fitstrct.sigcoeff = sigcoeff
     fitstrct.sigdatanorm = signorm
  endif 

  ;; Will just overwrite the error structure
  errstrct = {coeff_grid:coeff_grid,alpha_grid:alpha_grid,norm_grid:norm_grid,$
              dlogl_lim:dlogl_cl_lim,cl:cl,logL:logL,surf:errsurf}
  
  if size(fit_fil,/type) eq 7 and keyword_set(clobber) then begin
     ;; Write-out (duplicate of sdss_fitfxw call)
     mwrfits,fitstrct,fit_fil,/create,/silent
     mwrfits,errstrct,fit_fil,/silent 
     spawn,'gzip -f '+fit_fil
     print,'sdss_fitfxw_calcfiterr(): created ',fit_fil
  endif 

  ofitstr = fitstrct
  oerrstr = errstrct

  return, dlogl_cl_lim
end                             ; sdss_fitfxw_calcfiterr()


pro sdss_fitfxw_prnt, fitstr_lst, list=list
  ;; print out a formatted table of useful information
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fitfxw_prnt, fitstr_lst, [/list]'
     return
  endif 

  if keyword_set(list) then begin
     readcol,fitstr_lst,fit_fil,format='a'
  endif else fit_fil = fitstr_lst
  nfil = (size(fit_fil,/dim))[0] > 1 ; foil singularity

  
  print,'Data','Nsys',$
        '<z>','zmin','zmax',$
        'minDat','maxDat','maxLim',$
        'k','Errlo','Errhi',$
        'alpha','Errlo','Errhi',$
        'dlogL','logLmax',$
        'D_KS','P_KS',$
        'Chi^2','Prob',$
        format='(a30,2x,a6,2x,3(a6,1x),1x,3(a6,1x),1x,3(a7,1x),1x,3(a7,1x),1x,a6,1x,a9,2x,a6,1x,a8,2x,a6,1x,a8)'
  for ff=0,nfil-1 do begin
     fitstr = xmrdfits(fit_fil[ff],1,/silent)
     errstr = xmrdfits(fit_fil[ff],2,/silent)

     prs = strsplit(fit_fil[ff],'/',/extract,count=nprs)
     fil = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'.',/reverse_search))
     
     print,fil,(size(fitstr.data,/dim))[0],$
           fitstr.zrng[1],fitstr.zrng[0],fitstr.zrng[2],$
           fitstr.datalim[0],fitstr.datalim[1],fitstr.intlim[1],$
           fitstr.coeff,fitstr.sigcoeff[0],fitstr.sigcoeff[1],$
           fitstr.alpha,fitstr.sigalpha[0],fitstr.sigalpha[1],$
           errstr.dlogl_lim,fitstr.logl,$
           fitstr.d_ks,fitstr.prob_ks,$
           fitstr.chi_sqr[0]/fitstr.chi_sqr[1],fitstr.chi_sqr[2],$
           format='(a30,2x,i6,2x,3(f6.4,1x),1x,3(f6.2,1x),1x,3(f7.4,1x),1x,3(f7.4,1x),1x,f6.3,1x,f9.2,2x,f6.2,1x,e8.2,2x,f6.2,1x,e8.2)'
  endfor 

end                             ; sdss_fitfxw_prnt


pro sdss_fitfxw_prntkstwo, kstwo_fil, out_fil=out_fil
  ;; Print carefully formmated file
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fitfxw_prntkstwo, kstwo_fil, [out_fil=]'
     return
  endif 
  
  if size(kstwo_fil,/type) eq 8 then $
     ks2str = kstwo_fil $
  else ks2str = xmrdfits(kstwo_fil,1,/silent)

  rowend = uniq(ks2str.file1)      ; should be sorted
  nrow = (size(rowend,/dim))[0] > 1


  if keyword_set(out_fil) then begin
     openw,1,out_fil
     if size(kstwo_fil,/type) eq 7 then $
        printf,1,kstwo_fil
  endif
  sfmt = '(4x,a21,$)'           ; header
  fmt0 = '(4x,f7.5,1x,i4,1x,f8.6,$)' ; data if Prob >= 1e-5 (original)
  fmt0_bd = '(4x,f7.5,1x,i4,1x,f7.5,"*",$)' ; if Prob < 0.05, consider bad
  afmt = '(4x,f7.5,1x,i4,1x,e7.1,"*",$)' ; data if prob < 1e-5 (alternate)

  for rr=0,nrow-1 do begin
     if rr eq 0 then begin
        ncol = rowend[rr] + 1
        rng = lindgen(ncol)
        ;; Print header
        ;; Prepend
        if keyword_set(out_fil) then $
           printf,1,' ',format=sfmt $
        else print,' ',format=sfmt 
        for cc=0,ncol-1 do begin
           if keyword_set(out_fil) then $
              printf,1,ks2str.file2_short[cc],format=sfmt $
           else print,ks2str.file2_short[cc],format=sfmt 
        endfor                  ; loop cc=ncol
        ;; End
        if keyword_set(out_fil) then $
           printf,1,' ' $
        else print,' ' 
     endif else begin
        ncol = rowend[rr] - rowend[rr-1]
        rng = rowend[rr-1] + 1 + lindgen(ncol)
     endelse 
     
     for cc=0,ncol-1 do begin
        if cc eq 0 then begin
           ;; Column of names
           if keyword_set(out_fil) then $
              printf,1,ks2str.file1_short[rng[cc]],format=sfmt $
           else print,ks2str.file1_short[rng[cc]],format=sfmt 
        endif 
        if ks2str.prob_ks2[rng[cc]] lt 1.e-5 then fmt = afmt $
        else if ks2str.prob_ks2[rng[cc]] lt 0.05 then $
           fmt = fmt0_bd else fmt = fmt0

        if keyword_set(out_fil) then $
           printf,1,ks2str.d_ks2[rng[cc]],ks2str.n_eff[rng[cc]],$
                  ks2str.prob_ks2[rng[cc]],format=fmt $
        else print,ks2str.d_ks2[rng[cc]],ks2str.n_eff[rng[cc]],$
                   ks2str.prob_ks2[rng[cc]],format=fmt 
     endfor                     ; loop cc=ncol

     ;; End
     if keyword_set(out_fil) then $
        printf,1,' ' $
     else print,' ' 
     
  endfor                        ; loop rr=nrow
  
  if keyword_set(out_fil) then begin
     close,1
     print,'sdss_fitfxw_prntkstwo: created ',out_fil
  endif 

end                             ; sdss_fitfxw_prntkstwo


pro sdss_fitfxw_kstwo, fit_lst, out_fil, psfil=psfil, plot=plot
  ;; compare everything in the list to everything else
  if n_params() ne 2 then begin
     print,'Syntax - sdss_Fitfxw_kstwo, fit_lst, out_fil, [psfil=, /plot]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  readcol,fit_lst,fit_fil,format='a',/silent
  nfil = (size(fit_fil,/dim))[0]

  ncompare = (nfil^2 - nfil) / 2 ; think of the box
  ostrct = {file1:strarr(ncompare), $
            file1_short:strarr(ncompare), $
            file2:strarr(ncompare), $
            file2_short:strarr(ncompare), $
            title:strarr(ncompare), $
            d_ks2:fltarr(ncompare,/nozero), $
            n_eff:lonarr(ncompare,/nozero), $
            prob_ks2:fltarr(ncompare,/nozero) $
           }

  if keyword_set(psfil) then begin
     plot = 1
     x_psopen,psfil,/maxs,/portrait 
     !p.multi = [1,1,1]
     !x.margin = [8.7,1.5]      ; left and right border
     !y.margin = [3.2,2]        ; bottom and top border
  endif


  count = 0
  for ff=0,nfil-2 do begin      ; don't have to go all way up; ff for first
     fitstr1 = xmrdfits(sdssdir+fit_fil[ff],1,/silent)
     prs = strsplit(fit_fil[ff],'/',/extract,count=nprs)
     nam = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'.',/reverse_search))
     ostrct.file1_short[count] = strmid(nam,strpos(nam,'_',/reverse_search)+1)
     ostrct.file1[count] = fit_fil[ff]

     for ss=ff+1,nfil-1 do begin ; ss for second
        fitstr2 = xmrdfits(sdssdir+fit_fil[ss],1,/silent)
        if ss gt ff+1 then begin
           ostrct.file1[count] = ostrct.file1[count-1]
           ostrct.file1_short[count] = ostrct.file1_short[count-1]
        endif 
        prs = strsplit(fit_fil[ss],'/',/extract,count=nprs)
        nam = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'.',/reverse_search))
        ostrct.file2_short[count] = strmid(nam,strpos(nam,'_',/reverse_search)+1)
        ostrct.file2[count] = fit_fil[ss]
        ostrct.title[count] = ostrct.file1_short[count] + ':' + $
                              ostrct.file2_short[count] 
        
        kstwo_wght, fitstr1.data, fitstr2.data, d, prob, $
                    n_eff=n_eff, wght1=fitstr1.xdata, wght2=fitstr2.xdata, $
                    plot=plot, title=ostrct.title[count]

        ;; Save result
        ostrct.d_ks2[count] = d
        ostrct.prob_ks2[count] = prob
        ostrct.n_eff[count] = n_eff ; = n1*n2 / float(n1+n2)

        ;; Increment
        count++ 
     endfor                     ; loop ss=nfil
     
  endfor                        ; loop ff=nfil-1

  if keyword_set(psfil) then begin
     x_psclose
     print,'sdss_fitfxw_kstwo: created ',psfil
  endif 

  
  ;; Write out
  mwrfits,ostrct,out_fil,/create,/silent
  spawn,'gzip -f '+out_fil
  print,'sdss_fitfxw_kstwo: created ',out_fil

  ;; Table too
  out_tab = strmid(out_fil,0,strpos(out_fil,'.',/reverse_search))+'.tab'
  sdss_fitfxw_prntkstwo, ostrct, out_fil=out_tab

end                             ;  sdss_fitfxw_kstwo



function sdss_fitfxw_int_freqdistr, var
  ;; Integral of (var/var_0)^px[1] d(var) 
  ;; or 
  ;;             exp(var/var_0 * px[1]) d(var) 
  ;; or
  ;;             (var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) d(var)
  common sdss_fitfxw_cmmn, ffit_intlim, ffit_data, ffit_sigdata, ffit_zdata, $
     ffit_xdata, ffit_zmed, ffit_cumx, ffit_ewseam, ffit_coeff, ffit_alpha, $
     ffit_dblt, ffit_datatype, ffit_flg, ffit_datanorm, ffit_sigdatanorm, $
     qromb_par, ffit_z, ffit_ewmax, ffit_form, ffit_ndim, ffit_debug, $
     ffit_civstrct, ffit_sat, ffit_unsat, ffit_nsat, ffit_nunsat

  if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
  else dat_nrm = ffit_datanorm

  rslt = sdss_calcfxwfit(var,option=ffit_form,coeff=1.,alpha=ffit_alpha,$
                         datanorm=dat_nrm)
  if keyword_set(ffit_debug) then begin
     test = where(finite(rslt) eq 0)
     if test[0] ne -1 then stop,'sdss_fitfxw_int_freqdistr() debug stop: non-finite result'
  endif 

  return, rslt
  
end                             ; sdss_fitfxw_int_freqdistr()


function sdss_fitfxw_getdxw, var0

  common sdss_fitfxw_cmmn

  if ffit_datatype eq 'NCOLM' then var = alog10(var0) else var = var0
  
  ;; Using ffit_zmed here because need some canonical value for the
  ;; integration 
  x_var = sdss_getdxw(ffit_cumx,ffit_zmed,var,dz=ffit_z,$
                      ewmax=ffit_ewmax, /silent) 
  return, x_var

end                             ; sdss_fitfxw_getdxw()


function sdss_fitfxw_int_freqdistrdx, var
  ;; Integral of (var/var_0)^px[1] * X(var) d(var) 
  ;; or 
  ;;             exp(var/var_0 * px[1]) * X(var) d(var) 
  ;; or
  ;;             (var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) * X(var) d(var)
  common sdss_fitfxw_cmmn

  x_var = sdss_fitfxw_getdxw(var)

  rslt = sdss_fitfxw_int_freqdistr(var)
;  if keyword_set(ffit_debug) then printcol,var,rslt,x_var,rslt*x_var
  return, rslt * x_var
  
end                             ; sdss_fitfxw_int_freqdistrdx()


function sdss_fitfxw_int_ddalpha_freqdistr, var
  ;; d(ln L)/d(px[1]) 
  ;; Integral of ln(var/var_0) * (var/var_0)^px[1]
  ;; or
  ;;             var/var_0 * exp(var/var_0*px[1]) 
  ;; or
  ;;             ln(var/exp(px[2])) * (var/exp(px[2]))^px[1] * exp(-var/exp(px[2]))
  common sdss_fitfxw_cmmn

  case ffit_form of
     'pow': rslt = alog(var/ffit_datanorm) * $
              (var/ffit_datanorm)^ffit_alpha 
     'exp': rslt = var/ffit_datanorm * $
                   exp(var/ffit_datanorm * ffit_alpha) 
     'sch': rslt = alog(var/exp(ffit_datanorm)) * $
                   (var/exp(ffit_datanorm))^ffit_alpha * $
                   exp(-var/exp(ffit_datanorm))
  endcase

  return,  rslt 

end                             ; sdss_fitfxw_int_ddalpha_freqdistr()


function sdss_fitfxw_int_ddalpha_freqdistrdx, var
  ;; d(ln L)/d(px[1]) 
  ;; Integral of ln(var/var_0) * (var/var_0)^px[1] * X(var) 
  ;; or
  ;;             var/var_0 * exp(var/var_0*px[1]) * X(var) 
  ;; or
  ;;             ln(var/exp(px[2])) * (var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) * X(var)
  common sdss_fitfxw_cmmn

  x_var = sdss_fitfxw_getdxw(var)
  
  rslt = sdss_fitfxw_int_ddalpha_freqdistr(var)

  return,  rslt * x_var

end                             ; sdss_fitfxw_int_ddalpha_freqdistrdx()


function sdss_fitfxw_int_ddnorm_freqdistr, var
  ;; d(ln L)/d(px[2]) 
  ;; Only for Schechter function (or any future 3-parameter fits)
  ;; Integral of (var/px[2])^(px[1]+2) * (-px[1]*px[2] + var)/var^2
  ;;             * exp(-var/px[2]) 
  common sdss_fitfxw_cmmn
  
  case ffit_form of
     'pow': rslt = 0
     'exp': rslt = 0
     'sch': $
        rslt = (var/exp(ffit_datanorm))^(ffit_alpha+2) * $
               (-ffit_alpha*exp(ffit_datanorm) + var)/var^2 * $
               exp(-var/exp(ffit_datanorm)) * $
               exp(ffit_datanorm) ; this converts d ln(L)/d N*_0 to d ln(L)/d N*
                                ; where N*_0 = exp(N*) 
  endcase

  return,  rslt 

end                             ; sdss_fitfxw_int_ddnorm_freqdistr()


function sdss_fitfxw_int_ddnorm_freqdistrdx, var
  ;; d(ln L)/d(px[2]) 
  ;; Only for Schechter function (or any future 3-parameter fits)
  ;; Integral of (var/px[2])^(px[1]+2) * (-px[1]*px[2] + var)/var^2
  ;;             * exp(-var/px[2]) * X(var)
  common sdss_fitfxw_cmmn

  x_var = sdss_fitfxw_getdxw(var)
  
  rslt = sdss_fitfxw_int_ddnorm_freqdistr(var)

  return,  rslt * x_var

end                             ; sdss_fitfxw_int_ddnorm_freqdistrdx()


function sdss_fitfxw_logL, px, gradient
  ;; px = [coeff, alpha]
  ;; 1 ln L = -exp(px[0])*INT((var/var_0)^px[1] * 
  ;; 2                       X(var) * d(var), var_min, var_sat)
  ;; 3       -exp(px[0])*INT((var/var_0)^px[1] * 
  ;; 4                       X(var) * d(var), var_sat, var_max) 
  ;; 5       + px[0] * (n_elements(data[unsat]) + n_elements(data[sat]))
  ;; 6       + px[1] * SUM(ln(data[unsat]/var_0))
  ;; 7       + SUM(ln(X(data[unsat])) 
  ;; 8       + n_elements(data[sat]) * 
  ;; 9         ln( INT((var/var_0)^px[1] * 
  ;;10                 X(var) * d(var), var_sat, var_max) )
  ;; if fitting exponential:
  ;;         (var/var_0)^px[1] terms --> exp(var/var_0*px[1])
  ;; 6       + px[1] * SUM(data[unsat]/var_0)
  ;; if fitting Schechter Function
  ;;   ln L = -exp(px[0]*INT((var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) * 
  ;;                         X(var) * d(var), var_min, var_sat) 
  ;;          -exp(px[0]*INT((var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) * 
  ;;                         X(var) * d(var), var_sat, var_max) 
  ;;          + SUM( px[1] * alog(data[unsat]/exp(px[2])) -
  ;;                 data[unsat]/exp(px[2]) + alog(X(data[unsat])) )
  ;;          + px[0] * (n_elements(data[unsat]) + n_elements(data[sat]))
  ;;          + n_elements(data[sat]) * 
  ;;            ln( INT((var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) * 
  ;;                    X(var) * d(var), var_sat, var_max) )
  
  common sdss_fitfxw_cmmn

  ;; Be smarter for grid searches
  if (size(px,/n_dim))[0] eq 1 then nelem = 1 $
  else nelem = (size(px,/dim))[0]

  ;; Faster not to input values
  logL = dblarr(nelem,/nozero)
  int_freqdistrdx = dblarr(nelem,2,/nozero) ; [min,sat] and [sat,max]
  ;; for ln L, d(lnL)/d(px[1]), and d(lnL)/d(px[2])
  sum_freqdistr_all = dblarr(nelem,4,/nozero) ; may not use 4th
  if arg_present(gradient) then begin
     ;; derivative is min to sat and sat to max
     int_ddalpha_freqdistrdx = dblarr(nelem,2,/nozero) 
     gradient = dblarr(nelem,ffit_ndim,/nozero)
     if ffit_form eq 'sch' then $
        int_ddnorm_freqdistrdx = dblarr(nelem,2,/nozero) 
  endif 

  if nelem eq 1 then begin
     coeff_arr = px[0]
     alpha_arr = px[1]
     if ffit_form eq 'sch' then datanorm_arr = px[2] $
     else datanorm_arr = ffit_datanorm
  endif else begin
     coeff_arr = px[*,0]
     alpha_arr = px[*,1] 
     if ffit_form eq 'sch' then datanorm_arr = px[*,2] $
     else datanorm_arr = replicate(ffit_datanorm,nelem)
  endelse 


  for ii=0L,nelem-1 do begin
     ;; Set for frequency distribution integrals
     ffit_alpha = alpha_arr[ii]
     ffit_datanorm = datanorm_arr[ii] ; only matters for Schechter
     if ffit_form eq 'sch' and exp(ffit_datanorm) gt ffit_intlim[2] then begin
        print,'sdss_fitfxw_logL(): ffit_datanorm = ',exp(ffit_datanorm)
        if arg_present(gradient) then $
           stop,'sdss_fitfxw_logL() stop: WARNING! cannot compute gradient!'
        return,0.               ; large logL
     endif 

     if ii gt 0 then begin
        if ffit_alpha eq alpha_arr[ii-1] and $
           ffit_datanorm eq datanorm_arr[ii-1] then begin ; matters for Schechter
           ;; Save computational time and don't repeat integration
           int_freqdistrdx[ii,*] = int_freqdistrdx[ii-1,*]
           sum_freqdistr_all[ii,*] = sum_freqdistr_all[ii-1,*]

           if keyword_set(gradient) then begin
              int_ddalpha_freqdistrdx[ii,*] = int_ddalpha_freqdistrdx[ii-1,*]
              if ffit_form eq 'sch' then $
                 int_ddnorm_freqdistrdx[ii,*] = int_ddnorm_freqdistrdx[ii-1,*]
           endif 

           continue             ; SKIPPING
        endif  
     endif 

     int_freqdistrdx[ii,0] = $
        qromb('sdss_fitfxw_int_freqdistrdx',ffit_intlim[0],ffit_intlim[1],$
              eps=qromb_par[0], /double, jmax=qromb_par[1]) 
     
     ;; for ln L, d(lnL)/d(px[1]), and d(lnL)/d(px[2])
     if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
     else dat_nrm = ffit_datanorm
     sum_freqdistr_all[ii,0] = total(alog(ffit_data[ffit_unsat]/dat_nrm)) ; calc once
     sum_freqdistr_all[ii,1] = total(ffit_data[ffit_unsat]/dat_nrm)       ; calc once
     sum_freqdistr_all[ii,2] = $                                           ; SUM(dX[unsat]) included below
        ffit_alpha * sum_freqdistr_all[ii,0] - sum_freqdistr_all[ii,1]
     sum_freqdistr_all[ii,3] = $
        total((-ffit_alpha*dat_nrm + ffit_data[ffit_unsat])/dat_nrm^2)

     ;; Save time
     if ffit_nsat gt 0 then $
        int_freqdistrdx[ii,1] = $
        qromb('sdss_fitfxw_int_freqdistrdx',ffit_intlim[1],ffit_intlim[2],$
              eps=qromb_par[0], /double, jmax=qromb_par[1]) $
     else int_freqdistrdx[ii,1] = 0


     ;; ;;;;;;;
     if arg_present(gradient) then begin
        ;; d(ln L)/d(px[1]) 
        int_ddalpha_freqdistrdx[ii,0] = $
           qromb('sdss_fitfxw_int_ddalpha_freqdistrdx',$
                 ffit_intlim[0],ffit_intlim[1], eps=qromb_par[0], $
                 /double, jmax=qromb_par[1])

        ;; Save time
        if ffit_nsat gt 0 then $
           int_ddalpha_freqdistrdx[ii,1] = $
           qromb('sdss_fitfxw_int_ddalpha_freqdistrdx',$
                 ffit_intlim[1],ffit_intlim[2], eps=qromb_par[0], $
                 /double, jmax=qromb_par[1]) $
        else int_ddalpha_freqdistrdx[ii,1] = 0

        if ffit_form eq 'sch' then begin
           ;; d(ln L)/d(px[2]) 
           int_ddnorm_freqdistrdx[ii,0] = $
              qromb('sdss_fitfxw_int_ddnorm_freqdistrdx',$
                    ffit_intlim[0],ffit_intlim[1], eps=qromb_par[0], $
                    /double, jmax=qromb_par[1])

           ;; Save time
           if ffit_nsat gt 0 then $
              int_ddnorm_freqdistrdx[ii,1] = $
              qromb('sdss_fitfxw_int_ddnorm_freqdistrdx',$
                    ffit_intlim[1],ffit_intlim[2], eps=qromb_par[0], $
                    /double, jmax=qromb_par[1]) $
           else int_ddnorm_freqdistrdx[ii,1] = 0
        endif 

     endif                      ; gradient set 
  endfor                        ; loop ii = nelem

  case ffit_form of
     'pow': sum_freqdistr = alpha_arr * sum_freqdistr_all[*,0]
     'exp': sum_freqdistr = alpha_arr * sum_freqdistr_all[*,1]
     'sch': sum_freqdistr = sum_freqdistr_all[*,2] ; already set; changes for grad
  endcase

  logL = -exp(coeff_arr) * (int_freqdistrdx[*,0] + $ ; [min,sat]
                            int_freqdistrdx[*,1] ) + $ ; [sat,max]
         coeff_arr * double(n_elements(ffit_data)) + $
         sum_freqdistr + total(alog(ffit_xdata[ffit_unsat]))
  if ffit_nsat gt 0 then logL = logL + $
                           ffit_nsat * alog( int_freqdistrdx[*,1] ) 


  if keyword_set(ffit_debug) then begin
;     print,'px ',exp(px[0]),px[1],exp(px[2])
     print,'ln L  ',format='(a,$)'
     if ffit_datatype eq 'NCOLM' then fmt = '(e11.3,2x,2(e11.3,2x),3(e10.2,2x))' $
     else fmt = '(f6.2,2x,2(f9.2,2x),3(f9.1,2x))'
     if nelem eq 1 then begin
        print,-exp(coeff_arr),int_freqdistrdx[0,0],int_freqdistrdx[0,1],$
              coeff_arr*n_elements(ffit_data),sum_freqdistr[0],$
              total(alog(ffit_xdata[ffit_unsat])), format=fmt
     endif else begin
        printcol,-exp(coeff_arr),int_freqdistrdx[*,0],int_freqdistrdx[*,1],$
                 coeff_arr*n_elements(ffit_data),sum_freqdistr,format=fmt
     endelse
     test = where(finite(logl) eq 0)
     if test[0] ne -1 then stop,'sdss_fitfxw_logL() debug stop: non-finite ln L'
;     if not keyword_set(gradient) then stop
  endif

  if keyword_set(gradient) then begin
     ;; gradient = [d(ln L)/d(px[0]), d(ln L)/d(px[1])]
     ;; 1 d(ln L)/d(px[0]) = -exp(px[0]) * INT((var/var_0)^px[1] * 
     ;; 2                                      X(var) * d(var), var_min, var_sat) 
     ;; 3                    -exp(px[0]) * INT((var/var_0)^px[1] * 
     ;; 4                                      X(var) * d(var), var_sat, var_max) 
     ;; 5                    + n_elements(data)
     ;; 6 d(ln L)/d(px[1]) = -exp(px[1]) * $
     ;; 7                     INT(ln(var/var_0)*(var/var_0)^px[1] * 
     ;; 8                          X(var) * d(var), var_min, var_max)
     ;; 9                    -exp(px[1]) * $
     ;;10                     INT(ln(var/var_0)*(var/var_0)^px[1] * 
     ;;11                          X(var) * d(var), var_sat, var_max)
     ;;12                    + SUM(ln(data[unsat]/var_0)) 
     ;;13                    + n_elements(data[sat]) * 
     ;;14                      INT(ln(var/var_0)*(var/var_0)^px[1] * 
     ;;15                          X(var) * d(var), var_sat, var_max) / 
     ;;16                      INT((var/var_0)^px[1] * X(var) * d(var), 
     ;;17                          var_sat, var_max)
     ;; if fitting exponential:
     ;;         (var/var_0)^px[1] terms --> exp(var/var_0*px[1])
     ;;         INT(ln(var/var_0)*(var/var_0)^px[1]) -->
     ;;                     INT(var/var_0*exp(var/var_0*px[1]))
     ;; if fitting Schechter function:
     ;;   d(ln L)/d(px[0]) = -exp(px[0]) * INT((var/exp(px[2]))^px[1] *
     ;;                                        exp(-var/exp(px[2])) * 
     ;;                                        X(var) * d(var), var_sat, var_max) 
     ;;   d(ln L)/d(px[1]) = -exp(px[0]) * 
     ;;                      INT(ln(var/exp(px[2])) * (var/exp(px[2]))^px[1] *
     ;;                          exp(-var/exp(px[2])) * X(var) * d(var),
     ;;                          var_min, var_sat)  - exp(px[0]) * 
     ;;                      INT(ln(var/exp(px[2])) * (var/exp(px[2]))^px[1] *
     ;;                          exp(-var/exp(px[2])) * X(var) * d(var),
     ;;                          var_sat, var_max) 
     ;;                      + SUM(alog(var[unsat]/exp(px[2])))
     ;;                      + n_elements(data[sat]) * 
     ;;                      INT(ln(var/exp(px[2])) * (var/exp(px[2]))^px[1] *
     ;;                          exp(-var/exp(px[2])) * X(var) * d(var),
     ;;                          var_sat, var_max) / 
     ;;                          INT((var/exp(px[2]))^px[1] *
     ;;                          exp(-var/exp(px[2])) * X(var) * d(var),
     ;;                          var_sat, var_max) 
     ;;   d(ln L)/d(px[2]) = ( -exp(px[0]) * 
     ;;                        INT((var/exp(px[2]))^(alpha+2) * (-px[1]*exp(px[2])
     ;;                            + var)/var^2 * exp(-var/exp(px[2])) * X(Var) *
     ;;                            d(var), var_min, var_sat) + -exp(px[0]) * 
     ;;                        INT((var/exp(px[2]))^(alpha+2) * (-px[1]*exp(px[2])
     ;;                            + var)/var^2 * exp(-var/exp(px[2])) * X(Var) *
     ;;                            d(var), var_sat, var_max) 
     ;;                        + SUM( (-px[1]*exp(px[2]) + var)/exp(px[2])^2 )
     ;;                        + n_elements(data[sat]) *
     ;;                        INT((var/exp(px[2]))^(alpha+2) * (-px[1]*exp(px[2])
     ;;                            + var)/var^2 * exp(-var/exp(px[2])) * X(Var) *
     ;;                            d(var), var_sat, var_max) /
     ;;                        INT((var/exp(px[2]))^px[1] * exp(-var/exp(px[2])) *
     ;;                            X(var) * d(var), var_sat, var_max)
     ;;                       ) * exp(px[2])

     gradient[*,0] = -exp(coeff_arr) * int_freqdistrdx[*,0] $
                     -exp(coeff_arr) * int_freqdistrdx[*,1] + $
                     double(n_elements(ffit_data))

     if ffit_form eq 'sch' then begin
        gradient[*,2] = -exp(coeff_arr)*int_ddnorm_freqdistrdx[*,0] $
                        -exp(coeff_arr)*int_ddnorm_freqdistrdx[*,1] $
                        + sum_freqdistr_all[*,3] * exp(datanorm_arr) ; this converts d ln(L)/d N*_0 to d ln(L)/d N*
                                ; where N*_0 = exp(N*) 
        if ffit_nsat gt 0 then $
           gradient[*,2] = gradient[*,2] + ffit_nsat * int_ddnorm_freqdistrdx[*,1] / $
                           int_freqdistrdx[*,1] 
        
        sum_freqdistr = sum_freqdistr_all[*,0] ; d/d(px[1])
     endif else sum_freqdistr = sum_freqdistr / ffit_alpha

     gradient[*,1] = -exp(coeff_arr)*int_ddalpha_freqdistrdx[*,0] $
                     -exp(coeff_arr)*int_ddalpha_freqdistrdx[*,1] $
                     + sum_freqdistr 
     if ffit_nsat gt 0 then $
        gradient[*,1] = gradient[*,1] + ffit_nsat * int_ddalpha_freqdistrdx[*,1] / $
                        int_freqdistrdx[*,1]

     if nelem eq 1 then $
        gradient = -gradient[0,*] $
     else gradient = -gradient  ; minimize (so negative)     
     if keyword_set(ffit_debug) then begin
        ;; only get here if only one element
        print,'d(ln L)  ',format='(a,$)'
        fmt = '(f9.1,2x,2(f7.2,2x,f9.2),$)'
        if ffit_datatype eq 'NCOLM' then $
           print,-logL[0],alog10(exp(px[0])),gradient[0],px[1],gradient[1],$
                 format=fmt $ ;'(f9.1,2x,e12.3,2x,f9.2,2x,f7.2,2x,f9.2,$)' $
        else $
           print,-logL[0],exp(px[0]),gradient[0],px[1],gradient[1],$
                 format=fmt
        if ffit_form eq 'sch' then begin
           fmt = '(2x,2(f9.2,2x))' 
           if ffit_datatype eq 'NCOLM' then $
              print,alog10(exp(px[2])),gradient[2],format=fmt $ ;'(2x,e11.3,2x,f9.2)' $
           else print,exp(px[2]),gradient[2],format=fmt
        endif else print,''           ; close line
     endif 
  endif                         ; gradient set

  if nelem eq 1 then begin
     return, -logL[0]           ; minimize
  endif else begin
     return, -logL
  endelse 
end                             ; sdss_fitfxw_logL()



function sdss_fitfxw_kscumf, var
  ;; CDF = INT( f(var) * d(var), var_min, var ) 
  ;;     = exp(coeff_norm)/(1 + alpha) * var_0^(-alpha) *
  ;;       (var^(1+alpha) - var_min^(1+alpha))
  ;; or
  ;;     = exp(coeff_norm) * var_0/alpha * (exp(alpha * var / var_0) 
  ;;       - exp(alpha * var_min/var_0))
  ;; or 
  ;;     = -exp(coeff_norm)*var_0*(igamma(1+alpha, var/exp(var_0)) - 
  ;;                               igamma(1+alpha, var_min/exp(var_0)))
  ;; coeff_norm = alog((1+alpha)*var_0^alpha / 
  ;;                   (var_sat^(1+alpha) - var_min^(1+alpha)))
  ;; or
  ;;            = alog(alpha/var_0 / (exp(alpha*var_sat/var_0) -
  ;;              exp(alpha * var_min/var_0))
  ;; or 
  ;;            = -1/var_0 * (igamma(1+alpha, var_sat/exp(var_0)) - 
  ;;                          igamma(1+alpha, var_min/exp(var_0)))^(-1)
  common sdss_fitfxw_cmmn

  nvar = (size(var,/dim))[0] > 1 ; foil singularity
  rng = dblarr(nvar,2,/nozero)
  rng[*,0] = ffit_intlim[0]
  rng[*,1] = var
  if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
  else dat_nrm = ffit_datanorm
  rslt = sdss_calcdndxfit(rng, option=ffit_form, coeff=1., alpha=ffit_alpha, $
                          datanorm=dat_nrm)

  return, exp(ffit_coeff) * rslt
end                             ; sdss_fitfxw_kscumf()


;function sdss_fitfxw_omciv, px, intlim=intlim, _extra=extra
;  ;; CIV mass density 
;  ;; Analytic formula
;  common sdss_fitfxw_cmmn
;
;  if not keyword_set(intlim) then intlim = [10.^13,10.^15]
;
;  ;; Be smarter for grid searches
;  if (size(px,/n_dim))[0] eq 1 then nelem = 1 $
;  else nelem = (size(px,/n_elem))
;
;  logL = dblarr(nelem)
;  int_freqdistrdx = dblarr(nelem,2)
;  if arg_present(gradient) then begin
;     int_ddalpha_freqdistrdx = dblarr(nelem,2)
;     gradient = dblarr(nelem,2)
;  endif 
;
;  if nelem eq 1 then begin
;     coeff_arr = px[0]
;     alpha_arr = px[1]
;  endif else begin
;     coeff_arr = px[*,0]
;     alpha_arr = px[*,1] 
;  endelse 
;
;;  omciv = dblarr(nelem)
;;  for ii=0L,nelem-1 do begin
;;     omciv[ii] = civ_omciv_fit(intlim, exp(coeff_arr[ii]), $
;;                               alpha_arr[ii], ffit_datanorm, $
;;                               dblt_name=ffit_dblt, _extra=extra)
;;  endfor 
;  
;  return, omciv
;end                             ; sdss_fitfxw_omciv()


pro sdss_fitfxw_fitall, cumcmplt_list, final=final, exp=exp, sch=sch, $
                        intlim=intlim, ncolm=ncolm, dz=dz, odir=odir, $
                        float=float,excl_hi=excl_hi,civstrct_fil=civstrct_fil,$
                        _extra=extra
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fitfxw_fitall, cumcmplt_list, [/final, /exp, /sch, '
     print,'                        intlim=, /ncolm, /dz, odir=, '
     print,'                        /float,excl_hi=,civstrct_fil=,_extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  if keyword_set(dz) then xstr = 'z' else xstr = 'x'
  if keyword_set(sch) then suffix = 'sch' $
  else if keyword_set(exp) then suffix = 'exp' $
  else suffix = 'pow' 

  prefix = 'f'+xstr
  if keyword_set(ncolm) then prefix = prefix+'n' $
  else prefix = prefix+'w'
  if not keyword_set(intlim) then intlim = [0.6,3.5,3.5]

  readcol, cumcmplt_list, cmplt_fil, format='(a)', /silent
  ncmplt = (size(cmplt_fil,/dim))[0] > 1

  ;; _extra includes dvqso= (override), rating= (override), ewlim=,
  ;; zlim=, dvem=
  civstr0 = sdss_getcivstrct(civstrct_fil,/default,final=final,_extra=extra)
  tags = tag_names(civstr0)
  if keyword_set(final) then begin
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
  endif else begin
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
  endelse 

  for cc=0,ncmplt-1 do begin
     cmpltstr = xmrdfits(sdssdir+cmplt_fil[cc],1,/silent)

     sub = where(civstr0.(ztag)[0] ge cmpltstr.zlim[0] and $
                 civstr0.(ztag)[0] lt cmpltstr.zlim[1],nsub)
     if nsub eq 0 then begin
        print,'sdss_fitfxw_fitall: no doublets in zlim = ',cmpltstr.zlim
        continue
     endif
     if keyword_set(excl_hi) then begin
        ;; But only want to exclude fraction of that which will be fit
        sub2 = where(civstr0[sub].(ewtag)[0] ge intlim[0],nsub2,$
                    complement=low)
        srt = sort(civstr0[sub[sub2]].(ewtag)[0])
        sub2 = sub[sub2[srt]]
        if excl_hi lt 1. then $
           nexcl = round(excl_hi * nsub2) $ ; Fractional difference
        else nexcl = excl_hi
        nsub2 = nsub2 - nexcl
        nsub = nsub - nexcl
        sub = [sub[low],sub2[0:nsub2-1]] ; trim off last bit
     endif 
     prs = strsplit(cmplt_fil[cc],'/',/extract,count=nprs)
     froot = odir+prefix + strmid(prs[nprs-1],5)
     froot = strmid(froot,0,strpos(froot,'.',/reverse_search))+suffix

     ;; _extra includes /ewmax, /plot, dblt_name=, data_norm=,
     ;; conv_factor_lim=, nbin_alpha=, nbin_coeff=, dcoeff=, dalpha=,
     ;; int_param=, /tcpu, /clobber, /afp_flg, biasuser_fil= [_extra=]
     print,''
     if keyword_set(float) then $
        intlim[1:2] = max(civstr0[sub].(ewtag)[0]+civstr0[sub].(sigewtag)[0])
     sdss_fitfxw, cmplt_fil[cc], $
                  intlim=intlim, exp=exp, sch=sch, civstrct_fil=civstr0[sub], $
                  ostrct_fil=froot+'.fit', psfil=froot+'.ps',/plot,$
                  dz=dz, ncolm=ncolm, _extra=extra
  endfor                        ; loop cc=ncmplt

  print,'sdss_fitfxw_fitall: done!'

end                             ; sdss_fitfxw_fitall


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_fitfxw, sens_fil, guess=guess, civstrct_fil=civstrct_fil, $
                 exp=exp, sch=sch, dnorm=dnorm, nbin_norm=nbin_norm, $
                 dblt_name=dblt_name, data_norm=data_norm, ewmax=ewmax, $
                 intlim=intlim, ncolm=ncolm, ostrct_fil=ostrct_fil, $
                 conv_factor_lim=conv_factor_lim, infit_fil=infit_fil,fitonly=fitonly,$
                 erronly=erronly,ksonly=ksonly,plot=plot,resmplerr=resmplerr,$
                 nbin_alpha=nbin_alpha,nbin_coeff=nbin_coeff,psfil=psfil,$
                 dcoeff=dcoeff,dalpha=dalpha,int_param=int_param,$
                 tcpu=tcpu,dz=dz,final=final,clobber=clobber,$
                 afp_flg=afp_flg, biasuser_fil=biasuser_fil, $
                 status=status, debug=debug,_extra=extra

  ;; Instantiate common block
  ;; (ffit_intlim, ffit_data, ffit_sigdata, ffit_xdata, 
  ;;  ffit_cumx, ffit_ewseam, ffit_coeff, ffit_alpha, ffit_dblt, ffit_datatype,
  ;;  ffit_flg, ffit_datanorm,  qromb_par, ffit_z, ffit_ewmax,
  ;;  ffit_form, ffit_ndim, ffit_debug, ffit_civstrct, ffit_sat,
  ;;  ffit_unsat, ffit_nsat, ffit_nunsat)
  common sdss_fitfxw_cmmn 

  ffit_debug = keyword_set(debug)
  if size(sens_fil,/type) eq 7 then $
     ffit_cumx = xmrdfits(sens_fil,1,/silent) $
  else ffit_cumx = sens_fil
  tags = tag_names(ffit_cumx)
  mtch = where(tags eq 'EWSEAM')
  if mtch[0] then ffit_ewseam = 0 else begin
     ffit_ewseam = ffit_cumx.(mtch[0])
     if keyword_set(ncolm) then ffit_ewseam = 10.^ffit_ewseam
  endelse 
  if keyword_set(int_param) then qromb_par = int_param $
  else qromb_par = [1.e-10, 20.]  ; eps=, jmax=
  ffit_z = keyword_set(dz)        ; f(N) is over path length
  ffit_ewmax = keyword_set(ewmax) ; force highest value automatically
  if keyword_set(exp) then begin
     ffit_form = 'exp'          ; fit exponential
     ffit_ndim = 2
  endif else begin
     if keyword_set(sch) then begin
        ffit_form = 'sch'       ; fit Schechter function
        ffit_ndim = 3
     endif else begin
        ffit_form = 'pow'       ; fit default power law
        ffit_ndim = 2
     endelse
  endelse 

  if not keyword_set(dblt_name) then dblt_name = 'CIV' 
  if size(dblt_name,/type) eq 7 then ffit_dblt = dblt_retrieve(dblt_name) $
  else ffit_dblt = dblt_name

  if keyword_set(erronly) or keyword_set(ksonly) then begin
     ;; Read in previous results
     fitstrct = xmrdfits(infit_fil,1,/silent)
;     civcand = xmrdfits(infit_fil,2,/silent)
     
     ;; Set up common values
     ffit_dblt = dblt_retrieve(fitstrct.dblt_name)
     ffit_datatype = strtrim(fitstrct.datatype,2)
     if ffit_form eq 'sch' then ffit_datanorm = alog(fitstrct.datanorm) $
     else ffit_datanorm = fitstrct.datanorm
     ffit_data = fitstrct.data
     ffit_sigdata = fitstrct.sigdata
     ffit_flg = fitstrct.flag
     ffit_xdata = fitstrct.xdata
     ffit_intlim = fitstrct.intlim
     ffit_coeff = alog(fitstrct.coeff)
     ffit_alpha = fitstrct.alpha

     if keyword_set(ksonly) and not keyword_set(erronly) then $
        errstrct = xmrdfits(infit_fil,3,/silet)
  endif else begin
     ;;;;;;;;;;;;;;
     ;; Solve for coeff and alpha
     ;;;;;;;;;;;;;;
     ;; _extra= includes rating=, zlim=, ewlim=, nlim=, dvgal=,
     ;; dvqso=, /noBAL, /unblend, /default, /dropbox
     if size(civstrct_fil,/type) eq 8 then begin
        ffit_civstrct = civstrct_fil ; assume trimmed right
        nciv = (size(ffit_civstrct,/dim))[0]
     endif else $
        ffit_civstrct = sdss_getcivstrct(civstrct_fil,count=nciv,final=final,$
                                         _extra=extra)
     tags = tag_names(ffit_civstrct)
     if keyword_set(final) then begin
        ztag = (where(tags eq 'ZABS_FINAL'))[0]
        ewtag = (where(tags eq 'EW_FINAL'))[0]
        sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
     endif else begin
        ztag = (where(tags eq 'ZABS_ORIG'))[0]
        ewtag = (where(tags eq 'EW_ORIG'))[0]
        sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
     endelse 

     if nciv eq 0 then $
        stop,'sdss_fitfxw: no '+ffit_dblt.ion


     ;; Data must be in linear space
     ;; Sort for kstest
     if keyword_set(ncolm) then begin
        ffit_datatype = 'NCOLM'
        if keyword_set(data_norm) then ffit_datanorm = data_norm $
        else begin
           ffit_datanorm = 10.d^14. 
           if ffit_form eq 'sch' then ffit_datanorm = alog(ffit_datanorm)
        endelse 
        ;; _extra includes /estflg
        ncolm  = sdss_calcncolmdblt(ffit_civstrct,ffit_dblt,signcolm=signcolm,$
                                    ncolmflg=ffit_flg,/silent,final=final,$
                                    _extra=extra)
        if keyword_set(resmplerr) then $
           ncolm = ncolm + randomn(seed,nciv)*signcolm
        srt = sort(ncolm)
        ffit_civstrct = ffit_civstrct[srt]
        ffit_data = ncolm[srt]  ; error-weighted
        ffit_sigdata = signcolm[srt]
        ffit_zmed = median(ffit_civstrct.(ztag)[0],/even)
        ffit_zdata = double(ffit_civstrct.(ztag)[0])
        ffit_flg = ffit_flg[srt]

        ffit_xdata = sdss_getdxw(ffit_cumx,ffit_zdata,alog10(ffit_data),$
                                 dz=ffit_z, ewmax=ffit_ewmax, /silent)

        if not keyword_set(ostrct_fil) then $
           ostrct_fil = 'fn_fit.fit'
        scl = 1.e2              ; for setting upper limit "infinity"
     endif else begin
        ;; Default
        ffit_datatype = 'EW'
        if keyword_set(data_norm) then ffit_datanorm = data_norm $
        else begin 
           ffit_datanorm = 1.d  ; Ang
           if ffit_form eq 'sch' then ffit_datanorm = alog(ffit_datanorm)
        endelse 
        if keyword_set(resmplerr) then $
           ffit_civstrct.(ewtag)[0] = ffit_civstrct.(ewtag)[0] + $
                                      randomn(seed,nciv)*ffit_civstrct.(sigewtag)[0]
        srt = sort(ffit_civstrct.(ewtag)[0])
        ffit_civstrct = ffit_civstrct[srt]
        ffit_data = double(ffit_civstrct.(ewtag)[0])
        ffit_sigdata = double(ffit_civstrct.(sigewtag)[0])
        ffit_zmed = median(ffit_civstrct.(ztag)[0],/even)
        ffit_zdata = double(ffit_civstrct.(ztag)[0])
        ffit_flg = replicate(1l,nciv) ; dummy
        ffit_xdata = sdss_getdxw(ffit_cumx,ffit_zdata,ffit_data,$
                                 dz=ffit_z,ewmax=ffit_ewmax, /silent)

        if not keyword_set(ostrct_fil) then $
           ostrct_fil = 'few_fit.fit'
        scl = 1.5               ; for setting upper limit "infinity"
     endelse

     ;; ;;;;;;;
     test = file_search(ostrct_fil+'*',count=ntest)
     if ntest ne 0 and not keyword_set(clobber) then begin
        print,'sdss_fitfxw: will not clobber file ',ostrct_fil
        return
     endif


     ;; Sub-sample, can be used to determine intlim[1] (saturated limit)
     unsat = where((ffit_flg and sdss_getlimflg(/lower)) eq 0,$
                   nunsat,complement=sat,ncomplement=nsat)

     ;; Set integration limits
     ffit_intlim = dblarr(3)
     if not keyword_set(intlim) then begin
        ;; Accept all data
        ffit_intlim[0] = min(ffit_data,max=mx)
        if nsat eq 0 then ffit_intlim[1] = scl * mx $ ; something large
        else ffit_intlim[1] = min(ffit_data[sat])
        ffit_intlim[2] = scl * mx ; something large 
     endif else begin
        if n_elements(intlim) ne 3 then $
           stop,'sdss_fitfxw: integration limit = 3-element array'
        ffit_intlim = intlim
     endelse 
     civcand0 = ffit_civstrct

     ;; Trim saturated absorbers from unsaturated limits
     ;; and absorbers below lowest integration limit
     bd = where((ffit_flg and sdss_getlimflg(/lower)) ge sdss_getlimflg(/lower)$
                and ffit_data lt ffit_intlim[1],nbd,complement=gd)
     if nbd ne 0 then begin
        print,'sdss_fitfxw: Number of saturated absorbers below intlim[1]',nbd
        ffit_data = ffit_data[gd]
        ffit_sigdata = ffit_sigdata[gd]
        ffit_flg = ffit_flg[gd]
        ffit_xdata = ffit_xdata[gd]
        ffit_zdata = ffit_zdata[gd]
        ffit_civstrct = ffit_civstrct[gd]
     endif 

     ;; Trim absorbers outside the limits of interest
     bd = where(ffit_data lt ffit_intlim[0] or $
                ffit_data gt ffit_intlim[2],nbd,complement=gd)
     if nbd ne 0 then begin
        print,'sdss_fitfxw: Number of absorbers outside intlim',nbd
        ffit_data = ffit_data[gd]
        ffit_sigdata = ffit_sigdata[gd]
        ffit_flg = ffit_flg[gd]
        ffit_xdata = ffit_xdata[gd]
        ffit_zdata = ffit_zdata[gd]
        ffit_civstrct = ffit_civstrct[gd]
     endif 
     ;; Saturation will be defined by >= ffit_intlim[1]
;     ;; Set unsaturated absorber flags above saturation limits to
;     ;; saturated
;     bd = where(ffit_data ge ffit_intlim[1] and $
;                (ffit_flg and sdss_getlimflg(/lower)) eq 0,nbd)
;     if nbd ne 0 and ffit_datatype eq 'NCOLM' then begin
;        print,'sdss_fitfxw: Number of absorbers that should be saturated',nbd
;        ffit_flg[bd] = sdss_setlimflg(ffit_flg[bd],/lower) ; saturated
;     endif 
     ;; Set saturated arrays from the get go, now that done fiddling
     ;; with them 
     ffit_unsat = where(ffit_data lt ffit_intlim[1],ffit_nunsat,$
                        complement=ffit_sat,ncomplement=ffit_nsat)


     ;; Check that there's still enough data!
     nciv = (size(ffit_data,/dim))[0] > 1 ; foil singularity
     if nciv le 1 then begin
        ;; Terrible!
        print,'sdss_fitfxw: insufficient data to fit',nciv
        print,'sdss_fitfxw: exiting without creating ',ostrct_fil
        if keyword_set(guess0) then guess = guess0 ; restore original
  
        return
     endif 


     if not keyword_set(guess) then begin
        ;; guess[0] is ln(normalization) in order to make its
        ;; magnitude (and gradient) approximately OOM with guess[1] (alpha)
        if keyword_set(infit_fil) then begin
           fitstrct = xmrdfits(infit_fil,1,/silent)
           if ffit_form eq 'sch' then $
              guess = [alog(fitstrct.coeff), fitstrct.alpha, $
                       alog(fitstrct.datanorm)] $
           else  guess = [alog(fitstrct.coeff), fitstrct.alpha]
           if not keyword_set(conv_factor_lim) then $
              conv_factor_lim = 0.5*fitstrct.conv_factor
        endif else begin
           ;; Use prior knowledge to set these 
           if keyword_set(ncolm) then begin
              ;; Column density
              case ffit_form of 
                 'pow': guess = [alog(4.48e-14/qffit_datanorm^(-1.55)*$
                                      ffit_datanorm^(-1.5)), -1.5d] 
                 'exp': guess = [alog(4.48e-14) + 2.7, -2.7d]
                 'sch': guess = [alog(4.48e-14) + 2.7, -2.7d, 1.]
              endcase 
           endif else begin
              ;; Equivalent width
              case ffit_form of 
                 'pow': guess = [alog(0.2554), -1.8d]
                 'exp': guess = [alog(3.71), -2.7d] 
                 'sch': guess = [alog(2.2), -1.8d, 1.d] 
              endcase 
           endelse 
        endelse                 ; no input
     endif else guess0 = guess  ; guess not set
     if ffit_form eq 'sch' then ffit_datanorm = guess[2]

     ;; ;;;;;;;
     ;; Initialize
     ;; Frequency distribution:
     ;; f(data) = exp(px[0]) * (data/data_0)^px[1]
     ;; or
     ;; f(data) = exp(px[0]) * exp(data/data_0 * px[1])
     ;; or 
     ;; f(data) = exp(px[0]) * (data/px[2])^px[1] * exp(-data/px[2])
     ;; where data is either column density or equivalent width (linear).
     ;; Cast with exp(px[0]) to make gradient of maximum likelihood
     ;; function (ln L) with respect to px[0] on order of gradient with
     ;; respect to px[1] (and px[2]).
     ;; Keep data within reasonable range of values by dividing by
     ;; "typical" (minimum, in this case) value. Especially useful for
     ;; f(N).
     if keyword_set(tcpu) then $
        print,'sdss_fitfxw: starting minf_conj_grad ',systime()
     rslt = dblarr(10L,ffit_ndim+1) 

     if keyword_set(plot) then begin
        angstrom = STRING("305B)   
        title = strarr(2)
        if ffit_datatype eq 'NCOLM' then title[0] = '!8N!X' $
        else title[0] = '!8W!X ('+angstrom+')'
        title[1] = '!8f!X('+title[0]+')'
        
        ;; Check that the setup correct by plotting initial guess 
        ;; _extra= includes ewbinsize=, ewlim=, silent=, bigewbin=,
        ;; /c_ech, _extra=extra
        ewclim = sdss_getewclim(ffit_cumx)
        fxw0 = sdss_calcfxw(civcand0, ffit_cumx, reform(ffit_cumx.zlim), $
                            dz=dz, final=final, /silent, n_per_bin=100, $
                            afp_flg=afp_flg, biasuser_fil=biasuser_fil, $
                            /noBAL, czn=(ffit_datatype eq 'NCOLM'))
        gdfxw = where(fxw0.numtot0 ne 0.)

        if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
        else dat_nrm = ffit_datanorm

        if ffit_datatype eq 'NCOLM' then begin
           yguess = sdss_calcfxwfit(10.^fxw0.ewcenter, option=ffit_form, $
                                    coeff=exp(guess[0]), alpha=guess[1],$
                                    datanorm=dat_nrm)
           
           xrng = [fxw0.ewcenter[0]-alog10(1-alog(10)*fxw0.sigewcenter[0,0]),$
                   alog10(ffit_intlim[2])] 
           if not finite(xrng[0]) then xrng[0] = 0.8*fxw0.ewcenter[0] ; 80%
        endif else begin
           yguess = sdss_calcfxwfit(fxw0.ewcenter, option=ffit_form, $
                                    coeff=exp(guess[0]), alpha=guess[1],$
                                    datanorm=dat_nrm)

           xrng = [fxw0.ewcenter[0]-fxw0.sigewcenter[0,0]*(ffit_form eq 'exp'),$
                   ffit_intlim[2]] 
        endelse 

        yrng = [min((fxw0.fxw[gdfxw]-fxw0.sigfxw[gdfxw,0])) > 1.e-30,$
                max(fxw0.fxw[gdfxw]+fxw0.sigfxw[gdfxw,1])]
        clr = getcolor(/load)
        if not keyword_set(psfil) then begin
           window, 1
           wset, 1
           plot, fxw0.ewcenter, fxw0.fxw, psym=3, symsize=2, color=clr.black,$
                 background=clr.white, xtitle=title[0], ytitle=title[1], $
                 title='50% Complete (green), Guess (red), Fit (blue)',/ylog,$
                 xlog=(ffit_form ne 'exp'),xrange=xrng,yrange=yrng,$
                 /xstyle,/ystyle,charsize=2
           oploterror,fxw0.ewcenter, fxw0.fxw, fxw0.sigewcenter[*,0], $
                      fxw0.sigfxw[*,0],$
                      /lobar,errcolor=clr.black,psym=3,thick=2
           oploterror,fxw0.ewcenter, fxw0.fxw, fxw0.sigewcenter[*,1], $
                      fxw0.sigfxw[*,1],$
                      /hibar,errcolor=clr.black,psym=3,thick=2
           oplot,[ewclim,ewclim],[1.e-10,1000],color=clr.limegreen,linestyle=3,$
                 thick=2                                                   ; 50% complete
           oplot,fxw0.ewcenter, yguess, color=clr.red, thick=2,linestyle=2 ; Guess
        endif                                                              ; psfil not set
;        stop
     endif                      ; /plot

     ;; _extra includes tolerance=, /use_deriv, /quadratic
     if keyword_set(ffit_debug) then begin
        print,''
        print,'entering initialization call'
     endif 
     minf_conj_grad, guess, f_min, conv_factor, func_name='sdss_fitfxw_logL',$
                     /initialize,_extra=extra
     if ffit_form eq 'sch' and exp(ffit_datanorm) gt ffit_intlim[2] then begin
        status = 1              ; badness
        goto, fail_exit
     endif else status = 0      ; default
     if keyword_set(ffit_debug) then print,'exiting initialization call'

     count = 0
     if ffit_ndim gt 2 then $
        rslt[count,*] = [exp(guess[0]),guess[1],guess[2],f_min] $
     else rslt[count,*] = [exp(guess[0]),guess[1],f_min]

     if not keyword_set(conv_factor_lim) then conv_factor_lim = 1.e-6


     ;; ;;;;;;;
     ;; Iteratively fit
     while (conv_factor gt conv_factor_lim) do begin
        if (count mod 10) eq 0 and keyword_set(tcpu) then $
           print,'sdss_fitfxw: minf_conj_grad loop = '+$
                 string(count,format='(i4)'),systime()
        if keyword_set(ffit_debug) then begin
           print,''
           print,'entering call',count
        endif 
        minf_conj_grad, guess, f_min, conv_factor, $
                        func_name='sdss_fitfxw_logL',$
                        _extra=extra ; /use_deriv, tolerance=[sqrt(1.e-y)]
        if ffit_form eq 'sch' and exp(ffit_datanorm) gt ffit_intlim[2] then begin
           status = count+1           ; badness
           goto, fail_exit
        endif 
        if keyword_set(ffit_debug) then print,'exiting call',count
        
        count = count + 1
        if ffit_ndim gt 2 then $
           rslt[count,*] = [exp(guess[0]),guess[1],guess[2],f_min] $
        else rslt[count,*] = [exp(guess[0]),guess[1],f_min]

        ihi = n_elements(rslt[*,0])
        if ihi-1 eq count+1 then begin
           ;; Expand
           tmp = dblarr(2*ihi,ffit_ndim+1)
           tmp[0:ihi-1,0] = rslt[*,0]
           tmp[0:ihi-1,1] = rslt[*,1]
           tmp[0:ihi-1,2] = rslt[*,2]
           if ffit_ndim gt 2 then tmp[0:ihi-1,3] = rslt[*,3]
           rslt = tmp
        endif 

;        stop
     endwhile 


     ;; ;;;;;;;
     ;; Save the normal coeff * N^alpha formalism
     fail_exit:
     ;; Set to store and use later
     ffit_coeff = guess[0]
     ffit_alpha = guess[1]
     if ffit_form eq 'sch' then ffit_datanorm = guess[2]
     ffit_sigdatanorm = replicate(!values.d_nan,2)

     fitstrct = {dblt_name:ffit_dblt.ion, $
                 datatype:ffit_datatype, $ ; store following to use as infit_fil
                 fittype:ffit_form, $      
                 data:ffit_data, $
                 sigdata:ffit_sigdata, $
                 flag:ffit_flg, $
                 afp_flg:keyword_set(afp_flg), $
                 xdata:ffit_xdata, $
                 datanorm:ffit_datanorm, $ ; fit param for Schechter
                 sigdatanorm:ffit_sigdatanorm, $
                 intlim:ffit_intlim, $ ; [lo, "saturated", hi] 
                 zrng:[min(ffit_civstrct.(ztag)[0]),$
                       median(ffit_civstrct.(ztag)[0]),$
                       max(ffit_civstrct.(ztag)[0])], $ ; after trim
                 datalim:[min(ffit_data,max=mx),mx], $
                 chi_sqr:fltarr(3), $
                 coeff: exp(ffit_coeff),$ ; best-fit results and error
                 sigcoeff: dblarr(2), $
                 alpha: ffit_alpha, $
                 sigalpha: dblarr(2), $
                 dndx:0., $     ; intlim[[0,2]]
                 sigdndx:fltarr(2,/nozero), $
;                 omciv: sdss_fitfxw_omciv(guess,dblt_name=ffit_dblt), $ 
;                 sigomciv: dblarr(2), $ ; defined by following locations
;                 sigcoeff_omciv:dblarr(2), $
;                 sigalpha_omciv:dblarr(2), $
                 logL:-f_min, $ ; maximum likelihood
                 d_ks:0.d, $    ; K-S test results
                 prob_ks:0.d, $
                 conv_factor:conv_factor, $ ; convergence factor
                 rslt_iter:rslt, $           ; info about iterations
                 status:status $
                }
     if not keyword_set(status) then begin
        if ffit_form eq 'sch' then fitstrct.datanorm = exp(fitstrct.datanorm)
        fitstrct.dndx = sdss_calcdndxfit(fitstrct.intlim[[0,2]], fitstrct_fil=fitstrct)
        ;; error later
        
        ;; ;;;;;;;
        ;; Chi^2 for just the data used in the fit
        fitstrct.chi_sqr = sdss_calcfxwchisq(ffit_civstrct,ffit_cumx,$
                                             dz=dz,final=final, $
                                             fitstrct_fil=fitstrct, /silent,$
                                             n_per_bin=100L) ; roughly
     endif
     
     ;; Write-out (safety)
     mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil
     if keyword_set(status) then begin
        spawn,'gzip -f '+ostrct_fil
        print,'sdss_fitfxw: cannot converge; exit status = ',strtrim(status,2)
        if keyword_set(guess0) then guess = guess0 ; restore original

        return
     endif
     
     ;; Set next steps
     ksonly = 1
     if not keyword_set(fitonly) then erronly = 1
  endelse                       ; max-L


  if keyword_set(psfil) then begin
     x_psopen,psfil,/maxs,/portrait 
     !p.multi = [1,1,1]
     !x.margin = [8.7,1.5]      ; left and right border
     !y.margin = [3.2,2]        ; bottom and top border
  endif


  if keyword_set(ksonly) then begin
     ;; ;;;;;;;
     ;; KS Test
     ;; Normalize function:
     ;; px[0] = ln(1 + px[1]) - ln( (var_max^(1+px[1])/var_min^px[1]) -
     ;;                         var_min)
     ;; or 
     ;; px[0] = ln(var_0 / px[1] * 1/(exp(px[1]*var_max) - exp(px[1]*var_min))
     ;; Only for unsaturated features
     if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
     else dat_nrm = ffit_datanorm

     ffit_coeff = sdss_calcdndxfit(ffit_intlim[0:1], option=ffit_form, coeff=1., $
                                   alpha=ffit_alpha, datanorm=dat_nrm)
     ffit_coeff = alog(1./ffit_coeff)
     prs = strsplit(ostrct_fil,'/',/extract,count=nprs)
     ttl = strmid(prs[nprs-1],0,strpos(prs[nprs-1],'.',/reverse_search))
     wgt = 1./ffit_xdata[ffit_unsat]
     if keyword_set(ffit_ewseam) and ffit_ewseam lt ffit_intlim[1] then begin
        ;; Need to actually normalize the weights
        hi = where(ffit_data[ffit_unsat] ge ffit_ewseam, complement=lo)
        if hi[0] ne -1 and lo[0] ne -1 then begin
           mxhi = max(ffit_xdata[ffit_unsat[hi]])
           mxlo = max(ffit_xdata[ffit_unsat[lo]])
           wgt[lo] = wgt[lo] * mxlo/mxhi
        endif
;        stop
     endif 
     ksone_wght, ffit_data[ffit_unsat], 'sdss_fitfxw_kscumf', d, ksprob, $
                 wght = wgt, plot=plot, title=ttl, xlog=(ffit_form ne 'exp'), _extra=extra

     ;; Save results
     fitstrct.d_ks = d
     fitstrct.prob_ks = ksprob 

     ;; Restore common values 
     ffit_coeff = alog(fitstrct.coeff)

     ;; Write-out again (safety)
     mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil

     if keyword_set(tcpu) then $
        print,'sdss_fitfxw: finished KS Test ',systime()
  endif                         ; /ksonly


  if keyword_set(plot) then begin
     ;; To plot:
     ;; surface,logL_grid,coeff_grid,alpha_grid,ax=50.,charsize=2
     if not keyword_set(psfil) then wset, 1
     clr = getcolor(/load)
     plot, fxw0.ewcenter, fxw0.fxw, psym=3, symsize=2, color=clr.black,$
           background=clr.white, xtitle=title[0], ytitle=title[1], $
           title='50% Complete (green), Guess (red), Fit (blue)',/ylog,$
           xlog= (ffit_form ne 'exp'),xrange=xrng,yrange=yrng,$
           /xstyle,/ystyle,charsize=2
     oploterror,fxw0.ewcenter, fxw0.fxw, fxw0.sigewcenter[*,0], $
                fxw0.sigfxw[*,0],$
                /lobar,errcolor=clr.black,psym=3,thick=2
     oploterror,fxw0.ewcenter, fxw0.fxw, fxw0.sigewcenter[*,1], $
                fxw0.sigfxw[*,1],$
                /hibar,errcolor=clr.black,psym=3,thick=2
     oplot,fxw0.ewcenter, yguess, color=clr.red, thick=2,linestyle=2 ; Guess
     oplot,[ewclim,ewclim],[1.e-10,1000],color=clr.limegreen,linestyle=3,$
           thick=2              ; 50% complete

     if ffit_datatype eq 'NCOLM' then $
        lbl = string(fitstrct.chi_sqr[0]/fitstrct.chi_sqr[1],$
                     fitstrct.chi_sqr[2],alog10(fitstrct.coeff), fitstrct.alpha, $
                     format='("!9C!X!E2!N = ",f6.3,"; !8P!X!D!9C!X!E2!N!N = ",f6.3,"; log !8k!X = ",f7.3,"; !9a!X = ",f6.3)') $
     else $
        lbl = string(fitstrct.chi_sqr[0]/fitstrct.chi_sqr[1],$
                     fitstrct.chi_sqr[2],fitstrct.coeff, fitstrct.alpha, $
                     format='("!9C!X!E2!N = ",f6.3,"; !8P!X!D!9C!X!E2!N!N = ",f6.3,"; !8k!X = ",f5.3,"; !9a!X = ",f6.3)')
     if ffit_Form eq 'sch' then begin
        if ffit_datatype eq 'NCOLM' then $
           lbl = lbl + string(alog10(fitstrct.datanorm),format='("; log N* = ",f6.3)') $
        else lbl = lbl + string(fitstrct.datanorm,format='("; W* = ",f6.3)')
     endif 

     if ffit_datatype eq 'NCOLM' then $
        yfit0 = sdss_calcfxwfit(10.^fxw0.ewcenter,fitstrct=fitstrct) $
     else yfit0 = sdss_calcfxwfit(fxw0.ewcenter,fitstrct=fitstrct)
     
     oplot, fxw0.ewcenter, yfit0, color=clr.blue, thick=2
     xyouts, 0.5, 0.15, lbl, alignment=0.5,/normal,charsize=1.5
  endif                         ; /plot



  if keyword_set(erronly) then begin
     ;; ;;;;;;;
     ;; Errors
     ;; Grid of coeff and alpha
     alpha_range = 1.
     coeff_range = 1.
     if not keyword_set(nbin_alpha) then nbin_alpha = 1000L
     if not keyword_set(nbin_coeff) then nbin_coeff = 1000L
     if not keyword_set(dcoeff) then dcoeff = 5.e-4 ;coeff_range/nbin_coeff ;
     if not keyword_set(dalpha) then dalpha = 5.e-4 ;alpha_range/nbin_alpha ;

     ;; Ensure grids have max values "exactly" (compiled in log space)
     ;; Also know that likelihood turns over faster at high ceff 
     ;; and high alpha
     ;; I'm not so convinced that uniform spacing in log space
     ;; is so good... perhaps why the error est. takes so long
     coeff_grid = ffit_coeff + dcoeff * (lindgen(nbin_coeff) - nbin_coeff/2)
     alpha_grid = ffit_alpha + dalpha * (lindgen(nbin_alpha) - nbin_alpha/2)

     ;; Set array of [coeff,alpha] to find likelihood function value
     px = dblarr(nbin_coeff*nbin_alpha,ffit_ndim,/nozero)
     for ii=0L,nbin_alpha-1 do begin
        ;; Keep grouped by alpha (and datanorm) for quicker integration
        indx = lindgen(nbin_coeff)+ii*nbin_coeff
        px[indx,0] = coeff_grid
        px[indx,1] = alpha_grid[ii]
     endfor                     ; loop ii=nbin_alpha
     if ffit_form eq 'sch' then begin
        ;; Now need to dupliate the whole current px array and set
        ;; datanorm 
        if not keyword_set(nbin_norm) then nbin_norm = 20L
        if not keyword_set(dnorm) then dnorm = 1.e-3
        norm_grid = ffit_datanorm + dnorm * (lindgen(nbin_norm) - nbin_norm/2)
        px0 = px
        px = dblarr(nbin_coeff*nbin_alpha*nbin_norm,ffit_ndim,/nozero)
        indx0 = lindgen(nbin_coeff*nbin_alpha) ; px0 index
        for ii=0L,nbin_norm-1 do begin
           indx = ii*(nbin_coeff*nbin_alpha) + indx0
           px[indx,0] = px0[*,0]
           px[indx,1] = px0[*,1]
           px[indx,2] = norm_grid[ii]
        endfor                  ; loop ii=nbin_norm

        if ffit_form eq 'sch' then dat_nrm = exp(px[*,2]) $
        else dat_nrm = px[*,2]
     endif else begin
        norm_grid = 0
        nbin_norm = 0
     endelse 

     if ffit_form eq 'sch' then dat_nrm = exp(ffit_datanorm) $
     else dat_nrm = ffit_datanorm
     dndx_grid = sdss_calcdndxfit(fitstrct.intlim[[0,2]],$
                                  option=ffit_form,$
                                  coeff=exp(px[*,0]),alpha=px[*,1],$ ; <--- exp(k)? 
                                  datanorm=dat_nrm,ndim=ffit_ndim)
     
     ;; Evaluate (could be expensive)
     if keyword_set(tcpu) then $
        print,'sdss_fitfxw: calculating grid of log L ',systime()
     logL = -sdss_fitfxw_logL(px) ; want max

     
     ;; Does all the error analysis
     ;; _extra= includes sigma=, delta=, mindlogl=, maxdlogl=
     dlogl_cl = sdss_fitfxw_calcfiterr(fitstrct, cl, logl=logl,$
                                       nalpha_bin=nbin_alpha,$
                                       ncoeff_bin=nbin_coeff,$
                                       nnorm_bin=nbin_norm,$
                                       alpha_grid=alpha_grid,$
                                       coeff_grid=exp(coeff_grid),$ ; linear
                                       norm_grid=norm_grid,$
                                       ofitstr=ofitstr,oerrstr=errstrct,$
                                       debug=ffit_debug,$
                                       _extra=extra)
     fitstrct = ofitstr
     
     ;; include dN/dX error estimate by just the stddev() within the
     ;; 1-sigma ellipse which isn't quite right b/c the surface
     ;; is skewed but don't want to bin and fit a Gaussian or something
     errstrct = create_struct(errstrct,'dndx_grid',dndx_grid)
     fitstrct.sigdndx = stddev(dndx_grid[errstrct.surf]) 

     ;; Reset common block values
     ffit_coeff = alog(fitstrct.coeff)
     ffit_alpha = fitstrct.alpha
     if ffit_form eq 'sch' then ffit_datanorm = alog(fitstrct.datanorm) $
     else ffit_datanorm = fitstrct.datanorm


     ;; Accepted false-positive adjustment
     if keyword_set(afp_flg) then begin
        newfitstrct = sdss_adjafpfxwfit(fitstrct, errstrct, newerrstrct=newerrstrct,$
                                       biasuser_fil=biasuser_fil)

        ;; Adjust the dN/dX
        newfitstrct.dndx = sdss_calcdndxfit(newfitstrct.intlim[[0,2]], $
                                            fitstrct_fil=newfitstrct)
        ;; since error was estimated as a stddev() within the 1-sigma
        ;; surface, just leave it (smaller dN/dX best but same error
        ;; is about right)
        fitstrct = newfitstrct
        errstrct = newerrstrct
     endif 

     if keyword_set(plot) then begin
        ;; Error contours
        if not keyword_set(psfil) then begin
           window, 2
           wset, 2              ; new window
        endif
        logL_grid = reform(logL,nbin_coeff,nbin_alpha)

        contour,logL_grid-fitstrct.logL,errstrct.coeff_grid,$
                errstrct.alpha_grid,$
                levels=[errstrct.dlogl_lim],background=clr.white,$
                c_colors=clr.black,color=clr.black, xtitle='!8k!X',$
                ytitle='!9a!X',charsize=2

        oplot,[fitstrct.coeff],[fitstrct.alpha],psym=4
        oplot,fitstrct.coeff-[fitstrct.sigcoeff[0],fitstrct.sigcoeff[0]],$
              [errstrct.alpha_grid[0],errstrct.alpha_grid[nbin_alpha-1]],$
              color=clr.black,linestyle=1
        oplot,fitstrct.coeff+[fitstrct.sigcoeff[1],fitstrct.sigcoeff[1]],$
              [errstrct.alpha_grid[0],errstrct.alpha_grid[nbin_alpha-1]],$
              color=clr.black,linestyle=1
        oplot,[errstrct.coeff_grid[0],errstrct.coeff_grid[nbin_coeff-1]],$
              fitstrct.alpha-[fitstrct.sigalpha[0],fitstrct.sigalpha[0]],$
              color=clr.black,linestyle=1
        oplot,[errstrct.coeff_grid[0],errstrct.coeff_grid[nbin_coeff-1]],$
              fitstrct.alpha+[fitstrct.sigalpha[1],fitstrct.sigalpha[1]],$
              color=clr.black,linestyle=1
        
     endif                      ; /plot

  endif                         ; /erronly

  
  if keyword_set(psfil) then begin
     x_psclose
     print,'sdss_fitfxw: created ',psfil
  endif 


  ;; Write-out (for real); duplicate of sdss_fitfxw_calcfiterr()
  mwrfits,fitstrct,ostrct_fil,/create,/silent ; fit for infit_fil
  mwrfits,errstrct,ostrct_fil,/silent 
  spawn,'gzip -f '+ostrct_fil
  print,'sdss_fitfxw: created ',ostrct_fil

  if keyword_set(guess0) then guess = guess0 ; restore original

  if keyword_set(tcpu) then $
     print,'sdss_fitfxw: finished ',systime()
end
