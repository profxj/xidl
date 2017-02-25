;+
; NAME:
;   first_mcmc
; PURPOSE:
;   Make a Markov Chain Monte Carlo chain.
;
; INPUTS:
;   seed       - random number seed for randomu()
;   step_func  - function that takes a step in parameter space
;   like_func  - function that computes the likelihood
;   nstep      - number of links
;   pars       - initial parameters (can be an array or structure)
; KEYWORDS:
;   log        - like_func returns log_e(likelihood), not straight
;                likelihood
; OUTPUTS:
;   pars       - array of parameters, sorted from most to least likely
;   like       - array of likelihoods for the pars
; BUGS:
; REVISION HISTORY:
;   2013-06-23  -- JXP
;;   Based (copied) on hogg_mcmc
;;   Applied to the f(N) model of Omeara+13 as a test
;-
pro x_first_mcmc, seed,step_func,like_func,nstep,pars,like,log=log
  compile_opt strictarr

npars= n_elements(pars)
oldpars= pars
pars= replicate(pars[0],npars,nstep)
like= dblarr(nstep)
for ii=0L,nstep-1L do begin
    x_first_mcmc_step, seed,oldpars,oldlike,step_func,like_func,newpars,newlike, $
      log=log
    pars[*,ii]= newpars
    like[ii]= newlike
    oldpars= newpars
    oldlike= newlike
    if (ii MOD 100) then begin
       junk = check_math()
       print, $
          strjoin(strarr(21)+string(byte(8)),''), $
          'JXP_MCMC: ',100L*ii/nstep,' percent', $
          format= '($,A,A,I2.2,A)'
    endif
endfor
print, strjoin(strarr(21)+string(byte(8)),'')+'JXP_MCMC: done      '
sindx= reverse(sort(like))
pars= pars[*,sindx]
like= like[sindx]
return
end
