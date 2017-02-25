;+
; NAME:
;   x_mcmc_mhastings
; PURPOSE:
;   Make a Markov Chain Monte Carlo chain using the Metropolis
;   Hastings algorithm
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
;   2013-08-28  -- JXP
;;   Based (copied) on hogg_mcmc
;-
pro x_mcmc_mhastings, seed,step_func,like_func,nstep,pars,like,log=log, $
                      VERBOSE=verbose, OUTFIL=outfil
  compile_opt strictarr

  junk = CHECK_MATH(0, 1)  ;; suppress error messages

  npars= n_elements(pars)
  oldpars= pars

  ;; Setup the Chain
  pars= replicate(pars[0],npars,nstep)
  like= dblarr(nstep)

  ;; Main LOOP
  for ii=0L,nstep-1L do begin
     x_mcmc_mhastings_step, seed,oldpars,oldlike,step_func,like_func,newpars,newlike, $
                        log=log
     pars[*,ii]= newpars
     like[ii]= newlike
     oldpars= newpars
     oldlike= newlike
     if (ii MOD 100) and keyword_set(VERBOSE) then begin
        print, $
           strjoin(strarr(21)+string(byte(8)),''), $
           'JXP_MCMC: ',100L*ii/nstep,' percent', $
           format= '($,A,A,I2.2,A)'
     endif
  endfor
  if keyword_set(VERBOSE) then print, strjoin(strarr(21)+string(byte(8)),'')+'JXP_MCMC: done      '

  ;; Write?
  if keyword_set(OUTFIL) then begin
     mwrfits, pars, outfil, /create
     mwrfits, like, outfil
  endif

  junk = CHECK_MATH(0, 0)  ;; suppress error messages
  return
end
