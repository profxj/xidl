;+
; NAME:
;   x_mcmc_mhastings_parallel
; PURPOSE:
;   Launch a series of NPROC calls to x_mcmc_mhastings
; INPUTS:
;   seed       - random number seed for randomu()
;   step_func  - function that takes a step in parameter space
;   like_func  - function that computes the likelihood
;   nstep      - number of links per Chain
;   pars       - initial parameters (can be an array or structure)
; KEYWORDS:
;   log        - like_func returns log_e(likelihood), not straight
;                likelihood
;   NPROC=    Number of processors to run on
;   NCHAIN=    Number of chains to make
; OUTPUTS:
;   newpars    - new parameters
;   newlike    - new likelihood
; BUGS:
; REVISION HISTORY:
;   2013-08-28  JXP
;-
pro x_mcmc_mhastings_parallel, seed,step_func,like_func,nstep,pars,like, $
                               log=log, NPROC=nproc, NCHAIN=nchain, INIT_FIL=init_fil, $
                               REWRITE=rewrite, INIT_FLG=init_flg, INIT_DIR=init_dir
  compile_opt strictarr

  if not keyword_set(NPROC) then nproc = 2L
  if not keyword_set(NCHAIN) then nchain = nproc

  npars= n_elements(pars)
  oldpars= pars


  close, /all
  count = 0
  tot_count = 0
  if not keyword_set(REWRITE) then begin
     for qq=0L,nchain-1 do begin
        ;; Counters (for NPROC)
        tot_count = tot_count + 1
        count=count+1
        
        ;; Init (dum +/-1 offset)
        initp = pars + (1.-randomu(seed, npars)*2)
        
        ;; Generate the dummy script
        lbl = x_padstr(strtrim(qq,2),2,'0',/rev)
        dumfil = 'mcmc_mhastings_'+lbl+'.tmp'
        outfil = 'mcmc_mhastings_'+lbl+'.fits'
        ifil = qq+33
        openw, ifil, dumfil
        dums = -1*fix(10000L*randomu(seed))
        printf, ifil, 'seed = '+strtrim(dums,2)+'L'
        printf, ifil, 'nstep = '+strtrim(nstep,2)+'L'
        printf, ifil, 'outfil = '+''''+outfil+''''
        printf, ifil, 'step_func = '+''''+step_func+''''
        printf, ifil, 'like_func = '+''''+like_func+''''
        printf, ifil, 'log = ', keyword_set(LOG)
        spars = 'pars=['
        for jj=0L,npars-2 do spars=spars+string(initp[jj],format='(e12.5)')+','
        spars=spars+string(initp[npars-1],format='(e12.5)')+']'
        printf, ifil, spars
        ;stop
        ;; INIT
        if keyword_set(INIT_FIL) then begin
           if keyword_set(INIT_DIR) then printf, ifil, 'cd, '+init_dir+', curr=curr'
           if keyword_set(INIT_FLG) then init_s = init_fil+', '+init_flg else init_s=init_fil
           printf, ifil, init_s
           if keyword_set(INIT_DIR) then printf, ifil, 'cd, curr'
        endif
        ;; Main Call
        printf, ifil, 'x_mcmc_mhastings, seed, step_func, like_func, nstep, pars, OUTFIL=outfil, LOG=log'
        close, ifil
        ;; Launch the process
        if count LT NPROC AND TOT_COUNT NE nchain then sbkg = '&' else begin
           sbkg = ''
           count = 0
        endelse
        spawn, 'idl < '+dumfil+sbkg
     endfor
  endif

  ;; Wait for completion
  flg = 0
  while (not flg) do begin
     print, 'x_mcmc_mhastings_parallel: Waiting for completion..'
     wait, 100
     flg = 1
     for qq=0L,nchain-1 do begin
        lbl = x_padstr(strtrim(qq,2),2,'0',/rev)
        outfil = 'mcmc_mhastings_'+lbl+'.fits'
        flg = flg * x_chkfil(outfil)
     endfor
  endwhile
  
  ;; Repack
  pars= replicate(pars[0],npars,nstep, nchain)
  like= dblarr(nstep, nchain)
  for qq=0L,nchain-1 do begin
     lbl = x_padstr(strtrim(qq,2),2,'0',/rev)
     outfil = 'mcmc_mhastings_'+lbl+'.fits'
     ;; Read
     tpars = xmrdfits(outfil,/silen)
     tlike = xmrdfits(outfil,1,/silen)
     ;; Push
     pars[*,*,qq] = tpars
     like[*,qq] = tlike
     ;; Remove (for later runs)
     spawn, '\rm '+outfil
  endfor

  return
end
