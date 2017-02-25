;+
; NAME:
;   x_mcmc_chain_stats
; PURPOSE:
;   Analyze a Markov Chain Monte Carlo chain using simple stats
;     -- Best values correspond to the greatest Likelihood
;
; INPUTS:
;   chain_file       - random number seed for randomu()
;
; KEYWORDS:
;   BURN_FRAC= -- Burn fraction (e.g 0.25) to implement.  Default=0.3
; OUTPUTS:
;   strct       - Structure of stats
; BUGS:
; REVISION HISTORY:
;   2013-09-16  -- JXP
;-
function x_mcmc_chain_stats, chain_file, BURN_FRAC=burn_frac

  ;; Read
  chain = xmrdfits(chain_file)
  like = xmrdfits(chain_file,1)

  sz = size(chain, /dim)
  nparm = sz[0]
  if n_elements(sz) LT 2 then nchain = 1 else nchain = sz[2]
  
  strct = { $
          fil: chain_file, $
          nparm: nparm, $
          best_p: fltarr(nparm), $
          sig: fltarr(nparm,2) $
          }

  ;; Burn
  if not keyword_set(BURN_FRAC) then burn_frac = 0.3 
  burn = round(burn_frac*sz[1]) 
  chain = chain[*,burn:*,*]    ;; Burn
  like = like[burn:*,*]    ;; Burn
  sz = size(chain, /dim)

  ;; Analyze the chain
  if n_elements(sz) EQ 3 then begin
     chain = reform(chain, sz[0], sz[1]*sz[2])
     like = reform(like, sz[1]*sz[2])
  endif
  sz = size(chain, /dim)

  ;; Maximize
  mx = max(like, imx)
  strct.best_p = chain[*,imx]

  ;; Confidence limits (68%)
  if not keyword_set(CL) then cl = 0.683
  for qq=0L,nparm-1 do begin
     ;; Sort
     srt = sort(chain[qq,*])
     ;; Simple CDF
     idx1 = srt[round(sz[1]*(1-CL)/2.)]
     lowv = (chain[qq,*])[idx1]
     strct.sig[qq,0] = abs(strct.best_p[qq]-lowv)
     idx2 = srt[round(sz[1]*(1.-((1-CL)/2.)))]
     hiv = (chain[qq,*])[idx2]
     strct.sig[qq,1] = abs(hiv-strct.best_p[qq])
  endfor

  return, strct
end
