;+
; NAME:
;   x_mcmc_chk_chain
; PURPOSE:
;   Generates a series of plots to explore an MCMC chain
;
; INPUTS:
;   chain_file  
;
; KEYWORDS:
; OUTPUTS:
; BUGS:
; REVISION HISTORY:
;   2013-09-17  -- JXP
;-
pro x_mcmc_chk_chain, chain_file, BURN=burn, PSFILE=psfile

  compile_opt strictarr

  if not keyword_set(CHAIN_FILE) then return

  chain = xmrdfits(chain_file)
  sz = size(chain, /dim)
  if not keyword_set(BURN) then burn = 4000L
  if n_elements(sz) LT 2 then nchain = 0 else nchain = sz[2]
  print, 'Nchain = ', nchain

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; PLOT
  if not keyword_set(PSFILE) then psfile = 'x_mcmc_chk_chain.ps'
  x_psopen, psfile, /maxs
  !p.multi=[0,2,2]
  clr = getcolor(/load)
  clrs = [clr.black, clr.blue, clr.green, clr.gray, clr.purple, clr.brown]

  ;; ;;;;;;;;;;;;;;;;;;;;;;
  ;;   CHAIN WALKING
  xmrg = [8,1]
  ymrg = [4.0,0]

  nplt = sz[0]/2
  for qq=0L,nplt-1 do begin

     ;; Axes
     x_idx = qq
     y_idx = qq + nplt
     xtit = x_padstr(strtrim(x_idx,2), 2, '0', /rev)
     ytit = x_padstr(strtrim(y_idx,2), 2, '0', /rev)

     xrng = [min(chain[x_idx,*,*], max=mx), mx]
     yrng = [min(chain[y_idx,*,*], max=mx), mx]
  
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           xthi=thk, ythi=thk

     ;; Loop on chains
     for kk=0L,nchain-1 do begin
        ;; POST BURN
        oplot, chain[x_idx,burn:*,kk], chain[y_idx,burn:*,kk], color=clrs[kk], psym=-1, symsiz=0.2, $
               thick=1
        ;; BURN
        oplot, chain[x_idx,0:burn-1,kk], chain[y_idx,0:burn-1,kk], color=clr.red, psym=-1, symsiz=0.2, $
               thick=5
     endfor

  endfor

  ;; ;;;;;;;;;;;;;;;;;;;;;;
  ;;   CHAIN WALKING

  nplt = sz[0]
  xtit = 'Chain#'
  xrng = [0, sz[1]-1]
  for qq=0L,nplt-1 do begin

     ;; Axes
     ytit = x_padstr(strtrim(qq,2), 2, '0', /rev)
     yrng = [min(chain[qq,*,*], max=mx), mx]
  
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, xtitle=xtit, $
           ytitle=ytit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           xthi=thk, ythi=thk

     for kk=0L,nchain-1 do begin
        ;; BURN
        oplot, lindgen(burn), chain[qq,0:burn-1,kk], color=clr.red,  thick=5
        ;; POST BURN
        oplot, burn+lindgen(sz[1]-burn), chain[qq,burn:*,kk], color=clrs[kk], thick=3
     endfor

  endfor

  ;; Close Ps
  if keyword_set( PSFILE ) and not keyword_set(NOPLOT) then x_psclose
  print, 'tstmcmc_fn_z25:  All done!'
  !p.multi=[0,1,1]

;  stop
       
  return
end
      
