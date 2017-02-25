pro hires_qcomb, indvfil, outfil

 ;; Loop on indvfil
 sz = size(indvfil,/dimensions)

 for qq=0L,sz[1]-1 do begin
     ;; Open all files
     spos = strpos(indvfil[0,qq], '_F')
     errfil = strmid(indvfil[0,qq],0,spos)+'_E.fits'
     d1 = x_readspec(indvfil[0,qq],head=head,wav=w1,$
                     FIL_SIG=errfil, SIG=e1, npix=npix)

     ;; Final spec
     if qq EQ 0 then begin
         fin_fx = fltarr(npix)
         fin_sig = fltarr(npix)
     endif

     ;; Fill up the big array
     fx = fltarr(npix,sz[0])
     sig = fltarr(npix,sz[0])
     scale = replicate(1., sz[0])
     snr = fltarr(sz[0])
     fx[*,0] = d1
     sig[*,0] = e1

     ;; Scale
     for jj=1L,sz[0]-1 do begin
         spos = strpos(indvfil[jj,qq], '_F')
         errfil = strmid(indvfil[jj,qq],0,spos)+'_E.fits'
         dat = x_readspec(indvfil[jj,qq],FIL_SIG=errfil, SIG=dsig, npix=npix)
         a = where(d1 GT 0. AND dat GT 0.)
         fitstr = x_setfitstrct(FUNC='LEGEND', NORD=7L, $
                                LSIG=2.5, HSIG=2.5, NITER=4,$
                                /flgrej, MINPT=1000L)
         fit = x_fitrej(w1[a], dat[a]/d1[a], fitstr=fitstr)
         fx[*,jj] = dat / x_calcfit(w1, fitstr=fitstr)
         sig[*,jj] = dsig / x_calcfit(w1, fitstr=fitstr)
     endfor

     ;; SNR
     for jj=0L,sz[0]-1 do begin
         gd = where(sig[*,jj] GT 0.)
         snr[jj] = x_mode(fx[gd,jj]/sig[gd,jj])
     endfor
     print, 'SNR = ', snr

     ;; Combine
     x_combspec, fx, (sig>0.)^2, ffx, ffv, SCALE=scale, SNR=snr
     
     ;; Write into fin_fx
     gd = where(ffv GT 0.)
     fin_fx[gd] = ffx[gd]
     fin_sig[gd] = sqrt(ffv[gd])
 endfor

 ;; Output
 print, 'comb_spec: Creating -- ', outfil
 mwrfits, fin_fx, outfil, head, /create
 spos = strpos(outfil, '_F')
 errfil = strmid(outfil,0,spos)+'_E.fits'
 mwrfits, fin_sig, errfil, head, /create

 print, 'comb_spec: All done!'
 return
end
