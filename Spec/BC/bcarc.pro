	function bcarc, image, invvar, hdr, guess=guess, tset=t3set, $
    goodrows=goodrows, linefile=linefile, func=func

   ndim = size(image,/n_dim)
   ncol = (size(image,/dimen))[0]
   nrow = (size(image,/dimen))[1]   

; blah

   cleanimage = djs_maskinterp(image, invvar LE 0, iaxis=1)
   sum = total(cleanimage,1)
   s   = sort(sum)
   s90 = s[nrow*9/10]
   top = sum[s90]

   if NOT keyword_set(goodrows) then goodrows = where(sum GT 0.5 * top)
   ngoodrows = n_elements(goodrows)
   if goodrows[0] EQ -1 then return, 0
   
   tracemask           = lonarr(nrow) 
   tracemask[goodrows] = 1
   traceivar           = invvar * (replicate(1,ncol) # tracemask)

   thresh = 2.0 * median(image[*,goodrows])
   xcen = trace_crude(image, traceivar, yset=ycen, thresh=thresh, $
            radius=1., ystart=s90) 
   xy2traceset, ycen[goodrows,*], xcen[goodrows,*], xset, ncoeff=3, $
         xmin=0, xmax=nrow-1.0, yfit=xfit, /silent

 
;      help,where(strpos(hdr,'300') GT -1)
;      help,where(strpos(hdr,'600') GT -1)
;      help,where(strpos(hdr,'1200') GT -1)


;;;;;;;;;;;
;   Parse the Linelist
;;;;;;;;;;;

   if NOT keyword_set(linefile) then begin
       linefile='bclowdisp.dat'
       ; Check for linefile
       a = findfile(linefile, count=count)
       if count EQ 0 then begin
           print, linefile, ' does not exist!'
           return, -1
       endif
   endif


   readcol, linefile, lwave, linten, lquality, $
                            lion, lion2, format='D,F,A,A,A' 

   good = where(lion EQ 'Ne' OR lion EQ 'Ar' OR lion EQ 'He', nlamps)
   if nlamps EQ 0 then return, 0
   linelist = replicate({ wave : 0.0, fit : 0.0, $
                         inten : 0.0, q : '', ion : '' }, nlamps)
   
   linelist.wave    = lwave[good]
   linelist.inten   = linten[good]
   linelist.q       = lquality[good]
   linelist.ion     = lion[good]

   if ndim EQ 2 then comp = djs_median(image[*,goodrows],2) else comp = image

   gratings = [300,600,1200]
   centraldisp = [5.44, 2.78, 1.42] 
   best        = fltarr(3)
   place       = lonarr(3)
   bestdisp    = fltarr(3)
   maxwave  = max(linelist.wave) < 10000.0
   minwave  = min(linelist.wave) > 3000.0

;;;;;;;;;;;
;  Determine grating from arc-lines
;;;;;;;;;;;

   for i=0,n_elements(gratings)-1 do begin
    bestcorr = 0.0
    for disp=centraldisp[i]-0.05, centraldisp[i]+0.05, 0.01 do begin

       nmodel = (maxwave - minwave)/disp
     
       ; spectrum reds red to blue
       modelwave = maxwave - findgen(nmodel)*disp
       modelflux = modelwave*0.0

       ; add simple lines

       modelres1 = ((maxwave - linelist.wave)/disp) 
       modelfrac = modelres1 - long(modelres1)

       modelres2 = [[modelres1-1],[modelres1],[modelres1+1]]
  
       modelflux[modelres2] = $
       modelflux[modelres2] + linelist.inten # replicate(1,3)

;       modelflux[modelres1+1] = modelflux[modelres1+1] + $
;           linelist.inten*(1.0-modelfrac)

       nmodel = n_elements(modelwave) 
       pad = nmodel - ncol
       lag = lindgen(nmodel/2)* 2

       comp1 = comp

       if pad LE 0 then begin
         lag = -1 * lag
         if pad LT 0 then modelflux = [modelflux, fltarr(-pad)]
       endif else comp1 = [comp1,fltarr(pad)]

       hmm = c_correlate(comp1,modelflux, lag)
;       hmm_fft=float(fft(fft(comp1)*fft(modelflux),1))
       
       thisbest = max(hmm,pl)
       if thisbest GT bestcorr then begin
         bestcorr = thisbest
         best[i] = bestcorr
         place[i] = lag[pl]
         bestdisp[i] = disp
       endif
    endfor
    print, gratings[i], bestdisp[i], best[i], place[i]
  endfor

   bestcorr = max(best, winner)

   if bestcorr LT 0.4 then begin
     print, 'Did not find a strong correlation, try another line list'
     return, 0
   endif

;  iterate the best fit
   disp = bestdisp[winner]
   print, 'Using a dispersion of ', disp ,' Angstroms'
   print, '  matched with grating: ', gratings[winner]

   good = where(linelist.q EQ 'GOOD',ngood)
   if good[0] EQ -1 then return, 0
   gl = linelist[good]

   guess = (maxwave - gl.wave)/disp - place[winner]
   match = lonarr(ngood) - 1
   mtrace = djs_median(xfit, 1)

   for i=0,ngood - 1 do begin
     diff = abs(guess[i] - mtrace)
     s = sort(diff)
     close = diff[s[0]]
     if close LT disp*1.2 then match[i] = s[0]
   endfor

   mm = where(match GE 0, nfit)
   if mm[0] EQ -1 then return, 0

   mlist = match[mm]  
   finallist = gl[mm] 

   print, 'We have traces near ', nfit, ' lines'

   traceset2xy, xset, yfit, xfit
   xx = xfit[*,mlist]
   ywave = finallist.wave # replicate(1,nrow)

   xy2traceset, transpose(xx), ywave, tset, ncoeff=3, $
      maxdev=disp/5.0, xmin=0.0, xmax=ncol-1., yfit=yyfit, /silent

   pix1 = traceset2pix(tset, gl.wave)
   further = where(total(pix1 GT 0,2) EQ nrow AND $
                   total(pix1 LT ncol,2) EQ nrow, nfurther)

   print, nfurther, ' lines may fall in the CCD region'

   if nfurther LT 3 then return,0 

   yy = replicate(1,nfurther) # findgen(nrow)
   mask = yy * 0 
   mask[*,goodrows] = 1

   pix2 = trace_fweight(image, pix1[further,*], yy, radius=1.5)
   pix2 = trace_fweight(image, pix2, yy, radius=1.5)
   pix2 = trace_fweight(image, pix2, yy, radius=1.5)
   pix3 = trace_fweight(image, pix2, yy, radius=1.5)

   xy2traceset, transpose(yy), transpose(pix3), finalset, ncoeff=3, $
         invvar=transpose(mask), yfit=xfinal, maxdev=0.1, $
         xmin=0, xmax=nrow-1, /silent

;
;  Now do final 2d wavelength fit
;

   wavefinal = gl[further].wave # replicate(1,nrow)
   xy2traceset, transpose(xfinal), wavefinal, t2set, ncoeff=4, $
          func=func, maxdev=disp/5.0, /silent, $
          xmin=0.0, xmax=ncol-1., yfit=yfinal2, outmask=outmask

   finally = where(total(outmask,2) GT nrow*0.8, nfinal3)
   print, nfinal3, ' lines remain after rejection'

   mask3 = (total(outmask,2) GT nrow*0.8) # replicate(1,nrow)
   xy2traceset, transpose(xfinal), wavefinal, t3set, invvar=mask3, $
       ncoeff=4, func=func, xmin=0.0, xmax=ncol-1., $
       yfit=yfinal3, outmask=outmask, /silent

   traceset2xy, t3set, pixnorm, waveimage

   ; Now let's output some diagnostics

   xdif = wavefinal[finally,*]-yfinal3[finally,*]
   print, 'Arcline  Wave(Ang)  Mean pixel    median dev     sigma dev '
   for k=0, nfinal3-1 do begin
     djs_iterstat, xdif[*,k], median=md, sigma=sg
     print, k, wavefinal[finally[k],0], mean(xfinal[*,finally[k]]), md, sg, $
     format='(I3,5x,F8.2,5x,F7.2,7x,F6.2,8x,F6.2)'
   endfor
   
   plot, wavefinal[finally,*], xdif,ps=3
   return, waveimage
end
   


