;+
; NAME:
;   deimos_skytweak
;
; PURPOSE:
;   adjust wavelength solution row by row to observed skylines
;
; CALLING SEQUENCE:
; 
; INPUTS:
;   slit structure (from spSlit file)
;
; OPTIONAL INPUTS:
;   slitwid  - slit width [pix]
;   skyline  - skyline list to use
;	
; KEYWORDS:
;
; OUTPUTS:
;   slit structure, with .dlam field modified
;
; OPTIONAL KEYWORD OUTPUTS:
;   dlamsky=dlamsky  - array[nrow] amount of shift [Ang]
;   skymodel=skymodel - evaluated sky model (2D array)
;   lineparams=lineparams - gaussian fits to sky lines, array of
;                           structures
;   skytweak=skytweak - structure containing overall and line-by-line
;       tweaks, evaluated only on good pixels
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   2002-Jul-14 - Written by D. Finkbeiner, Princeton
;   2002-Jul-19 - Extensive revisions
;   2002-Jul-25 - call gfit_line - DPF
;
;----------------------------------------------------------------------
pro deimos_skytweak, slit, slitwid=slitwid, dlamsky=dlamsky, skyline=skyline, $
          skymodel=skymodel, lineparams=g, plot=plot, skytweak=skytweak, $
          crmask=crmask


    zerotweak = {medshift: 0., $ ; median shift overall
                 mnshift: 0., $   ; mean shift overall
                 sigshift: 0., $ ; RMS about that shift
                 stdshift: 0., $ ; std. error in the shift
                 medsigma: 0.} ; median Gaussian sigma of line

  if n_elements(crmask) eq 0  then $
    flag_cr, slit, tempivar, crmask $
    else if total(crmask ne 0) eq 0. then $
       flag_cr, slit, tempivar, crmask

  spec     = djs_maskinterp(slit.flux, crmask, iaxis=1)
  specivar = slit.ivar*(slit.mask EQ 0)*(crmask EQ 0)

  if NOT keyword_set(slitwid) then slitwid = 3 ; pixels
  npix   = (size(spec, /dimens))[0]
  nrow   = (size(spec, /dimens))[1]
  ymid   = nrow/2
  nbuf = 6
  lambda = lambda_eval(slit,/double) ; first evaluation

	specivar[*,0]=0.
	specivar[*,nrow-1]=0.


; -------- read sky line list
  if NOT keyword_set(skyline) then begin 
     dir = getenv('DEEP_DIR')+'/'
     if dir eq '/' then message, 'you must set $DEEP_DIR'
     fname = dir+'spec2d/etc/deep_skylines.dat'
;     fname = dir+'spec2d/etc/yuko_skylines.dat'
; should probably make this fail if the file is not there!
     flist = findfile(fname)
     if flist[0] ne fname then fname = dir+'spec2d/etc/sdss_skylines.dat'
     skyline = read_lampfile(fname)
     skyline = skyline[where(skyline.good)]
  endif 

; -------- trim line list to lines fully on chip
  lamb = lambda[*, 0]
  lamt = lambda[*, nrow-1]
  xstartb = interpol(findgen(npix), lamb, skyline.lambda)
  xstartt = interpol(findgen(npix), lamt, skyline.lambda)

	; two checks for vignetting: at central row and using the whole slit
        ; take the most generous of the two

   vigpix = (slit.mask[*, ymid] AND 8B) EQ 8B OR slit.ivar[*, ymid] EQ 0
   wheregood = where(vigpix eq 0)
   vigpix2 = total((slit.mask AND 8B) EQ 8B, 2) gt (0.2*nrow > 10)
   wheregood2=where(vigpix2 eq 0) 

   wheregood=[wheregood,wheregood2]	
   if max(wheregood) ge 1 then begin
	   goodrange = minmax(wheregood[where(wheregood ge 0)])
   endif else begin
     print, 'Sky tweak found no good columns'
     dlamsky = 0
     skytweak=zerotweak
     return
   endelse

  vigrow = total(slit.mask NE 0 OR slit.ivar EQ 0, 1) gt 0.7*npix

   goodrow = indgen(nrow) ge nbuf AND indgen(nrow) lt nrow-nbuf-1 AND NOT vigrow
   if total(goodrow) ge 1 then goodind = where(goodrow)


   w = where( ((xstartb < xstartt) GT (goodrange[0] > 11)) and ((xstartb > xstartt) LT (goodrange[1] < (npix-12))), nw) 


  if nw lt 1 or total(goodrow) lt 1 then begin 
     print, 'Sky tweak found no good lines! (or only 1)'
     dlamsky = 0
     skytweak=zerotweak
     return
  endif 

  skylinew = skyline[w]
  xstart = (xstartb[w]+xstartt[w])/2.


; -------- trace sky lines
  rad =  fix(4+slitwid/2)  ; make sure this is an int
  xcen1 = trace_crude(spec, specivar, xstart=xstart, ystart=ymid, $
            radius=rad, yset=ycen, nave=1, nmed=1, $
            maxerr=0.5, maxshifte=.9, maxshift0=2., xerr=xerr)

; ---- THIS STILL ISN'T THE RIGHT WAY TO DO THIS...
  if slitwid LT 15 then begin 
;     specmed = (median(spec[*, ymid], 99)#(fltarr(nrow)+1))
     ;can't do this with sky spectrum; features too dense
     specmed = median(spec)
  endif else begin 
     specmed = median(spec)
  endelse
  specsub = spec-specmed

  xcen2 = trace_fweight(specsub, xcen1, ycen, $
                        radius=rad, xerr=xerrx,  invvar=specivar) 
  
; smaller radius on second pass
  xcen = trace_fweight(specsub, xcen2, ycen, $
                       radius=rad-2, xerr=xerr,  invvar=specivar) 

  neterr=1/sqrt(total(1./xerr^2,1))

  fracbad = total(xerr GT 5, 1)/nrow

  interior = total((xcen LT rad) AND (xcen GT npix-rad-1), 1) EQ 0
  keep = where(interior AND (neterr LT 0.04), nkeep)

  
  if nkeep le 1 then begin 
     print, 'Sky tweak found no good lines! (or only 1)'
     dlamsky = 0
     skytweak=zerotweak
     return
  endif 
; toss bad lines
  xcen = xcen[*, keep]
  xerr = xerr[*, keep]
  ycen = ycen[*, keep]
  skylinew = skylinew[keep]  ; ugly!
  nw = n_elements(keep)

; --------  gauss fit
; Fit parameters: A0 = peak height, A1=cen, A2=sigma, A3=baseline

  if slitwid lt 10 then begin 
     g = gfit_line(xcen, rad, specsub, specivar)
     
     goodmask = g.mask*0b
     goodmask[goodind, *] = 1b

; correct slit structure for new lambda offsets
     lamobs  = interpolate(lambda, xcen, ycen)
     dlambdadx=interpolate(lambda, xcen+1., ycen)-lamobs
     dlamsky = lamobs-transpose(skylinew.lambda#(fltarr(nrow)+1))
     wgood = where(g.mask ne 0 AND goodmask ne 0 AND xerr lt 5, ngood)

; if slit is no good, do nothing to it
     if ngood lt 5 then begin
        print, 'no information in skytweak!'
        message, 'no information in skytweak!', /INFO
        dlamsky =  0
	skytweak=zerotweak
        return
     endif
;

     medshift = median(dlamsky[wgood])
     mnshift = djs_mean(dlamsky[wgood])
     sigshift = stdev(dlamsky[wgood])
     stdshift = stdev(dlamsky[wgood])/sqrt(ngood)
     medsigma=median(g[wgood].sig*(dlambdadx[wgood]))

     medbyline = fltarr(nw)
     mnbyline = medbyline
     sigbyline = medbyline
     stdbyline = medbyline
     lambyline = skylinew.lambda

     for i=0, nw -1 do begin
        rightline = g.mask*0
        rightline[*, i] = 1
        wgoodl = where(g.mask ne 0 and goodmask and rightline $
                       and xerr lt 5, ngoodl)
        if ngoodl gt 1 then medbyline[i] = median(dlamsky[wgoodl])
        if ngoodl gt 1 then mnbyline[i]  = djs_mean(dlamsky[wgoodl])
        if ngoodl gt 1 then sigbyline[i] = stdev(dlamsky[wgoodl])
        if ngoodl gt 1 then stdbyline[i] = stdev(dlamsky[wgoodl])/sqrt(ngoodl)
     endfor


     skytweak = {medshift: medshift, $ ; median shift overall
                 mnshift: mnshift, $   ; mean shift overall
                 sigshift: sigshift, $ ; RMS about that shift
                 stdshift: stdshift, $ ; std. error in the shift
                 linelambda: lambyline, $; wavelengths of lines
                 linemedshift: medbyline, $; median shift for each line
                 linemnshift: mnbyline, $ ; mean shift for each line
                 linesigshift: sigbyline, $; RMS about the shift for each line
                 linestdshift: stdbyline, $ ; std error for each line
                 medsigma: medsigma} ; median Gaussian sigma of line


;  slit.dlam = slit.dlam-total(dlamsky*g.mask,2)/total(g.mask,2)
     
; -------- force dlam change to be a constant
     slit.dlam = slit.dlam-medshift
     lambda = lambda_eval(slit,/double) ; evaluate again with new dlam
  endif else begin 
     medshift = 0
     dlamsky = fltarr(nrow)
  endelse 

  if keyword_set(plot) and medshift ne 0 then begin 
     
     plot, lamobs, ycen, ps=6, syms=.25, xrange=minmax(lamb)+[-200, 200], $
       /xst, /nodata, title='Sky sanity check', $
       xtitle='wavelength [Ang]', ytit='row', chars=1.5
     
     for i=0, nw-1 do begin 
         if sigbyline[i] LT 0.06 then begin
             oplot, lamobs[*, i], ycen[*, i], ps=6, syms=.25
             oplot, lamobs[*, i]+dlamsky[*, i]*250, ycen[*, i], ps=-1
        endif
        ;mn = djs_mean(dlamsky[*, i])
        ;sig = djsig(dlamsky[*, i])
        print, 'Line:', i, '  lambda:', skylinew[i].lambda, '  offset:', mnbyline[i], $
          ' +/-', sigbyline[i], ' Ang', format='(A,I2,A,F8.2,A,F7.2,A,F6.2,A)'
        
     endfor 

     print, 'median shift: ', medshift
     print,'mean shift: ',mnshift
;     print, 'median inc. edges: ', medshifto
     print, 'RMS scatter: ', stdev(dlamsky[wgood])
     print, 'std. error: ', stdev(dlamsky[wgood])/sqrt(ngood)

  endif 
  












;===================== presently unused
  if NOT arg_present(skymodel) then return

;  skyind = slit.skyind
   skyind=where(slit.skyrow)
; sky fit (again!)
  skywave = lambda[*, skyind]
  skyflux = slit.flux[*, skyind]
  skyivar = slit.ivar[*, skyind]


  everyn = 2.*n_elements(skyind) /3.
  outmask=1
  sind = sort(skywave)
  ind = sind[where(skyivar[sind] NE 0)]

  sset = bspline_iterfit(skywave[ind], skyflux[ind], $
               invvar=skyivar[ind], upper=20, lower=20, $ 
               maxiter=3, everyn=everyn, $
               outmask=outmask1, /silent)

  outmask = byte(skywave*0)
  outmask[ind] = outmask1
  wavmask = byte(lambda*0)
  wavmask[*, skyind] = 1B - outmask

  lamrange = minmax(sset.fullbkpt)
  lambad = (lambda LT lamrange[0]) OR (lambda GT lamrange[1])
  slit.mask = slit.mask OR lambad*32B  ; mask pixels with no lam soln
  lambda = lambda > lamrange[0]
  lambda = lambda < lamrange[1]
  skymodel = bspline_valu(lambda , sset)     

  s = gfit_line(xcen[5, *], rad, skymodel[*, 5], skymodel[*, 5]*0+1)

; correct for observed slit function
  
  slitamp = total(g.amp*g.mask/((1+fltarr(nrow))#s.amp),2)/total(g.mask, 2)
  if total(finite(slitamp) eq 0) eq 0 then begin 
     wgood = where(total(g.mask, 2) NE 0)
     poly_iter, wgood, slitamp[wgood], 1, 3, coeff=coeff
     slitcorr = poly(findgen(nrow), coeff)

     if min(slitcorr) lt .8 OR max(slitcorr) gt 1.2 then print, 'Slitcorr too big!'
     slitcorr = (slitcorr < 1.2) > 0.8
     
     skysub = skymodel-specmed
     for i=0, nrow-1 do skymodel[*, i] = skysub[*, i]*slitcorr[i]
     skymodel = skymodel+specmed
  endif 

  if keyword_set(plot) then begin 
     atv, spec-skymodel, max=10, min=-10
     
     for i=0, nkeep-1 do $
       atvplot, xcen1[*, i], ycen[*, i], color=1
     for i=0, nkeep-1 do $
       atvplot, xcen[*, i], ycen[*, i], color=2
         
  endif 
  
  return
end
