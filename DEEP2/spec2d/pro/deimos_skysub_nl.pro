;+
; NAME:
; deimos_skysub_nl
;
; PURPOSE:
; Does nonlocal sky subtraction for DEEP2 slitlets.  The method is as
; follows.  (1) Find the nearest sky-only slitlets to the given object
; slitlet that maximize wavelength coverage.  (2) Do variable sky
; tweaking.  (3) Run a new bspline,  using all rows, on the object
; slitlet.  (4) Find the convolution kernels necessary to make each of
; the sky bsplines look like the object bsplines.  The default here is
; to use a kernel that varies with peak location along the spectrum,
; but it is also possible to use a constant kernel. (5) Convolve each
; sky spectrum with the appropriate kernel and subtract the result
; from the object bspline.  (6) Fit the subtracted spectrum with a
; bspline and evaluate for a 2d subtracted spectrum.
;
; CALLING SEQUENCE:
; deimos_skysub_nl, file, nhalf, subflux
;
; INPUTS:
;
; file     - File name for the spSlit.fits file to be considered.
; nframe   - Number of frame to be considered within the spSlit file
;            (1,2, or 3).
; nhalf    - Half-width of kernel to use for subtraction (kernel is
;            2*nhalf+1 elements wide).
;
; OPTIONAL INPUTS:
;	
; deweighting - Factor by which to deweight the variance of continuum
;               regions before fitting the kernel (default: 1.). 
; npoly       - polynomial order plus one for fitting a legendre
;               polynomial to find the wavelength-dependent kernel
;               (default: 2).
;
; KEYWORDS:
;
; const    - Causes the routine to fit a constant kernel to the entire
;            spectrum, rather than fitting to many windows and fitting
;            a smooth function to give a varying kernel.
; presmooth - Smooth the 2d object (and sky) data with a 1-2-1 kernel before
;             anything else.
; fft       - Compute kernel function via FFT technique
; skybsplines - Compute bsplines for all sky-only slits on the mask
;               and save this information in .sav files in the current
;               directory.  (This should be removed for v0_10 reductions)
; plot     - Plot the bspline subtraction for one of the sky slitlets:
;            Red=object spectrum; White=sky spectrum; Green=convolved
;            sky spectrum (model); Blue=subtraction residual.
;
; OUTPUTS:
;
; subflux   - 2d sky-subtracted spectrum, from non-local subtraction.
;
; OPTIONAL OUTPUTS:
;
; selfsubflux - 2d sky-subtracted spectrum, from LOCAL subtraction.
; skymodel   - 2d sky model used for 2d subtraction.
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; Calls the following deep routines:
; find_kernel_reg.pro
; getsky.pro
; get_slitpos.pro
; var_kern_fit.pro
; convol_var.pro
; deimos_skytweak.pro
; deimos_skytweak_var.pro
; find_skyslits.pro
;
; REVISION HISTORY:
;
; Written by BFG, Summer 2002
; Made presentable 2002Oct18.
; Revised by JAN Nov 02
; Revised by BFG 2002Nov11
;
;----------------------------------------------------------------------


 
pro deimos_skysub_nl,  file,  nframe, nhalf, subflux, $
          npoly=npoly, selfsubflux=selfsubflux, fft=fft, $
          skymodel=skymod2d, const=const, presmooth=presmooth, $
          skybsplines=skybsplines, deweighting=deweighting, plot=plot, $
;---The rest are all to give back 1d spectra and kernels for testing 
;---and debugging:
  wave1d1= fullwave1,  wave1d2=fullwave2, $
    sky1d1= fullsky1,  sky1d2= fullsky2, $
    obj1d1=fullobj1, obj1d2=fullobj2, $
    mod1d1=skymodel1, mod1d2=skymodel2, $
    sub1d1=sub1,  sub1d2=sub2, $
    kern1=kernel1, kern2=kernel2,  fitk1=fitkern1,  fitk2=fitkern2, $
    fitl1=fitlam1, fitl2=fitlam2,  medobj1=medobj1,  medobj2=medobj2, $
    medsky1=medsky1,  medsky2=medsky2,crmask=crmask

;---get info about slit and mask from filename
  splitfile =  strsplit(file, '.', /extract)
  masknum = fix(splitfile[1])
  slitnum = fix(strmid(splitfile[2], 0, 3))
  color =  strmid(splitfile[2], 3, 1)

;---get name of data directory (if specified)
  splitdir =  strsplit(file, 'spSlit', /regex, /extract)
  datadir = splitdir[0]

  
  objbspline = 1
  skybsplines = 0

;---set amount of buffering on slit ends.
  nbuf = 5

;---find position of slit on mask
  get_slitpos,  datadir,  slitnum, slitx, slity

;---find sky slits and compute distances from slit


  getsky,  masknum,  datadir, color,  sky_slitn,$
    sky_slitx=sky_slitx,  sky_slity=sky_slity, $
    sky_minmaxwave=sky_minmaxwave,  skyspfile=skyspfile, $
    nexp=1, bspline=skybsplines, $ 
    nbuf=nbuf, /silent

;---read in object slit  
  objslit = mrdfits(file, 2*nframe-1, objhdr,/SILENT)
  objsset = mrdfits(file, 2*nframe,/SILENT)

  hasinfo=total(tag_names(objslit) eq 'INFOMASK') gt 0
  
  objx =  n_elements(objslit.ivar[*, 0])
  objy =  n_elements(objslit.ivar[0, *])
  ncol=objx

  oldobjwave =  lambda_eval(objslit,/double)
  objflux =  objslit.flux

  if hasinfo then whgood = where((objslit.infomask AND 1b) eq 0b $
                                 AND (objslit.mask AND 22b) eq 0b,goodct) $ 
  else whgood = where(objslit.mask eq 0b,goodct)

  obj_minmaxwave = minmax(oldobjwave[whgood])
                                ;^throws out vignetted regions. 

; check to see if we have existing sky slits
  goodframes=where(sky_slitx lt 1E4,goodct)
  if goodct ge 2 then begin
     sky_slitx=sky_slitx(goodframes)
     sky_slity=sky_slity(goodframes)
     skyspfile=skyspfile(goodframes)
     sky_minmaxwave=sky_minmaxwave[goodframes,*]
     dist = sqrt((slitx-sky_slitx)^2 + (slity-sky_slity)^2)
     
;---find nearest sky slitlets 
     whposx=where(sky_slitx gt 0,posct)
     whnegx=where(sky_slitx le 0,negct)

     if posct gt 0 and negct gt 0 then begin
         isdposx = bsort(dist[whposx])

         isdnegx = bsort(dist[whnegx])


         ilt = whnegx[isdnegx[0]]
         igt = whposx[isdposx[0]]


     endif else begin
         isd = bsort(dist)
         ilt=isd[0]
         igt=isd[1]
     endelse

 
  endif

;-----do sky tweaking

;  new_obj_fit = deimos_skytweak_var(objslit, objsset,  color, chisq=objchi)
;  objslit.lambdax = new_obj_fit

; evaluate object bspline
  objwave =  oldobjwave         ; lambda_eval(objslit)

;----flag cosmic rays and interpolate over bad pixels


  if n_elements(crmask) eq 0  then $
    flag_cr, objslit, newobjivar, crmask $
    else if total(crmask ne 0) eq 0. then $
       flag_cr, objslit, newobjivar, crmask

  if n_elements(newobjivar) eq 0 then newobjivar=objslit.ivar*(crmask eq 0)

;  flag_cr, objslit,  newobjivar,crmask
  objmask = newobjivar eq 0
  objflux = djs_maskinterp(objflux,  objmask, iaxis=1)


;---presmooth object data if called for.
  if keyword_set(presmooth) then begin
     skern = fltarr(3, 1)
;     skern[*, 0] =  [1/4., 1/2., 1/4.]
     skern[*, 0] = [1., 1., 1.]/3.
     objflux =  convol(objflux, skern, /center)
     ivarmask=newobjivar eq 0.
     ivarmask=dilate(ivarmask,intarr(3,1)+1) OR ivarmask

     invivar=1./newobjivar
     wh=where(finite(invivar) eq 0B, badct)
     if badct ge 1 then invivar[wh]=5.
     newobjivar=1./(1./16.*shift(invivar,1,0)+ $
                    1./16.*shift(invivar,-1,0)+ $
                    1./4.*invivar) 
     newobjivar=newobjivar*(ivarmask eq 0B)

  endif



  

;---pad ends (hopefully this will improve!)

  fullobjflux=objflux
  fullwave=objwave

  if nbuf gt 0 then begin
     objflux =  objflux[*, nbuf:objy-1-nbuf]
     oldobjwave = oldobjwave[*, nbuf:objy-1-nbuf]
     newobjivar =  newobjivar[*, nbuf:objy-1-nbuf]

     objwave = objwave[*, nbuf:objy-1-nbuf]

  endif

  skyind=findgen(n_elements(objwave[0,*]))
  nsky=n_elements(skyind)
  skyx=lindgen(ncol,nsky) MOD ncol

  badsky=(newobjivar eq 0.)
  whok=where(badsky eq 0b,okct)
  if okct gt 0 then begin
       meanivar=mean(newobjivar[whok])
       whbad=where((badsky AND (skyx lt 7 OR skyx gt (ncol-8))) ,badct)
       if badct gt 0 then newobjivar[whbad]=meanivar/250.
  endif

  whcr=where(crmask[*,nbuf:objy-1-nbuf] ne 0,crct)
  if crct gt 0 then newobjivar[whcr]=0.

;---do bsplines with flagged CR regions

    everyn = 2.*n_elements(skyind)/3.
    bkspace = (max(objwave)-min(objwave))/(1.5*ncol)
    testrange=objwave[ncol/2,*]
    testrange=testrange[sort(testrange)]
 ; determine range of wavelengths in a given column
    deltalambda=abs(testrange[fix(0.9*nsky)]-testrange[fix(0.1*nsky)])


  bkspace =  0.2
  objsset_old = objsset

  if objbspline then begin
     isort =  sort(objwave) 
     sortwave =  objwave(isort)
     sortflux =  objflux(isort)
     sortivar =  newobjivar(isort)

     if deltalambda lt 0.6*bkspace then begin

              row=floor(nsky/2.)
              print,'Using 1 breakpoint/pixel'
              breakpoints = (objwave[*, row:row+1])
              breakpoints=total(breakpoints,2)/2.
              objsset =  bspline_iterfit(sortwave,sortflux,invvar=sortivar, $
                                upper=20,  lower=20,  maxiter=3,$
                                fullbkpt = breakpoints, /silent)
      endif else begin       

              objsset =  bspline_iterfit(sortwave,  sortflux,  invvar=sortivar, $
                                upper=20,  lower=20,  maxiter=3,$
                                bkspace= bkspace, /silent) 
      endelse        

  endif   

  if goodct lt 2 then begin
     obj2d=bspline_valu(fullwave,objsset)
     subflux=fullobjflux-obj2d
     selfsubflux=subflux

     print,'No sky slits found, doing 1-2-1 local sky instead'
     skymod2d=obj2d
     wave1d1= 0 &  wave1d2=0 & $
       sky1d1= 0 &  sky1d2=0 & $
       obj1d1= 0 & obj1d2= 0 & $
       mod1d1=0 & mod1d2= 0& $
       sub1d1= 0 &  sub1d2=0 & $
       kern1= 0 & kern2=0 & fitk1=0 &  fitk2=0 & $
       fitl1=0 & fitl2=0 &  medobj1=0 & medobj2=0 & $
       medsky1=0 &  medsky2=0

     return
  endif

;---Get sky bsplines either from spslit files or from files output by getsky
  skybfile = 'sky_sset.'+string(ilt,format='(i1)')+'.0'+color+'.sav'
  skybfile = findfile(skybfile, count=count)
  if count eq 1 then begin
     restore,  file=skybfile
     skysset1 =  sky_sset
  endif else begin

     skyfile1 = skyspfile[ilt]
     skyfile1 = skyfile1[0]
     skysset1 = mrdfits(skyfile1, 2*nframe,/SILENT)
  endelse

;---same as above for skyslit 2

  skybfile = 'sky_sset.'+string(igt,format='(i1)')+'.0'+color+'.sav'
  skybfile = findfile(skybfile, count=count)
  if count eq 1 then begin
     restore,  file=skybfile
     skysset2 = sky_sset
  endif else begin

     skyfile2 = skyspfile[igt] 
     skyfile2 = skyfile2[0]  
     skysset2 = mrdfits(skyfile2, 2*nframe,/SILENT)

  endelse

  skysplit = strsplit(skyfile1, '.',  /extract)
  skyn1 = string(skysplit[2])
  skysplit = strsplit(skyfile2, '.',  /extract)
  skyn2 = string(skysplit[2])

  print,  'Using sky slits '+skyn1+' and '+skyn2

;---determine common wavelength range; pad slightly.  
  bad_angst = 2. ; number of angstroms on each end about which we have bad angst
;                  (i.e., those we think might have bspline end-effect problems)
  minlam1 = max([sky_minmaxwave[ilt, 0], obj_minmaxwave[0]])+bad_angst
  maxlam1 = min([sky_minmaxwave[ilt, 1], obj_minmaxwave[1]])-bad_angst

  minlam2 = max([sky_minmaxwave[igt, 0], obj_minmaxwave[0]])+bad_angst
  maxlam2 = min([sky_minmaxwave[igt, 1], obj_minmaxwave[1]])-bad_angst


;---create regular wavelength grid
  dlam = 0.2
  npts1 = floor((maxlam1-minlam1)/dlam) > 1
  npts2 = floor((maxlam2-minlam2)/dlam) > 1
  fullwave1 = findgen(npts1)*dlam + minlam1 
  fullwave2 = findgen(npts2)*dlam + minlam2 

	badsky=(npts1 lt 3 OR npts2 lt 3)

;---evaluate bsplines on this grid.  
  fullsky1 = bspline_valu(fullwave1,  skysset1)
  fullsky2 = bspline_valu(fullwave2,  skysset2)



  fullobj1 = bspline_valu(fullwave1,  objsset)
  fullobj2 = bspline_valu(fullwave2,  objsset)


; THE FOLLOWING IS CURRENTLY THE RATE-LIMITING STEP!!!!!

;---median smooth (20 A box) to determine continuum level
;  medobj1 =  djs_median(fullobj1, width=1000/dlam, boundary='reflect')
;  medobj2 =  djs_median(fullobj2, width=1000/dlam, boundary='reflect')
;  medsky1 =  djs_median(fullsky1, width=1000/dlam, boundary='reflect')
;  medsky2 =  djs_median(fullsky2, width=1000/dlam, boundary='reflect')

;---model continuum by taking median of lowest 1/frac of spectrum
  frac = 1.5
  sfo1 = bsort(fullobj1)
  sfullobj1 = fullobj1[sfo1]
  sfo2 = bsort(fullobj2)
  sfullobj2 = fullobj2[sfo2]
  sfs1 = bsort(fullsky1)
  sfullsky1 = fullsky1[sfs1]
  sfs2 = bsort(fullsky2)
  sfullsky2 = fullsky2[sfs2]
  medobj1=median(sfullobj1[0:npts1/frac])
  medobj2=median(sfullobj2[0:npts2/frac])
  medsky1=median(sfullsky1[0:npts1/frac])
  medsky2=median(sfullsky2[0:npts2/frac])

  oldobj1 =  bspline_valu(fullwave1, objsset_old)
  soo1 = bsort(oldobj1)
  soldobj1 = oldobj1[soo1]
;  medskymodel1 =  djs_median(oldobj1, width=1000/dlam, boundary='reflect')
  medskymodel1=median(soldobj1[0:npts1/frac])

  oldobj2 =  bspline_valu(fullwave2, objsset_old)
  soo2 = bsort(oldobj2)
  soldobj2 = oldobj2[soo2]
;  medskymodel2 =  djs_median(oldobj2, width=1000/dlam, boundary='reflect')
  medskymodel2=median(soldobj2[0:npts2/frac])

;---estimate continuum error
  rn =  2.32                    ; read noise
  exposure_time = float(sxpar(objhdr, 'EXPTIME')) ;get exp. time in seconds
  flux_corr=3600./exposure_time ;flux renormalization
  
  sigcont1 =  sqrt((flux_corr*medobj1 + (flux_corr*rn)^2)/(objy-2*nbuf))
  
  sigcont2 =  sqrt((flux_corr*medobj2 + (flux_corr*rn)^2)/(objy-2*nbuf))

;---estimate error at each point in object-slit bspline
  
  measerr1 = sqrt(((flux_corr*fullobj1>5.)+(flux_corr*rn)^2)/(objy-2*nbuf)) 
                                ;measurement error
  measerr2 = sqrt(((flux_corr*fullobj2>5.)+(flux_corr*rn)^2)/(objy-2*nbuf)) 


  if keyword_set(presmooth) then begin
     correction=sqrt(3./8.)
     sigcont1=sigcont1*correction
     sigcont2=sigcont2*correction
     measerr1=measerr1*correction
     measerr2=measerr2*correction
  endif

;TEST!!
  measerr1 = measerr1*0. +1.
  measerr2 = measerr2*0. +1.

;---preliminaries for kernel fit:
  

  if not keyword_set(deweighting) then $
    deweighting = 1.            ;continuum variance increased by this factor.

  window = 100.                 ;length of region to fit kernel over (Angstroms)
  winstep = 25.                 ;amount by which to step window in variable kernel fitting 
  if not keyword_set(npoly) then $
     npoly = 2                     ;order of polynomial fit +1 (spectral direction)

   threshold1 =  medobj1+6.*sigcont1
   threshold2 =  medobj2+6.*sigcont2
;  threshold1 = -1
;  threshold2 = -1 ;disable thresholds for test of FFT method

  if keyword_set(const) then begin
;---fit a constant kernel and convolve

;---find regions with signal; broaden slightly to get wings
  
     fullsignal1 =  (fullobj1 gt threshold1)
     fullsignal2 =  (fullobj2 gt threshold2)

     fullsignal1=fullsignal1 OR dilate(fullsignal1,fltarr(3*nhalf)+1)
     fullsignal2=fullsignal1 OR dilate(fullsignal2,fltarr(3*nhalf)+1)

;---deweight continuum regions by increasing errors
     icont1 =  where(fullsignal1 eq 0)
     measerr1[icont1] = sqrt(deweighting)*measerr1[icont1]
     icont2 =  where(fullsignal2 eq 0)
     measerr2[icont2] = sqrt(deweighting)*measerr2[icont2]
;---subtract off continuum from both spectra     
     fullsky1 =  fullsky1-medsky1
     fullobj1 = fullobj1-medobj1
     fullsky2 =  fullsky2-medsky2
     fullobj2 = fullobj2-medobj2

;---find kernel to go from sky to object spectrum, leaving out ends
     ifit1 =  where((fullwave1 ge minlam1+20) and (fullwave1 le maxlam1-20),$
                    l1)
;---first set up weighting array to apodize spectra 
     h = hanning(4*nhalf+1)
     weight = fltarr(l1)+1.
     weight[0:2*nhalf] = h[0:2*nhalf]
     weight[l1-2*nhalf-2:l1-2] = h[2*nhalf:4*nhalf]
     weight[l1-1]=0.
;---apodize spectra
     fitsky1 = fullsky1[ifit1]*weight
     fitobj1 = fullobj1[ifit1]*weight

;---set up array of shifts
     shiftsky1 = fltarr(2*nhalf+1, l1)
     for i=-nhalf, nhalf do begin
        shiftsky1(i+nhalf, *) = shift(fitsky1, -i)  
               ;to keep the kernel defined as for the fft, we need this choice of sign
     endfor   

     kernel1=  find_kernel_reg(shiftsky1,  fitobj1,  nhalf, $
                               measerr1, fit1,  sigout1,  const1, fft=fft) 

;---same thing for the second spectrum
     ifit2 =  where((fullwave2 ge minlam2+50) and (fullwave2 le maxlam2-50),$
                    l2)

     h = hanning(4*nhalf+1)
     weight = fltarr(l2)+1
     weight[0:2*nhalf] = h[0:2*nhalf]
     weight[l2-2*nhalf-2:l2-2] = h[2*nhalf:4*nhalf]
     weight[l2-1]=0.

     fitsky2 = fullsky2[ifit2]*weight
     fitobj2 = fullobj2[ifit2]*weight
     shiftsky2 = fltarr(2*nhalf+1, l2)
     for i=-nhalf, nhalf do begin
        shiftsky2(i+nhalf, *) = shift(fitsky2, -i)  
              ;to keep the kernel defined as for the fft, we need this choice of sign
     endfor   
     kernel2=  find_kernel_reg(shiftsky2,  fitobj2,  nhalf, $
                               measerr2, fit2,  sigout2,  const2, fft=fft) 
     
;---convolve sky with kernel
   if NOT keyword_set(fft) then begin ;kernel is different in FFT case
       skymodel1 = convol(fullsky1,  kernel1, /center,  /edge_wrap)
       skymodel2 = convol(fullsky2,  kernel2, /center, /edge_wrap)
   endif else begin
       skymodel1 = convol(fullsky1/measerr1,  kernel1, /center, /edge_wrap)* $
          measerr1 ;kernel is for a prewhited spectrum, so undo this
       skymodel2 = convol(fullsky2/measerr2,  kernel2, /center, /edge_wrap)* $
          measerr2
   endelse

;---add back medians
     fullsky1 =  fullsky1+medsky1
     fullobj1 = fullobj1+medobj1
     fullsky2 =  fullsky2+medsky2
     fullobj2 = fullobj2+medobj2
     skymodel1 = skymodel1+medskymodel1
     skymodel2 = skymodel2+medskymodel2
;---end of constant kernel convolution.

  endif else begin
;---fit a variable kernel and convolve
     if badsky eq 0 then begin     

             minwav1 = min(fullwave1, max=maxwav1)
             noregions = (maxwav1-minwav1)/winstep lt 2.
         
	     if noregions eq 0 then $
                skymodel1= var_kern_fit(fullwave1, fullsky1, $
				fullobj1, measerr1, window, $
                             nhalf, npoly, winstep=winstep, contin1=medsky1, $
                             contin2=medobj1,  cont2prime=medskymodel1, $
                             conterr2=sigcont1, fft=fft,  $
                             deweight=deweighting,  fullsignal=fullsignal1, $
                             /contsub, $
                             kernels=kernel1,  xfit=fitlam1, $
                             fitkern=fitkern1,threshold=threshold1) $
                else skymodel1=-1.


             minwav2=min(fullwave2,max=maxwav2)
             noregions = (maxwav2-minwav2)/winstep lt 2.

	     if noregions eq 0 then $
               skymodel2= var_kern_fit(fullwave2, fullsky2, $
				fullobj2, measerr2, window, $
                             nhalf, npoly,  winstep=winstep, contin1=medsky2, $
                             contin2=medobj2, cont2prime=medskymodel2, $ 
                             conterr2=sigcont2, fft=fft,  $
                             deweight=deweighting,  fullsignal=fullsignal2, $
                             /contsub, $ 
                             kernels=kernel2,  xfit=fitlam2, $
                             fitkern=fitkern2,threshold=threshold2) $
                else skymodel2=-1.

	endif

;---end of variable kernel convolution.
  endelse

  bad1=(n_elements(skymodel1) eq 1) OR (npts1 lt 3)
  bad2=n_elements(skymodel2) eq 1 OR (npts2 lt 3)

  if bad1 OR bad2 then begin
      if bad1 then message,'Skyslit '+skyn1+' failed!',/INFO
      if bad2 then message,'Skyslit '+skyn2+' failed!',/INFO
      message,'for slit '+file+' frame '+string(nframe),/INFO

      if bad1 then print,'Skyslit '+skyn1+' failed!'
      if bad2 then print,'Skyslit '+skyn2+' failed!'
      print,'for slit '+file+' frame '+string(nframe)

      subflux=-1
      return
  endif

;---subtract off model from bsplines
  sub1 = fullobj1-skymodel1
  
  sub2 = fullobj2-skymodel2
  sub2[0:5]=-1.
  sub2[n_elements(fullobj2)-6:n_elements(fullobj2)-1]=-1.

;---plot modeling and residuals if desired.
  if keyword_set(plot) then begin
     splot,  fullwave2,  fullsky2,  yrange=[min(sub2)-5, max(fullobj2)+10]
     soplot,  fullwave2,  fullobj2,  color=1
     soplot,  fullwave2, skymodel2,  color=2 
     soplot,  fullwave2,  sub2,  color=3
     soplot,fullwave1,sub1,color=4
  endif


;---to do the subtraction, fit the skymodel with a bspline,
;---so it can be evaluated at the data pixels.
  
  skymodsset1 = bspline_iterfit(fullwave1,  skymodel1,  $   
                                maxiter=0, bkspace=0.21, /silent)

  skymodsset2 = bspline_iterfit(fullwave2,  skymodel2,  $   
                                maxiter=0, bkspace=0.21, /silent)

  skymod2d1 =  bspline_valu(fullwave,  skymodsset1)
  skymod2d2 =  bspline_valu(fullwave,  skymodsset2)
  obj2d =  bspline_valu(fullwave,  objsset)


  w1 = (fullwave gt minlam1+5) and (fullwave lt maxlam1-5)
  w2 = (fullwave gt minlam2+5) and (fullwave lt maxlam2-5)
  w3 = (w1 EQ 0 AND w2 EQ 0)



  if hasinfo then begin
     vigpix=where((objslit.infomask AND 1B) eq 1B,vigct)
     if total(vigct) gt 0 then begin
        w1[vigpix]=0
        w2[vigpix]=0
        w3[vigpix]=1
     endif
  endif

; average models where possible; else use self-subtraction.
  skymod2d =  (w1*skymod2d1 + w2*skymod2d2+w3*obj2d)/(w1+w2+w3) 

  subflux =  fullobjflux-skymod2d
  selfsubflux =  fullobjflux-obj2d  

; make them zero where we are vignetted
  vigpix=where((objslit.mask AND 8B) eq 8B,vigct)
  if total(vigct) gt 0 then begin
     subflux[vigpix]=0.
     selfsubflux[vigpix]=0.
  endif
end




