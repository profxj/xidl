
pro deimos_makeflat, flat_image, flat2d, flat1d, vigcorr = vigcorr, varslitfn = varslitfn, ivar=ivar, mask = mask, bitmask = bitmask, quick=quick
;+
; NAME:
;   deimos_makeflat
;
; PURPOSE:
;   Generate a 1-d slit function and 2-d flat field
;
;     Creates a 1-d slit function and 2-d flat field (fringe frame) 
;     for purposes of removing effects of variations in slit thickness
;     and fringing.
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   deimos_makeflat, flat_image, flat2d, flat1d,[vigcorr=vigcorr, $
;         varslitfn=varslitfn,ivar= ivar, bitmask=bitmask,mask=mask]
; 
; INPUTS:
;   flat_image -- rectified flat field image for a slitlet
;
; OPTIONAL INPUTS:
;   ivar     -- inverse variance of input flat image
;
; KEYWORD PARAMETERS:
;   mask --  (output) 1 for good data, 0 for bad
;   bitmask -- rectified bad pixel mask for this slit        
;
; BUGS:
;
; OUTPUTS:
;   flat2d  -- 2d-flat (fringe frame), MINUS ONE (to facilitate 
;              floatcompressing).  Defined in the sense that one
;              normally wants to divide by (flat2d+1) to remove fringes
;   flat1d  -- 1d-flat (slit function). Defined in the sense that one
;              normally wants to divide by this to correct for the
;              slit function.  Note that slit function information has
;              been removed from flat2d!
;   VIGCORR=vigcorr -- vignetting correction array.  Divide by this to roughly
;              remove vignetting and illumination from an arc spectrum
;   VARSLITFN=varslitfn -- 2d slit function array, varying in the
;                          spectral direction.  REBIN this to the full
;                          dimensions of the slit to get the slit
;                          function at a given pixel.
;
; MODIFICATION HISTORY:
;   C. Marinoni, Spring 2001
;   MD Jan-Apr-2002 
;   JAN 3 July 02
;   14-Jul-2002   Major renovation - Doug Finkbeiner
;   22-Oct-2002   Vignetting corrections and spectrally-varying slit
;   function implemented
;-

  if (NOT keyword_set(ivar)) then ivar=float(flat_image)*0.+1.
  if (n_elements(mask) eq 0) then mask=byte((flat_image)*0+1)
  if (n_elements(bitmask) eq 0) then bitmask=byte((flat_image)*0)
  if n_elements(quick) eq 0 then quick = 0
  
  fullvig = (bitmask AND 2b) eq 2b

  partvig = (bitmask AND 4b) eq 4b OR fullvig
  partvig = partvig OR dilate(partvig, intarr(101, 1)+1)
  fullvig = fullvig OR dilate(fullvig, intarr(101, 1)+1)


  nn = (size(flat_image))[2]    ;spatial range
  mm = (size(flat_image))[1]    ;spectral range

  flatinput = flat_image
  minval=min(flatinput)
  if minval LT 0 THEN flatinput = flatinput - minval + 1.

                                ; patch up holes in flat
  edge = (round(0.1*nn) > 5) < 10

  badpixels = (ivar eq 0) OR ((bitmask AND 1b) EQ 1b)

;  maybevig = (bitmask AND 2b) EQ 2
;  nvig = max(total(maybevig, 1))

;  maybevig = maybevig OR dilate(maybevig, intarr(long(nvig) > 1)+1)
;  maybevig = maybevig OR dilate(maybevig, intarr(long(nvig > 1)*3/4)+1)

; throw out blocked columns for normalizing maximum
  medcols=djs_median((flat_image*(ivar ne 0))[mm/2-mm/8:mm/2+mm/8,*],1)
  maxgood=max(medcols)
  lowcols=medcols/maxgood lt 0.05
  
                                ; These are DEIMOS columns: the spectral direction
  badcols=where((total(badpixels, 1) gt 0.5*mm) OR lowcols,badcolct)


; these are DEIMOS rows: the spatial direction
  badrows = where(total(badpixels[*, 0:edge], 2) gt edge OR $
                  total(badpixels[*, nn-edge-1:nn-1], 2) gt edge, badct)
  if badct gt 0 then badpixels[badrows, *] = 0
  
;  flat_processed = djs_maskinterp(float(flatinput), badpixels, iaxis=1, /const)
  flat_processed=flatinput
; need to define ncol in a way that will avoid vignetted or otherwise
; bad parts of the spectrum.


  sumwt = fltarr(nn)+1.
  sumwt[0:edge-1] = 1.e-20
  sumwt[nn-edge:nn-1] = 1.e-20
  sumwt = (fltarr(mm)+1) # sumwt
  if total(badpixels) gt 0 then sumwt[where(badpixels)] = 1.e-20


; first, define a totalled spectrum to allow us to correct for
; illumination gradients and test for low-level regions 
  summedspec = total(sumwt*flat_processed, 2)/total(sumwt, 2)

; we want to define the slit function where the data isn't too
; vignetted and is relatively commensurate with the peak.

  medlevel = median(summedspec[1200:mm-1201])

  

; do this on both sides - typically, one side will be problematic but
;                         not the other.
  goodspec1 = where(djs_median(summedspec[200:mm/2], width=5, boundary='reflect') gt 0.9*medlevel AND NOT partvig, goodct1)
  goodspec2 = where(djs_median(summedspec[mm/2+1:mm-201], width=5, boundary='reflect') gt 0.9*medlevel AND NOT partvig, goodct2)
  
  ncol = (goodct1 < goodct2)*0.8
  ncol = (ncol > 300) < 800



; now, define a correction for lambda-dependence of illumination.
; NEW: try to _ONLY_ take out the spectral shape from the central
; region

; for now, try a 4th order polynomial
  polydeg = 3
; LEGENDRE version:
  xnorm = (findgen(mm)-mm/2.)*2/mm

  medflat=median((flat_processed)[mm/2-ncol:mm/2+ncol-1, edge:nn-edge-1])
  flat_processed=flat_processed/medflat
  ivar=ivar*medflat^2

; -------- first determine slitfn & spatial illumination gradient - ignore mask


; do not take out any gradients in the _spectral_ direction before summing
  slitfntot = total((flat_processed)[mm/2-ncol:mm/2+ncol-1, *], 1)/(2.*ncol) ; UNWEIGHTED mean
  slitfnwt = total((ivar)[mm/2-ncol:mm/2+ncol-1, *], 1) ;but pretend it's weighted in this line, no biggie
  
  slitillum = slitfntot
  slitillumwt = slitfnwt
  
  slitillumwt[0:edge-1] = 1e-20 <  slitillumwt[0:edge-1]
  slitillumwt[nn-edge:nn-1] = 1e-20 < slitillumwt[nn-edge:nn-1]

; mask really bad pixels in the slitfn 
  slitfninterior = slitfntot[edge:nn-edge-1]
  medmiddle = median(slitfninterior)
  badslitfngrad = slitfntot*0
;  badslitfn[edge:nn-edge-1] = slitfninterior/medmiddle lt 0.6 $
;        OR slitfninterior/medmiddle gt 1.5
  
  badslitfngrad[0:edge]=1
  badslitfngrad[nn-edge-1:nn-1]=1

  if badcolct gt 0 then badslitfngrad[badcols] = 1
  badslitfnsave=badslitfngrad

  slitillumwt = slitillumwt > 1.E-20

  sample = findgen(1E4)*n_elements(slitilluwt)/(1E4+23)

  maxerr = max(1./sqrt(slitillumwt[sample]))

;  rollmed=djs_median(slitillum,width=9,boundary='reflect')
;  bigmed=djs_median(slitillum,width=(31 < n_elements(slitillum)-1),$
;                    boundary='reflect')
  


  niter=2
  i=0

  slitfnderiv=deriv(slitfntot)
  smderiv=smooth(deriv(slitfntot),3)
  slitfnder2=deriv(slitfnderiv)

  whgood=where(badslitfngrad eq 0 AND slitfnder2 lt 3E-3,goodct)
  whfull=whgood
  if goodct eq 0 then nogood = 1 else nogood = 0

  if nogood eq 0 then begin
   while goodct gt 0 and i lt niter do begin	
     djs_iterstat,((slitfntot[whgood])),med=medval,max=20
     djs_iterstat,(slitfnderiv[whgood]),med=medder, $
       sigma=sigder,max=20

     goodderiv=abs(slitfnderiv-medder) lt sigder $
       AND abs(smderiv) lt 0.01 AND slitfnder2 lt 2.5E-3

     srt=sort(slitfntot[whgood])

     maxslitfn=slitfntot[whgood[srt[(goodct-2) > 0]]]
     
     highslitfn=slitfntot gt maxslitfn-0.02 $
       AND slitfnder2 lt 2.5E-3
     
     whgood=where( (goodderiv OR highslitfn)$
                   AND badslitfngrad eq 0,goodct)

     if goodct eq 0 then whgood=whfull


     badslitfngradtmp=badslitfngrad $
       OR (goodderiv OR highslitfn) eq 0

     badslitfngrad=badslitfngradtmp
     whgood=where(badslitfngrad eq 0,goodct)

     if goodct eq 0 then whgood=whfull

                                ; if okct lt 5 then whgood=where(badslitfngradtmp EQ 0,okct)

     djs_iterstat,(slitfnderiv[whgood]),med=medder,sigm=sigder,max=20
     corrslitfn=slitfntot/(1+findgen(nn)*medder)
     djs_iterstat,corrslitfn[whgood],med=medval,max=20

     temp=slitillum/(1+findgen(nn)*medder)/medval

     limit=0.03*(1+0.5*(goodct lt 8))
     siglimit=1*(1+0.25*(goodct lt 8))
     tighten=1.5


     goodderiv=abs(slitfnderiv-medder) lt sigder*siglimit $
       AND abs(smderiv) lt 0.01 AND slitfnder2 lt 2.5E-3
     
     whgood2=where( ((abs(temp-1) lt limit and goodderiv) $
                     OR highslitfn)  $
                    AND badslitfnsave eq 0,goodct)

     if goodct eq 0 then whgood2=whfull
     
     for i=0,niter-1 do begin

        djs_iterstat,(slitfnderiv[whgood2]),med=medder,max=10
        corrslitfn=slitfntot/(1+findgen(nn)*medder)
        djs_iterstat,corrslitfn[whgood2],med=medval,max=20

        temp=slitillum/(1+findgen(nn)*medder)/medval
        limit=limit/tighten
        siglimit=siglimit*sqrt(tighten)

        goodderiv=abs(slitfnderiv-medder) lt sigder*siglimit $
          AND abs(smderiv) lt 0.01 AND slitfnder2 lt 2.5E-3

        whgood2=where( highslitfn $
                       AND badslitfnsave eq 0,goodct)
        if goodct eq 0 then whgood2=whfull

        if goodct lt 6 then $
          whgood2=where( ((abs(temp-1) lt limit and goodderiv) $
                          OR highslitfn)  $
                         AND badslitfnsave eq 0,goodct)
        
        if goodct lt 3 then whgood2=whfull
      endfor

      if n_elements(whgood2) le 1 then mostly1 = 1 $
          else $
     mostly1 = total(slitillum[whgood2] eq 1.)/goodct gt 0.6 $
       OR stdev(slitillum[whgood2]) gt 1E-4

     if goodct gt 4 AND (mostly1 eq 0) $
        then $
       poly_iter,whgood2,slitillum[whgood2],1,3.,yfit,coeff=coeff $
     else $
       if goodct gt 4 then $
       coeff=svdfit(whgood2,slitillum[whgood2],2,yfit=yfit) $
     else coeff=[1.,0.]
  
; fit gradient to top of slitfn
;  coeff = linfit(findgen(nn), slitillum, $
;                 measure_error = 1./sqrt(slitillumwt)+100*maxerr*badslitfngrad, $
;                 yfit=yfit)
;  illumarr = yfit
     illumarr=findgen(nn)*coeff[1]+coeff[0]

   endwhile

  if goodct lt 4 then illumarr=slitfntot*0+1.

  yfit=illumarr

; -------- Slit function
; this slitfn has a linear spatial gradient taken out.  Such gradients
; are usually illumination, not slit function. 

  slitfn = slitfntot/illumarr 
  slitfn=slitfn/max(slitfn[edge+2<(nn/2-2):(nn-edge-3)>(nn/2+2)]) > 0.1 
;  badslitfn=abs(slitfn-1) gt 0.03 OR badslitfnsave

  slitfninterior = slitfn[edge:nn-edge-1]
  badslitfn = slitfn*0
  badslitfn[edge:nn-edge-1] = slitfninterior lt 0.7 $
    OR slitfninterior gt 1.4

                                ; refined masking

;  whbad = where(badslitfn, badct)
  badslitfnarr =  (1+intarr(mm)) # badslitfn


  illumarr  = (1+fltarr(mm)) # yfit
  slitfnarr = (1+fltarr(mm)) # slitfn
  badslitfnarr = (1+intarr(mm)) # badslitfn


;NOW: DEAL WITH VARIABLE SLIT FUNCTION!!!!!!!!!!!!!!!!!!!!!
; As usual, the following code is pretty messy, but works even for
; badly vignetted slits

  endif

  if quick le 0 AND nogood eq 0 then begin

     fitdata = flat_processed/illumarr ;/LAMPSPECTROCORR
     fitivar = ivar*(illumarr)^2 ;*LAMPSPECTROCORR)^2

; number of slitfn bins  
     nslitfn = 64
     medslitfarr = fltarr( nslitfn, nn )
     badpix=medslitfarr
     oldslitfn=(1+fltarr(nslitfn)) # slitfn
     yfitarr = medslitfarr
     coeff0=fltarr(nslitfn)
     coeff1=coeff0

     for i = 0, nslitfn-1 do begin
; define the bin and calculate median flux
        minuse = (mm / nslitfn)*i
        maxuse = (minuse+mm/nslitfn) < (mm-1)
        nuse = float(maxuse-minuse)

                                ; normalize out the median in this box
        okpix=where(badpixels[minuse:maxuse,*] eq 0,okct)
        tmpdata=fitdata[minuse:maxuse,*]
        if okct gt 51 then mednorm=median(tmpdata[okpix])

; note that here we use a median slitfn, not mean.  This worked better.
        slitfntot = fltarr(nn)
        slitfnwt = total(fitivar[minuse:maxuse, *]*mednorm^2, 1) 

        for j = 0, nn-1 do begin
           whok=where(badpixels[minuse:maxuse,j] eq 0,okct)
           if okct gt 21 then slitfntot[j] = $
             median(fitdata[minuse+whok, j], /even)/mednorm else $
             slitfntot[j]=median(fitdata[minuse:maxuse,j])/mednorm
           if okct le 20 then begin
              slitfnwt[j]=1.e-20
              badpix[i,j]=1
           endif
        endfor

; pretend it's a weighted average in this line, no biggie  
        
;    whbad=where(slitfntot eq -1,badct)
;    if badct gt 0 then slitfnwt[whbad]=1.e-20
        
        slitillum = slitfntot
        slitillumwt = slitfnwt > 1.e-20
        edge2=edge+4

; deweight slit edges from gradient fit  
        slitillumwt[0:edge2-1] = 1e-20 <  slitillumwt[0:edge2-1]
        slitillumwt[nn-edge2:nn-1] = 1e-20 < slitillumwt[nn-edge2:nn-1]
        
        whbadslitfn=where(slitfn lt 0.8 or slitfn gt 1.2,badct)
        if badct gt 0 then $
          slitillumwt[whbadslitfn]=1.e-20<slitillumwt[whbadslitfn]

        badregion=(nn-total(slitillumwt gt 1.e-20)) lt 5

        


; now that we've thrown out bad pixels, fit gradient to top of slitfn
        if badregion eq 0 then $
                                ;coeff = linfit(findgen(nn), slitillum/slitfn, $
                                ;   measure_error = 1./sqrt(slitillumwt)*slitfn^2, $ ;+100*maxerr*badslitfngrad, $
                                ;   yfit=yfit) $
        coeff = svdfit(findgen(nn), slitillum/slitfn, 2, $
                       measure_error = 1./sqrt(slitillumwt)*slitfn^2, $ 
                       yfit=yfit) $
        else coeff=[1.,0.]

        medslitfarr[i, *] = slitfntot

        coeff0[i]=coeff[0]
        coeff1[i]=coeff[1]
        yfitarr[i, *] = yfit


     endfor


;  medslitfarr=djs_maskinterp(medslitfarr,medslitfarr eq -1,iaxis=0)

     badrows=total(badpix,1) gt nslitfn*.4 AND lindgen(nn) ge edge-2 $
       AND lindgen(nn) le (nn-edge+1)
     arebadrows=total(badrows) gt 0
     if arebadrows then begin
        badmask=(fltarr(nslitfn)+1) # (badrows OR dilate(badrows,bytarr(3)+1)) $
          AND badpix

        badmask=badmask OR (dilate(badmask,intarr(3,3)+1) AND badpix)

        medslitfarr=djs_maskinterp(medslitfarr,badmask,iaxis=1)

;atv,badmask
     endif 
     


; smooth out gradient fit array - to be used in vignetting correction
     for j = 0, nn-1 do yfitarr[*, j] = $
       smooth(djs_median(reform(yfitarr[*, j]), width = 3, $
                         boundary = 'reflect'), 3)

     coeff0r=coeff0
     coeff0s=smooth(djs_median(coeff0, width = 3, $
                               boundary = 'reflect'), 3)

     coeff1s=smooth(djs_median(coeff1, width = 3, $
                               boundary = 'reflect'), 3)

     yfitarr2=yfitarr


     slitfnvig= (rebin(partvig[*,nn/2],mm/nslitfn)) NE 0

     smoothcoeff= (slitfnvig OR dilate(slitfnvig,bytarr(5)+1)) EQ 0
     smoothcoeff[0]=0
     smoothcoeff[mm/nslitfn-1]=0

     whsmooth=where(smoothcoeff ne 0, smoothct)


     if smoothct gt 1 then begin
        coeff0[whsmooth]=coeff0s[whsmooth]
        coeff1[whsmooth]=coeff1s[whsmooth]
     endif



;  slitfnvig=rebin(fullvig,mm/nslitfn,nn)

     
     for i=0,nslitfn-1 do yfitarr2[i,*]=coeff0[i]+coeff1[i]*findgen(nn)


; flatten out before smoothing
     medslitfarr = medslitfarr/yfitarr2
     medslitfarr=(medslitfarr > 0.1) < 1.2


     for j = 0, nn-1 do medslitfarr[*, j] = $
       smooth(djs_median(reform(medslitfarr[*, j]), width = 3, $
                         boundary = 'reflect'),3)


; where we are vignetted, use the nearest unvignetted slit function
; for safety

     slitfnvig=rebin(partvig,mm/nslitfn,nn)
;  slitfnvig=rebin(fullvig,mm/nslitfn,nn)

     nvig=total(slitfnvig,2)
     badslitfn=(nvig gt (nn-15))
;  badslitfn=badslitfn OR dilate(badslitfn,intarr(5)+1)
     badcols=minmax(where(badslitfn))
     allbad=badcols[0] eq 0 and badcols[1] eq nslitfn-1

     if (badcols[0] ne -1) AND (allbad eq 0) then begin
        leftvig=max(badcols) lt (mm/nslitfn - 1)

        if leftvig then badcols[1]=(badcols[1]-1) > 0 else $
          badcols[0]=(badcols[0]+1) < (mm/nslitfn-1)

        if leftvig then medslitfarr[0:badcols[1],*] = $
          (fltarr(badcols[1]+1)+1) # reform(medslitfarr[badcols[1]+1,*]) $
        else medslitfarr[badcols[0]:(mm/nslitfn)-1,*] = $
          (fltarr(mm/nslitfn-badcols[0])+1) # reform(medslitfarr[badcols[0]-1,*])
     endif

     slitfnarr = rebin(medslitfarr, mm, nn)

     varslitfn = medslitfarr


; -------- Fringe pattern
; get things as flat as possible, then median filter
; NOTE THAT A CORRECTION IN THE SPECTRAL DIRECTION IS NOW APPLIED!!!

     flatflat = fitdata/slitfnarr

     flatinterp = djs_maskinterp(flatflat, (ivar eq 0), iaxis=1, /const)
     
; the following is _MUCH_ quicker than a median smooth  
     bgslope= rebin(yfitarr2, mm, nn)
     
     filtered = sigma_filter(flatinterp/bgslope, 11, /all)

; filter out large-scale issues
     background=boxcar2d(filtered,[75,15])

;  background=filter_image(filtered,smooth=31,/all)

; added the division to help deal with near-vignetted regions
     fringe = djs_median(flatinterp/background/bgslope, width=5, boundary='reflect') 

; put all the large-scale parts of the flat into the vigcorr array
; divide by it to correct arcs for vignetting

     vigcorr = bgslope*illumarr*background
     vigcorr = vigcorr/median(vigcorr)

     flat1d = slitfn
     flat2d = fringe
     
; Rescale to the median (just to make sure we're normalized to 1) and
; subtract 1 (to facilitate FLOAT_COMPRESS-ing
     flat2d = flat2d/median(flat2d[mm/2-ncol:mm/2+ncol-1, *])-1. 

; don't muck with bad data
;  if badct gt 0 then slitfnarr[where(badslitfn)] = 1.
     flat2d[*,0]=0.
     flat2d[*,nn-1]=0.
     flat2d[0:16,*]=0.
     flat2d[mm-17:mm-1,*]=0.

; don't throw out data at edges yet!
;  badslitfnarr=badslitfnarr OR slitfnarr LT 0.8 OR slitfnarr GT 1.2

     mask = (fringe GT 0.8) AND (fringe LT 1.2) AND (badslitfnarr EQ 0) 
     mask[0:3, *] = 0B
     mask[mm-4:mm-1, *] = 0B
     if arebadrows then mask=(mask AND (ceil(rebin(float(badmask),mm,nn)) eq 0))

; bad pixels/columns are fixable in a fringe frame - do not mask
     outmask=mask
     mask=mask AND (ivar NE 0)


; where things are bad, simply apply no correction
     wherebad=where(mask eq 0,badct)
     if badct gt 0 then flat2d[wherebad] = 0.
     mask=outmask

   endif else begin

;  badslitfnarr =  (1+intarr(mm)) # badslitfn


;  illumarr  = (1+fltarr(mm)) # yfit
;  slitfnarr = (1+fltarr(mm)) # slitfn
;  badslitfnarr = (1+intarr(mm)) # badslitfn
      if n_elements(slitfn) gt 0 then flat1d = slitfn $
        else flat1d = badslitfngrad*0+1.
      flat2d = ivar*0.

      if n_elements(badslitfnarr) gt 0 then $
        mask = (badslitfnarr eq 0) AND (ivar ne 0) $
        else mask = ((1+intarr(mm)) # badslitfngrad) AND (ivar ne 0) 
      vigcorr = 1.
      varslitfn = flat1d
   endelse


   
   return
  end













