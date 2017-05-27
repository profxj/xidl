;+
;
; NAME
;      extract1d.pro
;
; PURPOSE
;      To extract a 1d spectrum from a sky-subtracted 2-d slit file.  
;      The function extract1d takes a slitfile as its only required
;      argument and returns a structure which contains the
;      1-dimensional spectrum for the given slit. The structure
;      contains the flux, wavelength, and inverse variance along the
;      slit. The 1-d spectrum is extracted according to a
;      user-specified algorithm (either a boxcar inverse
;      variance weighted extraction or an optimal extraction
;      algorithm) with presently the inverse variance extraction set
;      as the default extraction setting. 
;
;
; SYNTAX
;      ss1d = extract1d(slitfile, pos, fwhm, [hdr=hdr, $
;                       /optimal, /boxcar, /horne, /boxsprof,
;                       /nonlocal, nsigma=nsigma])
;
; INPUTS
;      slitfile = a structure containing the sky-subtracted,
;                 cosmic-ray rejected 2-d flux values, the wavelength
;                 solution, and the inverse variance
;                 information. Also, recall that the 2-d spectrum from
;                 a given slit is broken into two slit files: red
;                 chip/blue chip. That is, a single slitfile only
;                 contains the blue portion of the spectrum or the red
;                 end of the spectrum. Can also be a file containing
;                 the structure (a slit file).
;      pos = the position of the object in the slit. pos gives the
;            pixel column about which to extract the object.
;      fwhm = a parameter specifying the width in the spatial
;             direction over which to extract the 1d spectrum. In the
;             optimal extraction algorithm this specifies the
;             full-width at half-maximum (fwhm) for the Gaussian
;             profile.
;      hdr = this optional parameter is provided so that if the user
;            passes a slit structure (rather than file), the header
;            for the slit FITS file can be passed to the routine. This
;            header is employed to locate the slitno and then the
;            corresponding spSlit file.
;
; KEYWORDS
;      optimal = if this keyword is set, then the spectrum is
;                extracted using the optimal extraction
;                algorithm. Note that in this approach we assume that
;                the spatial profile is described by a Guassian. And
;                also not that this is a modified version of the
;                optimal extraction described by K. Horne (see /horne
;                keyword). 
;                Recall that if the /optimal keyword is not set, then
;                the spectrum is extracted according to a boxcar
;                inverse variance weighted technique meaning that
;                each pixel within a given spatial width is weighted
;                according to its inverse variance. That is, outside a
;                given width the weight is 0 and inside the defined
;                spatial range the weight is given by the inverse
;                variance.
;      horne = if this keyword is set, then the spectrum is extracted
;              using the optimal extraction algorithm as detailed by
;              K. Horne (1986, PASP, 98, 609).
;      boxsprof = if this keyword is set, then the spectrum is
;                 extracted using an inverse variance weighted tophat
;                 extraction algorithm. This algorithm differs from
;                 the /boxcar extraction in how it handles
;                 bad-pixels. Here bad-pixels are masked and the flux
;                 missed within these pixels is determined using the
;                 object's spatial profile. The final spectrum is
;                 scaled so as to account for this missed flux.
;      boxcar = this is the default setting. If this keyword, is set
;               then an inverse-variance weighted tophat extraction is
;               employed. For the /boxcar extraction, all bad-pixels
;               will be linearly interplated over and then the
;               extraction is performed.
;      nonlocal = if this keyword is set, then the extracted 1-d
;                 spectrum is extracted from the non-local
;                 sky-subtracted 2-d slit spectrum.
;
; OUTPUTS
;      ss1d = a structure containing the 1-d extracted spectrum:
;             spec.spec = a vector containing the flux values as a
;                         function of wavelength. Note that the
;                         spec.flux array has curious dimensions. 
;             spec.ivar = a vector containing the inverse variance of
;                         the extracted 1-d spectrum. 
;             spec.lambda = a vector containfng the wavelength values
;                           corresponding to the the flux and inverse
;                           variance values given in spec.flux and
;                           spec.ivar. 
;             spec.ivarfudge = an empirical factor by which to
;             multiply spec.ivar to make it match observed
;             fluctuations (in a median sense).  multiplying by
;             (spec.ivarfudge < 1) may be advisable in many situations.
;
; PROCEDURES CALLED 
;      lambda_eval
;      find_object
;      mcc_polyfit
;
; EXAMPLES
;      None.
;
; HISTORY
;      Created June 27, 2002 by mcc.
;      Revised July 8, 2002 by mcc - incorporated additional code
;         (courtesy of MD) to calculate the wavelength solution from
;         the slitfile information.  
;      Revised July 10, 2002 by mcc - added lines to check if proper
;         number of parameters were passed to function; added optimal
;         extraction algorithm. I also cleaned up the default boxcar
;         extraction and the inverse variance weighted extraction so
;         that iteration (via a DO loop) is no longer employed. Note
;         that presently the optimal extraction algorithm is
;         hard-wired to exclude the first 5 and last 5 pixels in the
;         spatial direction when fitting to a Gaussian profile.
;      Revised July 16, 2002 by mcc - added keyword /sp so that the
;         user has the option of not using or using the spSlit file to
;         determine the location of the object in the slit. If the /sp
;         keyword is set, then extract1d.pro will search for the
;         corresponding spSlit file and use it. Eventually this
;         feature will be removed...eventually being when the sky
;         subtraction if good enough. Also, I removed the keyword
;         ivar_weight and made the inverse variance weighted
;         extraction the default algorithm. 
;      Revised July 19, 2002 by mcc - revised code so that the tilt in
;         the 2d spectrum is appropriately accounted for in the iverse
;         variance wieghted extraction algorithm. Note that this is
;         not yet properly handled in the optimal extraction
;         technique. Also, the position of the object and boxcar width
;         are now pulled from the design specifications of the slit
;         mask. 
;      Revised July 22, 2002 by mcc - made changes to the optimal
;         extraction algorithm in an attempt to account for the tilt
;         of the 2d spectra. 
;      Revised July 26, 2002 by mcc - entire format of function was
;         changed. It now takes the fwhm and pos positions and only
;         does the extraction. This simplified the routine
;         greatly. Also, it no longer takes a slitfile as its first
;         argument, but instead takes a slit structure.
;      Revised August 5, 2002 by mcc - now routine takes a slitfile or
;         a slit structure as first argument. Also, due to changes in
;         the format of the structures stored in slit files, the
;         manner in which tilts in the spectra are corrected has been
;         revised. 
;      Revised August 12, 2002 by mcc - omproved optimal extraction
;         technique so that it now properly handles the tilt in 2d
;         spectra. In the future, this may be revised to use a
;         interpolation algorithm rather than the simply SHIFT
;         technique presently employed.
;      Revised Sep. 18, 2002 by JAN - properly uses old or new
;         slit.MASK/slit.BITMASK
;      Revised Oct 17, 2002 by MD - interpret over mask for missing
;         data.
;      Revised March 15, 2003 by mcc & JAN - added two new extraction
;         techniques (horne and boxsprof). Also added code to
;         incorporate the uncertainties in the b-spline sky model.
;
;-

; estimate how far off our ivar estimate is, by comparing deviations
; from median-smooth version of spectrum to expected deviations
function ivarfudge, spectrum, ivar

  nspec=n_elements(spectrum)
  nbins=8

  fudgearr=fltarr(nbins)
  binarr=fudgearr

  for i=0,nbins-1 do begin
      whbin = where(lindgen(nspec) ge (nspec/8.*i>50) and $
                    lindgen(nspec) le  (nspec/8.*(i+1)<(nspec-50)) and $
                    ivar ne 0. and spectrum ne 0.,binct)
      binarr[i]=binct
      if binct gt 100 then begin
          medlevel = median(spectrum,50)
          fudgearr[i] = djsig((spectrum[whbin]-medlevel[whbin]) * $
                              sqrt(ivar[whbin]))
      endif else fudgearr[i] = 0.
  endfor

  whok=where(binarr gt 200,okct)

  if okct gt 0 then $
    fudgefactor=1/(median(fudgearr[whok],/even))^2 $
  else fudgefactor=1.
;  if okct gt 0 then print,median(fudgearr[whok],/even)

  return,fudgefactor
end



function extract1d, slitfile, pos, fwhm, nsigma=nsigma, $
                    optimal=optimal, horne=horne, boxsprof=boxsprof, $
                    nonlocal=nonlocal, hdr=hdr
; check that the proper number of parameters was passed.
  if n_params() lt 3 then begin
      print, 'CALLING SEQUENCE: ss1d = extract1d(slitfile, pos, fwhm, '+ $
        'nsigma=nsigma, /optimal, /horne, /boxsprof, /nonlocal)'
      retall
  endif

; check the various keyword settings.
  if n_elements(optimal) gt 0 then $
    optimal = optimal[0] ge 1 else optimal = 0
  if n_elements(horne) gt 0 then $
    horne = horne[0] ge 1 else horne = 0
  if n_elements(boxsprof) gt 0 then $
    boxsprof = boxsprof[0] ge 1 else boxsprof = 0
  if n_elements(nonlocal) gt 0 then $
    nonlocal = nonlocal[0] ge 1 else nonlocal = 0
  if (optimal eq 0) and (horne eq 0) and $
    (boxsprof eq 0) then boxcar = 1 else boxcar = 0

; check if the argument slitfile is a FITS file or an IDL structure
; and then define slit accordingly.
  if size(slitfile, /tname) eq 'STRING' then begin
      if keyword_set(nonlocal) then begin
          fits_info, slitfile, /silent, n_ext=n_ext
          if n_ext gt 2 then begin
              slit = mrdfits(slitfile, 3, hdr, /silent)
          endif else return, 0
      endif else slit = mrdfits(slitfile, 1, hdr, /silent)
      slitno = strcompress(sxpar(hdr, 'SLITNO'), /rem)
  endif else begin
      slit = slitfile
      if keyword_set(hdr) then $
        slitno = strcompress(sxpar(hdr, 'SLITNO'), /rem) $
      else slitno = ''
  endelse

; construct an empty spectrum: filled with zeros. this will be
; returned as the result in case of various errors in the extraction.
  n = n_elements(slit.flux[*,0])

  skyspec=fltarr(n)
  usesky=0
  nsky=0
  readnoise = 2.32


; JAN first sky code here
  if slitno ne '' then begin
; find the spSlit file corresponding to the giv4en slit 7file/structure.
      endpos = strpos(slitfile, '.fits')
      spname = 'spSlit' + strmid(slitfile, 4, endpos-4) + '.fits*'
      spslitfile = findfile(spname, count=numsp)
      if numsp gt 0 then begin
          spslitfile = spslitfile[0]
          fits_info, spslitfile, /silent, n_ext=n_ext
          nbsplines=(n_ext / 2)
          usesky=1
          skyarr=fltarr(n,nbsplines)
          spslit=mrdfits(spslitfile, 1, /silent)
                                ; keep track of # of good sky pixels
                                ; in each column.  this should roughly
                                ; correspond to the extracted
                                ; pixels... INCORPORATE CRMASK!
          skydex = where(spslit.skyrow, skydexcnt)
          if skydexcnt gt 0 then $
            nsky=total((spslit.mask[*,where(spslit.skyrow)] AND 22b) eq 0b,2) $
          else begin
              print, 'No skyrows in slit!'
              nsky = 0
          endelse
      endif
      exptime=0.
  endif


; define a 1-d spectrum structure to return in case of error.
  zero1d = {spec:fltarr(n), lambda:fltarr(n), ivar:fltarr(n), $
            crmask:intarr(n), bitmask:intarr(n), ormask:intarr(n), $
            nbadpix:intarr(n), infomask:intarr(n), $
            objpos:float(pos), fwhm:float(fwhm), $
            nsigma:0.0, r1:0, r2:0, skyspec:fltarr(n), ivarfudge:1.} 

; check that the pos argument is valid. if not valid, then return a
; spectrum of all zeros and print an error message.
  nrows = n_elements(slit.flux[0,*])
  if pos lt 0 or pos gt nrows-1 or not(finite(pos)) then begin
      if size(slitfile, /tname) eq 'STRING' then begin
          hdr = headfits(slitfile, ext=1, /SILENT)
          slitn = string(sxpar(hdr, 'SLITNO'))
          print, '(extract1d.pro) ERROR: Invalid position (' + $
            strcompress(string(pos), /rem) + $
            ') for slit ' + strcompress(slitn, /rem) + '!!!'
      endif else print, '(extract1d.pro) ERROR: Invalid position (' + $
        strcompress(string(pos), /rem) + ')!!!'
      return, zero1d
  endif
; check that the fwhm argument is valid. if not valid, then set fwhm =
; 10 and print an error message.
  if fwhm le 0.0 or not(finite(fwhm)) then begin
      if size(slitfile, /tname) eq 'STRING' then begin
          hdr = headfits(slitfile, ext=1, /SILENT)
          slitn = string(FXPAR(hdr, 'SLITNO'))
          print, '(extract1d.pro) ERROR: Invalid fwhm (' + $
            strcompress(string(fwhm), /rem) + $
            ') for slit ' + strcompress(slitn, /rem) + '!!!'
      endif else print, '(extract1d.pro) ERROR: Invalid fwhm (' + $
        strcompress(string(fwhm), /rem) + ')!!!'
      fwhm = 10.
  endif

; get the 2-d lambda array.
  lambda2d = lambda_eval(slit)
; extract the wavelength values along the center of the object.
  cwave = lambda2d[*,pos]
  npix = n_elements(cwave)
; determine the dispersion level (angstrom/pixel) along the slit.
  dldx = (cwave[npix-1] - cwave[0]) / npix
  tiltx = (lambda2d[npix/2,nrows-1] - lambda2d[npix/2,0]) / nrows
  dxdp = tiltx/dldx ;slope of constant lambda/vert pix


; build up sky spectrum, for use in determining its inverse variance
; JAN - sky ivar code here
  if usesky gt 0 then begin
      for i=0,nbsplines-1 do begin
          bspline=mrdfits(spslitfile,i*2+2, /silent)
          skyarr[*,i]=bspline_valu(cwave,bspline)
          exptime=exptime + $
            sxpar(headfits(spslitfile,ext=i*2+1, /silent),'EXPTIME')
      endfor

; assume average spec. is ok for statistics, and that don't have to
; worry about bad columns, etc. in sky region.  
      avgexptime=exptime/nbsplines
      if nbsplines gt 1 then skyspec=total(skyarr,2)/nbsplines $
      else skyspec = skyarr
; note the bspline spectrum is an average, not a sum, of the values at
; each pixel that went into it.
      skyctsper2dpixperspec=skyspec*(avgexptime/3600.)
      skyivar= nsky*nbsplines/(skyctsper2dpixperspec+readnoise^2)
      skyivar=skyivar*(avgexptime/3600.)^2

; maybe add...
;      skyivar=skyivar*(2./3.) ; as typically are using ~nsky*2/3 pixel
        ; to define the sky value at a breakpoint. only slightly kludgey.

; this is the inverse-variance in the sky spectrum at each pixel of
; the extraction.  Note that this contributes, wholly covariantly, to
; EVERY pixel at a given wavelength.
  endif else begin
      print, 'No sky info found!!!'
      skyivar = fltarr(n) + 1.0
  endelse

; locate bad pixels in the flux array and interpolate over them.
; interpolate across CTE problems, bad spots in pixmap, in spatial
; direction. do not interpolate across vignetted regions.
  interpolate = (slit.mask and 22b) gt 0
; if the boxsprof keyword is set, then do NOT interpolate. just set
; the bad pixels to zero and we will compensate for them using the
; spatial profile. for the horne and optimal extractions, set the ivar
; for interpolated pixels to zero.
  if boxsprof then flux2d = slit.flux * (1 - interpolate) $
  else flux2d = djs_maskinterp(slit.flux, interpolate, iaxis=1, /const)
; interpolate the ivar values as well.
  if boxsprof then ivar2d = slit.ivar * (1 - interpolate) $
  else ivar2d = djs_maskinterp(slit.ivar, interpolate, iaxis=1, /const)
  interpix = where(interpolate, inter_cnt)
; increase the variance in the interpolated pixels, but exact factor
; depends on number of rows interpolated over.
  if boxcar then begin
      if inter_cnt gt 0 then ivar2d[interpix] = ivar2d[interpix]/4.
  endif else begin
      if inter_cnt gt 0 then ivar2d[interpix] = 0.
  endelse


; ----------------------------
; if keyword optimal or horne is passed then extract the 1-d spectrum
; according to one of the two variations of the optimal extraction
; algorithm. for more see Horne, K. 1986, /pasp, 98, 609.
  if optimal or horne then begin
; define range in the spatial direction [r1:r2] over which to do
; extraction. exclude the first and last few pixels due to slit edge
; effects.
      r1 = 4
      r2 = nrows-5
; estimate the spatial profile map (P) by modeling the spatial profile
; as a gaussian. let's define the width (sigma) of the gaussian
; according to the fwhm of the spatial profile as measure by the
; find_object routine. and take the center of the gaussian (x0) to be
; at the object position pos (as determined using find_object + peakinfo).
      sigma = fwhm / 2.35482
      x0 = pos
; now do a linear least-squares fit to solve for the continuum and
; amplitude of the gaussian y = a0 + a1*e^(-(x-x0)^2/(2*sigma)).
;      x = (findgen(nrows))[r1:r2] 
;      xdata = exp( -(x - x0)^2 / (2.*sigma^2) )
;      sprof = find_object(slitfile) ;, /CR, /BPM, /NOSUBTRACT, /USETILT)
;    mcc_polyfit, xdata, sprof[r1:r2], [0, 1], a=a
; finally, construct the spatial profile map (P) according to the
; gaussian profile as inferred from the spatial profile given by
; find_object. 
      nlams = n_elements(slit.flux[*,0])
      x = findgen(nrows)
      xdata = exp( -(x - x0)^2 / (2.*sigma^2) )
      P = (xdata) ## (fltarr(nlams)+1.) 
; and make sure that P is >= 0.
;      P = P > 0
      P = abs(P)
; then normalize P.
      Ptot = total(P, 2)
      P = P / (Ptot # (fltarr(nrows)+1.) )
;      for i=0,nlams-1 do P[i,*] = P[i,*] / Ptot[i]



; construct a bad pixel mask (a zero value in the slit.badmask
; corresponds to a good pixel).
;    M = ABS( (slit.mask - 1) < 0 )

; WE REALLY SHOULDN'T HAVE BEEN DISABLING THE BAD PIXEL MASK FOR THE
; OPTIMAL.  _REALLY_.
;      M = slit.mask*0. + 1.
      M = (slit.mask EQ 0b)
; now define some arrays of the proper length which we will fill with
; the extracted 1-d spectrum, ivar, and bad pixel masks.
      nbits = n_elements(slit.flux[*,0])
      spec_num = fltarr(nbits)
      spec_denom = fltarr(nbits)
      var_num = fltarr(nbits)
      sky_num=var_num
      var_denom = fltarr(nbits)
      crmask = intarr(nbits)
      bitmask = intarr(nbits)
      ormask = intarr(nbits)
      infomask = intarr(nbits)
      nbadpix = intarr(nbits)
      var2d = slit.ivar*0. 
      ivar = fltarr(nbits)
      spec = ivar
; extract the 1-d spectrum. iterate row by row and account for the
; tilt. define the range (in spatial direction) over which to do
; extraction. take only those pixel columns for which we are within
; nsigma*sigma from the object's position (pos).
; 1.5sigma ---> 0.12952 Bevington page 251.
; 2.0sigma ---> 0.05400
      if keyword_set(nsigma) then nsigma=nsigma[0] else nsigma=1.5
      Pgauss = 1.0 / sqrt(2.0 * !pi) / sigma * exp( (-1.0/2.0) * nsigma^2)
      incs = where(P[0,*] ge Pgauss, incnt)
      if incnt gt 4 then begin
          r1 = min(incs)
          r2 = max(incs)
      endif else begin
          print, 'ERROR: not enough rows selected in optimal extraction!!'
          return, zero1d
      endelse

; begin iteration...
      for i=r1,r2 do begin
; determine the shift needed to rectify the rows in wavelength.
        cshift = round((i-pos)*dxdp)
; find all of the points in the ith row which will be improperly
; wrapped when vector shift operations are applied.
        badpts = where(cwave[0] gt lambda2d[*,i] or $
                       cwave[nbits-1] lt lambda2d[*,i], bcnt)
; at all points where we improperly wrap, set the ivar to zero.
        if bcnt gt 0 then ivar2d[badpts,i] = 0.0
; fill the var2d array with ivar at the good points....exclude the
; points where ivar=0 since this will blow up upon inversion. the
; variance values at the bad pixels will be set at zero, but this
; isn't a concern since these variance values are not used anywhere.
        good = where(ivar2d gt 0, goodcnt)
        if goodcnt gt 0 then var2d[good] = 1./ivar2d[good]
; now perform the optimal extraction (truly a modified optimal
; extraction) or the horne extraction (which is a true optimal
; extraction). recall the these vectors (spec2d_num, spec2d_denom,
; etc.) are initially filled with zeros. 
        if optimal then begin
            spec_num = spec_num + $
              shift(M[*,i] * P[*,i] * flux2d[*,i], cshift) 
            spec_denom = spec_denom + $
              shift(M[*,i] * (P^2)[*,i], cshift) 
            var_num = var_num + $
              shift(M[*,i] * P[*,i]^2 * var2d[*,i], cshift) 
            var_denom = var_denom + $
              shift(M[*,i] * (P[*,i])^2, cshift) 
            sky_num = sky_num + $
              shift(M[*,i] * (P[*,i]), cshift) 
        endif
        if horne then begin
            spec_num = spec_num + $
              shift(M[*,i] * P[*,i] * flux2d[*,i] * ivar2d[*,i], cshift) 
            spec_denom = spec_denom + $
              shift(M[*,i] * (P^2)[*,i] * ivar2d[*,i], cshift) 
            var_num = var_num + $
              shift(M[*,i] * P[*,i]^2* ivar2d[*,i], cshift) 
            var_denom = var_denom + $
              shift(M[*,i] * (P[*,i])^2 * ivar2d[*,i], cshift)
            sky_num = sky_num + $
              shift(M[*,i] * (P[*,i]) * ivar2d[*,i], cshift)  
        endif
; determine the number of pixels at each wavelength which have a CR
; hit in them.
        crmask = crmask + shift(slit.crmask[*,i] gt 0, cshift)
; take the minimum bitmask value (along the row) at each wavelength.
        bitmask = bitmask and shift(slit.mask[*,i], cshift)
; also take the OR of the bitmask.
        ormask = ormask or shift(slit.mask[*,i], cshift)
; determine the number of badpixels at each wavelength.
        nbadpix = nbadpix + shift(slit.mask[*,i] ne 0b, cshift)
; take the OR of the infomask.
        infomask = infomask or shift(slit.infomask[*,i], cshift)
    endfor
; now calculate the final spectrum. only do this where things wont
; blow up, that is, when the denominator is non-zero. at places where
; we have a denominator equal to zero, the spectrum will be zero.
    gpix = where(spec_denom gt 0, goodcnt)
    if goodcnt gt 0 then $
      spec[gpix] = spec_num[gpix] / spec_denom[gpix] $ ;otherwise 0.
    else print, '(extract1d.pro) ERROR: No good pts in denominator ' + $
      'of spectrum for slit ' + slitno + '!!!'
; and calculate the final ivar array. only do this where things wont
; blow up (just like for spec).
    gpix = where(var_denom ne 0, goodcnt)
    if goodcnt gt 0 then begin
; I _think_ I got the propagation of errors right...
        if optimal then $
          ivar[gpix]= 1/ (var_num[gpix]/(var_denom[gpix])^2 + $
               sky_num[gpix]^2/var_denom[gpix]^2/skyivar[gpix])
        if horne then $
          ivar[gpix]= 1 / (var_num[gpix]/var_denom[gpix]^2 + $
               sky_num[gpix]^2/var_denom[gpix]^2/skyivar[gpix])
    endif else print, '(extract1d.pro) ERROR: No good pts in numerator ' + $
      'of variance or denominator of spectrum for slit ' + slitno + '!!!'
; construct the output structure and return it.
    ss1d =  {spec:spec, lambda:cwave, ivar:ivar, crmask:crmask, $
             bitmask:bitmask, ormask:ormask, nbadpix:nbadpix, $
             infomask:fix(infomask), $
             objpos:float(pos), fwhm:float(fwhm), $
             nsigma:nsigma, r1:r1, r2:r2, skyspec: skyspec, $
             ivarfudge:ivarfudge(spec,ivar)}
    return, ss1d

  endif else begin
; ----------------------------
; else extract the 1-d spectrum according to the boxcar or boxsprof
; extraction alogorithms.

; define some arrays of the proper length which we will fill with the
; extracted 1-d spectrum, ivar, and bad pixel masks.
    nbits = n_elements(slit.flux[*,0])
    spec = fltarr(nbits)
    count = spec
    ivar = fltarr(nbits)
    var = fltarr(nbits)
    var2d = slit.ivar*0.
    crmask = intarr(nbits)
    bitmask = intarr(nbits)
    ormask = intarr(nbits)
    infomask = intarr(nbits)
    nbadpix = intarr(nbits)
; if the keyword nwidth is set, then multiply the fwhm by nwidth. this
; will determine the extraxtion width spatially.
    if keyword_set(nsigma) then nsigma = nsigma[0] $
    else nsigma = 1.5
    width = fwhm * nsigma

; do a flux conserving boxcar extraction...no longer inverse variance
; weighted. iterate and column by column apply the appropriate shift
; so as to remove any tilt in the spectrum.
; determine the pixels (in the spatial direction) over which to
; extract.
    ext = (ceil(width/2.)) < pos < (nrows-1-pos)
    r1 = (round(pos-ext)) > 4
    r2 = (round(pos+ext)) < (nrows-5)
    if r2 lt r1 then begin
        print, '(extract1d.pro) ERROR: invalid extraction extrema defined!'
        print, '(r1, r2, pos):', r1, r2, pos
        return, zero1d
    endif
; but before we begin to iterate, let's calculate the scaling factor
; (from the spatial profile) to be used in the boxsprof extraction. 
    if keyword_set(boxsprof) then begin
; create a guassian spatial profile over the row limits r1:r2.
        sigma = fwhm / 2.35482
        x0 = pos
        x = findgen(r2-r1+1) + r1
        sprof1d = exp( -(x - x0)^2 / (2.*sigma^2) )
; determine at each wavelength the portion of this spatial profile
; that was designated as bad (set to be interpolated over in the
; standard boxcar extraction). calculate the factor by which we must
; multiply the extracted spectrum to account for this missing flux
; (the flux in the interpolated pixels was set to zero for the
; boxsprof algorithm...see above code).
        sprof2d = sprof1d ## (fltarr(nbits)+1.0)
        if n_elements(sprof1d) gt 1 then $
          intprof = total(sprof2d * interpolate[*,r1:r2], 2) $
        else intprof = sprof2d * interpolate[*,r1:r2]
        intscl = intprof / total(sprof1d)
; this failed at one end!!!!!!!!!!! added a trap for now.JAN
        intscl = 1.0 / (1.0 - intscl)
        whbad=where(finite(intscl) eq 0,badct)
        if badct gt 0 then intscl[whbad]=1.
    endif

    for i=r1,r2 do begin
; determine the shift needed to rectify the rows in wavelength.
        cshift = round((i-pos)*dxdp)
; find all of the points in the ith row which will be improperly
; wrapped when vector shift operations are applied to rectify the slit
; in wavelength.
        badpts = where(cwave[0] gt lambda2d[*,i] or $
                       cwave[nbits-1] lt lambda2d[*,i], bcnt)
; at all points where we improperly wrap, set the ivar to zero.
        if bcnt gt 0 then ivar2d[badpts,i] = 0.
; fill the var2d array with ivar at the good points....exclude the
; points where ivar=0 since this will blow up upon inversion. the
; variance values at the bad pixels will be set at zero, but this
; isn't a concern since these variance values are not used anywhere.
        good = where(ivar2d gt 0, gcnt)
        if gcnt gt 0 then var2d[good] = 1./ivar2d[good]
; add the flux in the spatial direction correcting for the shift.
        spec = spec + shift(flux2d[*,i], cshift)
; sum the variance across the rows, too...keeping track of how many
; good pixels are present at each wavelength. recall that var2d = 0 at
; all bad pixels since we initiated the var2d array as a fltarr filled
; with zeros and only filled it at the good pixels.
        var = var + shift(var2d[*,i], cshift)
        count = count + shift(var2d[*,i] gt 0, cshift)
; determine the number of pixels at each wavelength which have a CR
; hit in them.
        crmask = crmask + shift(slit.crmask[*,i] gt 0, cshift)
; take the minimum bitmask value (along the row) at each wavelength.
        bitmask = bitmask and shift(slit.mask[*,i], cshift)
; also take the OR of the bitmask.
        ormask = ormask or shift(slit.mask[*,i], cshift)
; determine the number of badpixels at each wavelength.
        nbadpix = nbadpix + shift(slit.mask[*,i] ne 0b, cshift)
; take the OR of the infomask.
        infomask = infomask or shift(slit.infomask[*,i], cshift)
    endfor

; now for the boxcar account for any missing rows by scaling the
; spectrum (at that wavelength) by the number of missed pixels (thus
; giving each pixel equal weight).
    if boxcar then begin
; ADD CONTRIBUTION FROM SKY
        var=var+count^2/skyivar

        gpix = where(count gt 0, goodcnt)
        if goodcnt gt 0 then $
          spec[gpix] = spec[gpix] / count[gpix] * (r2-r1+1)
; convert the variance into an inverse variance.
        if goodcnt gt 0 then $
          ivar[gpix] = 1.0 / var[gpix] * (count[gpix] / (r2-r1+1))^2
    endif

; and for the boxsprof account for any missing pixels by scaling the
; spectrum according to the spatial profile.
    if boxsprof then begin
        gpix = where(var gt 0, goodcnt)
        var=var+count^2/skyivar

; calculate the ivar values from the variance array.
        if goodcnt gt 0 then ivar[gpix] = 1.0 / var[gpix]
; finally, scale the 1-d spectrum by the scale factor.
        spec = spec * intscl
; acount for this scaling of the spectrum in the ivar array.
; THIS HAD BEEN /intscl - REMEMBER, IVAR GOES AS 1/factor^2 WHEN YOU
;                         MULTIPLY BY FACTOR
        ivar = ivar / intscl^2
    endif

; construct the output structure and return it.
    ss1d = {spec:spec, lambda:cwave, ivar:ivar, $
            crmask:crmask, bitmask:bitmask, $
            ormask:ormask, nbadpix:nbadpix, $
            infomask:fix(infomask), $
            objpos:float(pos), fwhm:float(fwhm), $
            nsigma:nsigma, r1:r1, r2:r2, $
            skyspec:skyspec, ivarfudge:ivarfudge(spec,ivar)}

    return, ss1d

  endelse

end








