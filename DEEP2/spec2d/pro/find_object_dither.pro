;+
;
; NAME
;       find_object
;
; PURPOSE
;       The function find_object takes a [sp]slit filename
;       or structure as its only
;       required argument and determines the spatial profile of the
;       light from the object (source under study) along the slit. If
;       optional keyword parameters are provided, then find_object
;       will mask bad pixels and/or do cosmic-ray rejection. 
;
;
; SYNTAX
;       sprof = find_object([slitfile, profivar=profivar, $
;           npix=npix, /CR, /BPM, /NOSUBTRACT, /MMM,/USETILT,$ 
;           /SKYPROFILE, pixrange=pixrange,modeout=modeout]
;
; INPUTS
;       slitfile = the (string) name of a slitfile or spSlitfile, or
;                  alternatively the name of a [sp]slit structure. 
;                  If it is a filename, slitfile
;                  will be read into a structure. From the structure,
;                  the flux, inverse variance, and wavelength
;                  information will be extracted. The flux information
;                  will be in the form of a two-dimensional array
;                  containing the flux along the slit where flux[*,i]
;                  is the flux as a function of wavelength along the
;                  ith column of the slit. The wavelength and inverse
;                  variance values corresponding to each pixel will be
;                  established in similar array structures.
;       profivar = a variable that will be set equal to the inverse
;                  variance of the derived spatial profile for the
;                  slit. 
;       npix = the number of pixels used in each row of the sum.
;              Useful for estimating S/N per pixel.
;
; KEYWORDS
;       NOTE: find_object will by default choose reasonable values for 
;       unset keywords based on whether the input structure is a slit or
;       spSlit structure.  Note that for an spslit structure,
;       NOSUBTRACT is not set by default.  Otherwise, there is
;       hopefully little need to use these keywords.
;
;       /CR = a keyword passed to the function. If this keyword is set,
;              then cosmic ray rejection will be preformed - this is
;              very beneficial when attempting to locate the object in
;              a spSlit file.
;       /BPM = if this keyword is passed, then the bad pixel mask in
;              the slitfile will be applied. 
;       /NOSUBTRACT = if this keyword is set, then a continuum level
;                     will not be subtracted from the extracted
;                     spatial profile. This is useful for running
;                     find_object on slit files vs. spSlit files. In
;                     the case of slit files we would want to pass
;                     this keyword since sky-subtraction has already
;                     been performed (and vice-versa for spSlit
;                     files).  
;       /MMM : set to subtract modal value instead of minimum for sky
;       /USETILT: set to use the wavelength solution to mask tilted
;       sky lines, rather than mask a larger region around all sky
;       lines (with the size of the masked region in that case
;       determined from the difference in wavelength at top and bottom
;       of slit) ; i.e., use a larger region with tilted edges,
;       instead of the default smaller region with straight edges.
;       /SKYPROFILE: set to make the profile on sky lines, instead of
;       off sky lines.  This might help us look at the slit function
;       issues.
;       PIXRANGE=PIXRANGE: a 2-element vector giving the min/max
;       columns (along the spectral direction) to sum over.  [0,4095] is
;       used by default.
;       MODEOUT= MODEOUT: returns the sky value subtracted off
;
; OUTPUTS
;       sprof = the spatial flux profile of the object in the
;               slit. sprof will be a vector with length equal to the
;               number of columns in the slit. 
;       additional outputs: profivar, npix (see INPUTS).
;
; PROCEDURES CALLED 
;       lambda_eval.pro
;
; EXAMPLES
;       IDL> file = $ 
;       IDL> '/deepscr2/marc/DEIMOS/1145/1145/spSlit.1145.145R.fits.gz'
;       IDL> sprof = find_object(slitfile=file, profivar=profivar, /CR)
;       IDL> PLOT, sprof[*,0], sprof[*,1], THICK=2, $
;       IDL> XTITLE="COLUMN NUMBER", YTITLE="FLUX", $ 
;       IDL> TITLE="SLIT PROFILE"
;
; HISTORY
;       Created June 14, 2002 by jnewman and mcc.
;       Revised July 9, 2002 by mcc. The revision changed the input
;          paramaters such that a slitfile was taken as the only
;          required argument/parameter. find_object was also changed
;          from a procedure to a function. 
;       Revised July 12, 2002 by mcc. 
;       Revised heavily Sep. 20, 2002 by jan
;
;-

FUNCTION find_object_dither, slitfile,  $
                      profivar=profivar, npix=npix, $
                      cr=cr, bpm=bpm, nosubtract=nosubtract, mmm=mmm, $
                      usetilt=usetilt, skyprofile=skyprofile, atv=atv,  $
                      pixrange = pixrange,ivarweight=ivarweight,$
                      modeout=modeout, align=align, quick=quick


;;; CHECK THAT ENOUGH PARAMETERS WERE SUPPLIED 
  IF N_ELEMENTS(slitfile) eq 0 THEN BEGIN
    print, 'ERROR: NOT ENOUGH PARAMETERS PASSED TO find_object!' 
    print, 'You must pass name of slitfile or a slit structure!'
    RETURN, 0
  ENDIF

  IF n_elements(usetilt) eq 0 then usetilt = 1
  if n_elements(skyprofile) eq 0 then skyprofile = 0
  if n_elements(interp) eq 0 then interp = 1
  if n_elements(atv) eq 0 then atv = 0
  if n_elements(pixrange) eq 0 then pixrange = [0, 4095]

  if n_elements(ivarweight) eq 0 then ivarweight = 0

  if n_elements(quick) eq 0 then quick = 0
  if quick then interp = 0

  
  mincol = min(pixrange) > 0
  maxcol = max(pixrange) <  4095
  if n_elements(align) eq 0 then align = 0

  bspline = 0
  on_error, 3

; decide how to subtract off BG


 

 tempname = slitfile

;;; READ-IN THE slitfile 
  IF n_elements(slitfile) NE 0 THEN BEGIN
     IF size(slitfile, /type) eq 7 then slit = MRDFITS(slitfile, 1, /SILENT) $
        else BEGIN 
          slit = slitfile
          tempname= ' '
        ENDELSE
     IF size(tempname, /type) eq 7 and skyprofile and strpos(tempname, 'spSlit') ge 0 then begin
        sset = mrdfits(slitfile, 2, /SILENT)
        bspline = 1
     endif
    flux = slit.flux
;;; ESTABLISH THE WAVELEGNTH INFORMATION 
    lambda = lambda_eval_reduced_dither(slit)
  ENDIF
;;; DETERMINE NUMBER OF COLUMNS IN SLIT.
  ncols = N_ELEMENTS(lambda[0,*])

  lambdarange=minmax(lambda)

; choose whether to use 600-line or 1200-line skybpm file
  use600 = lambdarange[0] lt 6300 OR lambdarange[1] gt 9500. $
    OR (lambdarange[1]-lambdarange[0]) gt 1500.

  if use600 then stem='g600' else stem='g1200'



;;; SET CALIBRATION DIRECTORY.
 calibdir = getenv('CALIB_DATA')
;;; IF 'CALIB_DATA' IS NOT DEFINED, LOOK IN 
;;; THE CURRENT DIRECTORY.
 IF strlen(calibdir) NE 0 THEN calibdir = calibdir+'/'

;;; READ-IN THE SKY BAD PIXEL MASK skybpm. 
   if skyprofile then skybpmfile = calibdir+'brightskybpm.' $
          +stem+'.fits.gz' $
   else skybpmfile = calibdir+'skyfreebpm.' $
          +stem+'.fits.gz'

   if findfile(skybpmfile) ne skybpmfile then skybpmfile= $
     calibdir+'skyfreebpm.fits.gz'
   
    skybpm=readfits(skybpmfile,/SILENT)

  if total(tag_names(slit) eq 'CRMASK') OR total(tag_names(slit) eq 'BITMASK') then begin
    ; settings for a slit file/structure
    if n_elements(nosubtract) eq 0 then nosubtract = 1
    if n_elements(mmm) eq 0 then mmm = 0
    if n_elements(bpm) eq 0 then bpm = 1
    if n_elements(cr) eq 0 then cr = 0
  endif else begin
    ; settings for an spSlit file/structure
    if n_elements(nosubtract) eq 0 then nosubtract = 0
    if n_elements(bpm) eq 0 then bpm = 1
    if n_elements(cr) eq 0 then cr = 1
    if n_elements(mmm) eq 0 AND NOT nosubtract then mmm = 1 else mmm = 0
  endelse

 if n_elements(nosubtract) eq 0 then nosubtract = 0
 if n_elements(mmm) eq 0 then mmm = 0

  bitmask = byte(slit.mask)
  ivar = slit.ivar 
  nearvig = (slit.infomask AND 1b) eq 1b
  crmaskin = (slit.infomask AND 8b) eq 8b

  if maxcol-mincol gt 100 then begin
    flux = flux[mincol:maxcol, *]
    lambda = lambda[mincol:maxcol, *]
    bitmask = bitmask[mincol:maxcol, *]
    ivar = ivar[mincol:maxcol, *]
    nearvig = nearvig[mincol:maxcol, *]
  endif else message, 'too small region selected-using whole slit', /INFO



;;; EXTRACT WAVELENGTH VALUES (wvlen) ALONG CENTER OF SLIT.
  wvlen = lambda[*,ncols/2]

;;; DETERMINE WAVELENGTH RANGE ALONG SLIT.
  lambdarange = MINMAX(lambda[mean(pixrange)-mincol, *])

  difflambda = lambdarange[1] - lambdarange[0]

;;; DEFINE A SMOOTHING SCALE AND BOXCAR SMOOTH 
;;; THE skybpm. THIS IS DONE SO THAT WE ARE SURE TO
;;; BLOCK OUT ALL AREAS OF THE SPECTRUM AFFECTED BY SKY
;;; LINES EVEN IN THE CASE OF TILTED SLITS.
;  Note that scale of skybpm=0.5 AA/pixel
  if usetilt eq 0 then smoothscale = 2*FLOOR(difflambda)+7+4*align $
     else smoothscale = 5+4*align
; IF LOOKING AT SKYLINES, DO NOT REDUCE SKYBPM REGION
  if skyprofile then skybpm[*,1] = $
      CEIL( SMOOTH(float(skybpm[*,1]), smoothscale) ) $
     else skybpm[*,1] = FLOOR( SMOOTH(float(skybpm[*,1]), smoothscale) )

;;;INTERPOLATE TO FIND INDICES IN skybpm THAT 
;;;CORRESPOND TO WAVELENGTHS SLONG SLIT.
  dex = FINDGEN( N_ELEMENTS(skybpm[*,0]) )
 
; KLUDGE!  Fix skybpm for 600-line data !!!!
  if skyprofile then missing=0. else missing=mean(lambdarange) lt 8000.

; two options: if usetilt = 0, we just keep our distance from sky
; lines, with the buffer size depending on the tilt of the line.
; if usetilt = 1, we use the actual tilt in each row to define the
; region used.
  if usetilt eq 0 then begin
     idex = INTERPOL(dex, skybpm[*,0], wvlen)

     ;;;INTERPOLATE TO FIND THE skybpm VALUE (0 OR 1)
     ;;;AT THE GIVEN INDICE (AT THE GIVEN WAVELENGTH).
     colmask = FLOOR( INTERPOLATE(skybpm[*,1], idex, MISSING=missing) )

     ;;;AND THEN MAKE A MASK FOR THE ENTIRE SLIT. 
     ;;;THAT IS, FOR ALL COLUMNS IN THE SLIT.
     skylinemask = colmask # (FLTARR(ncols)+1.)
  endif else begin
     idex = INTERPOL(dex, skybpm[*,0], lambda)

     ;;;INTERPOLATE TO FIND THE skybpm VALUE (0 OR 1)
     ;;;AT EACH PIXEL (I.E. AT THE GIVEN WAVELENGTH).
     skylinemask = FLOOR( INTERPOLATE(skybpm[*,1], idex, MISSING=missing) )
   endelse

; if we were so far off the standard wavelengths the skybpm is
; useless, just do an IVAR-weighted sum
   if total(skylinemask) lt 500. then begin
     skylinemask = skylinemask*0+1
     print, 'WARNING: skylinemask was useless, IVAR-weighted sum used'
   endif

; NOTE: right now, masks are 1 on good pixels, 0 on bad.
; This needs to be reversed at every djs_maskinterp .
   
   slitmask = skylinemask*0+1
   crmask = slitmask*0
;;; IF KEYWORD bpm WAS PASSED, THEN APPLY BAD PIXEL MASK.
  IF (bpm) THEN $
    bpmmask = ((bitmask AND 2b) EQ 0b AND (bitmask AND 1b) EQ 0b) $
      OR (bitmask AND 4b) eq 4b else bpmmask = 1

  pixmapmask = (bitmask AND 16b) EQ 0b

  satmask = pixmapmask*0
  satpix = where(flux GT 6E4, satct)
  if satct gt 0 then begin
    satmask(satpix) = 1
    satmask = satmask OR dilate(satmask, fltarr(5, 5)+1)
;    if quick le 0 then flux = djs_maskinterp(flux, satmask, iaxis = 1)
     if quick le 0 then flux = flux
  endif

;;; IF KEYWORD cr WAS PASSED, THEN FIND INDICES (whcr) OF COSMIC RAYS. 
  IF KEYWORD_SET(cr) THEN BEGIN

      if quick le 0 then $
        medflux = djs_MEDIAN(flux, width=7, boundary='reflect')

  ;  wh = where(ivar EQ 0, whct)
  ;  IF whct GT 0 THEN ivar[wh] = 1.d-20
  ;  IF whct GT 0 then slitmask[wh] = 0.
      usecrmaskin = total(crmaskin) gt 0

      if  usecrmaskin then $
         whcr = WHERE( crmaskin OR satmask, whct) $
      else $
        whcr = WHERE( flux-medflux GT 6./SQRT(ivar > 1.e-5) OR satmask, whct)


;;; SET crmask VALUE TO 0 WHEREVER COSMIC RAY OCCURS.
    IF whct NE 0 THEN begin
      crmask[whcr] = 1
      if usecrmaskin eq 0 then $
        crmask = dilate(crmask, fltarr(5, 5)+1) OR crmask
      if quick le 0 then  flux[whcr] = medflux[whcr]
    ENDIF
  ENDIF
  crmask = 1-crmask

; DEBUG CODE:
  if atv then  atv, (skylinemask)*flux, min=-10., max=20.
  if atv and skyprofile then plothist, flux[where(skylinemask gt 0)], $
    xrange = [0, 100]
  if atv and skyprofile then print, 'mean counts', $
    mean(flux[where(skylinemask gt 0)])

;;; CHECK FOR PLACES WHERE FLUX IS INFINITE. 
  whinf = WHERE(FINITE(flux) EQ 0, whct)
  IF whct NE 0 THEN slitmask[whinf] = 0
  IF whct NE 0 THEN flux[whinf] = 0.

; if making a profile on skylines, want to renormalize by sky brightness (so
; losing a few pixels here and there won't affect us)

  if bspline then begin
     skymodel = bspline_valu(lambda, sset)
     flux = flux/skymodel
     wh = where(finite(flux) eq 0, whct)
     if whct gt 0 then flux[wh] = 0
     if whct gt 0 then slitmask[wh] = 0
     
   endif

;;; DETERMINE THE NUMBER OF PIXELS IN EACH ROW.    

; begin interpolated scenario
   if interp then begin



; interpolate over bad pixels/columns in the spatial direction
     if bpm AND quick le 0 $
;       then flux = djs_maskinterp(flux, 1-bpmmask, iaxis = 1)
       then flux = flux



; interpolate over  non-finite pixels and skylines

;     tmpcr = (1-crmask) OR dilate(1-crmask, fltarr(3)+1)
;     flux = djs_maskinterp(flux, (1-slitmask) OR (1-pixmapmask), iaxis = 0)


; include interpolated pixels in ivar-weighted sum, with low weight
     whok = where(ivar ne 0., okct)
     whbad = where(1b-bpmmask, badct)

     if okct gt 0 and badct gt 0 then begin
       minivar =  min(ivar[whok])
       ivar[whbad] = minivar
     endif

     ivar = ivar*slitmask*pixmapmask
     if quick then ivar = ivar*bpmmask

; now, need to determine the useful wavelength range along each slit.
; We want to use the same number of pixels in each row!!!

; first, find the regime where no pixels in any column are vignetted,
; and we are also not on a skyline in some rows but not others

;     badregions =  ( (nearvig ne 0 AND (maxcol-mincol gt 1000))) OR ((bitmask AND 8b) EQ 8b)
      badregions =  (bitmask AND 8b) EQ 8b OR (nearvig ne 0)
     ivar = ivar*(badregions eq 0)
     minmaxes = fltarr(ncols, 2)

; The following seems to work fine. 

         for i = 0, ncols -1 do minmaxes(i, *) = $
           minmax(where(badregions(*, i) EQ  0 AND skylinemask[*,i]))

         whminok = where(minmaxes(*, 0) NE -1, minct)

;     whminok = where(minmaxes(*, 0) NE -1 AND minmaxes(*,0) lt 1200 $
;	AND minmaxes(*,1) gt 2900, minct)
     if minct gt 0 then minpix = max(minmaxes(whminok, 0)) ELSE minpix = 0
     if minct gt 0 then maxpix = min(minmaxes(whminok, 1)) ELSE maxpix = maxcol-mincol-1

     maxpix = maxpix-10
     minpix = minpix+10
     npix = maxpix-minpix+1

     if npix lt 10 then begin
       minpix = 20 ;1000 > mincol
       maxpix = maxcol-mincol-20 ; 3000 < maxcol
       print, 'warning: entire slit may be vignetted!'
     endif
   
; unweighted
     npix = total(skylinemask[minpix:maxpix, *]* $
                  (ivar[minpix:maxpix, *] NE 0), 1)

;     minivar = min(ivar(where(ivar GT 0)))
;     sprof = total(flux[minpix:maxpix, *]*skylinemask[minpix:maxpix, *], 1)/ $
;             npix
;     profivar = 1/ $
;                total( skylinemask[minpix:maxpix, *]/ $
;                       (ivar[minpix:maxpix, *] > minivar), 1)*npix^2

; ivar weighted

     if skyprofile  then ivarweight = 0

     weight = ivar[minpix:maxpix, *]

     if NOT ivarweight then weight = (badregions[minpix:maxpix, *] eq 0b)

     if skyprofile and bspline then weight = weight* $
       dilate(((skymodel[minpix:maxpix, *]-median(skymodel) gt 40)), $
        intarr(7, 7)+1)

     profivar = total(skylinemask[minpix:maxpix, *]*weight, 1)

     sprof = total(flux[minpix:maxpix, *]*skylinemask[minpix:maxpix, *] $
                   *weight, 1)/ profivar

     if not ivarweight then profivar = $
       total(skylinemask[minpix:maxpix, *]*ivar[minpix:maxpix,*], 1)

   endif else begin

; begin non-interpolated scenario


; note that slitmask is 1 on GOOD pixels
     slitmask = slitmask AND skylinemask AND crmask AND bpmmask $
                AND pixmapmask AND (ivar ne 0)

     npix = TOTAL(slitmask, 1)

;;; NOW MULTIPLY flux BY slitmask AND TOTAL 
;;; TO DETERMINE THE SPATIAL DISTRIBUTION OF 
;;; THE OBJECT FLUX ACROSS THE SLIT. DIVIDE 
;;; BY THE NUMBER OF PIXELS IN EACH ROW SO THAT
;;; THE FLUX VALUES IN sprof ARE AN AVERAGE OF 
;;; THE FLUX PER PIXEL ALONG THE SPATIAL DIRECTION
;;; OF THE OBJECT SPECTRUM. FINALLY, REMOVE THE
;;; FIRST 3 COLUMNS AND FINAL 3 COLUMNS TO REMOVE 
;;; ANY EDGE EFFECTS.
     sprof = TOTAL( (slitmask * flux), 1 ) / FLOAT(npix)

;;; MAKE SURE THAT THE SPATIAL PROFILE IS FINITE IN
;;; ALL PIXEL COLUMNS.
     wh = WHERE(FINITE(sprof) EQ 0, whct)
     IF whct GT 0 THEN sprof(wh) = 0.

;;; DETERMINE THE INVERSE VARIANCE FOR THE SPATIAL PROFILE.
     IF N_ELEMENTS(ivar) NE 0 THEN $
       profivar = FLOAT(npix)^2/(TOTAL(slitmask/(ivar > 1.E-10), 1)) $
       ELSE profivar = FLOAT(npix)/sprof

   ENDELSE

   whinf=where(finite(sprof) eq 0, infct)
	if infct gt 0 then sprof[whinf]=-1.


;;; CONSTRUCT ARRAY CALLED colnum AND SET THE BUFFER
;;; TO BE 5 COLUMNS (5 PIXELS ON EITHER END IN THE 
;;; SPATIAL DIRECTION).  
  colnum = FINDGEN(ncols)
  nbuf = 5
;;; EXLUDE THE FIRST AND LAST 5 COLUMNS IN THE SLIT 
;;; AND ANY COLUMNS FOR WHICH THE SPATIAL PROFILE IS 
;;; INFINITE. CHECK HOW MANY GOOD COLUMNS WE HAVE.
  whuse = WHERE( colnum GE nbuf AND colnum LE ncols-nbuf-1 AND $
                 FINITE(sprof) EQ 1, whct)
;;; IF THERE ARE NO USABLE COLUMNS THEN SET whuse=3. 
;;; This prevents the code from dying completely if there is no useful
;;; region of the slitlet.
  IF N_ELEMENTS(whuse) EQ 0 THEN whuse = 3

  errlabel = 0
;;; DEFINE THE COLUMNS TO USE. 
  use = WHERE(LINDGEN(ncols) GT nbuf AND LINDGEN(ncols) LT (ncols-nbuf-1))

  CATCH, errlabel
  IF errlabel LT 0 THEN BEGIN
    PRINT, 'Error index: ', Errlabel
    PRINT, 'Error message: ', !ERROR_STATE.MSG
    colmask = LINDGEN(ncols) GT nbuf AND LINDGEN(ncols) LT (ncols-nbuf-1)
    
    slitmask2 = (FLTARR(N_ELEMENTS(slitmask)/ncols)+1) # colmask

    WH = WHERE(abs(SLITMASK*slitmask2-1.) LT 0.01, whct)
    IF whct EQ 0 THEN wh = WHERE(abs(slitmask2-1.) LT 0.01)
    f = flux[wh]

     srt=sprof[use]
    srt=srt[sort(srt)]
    minval=srt[3 < (n_elements(srt)-1)]
;    djs_iterstat,f[where(f-minval lt 8.)],median=modeout
    modeout=minval

  ENDIF

  if n_elements(modeout) ne 0 then if finite(modeout) eq 0 then modeout =  0.

  if mmm AND errlabel eq 0 AND align eq 0 $
    then mmm, sprof[use], modeout, modesig

;print, modeout

  if align eq 0 then begin
     if mmm then sprof =  sprof-modeout*(whct ne 0) $
     ELSE if NOT nosubtract then begin
        modeout=MIN(sprof[whuse])
        sprof = sprof-modeout*(whct ne 0)
     endif
  endif else begin
     ddprof = deriv(deriv(sprof))
     nprof = n_elements(ddprof)
     left = min(where(ddprof lt 0 AND shift(ddprof, -1) gt 0))
     right = max(where(ddprof lt 0 AND shift(ddprof, 1) gt 0))
     left2 = left-ddprof[left]/(ddprof[left+1]-ddprof[left])
     if right eq nprof-1 then right2 = right $
       else right2 = right-ddprof[right]/(ddprof[right]-ddprof[right-1 ])

     sprof = sprof-mean(interpolate(sprof, [left2, right2]))

     maxout = sprof[0] > sprof[nprof-1]
     if maxout eq sprof[0] then begin
        align = max(where(sprof gt maxout))/2. 
     endif else begin
        align = (nprof-1+min(where(sprof gt maxout)))/2.
     endelse

;     align = [left2+right2]/2.
  endelse
;;;;;;;;;;;;;;;;;;

;;; RETURN THE SPATIAL PROFILE.
  RETURN, sprof

END





