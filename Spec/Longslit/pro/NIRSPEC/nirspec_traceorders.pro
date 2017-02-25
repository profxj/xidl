;+
; NAME:
;   gnirs_traceorders
;
; PURPOSE:
;   Determine the positions of the slits on the image, and return tracesets
;
; CALLING SEQUENCE:
;   gnirs_traceorders
;
; INPUTS:
;   filename   - Image for finding the slits, which would typically
;                be a flat-field image, an arc image, or a sum of those
;   outfile    - Output file name with slit mask positions
;
; OPTIONAL INPUTS:
;   y1         - Starting row number for smashing image to initially
;                identify slits; default to 0.40*NY
;   y2         - Ending row number for smashing image to initially
;                identify slits; default to 0.60*NY
;   nmed       - Width for median-filtering the image first in the wavelength
;                direction (to remove cosmics)
;   ksize      - Half-kernel size for sharpness filter; default to 5 pix
;   peakthresh - Flux threshhold for finding slits; the flux must be
;                at least this fraction of the brightest slit; default to 0.02
;   radius     - Keyword for TRACE_CRUDE; default to same value as KSIZE
;   nave       - Keyword for TRACE_CRUDE; default to 3
;   maxshifte  - Keyword for TRACE_CRUDE; default to 0.1
;   maxshift0  - Keyword for TRACE_CRUDE; default to 1.0
;   func       - Keyword for XY2TRACESET; default to 'legendre'
;   ncoeff     - Keyword for XY2TRACESET; default to 3
;   verbose    - Verbose if set
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
PRO nirspec_traceorders, flatfile, darkfile,outfile $
                              , y1 = y1, y2 = y2, nmed = nmed, ksize = ksize $
                              , peakthresh = peakthresh, radius = radius $
                              , nave = nave, maxshifte = maxshifte $
                              , maxshift0 = maxshift0, func = func $
                              , ncoeff = ncoeff $
                              , verbose = verbose $
                              , CRUDE = CRUDE

nirspec_proc,flatfile,flatimg
nirspec_proc,darkfile,darkimg
image = flatimg-darkimg
;crude=1  
t0 = systime(1)

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]

;----------
; Set defaults

IF NOT KEYWORD_SET(CRUDE) THEN archive_file= $
   getenv('XIDL_DIR') + '/Spec/Longslit/calib/flats/NIRSPEC/nirspec_archive_Korder_2009.fits'

if (n_elements(y1) EQ 0) then y1 = 450.0
if (n_elements(y2) EQ 0) then y2 = 460.0
if (y1 GT y2 OR y1 LT 0 OR y2 GT ny-1) then $
  message, 'Invalid values for Y1,Y2'
if (NOT keyword_set(ksize)) then ksize = 3
if (ksize LT 1 OR ksize GE nx/2-1) then $
  message, 'Invalid kernel size KSIZE'
if (NOT keyword_set(peakthresh)) then peakthresh = 0.01
if (NOT keyword_set(nmed)) then nmed = 5

;  Default values for TRACE_CRUDE
if (NOT keyword_set(radius)) then radius = ksize
if (NOT keyword_set(nave)) then nave = 5
if (NOT keyword_set(maxshifte)) then maxshifte = 5.0
if (NOT keyword_set(maxshift0)) then maxshift0 = 5.0
if (NOT keyword_set(maxerr)) then maxerr = 1.0

;  Default values for parametrizing the traces (XY2TRACESET)
if (NOT keyword_set(func)) then func = 'legendre'
if (NOT keyword_set(ncoeff)) then ncoeff = 3

;   Convolve the input image with the sharpness filter
kern = [-(findgen(ksize)+1), reverse(findgen(ksize)+1)] / (2.*ksize)
cimg = convol(image, kern, /center, /edge_truncate)

norders = 1
xstart = dblarr(norders)
xend   = dblarr(norders)
;   Smash the center of the sharpness-filtered image to a 1-dimensional vector
fsum = djs_median(cimg[*, y1:y2], 2)

xleft_left   = [490]
xright_left  = [760]

xleft_right  = [640]
xright_right = [870]

xstart = dblarr(norders)
xend   = dblarr(norders)
FOR iorder = 0L, norders-1L DO BEGIN
   ;; Find candidate starting positions
   fsum_iord = djs_median(fsum[xleft_left[iorder]:xright_left[iorder]] $
                          , width = 4, boundary = 'reflect')
   xstart1 = long_find_nminima(-fsum_iord, nfind = 1, minsep = 40L, width = 3)
   ystart = interpolate(fsum_iord, double(xstart1))
   xstart[iorder] = double(xstart1) + double(xleft_left[iorder])
   fsum_iord = djs_median(fsum[xleft_right[iorder]:xright_right[iorder]] $
                          , width = 4, boundary = 'reflect')
   xend1  = long_find_nminima(fsum_iord, nfind = 1, minsep = 40L, width = 3)
   xend[iorder]   = double(xend1) + double(xleft_right[iorder])
ENDFOR


;; Make sure that xstart is always before xend
FOR iorder = 0L, norders-1L DO BEGIN
   IF xend[iorder] LT xstart[iorder] THEN BEGIN
      xs1 = xstart[iorder]
      xe1 = xend[iorder]
      xend[iorder] = xs1
      xstart[iorder] = xe1
   ENDIF
ENDFOR

ycen = (y1+y2)/2.0
minval = -djsig(cimg)     ; Clip negative values on image when tracing

IF KEYWORD_SET(ARCHIVE_FILE) THEN BEGIN 
;; Use archived slit traces as crutch for tracing
    tset_guess = mrdfits(archive_file, 1) 
    traceset2xy, tset_guess[0], rows, left_edge
    traceset2xy, tset_guess[1], rows, right_edge
    left_crutch = 0*left_edge
    right_crutch = 0*right_edge
    FOR iorder = 0L, norders-1L DO BEGIN
        left_crutch[*, iorder] = left_edge[*, iorder] - $
          (interpol(left_edge[*, iorder], rows[*, iorder], ycen))  + xstart[iorder]
        right_crutch[*, iorder] = right_edge[*, iorder] - $
          (interpol(right_edge[*, iorder], rows[*, iorder], ycen))  + xend[iorder]
    ENDFOR
 ENDIF ELSE BEGIN
    ;; use trace crude to determine the crutch for slit tracing
    radius = 3
;   Trace the starting+end positions on the CCD.
    xpos1 = trace_crude(cimg > minval, xstart = xstart $
                        , ystart = replicate(ycen, n_elements(xstart)) $
                        , radius = radius, maxerr = maxerr  $
                        , nave = nave, maxshifte = maxshifte $
                        , maxshift0 = maxshift0, xerr = xerr1)
    xpos2 = trace_crude(-cimg > minval, xstart = xend $
                        , ystart = replicate(ycen, n_elements(xstart)) $
                        , radius = radius, maxerr = maxerr $
                        , nave = nave, maxshifte = maxshifte $
                        , maxshift0 = maxshift0, xerr = xerr2)
    ;;  Fit these traces by a polynomial function
    xy2traceset, rebin(findgen(ny), ny, norders), xpos1, sset1 $
                 , invvar = 1./xerr1^2, func = func, ncoeff = ncoeff
    xy2traceset, rebin(findgen(ny), ny, norders), xpos2, sset2 $
                 , invvar = 1./xerr2^2, func = func, ncoeff = ncoeff
    traceset2xy, sset1, rows, left1
    traceset2xy, sset2, rows, right1
    xerr1 = xpos1-left1
    xerr2 = xpos2-right1
    ivar1 = (abs(xerr1) GT 0.0 AND abs(xerr1) LE 5.0)/ $
            (xerr1^2 + (xerr1^2 LE 0.0)) <  1.0d4
    ivar2 = (abs(xerr2) GT 0.0 AND abs(xerr2) LE 5.0)/ $
            (xerr2^2 + (xerr2^2 LE 0.0)) <  1.0d4
    
    ;; refit with these new errors based on the fits
    xy2traceset, rebin(findgen(ny), ny, norders), xpos1, sset1 $
                 , invvar = ivar1, func = func, ncoeff = ncoeff
    xy2traceset, rebin(findgen(ny), ny, norders), xpos2, sset2 $
                 , invvar = ivar2, func = func, ncoeff = ncoeff
    tset_guess = [sset1, sset2]
    tset_guess = struct_addtags(tset_guess $
                                , replicate(create_struct('DIMS', dims) $
                                            , size(tset_guess, /dimens)))
    traceset2xy, tset_guess[0], rows, left_crutch
    traceset2xy, tset_guess[1], rows, right_crutch
 ENDELSE

niter = 12L
xfit1 = left_crutch
FOR i = 1L, niter DO BEGIN
    IF i LT niter/3 THEN rad = 3.0*ksize $
    ELSE IF (i GE niter/3) AND (i LT 2*niter/3) THEN rad = ksize $
    ELSE rad = ksize/2.0d
    xpos1 = trace_fweight(cimg > 0.0, xfit1, rows, radius = rad)
    xy2traceset, rows, xpos1, left_set1, ncoeff = 5, yfit = xfit1 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit2 = xfit1
FOR i = 1L, niter DO BEGIN
    xpos2 = trace_gweight(cimg > 0.0, xfit2, rows $
                          , sigma = KSIZE/3.0D/2.3548D)
    xy2traceset, rows, xpos2, left_set, ncoeff = 5, yfit = xfit2 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit1 = right_crutch
FOR i = 1L, niter DO BEGIN
    IF i LT niter/3 THEN rad = 3.0*ksize $
    ELSE IF (i GE niter/3) AND (i LT 2*niter/3) THEN rad = ksize $
    ELSE rad = ksize/2.0d
    xpos1 = trace_fweight(-cimg > 0.0, xfit1, rows, radius = rad)
    xy2traceset, rows, xpos1, right_set1, ncoeff = 5, yfit = xfit1 $
                 , maxdev = 1.0D, /silent
ENDFOR
xfit2 = xfit1
FOR i = 1L, niter DO BEGIN
    xpos2 = trace_gweight(-cimg > 0.0, xfit2, rows $
                          , sigma = KSIZE/3.0D/2.3548D)
    xy2traceset, rows, xpos2, right_set, ncoeff = 5, yfit = xfit2 $
                 , maxdev = 1.0D, /silent
ENDFOR
tset_slits = [left_set, right_set]
tset_slits = struct_addtags(tset_slits $
                            , replicate(create_struct('DIMS', dims) $
                                        , size(tset_slits, /dimens)))
;  Generate a mask of the order numbers
slitmask = long_slits2mask(tset_slits)

splog, 'Number of orders = ', norders
splog, 'Elapsed time = ', systime(1)-t0, ' sec'

mwrfits,slitmask,outfile,/create
mwrfits,tset_slits,outfile

end
;------------------------------------------------------------------------------
