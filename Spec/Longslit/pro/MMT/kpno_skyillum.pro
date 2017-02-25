;+
; NAME:
;   knpo_skyillum
;
; PURPOSE:
;   Algorithm to generate an illumination flat, ideally from a
;   twilight flat.  This routine is tuned to observations taken at
;   KPNO.
;
; CALLING SEQUENCE:
; kpno_skyillum, filenames, illumflatfile, pixflatfile $
;
; INPUTS:
;   filenames  --  List of sky flat frames
;   pixflatfile -- Filename of the pixel flat.  Presumably created
;                  from internal or dome exposures.
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;  illumflatfile -- Filename for the illumination flat
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and 
;                Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
PRO kpno_skyillum, filenames, illumflatfile, pixflatfile $
                   , biasfile = biasfile $
                   , nsmooth = nsmooth, NSAMPLE = NSAMPLE $
                   , WAVEFILE = WAVEFILE, VERBOSE = VERBOSE $
                   , INDIR = INDIR, TEMPDIR = TEMPDIR, WRITE_FLATS = WRITE_FLATS

nfiles = n_elements(filenames)

if (keyword_set(tempdir)) then $
    spawn, '\mkdir -p '+tempdir

IF NOT keyword_set(indir) then indir = './'
IF NOT KEYWORD_SET(TEMPDIR) THEN tempdir = './'
IF NOT KEYWORD_SET(pixflatfile) THEN $
  message, 'ERROR: pixel flat file must be included for illumination flats'

tempfile = djs_filepath('tempillum' + filenames[0], root_dir = tempdir)

IF (NOT keyword_set(sigrej)) then begin
    if (nfiles LE 2) then sigrej = 1.0 $ 
; Irrelevant for only 1 or 2 files
    else if (nfiles EQ 3) then sigrej = 1.1 $
    else if (nfiles EQ 4) then sigrej = 1.3 $
    else if (nfiles EQ 5) then sigrej = 1.6 $
    else if (nfiles EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0
ENDIF
if NOT keyword_set(xsamp) then xsamp = 3.0
IF NOT keyword_set(nsmooth) then nsmooth = 11

upper = 3.0
lower = 3.0
t0 = systime(1)
i = 0

long_proc, filenames[i], flat1, ivar1, hdr = hdr, biasfile = biasfile $
           , pixflatfile = pixflatfile, VERBOSE = VERBOSE
rnoise = double(sxpar(hdr, 'RDNOISE'))
dims = size(flat1, /dimen) 
nx = dims[0]
ny = dims[1]
flatfull = rebin(flat1[*], nx*ny, nfiles)
;   Read in wavelength solution structures
pixset  = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 1)
wset    = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 2)
tset_slits = mmt_slitset(nx, ny)
piximg = long_wpix2image(pixset, tset_slits, wset = wset, waveimg = waveimg)
wave_max = max(waveimg[WHERE(waveimg GT 0.0)])
wave_min = min(waveimg[WHERE(waveimg GT 0.0)])

; Read in all of the object exposures
FOR i = 0, nfiles-1 do begin 
    splog, 'Working on file ', filenames[i]
    IF i GT 0 then BEGIN
        long_proc, filenames[i], flat1, ivar1, biasfile = biasfile $
                   , pixflatfile = pixflatfile, VERBOSE = VERBOSE
    ENDIF
    flatfull[*, i] = flat1
ENDFOR
; Create a super sky flat from all the object exposures
if nfiles GT 1 then image = reform(djs_median(flatfull, 2), nx, ny) $
  else image = reform(flatfull, nx, ny)
IF KEYWORD_SET(WRITE_FLATS) THEN mwrfits, image, tempfile, /create

; This will serve as an approximate inverse variance image for the fits
ivar = 1.0/(abs(image - sqrt(2.0)*rnoise) + rnoise^2)

;smash = djs_median(image[*, 20:ncol-20], 2)
xmin = min(piximg)
xmax = max(piximg)
xrange = xmax - xmin
nxbkpt = long(xrange/xsamp) + 1
;       take the good pixels
good = where(ivar GT 0 AND waveimg GT 0.0, ngood, COMPLEMENT = bad $
             , NCOMPLEMENT = nbad)
mask = lonarr(nx, ny)+1L
IF nbad GT 0 THEN mask[bad] = 0L

pix_good = piximg[good]
log_image = alog(image[good] > 1)
log_ivar = 1.0*(image[good] GT 1 AND ivar[good] GT 0)
xs = sort(piximg[good])
splog, 'Spectral fit of flatfield', ngood, ' Pixels to fit'
spec_set = bspline_longslit(pix_good[xs], log_image[xs], log_ivar[xs] $
                            , xs*0.+1, everyn = ngood/nxbkpt $
                            , nord = 2, upper = 0.2, lower = 0.2 $
                            , maxrej = 5, /groupbadpix $
                            , outmask = outmask1, /silent) 
image_fit_log = bspline_valu(piximg, spec_set)
image_fit = exp(image_fit_log)
outmask = lonarr(nx, ny)
outmask[good[xs]] = outmask1
outmask = outmask*mask
qgood = image_fit GT 0
;smash_max = max(smash)
norm1 = mask*image/(image_fit + (qgood EQ 0))
norm_ivar1 = outmask*ivar*(image_fit^2)
norm = 0.0*norm1
norm_ivar = 0.0*norm_ivar1
;Now do a rolling median of the image in the spatial direction to remove 
;and remaining object traces in the residual image
FOR j = 0, ny-1L DO BEGIN
    norm[*, j] = djs_median(norm1[*, j], width = nsmooth)
    norm_ivar[*, j] = djs_median(norm_ivar1[*, j], width = nsmooth)
ENDFOR

; now polynomial fit each row in the spectral direction to removed
; the bumps and wiggles from the sky line residuals that didn't normalize out
xcol = findgen(ny)#replicate(1.0, nx)
xy2traceset, xcol, transpose(norm), tset, ncoeff = 5, yfit = normfit_t $
             , invvar = transpose(norm_ivar), /silent
; convolve the fit with a spatial boxcar to smooth things out over the 
; scale of a trace
normfit = transpose(normfit_t)
kernel = replicate(1.0, nsmooth)/float(nsmooth)
flatfinal = mask*convol(normfit, kernel, /EDGE_TRUNCATE)

mwrfits, flatfinal, illumflatfile, /create
splog, 'Elapsed time = ', systime(1)-t0, ' sec'

return
end


