PRO gmos_illumflat, filenames, illumflatfile, pixflatfile $
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

tempfiles = djs_filepath('tempillum-' + fileandpath(filenames) $
                         , root_dir = tempdir)

IF (NOT keyword_set(sigrej)) then begin
    if (nfiles LE 2) then sigrej = 1.0 $ 
; Irrelevant for only 1 or 2 files
    else if (nfiles EQ 3) then sigrej = 1.1 $
    else if (nfiles EQ 4) then sigrej = 1.3 $
    else if (nfiles EQ 5) then sigrej = 1.6 $
    else if (nfiles EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0
ENDIF

upper = 3.0
lower = 3.0
t0 = systime(1)
i = 0


lris_proc, filenames[i], flat_full, invvar_full, biasfile = biasfile $
           , pixflatfile = pixflatfile, VERBOSE = VERBOSE
; unmosaic files since we will fit each CCD separately
flat1 = gunmosaic(flat_full)
invvar1 = gunmosaic(invvar_full)
dims = size(flat1, /dimen) 
ncol = dims[0]
nrow = dims[1]
nccd = dims[2]
flatfull = rebin(flat1[*], ncol*nrow*nccd, nfiles)

waveimg_full = mrdfits(wavefile,  silent = (keyword_set(verbose) EQ 0), 0)
dims_full = size(waveimg_full)
nwave = dims_full[2]
wave_max = max(waveimg_full[WHERE(waveimg_full GT 0.0)])
wave_min = min(waveimg_full[WHERE(waveimg_full GT 0.0)])
dispersion =  (wave_max-wave_min)/double(nwave)
piximg_full = waveimg_full/dispersion
; unmosaic waveimg and piximg
waveimg_arr = gunmosaic(waveimg_full)
piximg_arr = gunmosaic(piximg_full)

if NOT keyword_set(xsamp) then xsamp = 3.0


IF NOT keyword_set(nsmooth) then nsmooth = 11
;IF NOT KEYWORD_SET(NSAMPLE) THEN NSAMPLE = 50

FOR i = 0, nfiles-1 do begin 
    splog, 'Working on file ', filenames[i]
    IF i GT 0 then BEGIN
        lris_proc, filenames[i], flat_full, invvar_full, biasfile = biasfile $
                   , pixflatfile = pixflatfile, VERBOSE = VERBOSE
        flat1 = gunmosaic(flat_full)
        invvar1 = gunmosaic(invvar_full)
    ENDIF
;   loop over the three CCDs    
    FOR j = 0, 2 DO BEGIN
        image = flat1[*, *, j]
        invvar = invvar1[*, *, j]
        waveimg = waveimg_arr[*, *, j]
        piximg  = piximg_arr[*, *, j]
        smash = djs_median(image[20:ncol-20, *], 1)
        smash_image = replicate(1, ncol) # smash
        smash_mask = convol(1.0*(smash_image LT 0.5*median(smash)), $
                            transpose(replicate(1, 5)), /edge_trun) EQ 0
;       take the good pixels
        good = where((smash_mask GT 0) AND (invvar GT 0) AND waveimg GT 0.0 $
                     , ngood)
        xmin = min(piximg)
        xmax = max(piximg)
        xrange = xmax - xmin
        nxbkpt = long(xrange/xsamp) + 1
;       under sample the sky pixels to speed up the b-spline fit
;        IF KEYWORD_SET(NSAMPLE) THEN BEGIN
;            pick = [0, (lindgen(nsample*ngood/ncol)+0.5)*ncol/nsample, ngood-1]
;            everyn = nsample/2L
;        ENDIF ELSE BEGIN
;            pick = lindgen(ngood)
;            everyn = ngood/2L
;        ENDELSE
        pix_good = piximg[good]
        log_image = alog(image[good] > 1)
        log_ivar = 1.0 * (image[good] GT 1 AND invvar[good] GT 0)
        xs = sort(piximg[good])
        splog, 'Spectral fit of flatfield', ngood, ' Pixels to fit'
        spec_set = bspline_longslit(pix_good[xs], log_image[xs], log_ivar[xs] $
                                     , xs*0.+1, everyn = ngood/nxbkpt $
                                     , nord = 2, upper = 0.2, lower = 0.2 $
                                     , maxrej = 5, /groupbadpix $
                                     , outmask = outmask1, /silent) 
        image_fit_log = bspline_valu(piximg, spec_set)
        image_fit = exp(image_fit_log)
        outmask = piximg*0.
        outmask[good[xs]] = outmask1
        outmask = outmask*smash_mask
        qgood = image_fit GT 0
        smash_max = max(smash)
        norm = smash_mask*image/(image_fit + (qgood EQ 0))
        norm_ivar = outmask*invvar*(image_fit^2)
;        ratio = im[good]/smash_image[good]
;        iv_ratio = iv[good]*(smash_image[good])^2
;       bspline fit the sky pixels 
        ;bset = bspline_iterfit(waveimg[good[ss[pick]]], ratio[ss[pick]] $
        ;                       , invvar = iv_ratio[ss[pick]], maxdev = 0.1 $
        ;                       , everyn = everyn, maxiter = 2)
;       evaluate the b-spline for the whole image
;        image_fit = bspline_valu(waveimg, bset)
;       normalize the illumination function
;        norm = im/image_fit/max(smash)
;       now polynomial fit each row in the spectral direction to removed
;       the bumps and wiggles from the sky line residuals that didn't normalize
;       out
        xcol = findgen(ncol)#replicate(1.0, nrow)
        xy2traceset, xcol, norm, tset, ncoeff = 5, yfit = normfit $
                     , invvar = norm_ivar, /silent
;       convolve the fit with a spatial boxcar to smooth things out over the 
;       scale of a trace
        kernel = transpose(replicate(1.0, nsmooth)/float(nsmooth))
        flat1[*, *, j] = smash_mask*convol(normfit, kernel)
    ENDFOR
    IF KEYWORD_SET(WRITE_FLATS) THEN BEGIN
        temp_array = gmosaic(flat1)
        mwrfits, temp_array, tempfiles[i], /create
    ENDIF
    flatfull[*, i] = flat1
ENDFOR

flatfinal = reform(djs_avsigclip(flatfull, 2, sigrej = sigrej) $
                   , ncol, nrow, nccd)

flat_mosaic = gmosaic(flatfinal)
mwrfits, flat_mosaic, illumflatfile, /create


splog, 'Elapsed time = ', systime(1)-t0, ' sec'

return
end


