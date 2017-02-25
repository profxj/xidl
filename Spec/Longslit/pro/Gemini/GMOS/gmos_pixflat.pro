PRO gmos_pixflat, filenames, superflatfile, SPECFIT = SPECFIT $
                , biasfile = biasfile, NSMOOTH = NSMOOTH, VERBOSE = VERBOSE $
                , INDIR = INDIR, TEMPDIR = TEMPDIR, WRITE_FLATS = WRITE_FLATS

nfiles = n_elements(filenames)

if (keyword_set(tempdir)) then $
    spawn, '\mkdir -p '+tempdir

IF NOT keyword_set(indir) then indir = './'
IF NOT KEYWORD_SET(TEMPDIR) THEN tempdir = './'

tempfiles = djs_filepath('temppixel-' + fileandpath(filenames) $
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

i = 0

lris_proc, filenames[i], flat_full, invvar_full $
           , biasfile = biasfile, VERBOSE = VERBOSE
; unmosaic files
flat1 = gunmosaic(flat_full)
invvar1 = gunmosaic(invvar_full)
dims = size(flat1, /dimen) 
ncol = dims[0]
nrow = dims[1]
nccd = dims[2]
flatfull = rebin(flat1[*], ncol*nrow*nccd, nfiles)

FORMAT   = '(%"Masked %d pixels for CCD # %d in flatfile # %d")'
FOR i = 0, nfiles-1 do begin 
    splog, 'Working on file ', filenames[i]
    IF i GT 0 then BEGIN
        lris_proc, filenames[i], flat_full, invvar_full, biasfile = biasfile $
                   , VERBOSE = VERBOSE
        flat1 = gunmosaic(flat_full)
        invvar1 = gunmosaic(invvar_full)
    ENDIF
;   loop over the three CCDs    
    FOR j = 0, 2 DO BEGIN
;       construct a mask in the slit 
        slit_med = djs_median(flat1[*, *, j], 1) > 0
        poly = poly_fit(findgen(nrow), slit_med, 2, yfit = slit_poly)
        slitmask = smooth(1.0*(slit_med LT 0.5*slit_poly), 3) EQ 0
;       divide out the median spatial profile
        flat1[*, *, j] = $
          flat1[*, *, j]/ $
          (replicate(1, ncol) # (slit_med + (slit_med EQ 0)))
;        invvar1[*, *, j] = $
;          invvar1[*, *, j]*fitmask*(replicate(1, ncol) # slit_med)
;       now do a polynomial fit column by column to the spatial residual 
;       construct a mask in the spec direction
        spec_med = djs_median(flat1[*, *, j], 2) > 0
        specmask = smooth(1.0*(spec_med LT 0.1), 3) EQ 0
        n_bad = total(specmask EQ 0)
        IF KEYWORD_SET(VERBOSE) THEN splog, format = format, n_bad, j, i
;       multiply to get a total mask
        mask =  specmask # slitmask
        specimage = spec_med # replicate(1.0, nrow)
        fitmask = (flat1[*, *, j] GT 0.8*specimage)*mask
        invvar1[*, *, j] = $
          invvar1[*, *, j]*fitmask*(replicate(1, ncol) # slit_med^2)
;       now do a polynomial fit column by column to the spatial residual 
        xrow = findgen(nrow)#replicate(1.0, ncol)
        xy2traceset, xrow, transpose(flat1[*, *, j]) $
                     , invvar = transpose(invvar1[*, *, j]) $
                     , tset, ncoeff = 5 $
                     , yfit = slitfit, upper = 3, lower = 3, /silent
;       smooth the polynomial coefficients
        smooth_set = tset
        FOR k = 0, n_elements(tset.COEFF[*, 0])-1L DO $
          smooth_set.COEFF[k, *] = median(transpose(tset.coeff[k, *]), 5)
        traceset2xy, smooth_set, xrow, slitsmooth
        slitsmooth = transpose(slitsmooth)
;       remove the spatial trend
        flat1[*, *, j] = $
          flat1[*, *, j] * mask/(slitsmooth + (slitsmooth EQ 0))
        test_dev = total(abs(flat1[*, *, j]-1.0) GT 0.1, 1)
        finalmask = (test_dev LT ncol/2) ## replicate(1.0, ncol)
        flat1[*, *, j] = flat1[*, *, j] * finalmask
    endfor
    IF KEYWORD_SET(WRITE_FLATS) THEN BEGIN
        temp_array = gmosaic(flat1)
        mwrfits, temp_array, tempfiles[i], /create
    ENDIF
    flatfull[*, i] = flat1 
endfor

flatfinal = reform(djs_avsigclip(flatfull, 2, sigrej = sigrej) $
                   , ncol, nrow, nccd)

flat_mosaic = gmosaic(flatfinal)
mwrfits, flat_mosaic, superflatfile, /create

return

END
