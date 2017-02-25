FUNCTION NIRI_SKYSUB, filename, skyfiles, flatfile, NOMASK = NOMASK1 $
                      , sciimg = sciimg, ivar = ivar $
                      , objstruct = obj_pos, WAVEIMG = WAVEIMG $
                      , TELLURIC = TELLURIC1, CHK = CHK, WVCHK = WVCHK $
                      , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                      , AVG_SKY = AVG_SKY

IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 50L

IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0

IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
    NOMASK = 1
ENDIF ELSE TELLURIC = 0

IF (N_params() LT 1) THEN BEGIN
    print, 'Syntax: niri_skysub,filenames,skyfiles,flatfile,[OBJFIND=, CHK=,VERBOSE=,]'
;    return, 0
endif

flat = mrdfits(flatfile, 0)
dims = size(flat, /dim)
nx = dims[0]
ny = dims[1]
tset_slits = niri_slitset(nx, ny)

slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)

; First pass average sky subtraction removes dark current and persistence
niri_skyproc, filename, skyfiles, flat, tset_slits $
              , sciimg = sciimg, ivar = ivar $
              , piximg = piximg, waveimg = waveimg, targdir = targdir $
              , TELLURIC = TELLURIC, AVG_SKY = AVG_SKY, hdr = hdr, WVCHK = WVCHK

;   Second-pass global bspline sky subtraction

IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Skip second-pass for standard stars and bright objects.
    print, "    Bright calibration source; skipping 2nd pass..."
    sky_model = avg_sky
    global_sky = 0.0*avg_sky
ENDIF ELSE BEGIN
;      Otherwise do 2nd pass
    print, "  Second-pass: Global sky subtraction..."
    
    y_img = replicate(1.0, nx)#findgen(ny)
    global_sky = fltarr(nx, ny)
    img_minsky = sciimg-avg_sky 
    buffer = 4
    inslit = where(slitmask NE 0 AND $
                   y_img GT 0 AND $ 
                   y_img LT 1024, nord)
    fitpix  = where(slitmask NE 0 AND $
                    finite(img_minsky) EQ 1 AND $
                    finite(ivar) EQ 1  AND $
                    ximg GT 0.05  AND $
                    ximg LT 0.95  AND $
                    abs(img_minsky) LE 5.0d4 AND $
                    ivar GT 0.0, npix)
    psort = sort(piximg[fitpix])
    sset = bspline_iterfit(piximg[fitpix[psort]], img_minsky[fitpix[psort]] $
                           , invvar = (ivar[fitpix[psort]] GT 0.0) $
                           , upper = 3, lower = 3 $
                           , bkspace = 1.1D, maxiter = 20, maxrej = 10)
    global_sky[inslit] = bspline_valu(piximg[inslit], sset)
    IF KEYWORD_SET(CHK) THEN BEGIN
        plotx = piximg[fitpix[psort]]
        ploty = img_minsky[fitpix[psort]]
        rms = sqrt(djs_median(ploty^2))
        x_splot, plotx, ploty, psym1 = 3, ymnx = [-30.0*rms, 30.0*rms] $
                 , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;, xrange = [min(plotx), max(plotx)] $

        wait, 1.5
    ENDIF
ENDELSE
; now break the image up into two quadrants and subtract out any systematic
; horizontal features for each chip
img_minsky = sciimg - avg_sky - global_sky
smash_left = djs_avsigclip(img_minsky[0:(nx/2-1L), *], 1, sigrej = 2.0)
corr_left = smash_left ## replicate(1.0, nx/2)
smash_right = djs_avsigclip(img_minsky[nx/2:*, *], 1, sigrej = 2.0)
corr_right = smash_right ## replicate(1.0, nx/2)
corr = fltarr(nx, ny)
corr[0:(nx/2-1), *] = corr_left
corr[nx/2:*, *] = corr_right
avg_sky = avg_sky + corr
img_minsky = sciimg - avg_sky - global_sky
;img_minsky = img_minsky - corr

print, "Finding objects in sky-subtracted image"
; NIRI plate scale is 0.1171 at f/6. Assuming 0.6" seeing FWHM = 5.12 pix
FWHM = 5.12
obj_pos = long_objfind(img_minsky, tset_slits = tset_slits $
                       , FWHM = FWHM, OBJMASK = OBJMASK_POS $
                       , NPERSLIT = KEYWORD_SET(TELLURIC))
IF NOT KEYWORD_SET(OBJ_POS) THEN BEGIN
    sky_model = avg_sky + global_sky 
    return,  sky_model
ENDIF
fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
; Search for negative objects that might be left over in the avg sky
obj_neg = long_objfind(-img_minsky, tset_slits = tset_slits $
                       , FWHM = fwhm, OBJMASK = OBJMASK_NEG)
IF KEYWORD_SET(obj_neg) THEN objtot = [obj_pos, obj_neg] $
ELSE objtot = obj_pos
OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)

nobj_tot = n_elements(objtot)
nobj = n_elements(obj_pos)
IF NOT KEYWORD_SET(NOMASK) THEN BEGIN
   img_minsky = sciimg-avg_sky 
    local_sky = global_sky
    print, "Doing third pass sky subtraction with object masking"
    FOR iobj = 0L, nobj-1L DO BEGIN
        skymask = lonarr(nx, ny)
        left  = obj_pos[iobj].xpos - SKYBUFFER 
        right = obj_pos[iobj].xpos + SKYBUFFER 
        FOR j = 0L, ny-1L DO BEGIN 
            xmin = floor(left[j]) >  0
            xmax = ceil(right[j]) < nx
            skymask[xmin:xmax, j] = 1L
            skymask[xmin:xmax, j] = 1L
        ENDFOR
        inobj = where(slitmask NE 0 AND   $
                      skymask  EQ 1 AND   $
                      y_img GT 0 AND      $ 
                      y_img LT 1024, nord)
        fitpix = where(slitmask NE 0 AND  $
                       skymask  EQ 1 AND  $
                       objmask  EQ 0 AND  $
                       finite(img_minsky) EQ 1 AND $
                       finite(ivar) EQ 1  AND $
                       ximg GT 0.05 AND $
                       ximg LT 0.95 AND $
                       abs(img_minsky) LE 5.0d4 AND $
                       ivar GT 0.0, npix)
        psort = sort(piximg[fitpix])
        ;; no weights on the sky-subtraction for now. 
        sset = bspline_iterfit(piximg[fitpix[psort]] $
                               , img_minsky[fitpix[psort]] $
                               , invvar = (ivar[fitpix[psort]] GT 0.0) $
                               , upper = 3, lower = 3 $
                               , bkspace = 1.1D, maxiter = 20, maxrej = 10)
        local_sky[inobj] = bspline_valu(piximg[inobj], sset)
        IF KEYWORD_SET(CHK) THEN BEGIN
            plotx = piximg[fitpix[psort]]
            ploty = img_minsky[fitpix[psort]]
            rms = sqrt(djs_median(ploty^2))
            x_splot, plotx, ploty, psym1 = 3 $
                     , ymnx = [-10.0*rms, 10.0*rms] $
                     , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;                     , xrange = [min(plotx), max(plotx)] $

        ENDIF
    ENDFOR
    sky_model = avg_sky + local_sky
ENDIF ELSE BEGIN
    sky_model = avg_sky + global_sky
ENDELSE

IF keyword_set(CHK) THEN BEGIN
    xatv, (sciimg-sky_model)*sqrt(ivar), wv = waveimg, /block
ENDIF

RETURN, sky_model
END
