PRO NIRI_SKYPROC, filename, skyfiles1, flat, tset_slits $
                  , img_minsky = img_minsky, ivar = ivar $
                  , piximg = piximg, waveimg = waveimg, targdir = targdir $
                  , TELLURIC = TELLURIC, hdr = hdr1, sigrej = sigrej1 $
                  , AVG_SKY = AVG_SKY, sciimg = sciimg, SIGLOWER = SIGLOWER $
                  , WVCHK = WVCHK

siglower = 1.0
IF KEYWORD_SET(TELLURIC) THEN skyfiles = telluric $
ELSE skyfiles = skyfiles1

if (n_elements(adderr) EQ 0) then adderr = 0.01

dims = size(flat, /dimen)
nx = dims[0]
ny = dims[1]
nfiles = n_elements(skyfiles)

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)

IF (NOT keyword_set(sigrej1)) then begin
    if (nfiles LE 2) then sigrej = 1.0 $ 
; Irrelevant for only 1 or 2 files
    else if (nfiles EQ 3) then sigrej = 1.1 $
    else if (nfiles EQ 4) then sigrej = 1.3 $
    else if (nfiles EQ 5) then sigrej = 1.6 $
    else if (nfiles EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0
ENDIF ELSE sigrej = sigrej1
IF NOT KEYWORD_SET(SIGLOWER) THEN SIGLOWER = 1.0D
sigrej = siglower*sigrej

skystack  = fltarr(nx*ny, nfiles)
ivarstack = fltarr(nx*ny, nfiles)
maskstack = lonarr(nx*ny, nfiles)

niri_proc, filename, sciimg, ivar, hdr = hdr1, flatimg = flat $
           , TELLURIC = KEYWORD_SET(TELLURIC)
mask = (ivar LE 0.0)            ; opposite convention for avsigclip

; Read in sky files (this is wavelength files for TELLURIC)
FOR j = 0L, nfiles-1L DO BEGIN
    niri_proc, skyfiles[j], imag1, ivar1, hdr = hdr, flatimg = flat
    skystack[*, j] = imag1
    ivarstack[*, j] = ivar1
    maskstack[*, j] = (ivar1 LE 0.0) ; opposite convention for avsigclip
ENDFOR

; Determine slit width from header
slit_hdr   =  strtrim(sxpar(hdr1, 'FPMASK'))
slit_str = strsplit(slit_hdr, '-.pix', /extr)
slit = double(slit_str[1])
pkwdth = slit
TOLER = slit/2.0D

FWHM = slit         
sigma_psf = FWHM/2.35482D
pkwdth = slit
TOLER = slit/2.0D

;   Take the PSF width to be that of the spectra direction (3 pixels) for 
;   the 3-pixel slit. This prevents the routine from rejecting sky lines
;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
;    
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
;    PSFVALS[0]          1.   PSFVALS[0]
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
ncr = 0
IF nfiles GT 1 THEN $
  avg_sky1 = djs_avsigclip(skystack, 2, sigrej = sigrej, inmask = maskstack) $
ELSE avg_sky1 = skystack
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
reject_cr, sciimg-avg_sky1, ivar, psfvals, rejects, nrejects = ncr $
           , ignoremask = mask
IF ncr NE 0 THEN mask[rejects] = 1L
;   Now overscan subtract masking possible CRs in the oscan
niri_oscan, sciimg, invvar = ivar, inmask = mask
FOR j = 0L, nfiles-1L DO BEGIN
    ncr = 0
    reject_cr, reform(skystack[*, j]-avg_sky1, nx, ny) $
               , reform(ivarstack[*, j], nx, ny), psfvals $
               , rejects, nrejects = ncr $
               , ignoremask = reform(maskstack[*, j], nx, ny)
    IF ncr NE 0 THEN maskstack[rejects, j] = 1L
    imag1 = reform(skystack[*, j], nx, ny)
    ivar1 = reform(ivarstack[*, j], nx, ny)
    mask1 = reform(maskstack[*, j], nx, ny)
    niri_oscan, imag1, inmask = mask1
    skystack[*, j] = imag1
ENDFOR
; Preliminary stacked sky image
IF nfiles GT 1 THEN BEGIN
    avg_sky_temp =  djs_avsigclip(skystack, 2, sigrej = sigrej $
                                  , inmask = maskstack $
                                  , outmask = outmask)
    smashmask = reform((total(outmask, 2) EQ nfiles), nx, ny) 
ENDIF ELSE BEGIN
    avg_sky_temp = skystack
    smashmask = reform(maskstack GT 0.0, nx, ny)
ENDELSE

IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN
    pixset = long_wavepix(reform(avg_sky_temp, nx, ny), tset_slits $
                          , FWHM = FWHM, pkwdth = pkwdth, toler = toler)
    piximg = long_wpix2image(pixset, tset_slits)
    inslit = where(slitmask NE 0)
    dsamp = 10

    IF nfiles GT 1 THEN BEGIN 
       FOR j = 0L, nfiles-1L DO BEGIN
          fitpix  = where(slitmask NE 0 AND $
                          finite(avg_sky_temp/skystack[*, j]) EQ 1 AND $
                          finite(ivarstack[*, j]) EQ 1 AND $
                          avg_sky_temp/skystack[*, j] GT 0.0 AND $
                          avg_sky_temp/skystack[*, j] LT 10.0 AND $
                          maskstack[*, j] EQ 0 AND $
                          smashmask EQ 0 AND $
                          skystack[*, j] GT 0.0 AND $
                          avg_sky_temp GT 0.0 AND $
                          ximg GT 0.1  AND $
                          ximg LT 0.9, npix)
          psort1 = sort(piximg[fitpix])
          isamp = lindgen(npix/dsamp)*dsamp
          psort = psort1[isamp]
          sset_j = bspline_iterfit(piximg[fitpix[psort]] $
                                   , avg_sky_temp[fitpix[psort]] $
                                   /skystack[fitpix[psort], j] $
                                   , nord = 3, upper = 3.0, lower = 3.0 $
                                   , bkspace = 1.2D, yfit = fitans)
          skystack[inslit, j] = $
             skystack[inslit, j]*bspline_valu(piximg[inslit], sset_j)
          img_minsky = reform(skystack[*, j]-avg_sky_temp, nx, ny)
;       Now go through and trace and mask bright objects
;       NIRI plate scale is 0.1171 at f/6. Assuming 0.6" seeing FWHM = 5.12 pix
          FWHM = 5.12
          obj = long_objfind(img_minsky, tset_slits = tset_slits $
                             , FWHM = FWHM, OBJMASK = OBJMASK $
                             , peakthresh = 0.1, /SILENT $
                             , OBJTHRESH = 0.01D) ;0.001D)
          IF KEYWORD_SET(obj) THEN BEGIN
             nfound = n_elements(obj)
             splog, 'Found nboj = ' + strcompress(string(nfound), /rem) $
                    + ' objects on skyimage: ' + skyfiles[j] +  ' for masking'
             objpix = WHERE(reform(objmask, nx*ny) EQ 1, nobj)
             IF nobj NE 0 THEN maskstack[objpix, j] = 1
          ENDIF
       ENDFOR
       avg_sky = djs_avsigclip(skystack, 2, sigrej = sigrej $
                               , inmask = maskstack $
                               , outmask = outmask)
       smashmask = reform((total(outmask, 2) EQ nfiles), nx, ny)
    ENDIF ELSE If nfiles EQ 1 THEN BEGIN
        avg_sky = skystack
        smashmask = reform(maskstack GT 0.0, nx, ny)
    ENDIF
    avg_sky = reform(avg_sky, nx, ny)
 ENDIF


;   if total = nfiles this means mask it was rejected for each image
;   switch mask convention here
finalmask = (mask EQ 0) AND (smashmask EQ 0)
wpix_img = reform(avg_sky_temp, nx, ny)
; since avg_sky has been scaled by a wavelength dependent scaling 
; possibly altering line shapes, use avg_sky_temp for wavelengths
hdr_wav = hdr

IF KEYWORD_SET(TELLURIC) THEN BEGIN
; Read in Telluric sky files
    nfiles1 = n_elements(skyfiles1)
    skystack  = fltarr(nx*ny, nfiles1)
    IF (NOT keyword_set(sigrej1)) then begin
        if (nfiles1 LE 2) then sigrej = 1.0 $ 
; Irrelevant for only 1 or 2 files
        else if (nfiles1 EQ 3) then sigrej = 1.1 $
        else if (nfiles1 EQ 4) then sigrej = 1.3 $
        else if (nfiles1 EQ 5) then sigrej = 1.6 $
        else if (nfiles1 EQ 6) then sigrej = 1.9 $
        else sigrej = 2.0
    ENDIF ELSE sigrej = sigrej1
    sigrej = siglower*sigrej
    FOR j = 0L, nfiles1-1L DO BEGIN
        niri_proc, skyfiles1[j], imag1, ivar1 $
                   , hdr = hdr, flatimg = flat, /TELLURIC
        skystack[*, j] = imag1
     ENDFOR
    IF nfiles1 GT 1 THEN avg_sky =  djs_median(skystack, 2) $
    ELSE avg_sky = skystack
    avg_sky = reform(avg_sky, nx, ny)
    finalmask = (ivar GT 0.0)
ENDIF

; Create wavelength QA file
qafile = targdir + '/waveQA-' + gnirs_fileprefix(filename[0]) + '.ps'

; Compute NIRI wavelengths
waveimg = niri_waveimg(wpix_img, hdr_wav, piximg = piximg, QAFILE = QAFILE $
                       , CHK = WVCHK)
;   First pass sky subtraction.  Also corrects for bias and dark counts
print, "  First pass sky subtraction"
;img_minsky = sciimg -avg_sky

; Mask bad pixels 
ivar = ivar*finalmask

RETURN
END

; RETURN
; END
; ENDIF ELSE IF nfiles GE 5 THEN BEGIN
;     avg_sky = djs_median(skystack, 2)
;     outmask = 0*maskstack + 1
; ;   swapped mask convention for djs_reject
;     inmask = maskstack EQ 0
;     FOR j = 0L, nfiles-1L DO BEGIN
;         qdone = djs_reject(skystack[*, j], avg_sky $
;                            , invvar = ivarstack[*, j] $
;                            , inmask = inmask[*, j] $
;                            , upper = sigrej, lower = sigrej $
;                            , outmask = outmask1)
;         outmask[*, j] = (outmask1 EQ 0)
;     ENDFOR
;     smashmask = reform((total(outmask, 2) EQ nfiles), nx, ny) 
