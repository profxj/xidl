; NAME:
;   long_superflat
;
; PURPOSE:
;   Generic routine to construct a unit normalized pixel flat for longslit
;   spectra. This is the same code as gsuperflat (which constructs 
;   pixelflats for GMOS) but without the gap masking which is specific to 
;   that detector. Also, things are transposed here relative to gsuperflat.
;
; CALLING SEQUENCE:
;      long_superflat, filenames, superflatfile, $
;                      [biasfile = , VERBOSE = ,INDIR = ,TEMPDIR = , $
;                       WRITE_FLATS = ]
;
; INPUTS:
;   filenames    - Input flat-field file names, either dome or twilight flats
;   superflatfile  - Output file name for super pixel flat
;
; OPTIONAL INPUTS:
;   biasfile     - Bias file to apply to raw images
;   verbose      - Verbose if set
;   tempdir      - Directory for temporary files when generating the flats;
;                  default to the current directory
;   write_flats  - If set write individual normalized flats to tempdir
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS: 
;
; COMMENTS:

; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   djs_avsigclip()
;   fileandpath()
;   lris_proc
;   mrdfits
;   mwrfits
;   splog
;
; REVISION HISTORY:
;   22-Aug-2005 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------

PRO niri_superflat, filenames, darkfiles, superflatfile $
                    , OBJFILE = OBJFILE, VERBOSE = VERBOSE $
                    , INDIR = INDIR, TEMPDIR = TEMPDIR $
                    , WRITE_FLATS = WRITE_FLATS, SLIT_MED = SLIT_MED

; Create the superdark
darkimg = niri_superdark(darkfiles)

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

i = 0
niri_proc, filenames[i], flat1, ivar1, darkimg = darkimg, VERBOSE = VERBOSE
dims = size(flat1, /dimen) 
nx = dims[0]
ny = dims[1]

oscan_left = 20
oscan_right = 975
; Allocate the flat stack
flatfull = fltarr(nx, ny, nfiles) + 1.0D

FOR i = 0, nfiles-1 do begin 
    splog, 'Working on file ', filenames[i]
    IF i GT 0 THEN  $
      niri_proc, filenames[i], flat1, ivar1, darkimg = darkimg $
      , VERBOSE = VERBOSE 
    subflat1 = flat1[oscan_left:oscan_right, *]
    subivar1 = ivar1[oscan_left:oscan_right, *]
    slit_med = djs_median(subflat1, 2) > 0
;   divide out the median spatial profile
    subflat1 = $
      subflat1/ $
      ((slit_med + (slit_med EQ 0)) # (replicate(1, ny)))
    subivar1 = subivar1*(slit_med^2 # replicate(1, ny))
;   divide out the median spectra profile
    spec_med = djs_median(subflat1, 1) > 0
    specimage = replicate(1.0, nx) # spec_med
    subflat1 = subflat1/(specimage + (specimage LT 0.0))
    subivar1 = subivar1*specimage^2
    nx_sub = n_elements(slit_med)
    fitmask = 0.0*subflat1
    FOR j = 0L, ny-1L DO BEGIN
        djs_iterstat, subflat1[*, j] $ ;, invvar = subivar1[*, j] $
                      , median = median, sigma = sigma, sigrej = 4.0 $
                      , mask = mask
        fitmask[*, j] = mask
;        IF total(mask) LE 10 THEN STOP
    ENDFOR
;    fitmask = (subflat1 GT 0.8 AND subflat1 LT 2.0)
;    invvar1 = invvar1*fitmask*(slit_med^2 # replicate(1, ny))
;   now do a polynomial fit column by column to the spatial residual 
    xrow = findgen(nx_sub)#replicate(1.0, ny)
    xy2traceset, xrow, subflat1, invvar = subivar1*fitmask $
                 , tset, ncoeff = 5 $
                 , yfit = slitfit, upper = 3, lower = 3, /silent $
                 , /sticky, /groupbadpix, grow = 1
;   smooth the polynomial coefficients
    smooth_set = tset
    FOR k = 0, n_elements(tset.COEFF[*, 0])-1L DO $
      smooth_set.COEFF[k, *] = median(transpose(tset.coeff[k, *]), 5)
    traceset2xy, smooth_set, xrow, slitsmooth
;   remove the spatial residuals
    subflat1 = subflat1/(slitsmooth + (slitsmooth EQ 0))
    subivar1 = subivar1*slitsmooth^2
;     fitmask2 = (subflat1 GT 0.0 AND subflat1 LT 2.0)
;     xrow = findgen(nx_sub)#replicate(1.0, ny)
;     xy2traceset, xrow, subflat1, invvar = subivar1*fitmask2 $
;                  , tset2, ncoeff = 5 $
;                  , yfit = slitfit2, upper = 3, lower = 3, /silent $
;                  , 
;     smooth_set2 = tset2
;     FOR k = 0, n_elements(tset2.COEFF[*, 0])-1L DO $
;       smooth_set2.COEFF[k, *] = median(transpose(tset2.coeff[k, *]), 3)
;     traceset2xy, smooth_set2, xrow, slitsmooth2
;     subflat1 = subflat1/(slitsmooth2 + (slitsmooth2 EQ 0))
    test_dev = total(abs(subflat1-1.0) GT 0.1, 1)
    finalmask = (test_dev LT nx_sub/2) ## replicate(1.0, nx_sub)
    subflat1 = subflat1*finalmask
    flatfull[oscan_left:oscan_right, *, i] = subflat1 
endfor

IF nfiles GT 1 THEN $
  flatfinal = reform(djs_avsigclip(flatfull, 3, sigrej = sigrej), nx, ny) $
ELSE flatfinal = reform(flatfull[*], nx, ny)
mwrfits, flatfinal, superflatfile, /create

return

END
