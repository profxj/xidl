; NAME:
;   kpno_superflat
;
; PURPOSE:
;   Generic routine to construct a unit normalized pixel flat for longslit
;   spectra. This is the same code as gsuperflat (which constructs 
;   pixelflats for GMOS) but without the gap masking which is specific to 
;   that detector. Also, things are transposed here relative to gsuperflat.
;   This routine is tuned to the instrument at KPNO
;
; CALLING SEQUENCE:
;      kpno_superflat, filenames, superflatfile, $
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
;   long_proc
;   mrdfits
;   mwrfits
;   splog
;
; REVISION HISTORY:
;   22-Aug-2005 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------

PRO kpno_superflat, filenames, superflatfile $
                    , biasfile = biasfile $
                    , VERBOSE = VERBOSE $
                    , INDIR = INDIR, TEMPDIR = TEMPDIR $
                    , WRITE_FLATS = WRITE_FLATS, SLIT_MED = SLIT_MED

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

long_proc, filenames[i], flat1, invvar1 $
           , biasfile = biasfile, VERBOSE = VERBOSE
dims = size(flat1, /dimen) 
nrow = dims[0]
ncol = dims[1]
flatfull = rebin(flat1[*], nrow*ncol, nfiles)

FORMAT   = '(%"Masked %d pixels for CCD # %d in flatfile # %d")'
FOR i = 0, nfiles-1 do begin 
    splog, 'Working on file ', filenames[i]
    IF i GT 0 THEN  $
      long_proc, filenames[i], flat1, invvar1, biasfile = biasfile $
      , VERBOSE = VERBOSE
    slit_med = djs_median(flat1, 2) > 0
;   divide out the median spatial profile
    flat1 = $
      flat1/ $
      ((slit_med + (slit_med EQ 0)) # (replicate(1, ncol)))
;   now do a polynomial fit column by column to the spatial residuals
    spec_med = djs_median(flat1, 1) > 0
    specimage = replicate(1.0, nrow) # spec_med
    fitmask = (flat1 GT 0.8*specimage)
    invvar1 = invvar1*fitmask*(slit_med^2 # replicate(1, ncol))
;       now do a polynomial fit column by column to the spatial residual 
    xrow = findgen(nrow)#replicate(1.0, ncol)
    xy2traceset, xrow, flat1, invvar = invvar1 $
                 , tset, ncoeff = 5 $
                 , yfit = slitfit, upper = 3, lower = 3, /silent
;       smooth the polynomial coefficients
    smooth_set = tset
    FOR k = 0, n_elements(tset.COEFF[*, 0])-1L DO $
      smooth_set.COEFF[k, *] = median(transpose(tset.coeff[k, *]), 5)
    traceset2xy, smooth_set, xrow, slitsmooth
;   remove the spatial residuals
    flat1 = flat1/(slitsmooth + (slitsmooth EQ 0))
    test_dev = total(abs(flat1-1.0) GT 0.1, 1)
    finalmask = (test_dev LT nrow/2) ## replicate(1.0, nrow)
    flat1 = flat1*finalmask
    IF KEYWORD_SET(WRITE_FLATS) THEN BEGIN
        mwrfits, flat1, tempfiles[i], /create
    ENDIF
    flatfull[*, i] = flat1 
endfor

IF nfiles GT 1 THEN $
  flatfinal = reform(djs_avsigclip(flatfull, 2, sigrej = sigrej), nrow, ncol) $
ELSE flatfinal = reform(flatfull[*], nrow, ncol)
mwrfits, flatfinal, superflatfile, /create

return

END
