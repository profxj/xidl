; NAME:
;   sinfoni_superflat
;
; PURPOSE:
;   Construct pixel and illumination flats for SINFONI
;
; CALLING SEQUENCE:
;   sinfoni_superflat, filenames, superpixflatfile, superillumflatfile, $
;    [ slitfile=, waveile=, biasfile=, $
;    npoly=, use_illum=, use_pixel=, tempdir=, $
;    sigrej=, maxiter=, /verbose, /write_flats ]
;
; INPUTS:
;   filenames    - Input flat-field file names, either dome or twilight flats
;   superpixflatfile  - Output file name for super pixel flat
;   superillumflatfile- Output file name for super illumination flat
;   wavefile     - File with wavelength solution
;
; OPTIONAL INPUTS:
;   slitfile     - File with trace sets describing the slit positions;
;                  if not set, then assume the entire image is a single slit.
;   biasfile     - Bias file to apply to raw images
;   pixflatfile  - Pixel flat to be used if an illumination 
;                  flat is being constructed
;     
;   tempdir      - Directory for temporary files when generating the flats;
;                  default to the current directory
;   npoly        - Number of polynomial terms for B-spline in the spatial
;                  direction; default to 3 for a cubic fit
;                  (But demand at least 10 pixels per row, on average,
;                  per degree of the polynomial.)
;   sigrej       - Rejection threshold in call to DJS_AVSIGCLIP()
;   maxiter      - Number of rejection iterations; default to 3.
;   use_illum    - Array the size of FILENAMES with 1s or 0s 
;                  which specifies which files are to be used for the 
;                  illumination function; default is to use all files
;   use_pixel    - Array the size of FILENAMES with 1s or 0s 
;                  which specifies which files are to be used for the 
;                  illumination function; default is to use all files
;   verbose      - Verbose if set
;   write_flats  - ???
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS: 
;
; COMMENTS:
;   Each slit of each input image has a model computed.  We divide each
;   input image by its corresponding model, and then take the sigma-clipped
;   mean of all those results.  The model is a B-spline in the wavelength
;   direction, and a low-order polynomial in the X direction.
;
;   The output header is copied from the first input file, with an
;   additional NEXP keyword (for the number of input files) and FILE*
;   keywords that list each input file name.
;  
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   djs_avsiglip()
;   fileandpath()
;   long_proc
;   long_slits2mask()
;   long_slits2x()
;   mrdfits
;   mwrfits
;   splog
;
; REVISION HISTORY:
;   10-Mar-2005 Written by J. Hennawi (UCB), D. Schlegel (LBL), S. Burles (MIT)
;   27-Feb-2012    Eliminated a 'bug' in the slitless flats with a
;   PIXIMG kludge [JXP]
;-
;------------------------------------------------------------------------------
PRO sinfoni_superflat, flatfiles, skyfiles, flatdarkfile, slitfile $
                       , pixflatfile, illumflatfile $
                       , darkfiles = darkfiles $
                       , tempdir = tempdir, verbose = verbose, chk = chk

  slitmask = mrdfits(slitfile, 0)
  tset_slits = mrdfits(slitfile, 1)

  waveimg = sinfoni_waveimg(skyfiles, tset_slits $
                            , darkfile = darkfile $
                            , piximg = piximg, CHK = WVCHK, QAFILE = QAFILE)
  nflat = n_elements(flatfiles)
  use_pixel = [lonarr(nflat) + 1L]
  use_illum = [lonarr(nflat) + 1L]
  
  long_superflat, flatfiles, pixflatfile, illumflatfile $
                  , waveimg = waveimg, piximg = piximg $
                  , slitfile = slitfile $
                  , darkfiles = flatdarkfile $
                  , use_illum = use_illum, use_pixel = use_pixel $
                  , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /SINFONI

  RETURN
  END


  

  
