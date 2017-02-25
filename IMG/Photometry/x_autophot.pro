;+
; NAME:
;      x_autophot
; PURPOSE:
;      Compute concentric aperture photometry (adapted from DAOPHOT) 
;    I have adopted this from the Goddard routines
; EXPLANATION:
;     APER can compute photometry in several user-specified aperture radii.  
;     A separate sky value is computed for each source using specified inner 
;     and outer sky radii.   
;
; CALLING SEQUENCE:
;     APER, image, xc, yc, [ mags, errap, sky, skyerr, phpadu, apr, skyrad, 
;                       badpix, /EXACT, /FLUX, PRINT = , /SILENT, SETSKYVAL = ]
; INPUTS:
;     IMAGE -  input image array
;     XC     - vector of x coordinates. 
;     YC     - vector of y coordinates
;
; OPTIONAL INPUTS:
;     PHPADU - Photons per Analog Digital Units, numeric scalar.  Converts
;               the data numbers in IMAGE to photon units.  (APER assumes
;               Poisson statistics.)  
;     APR    - Vector of up to 12 REAL photometry aperture radii.
;     SKYRAD - Two element vector giving the inner and outer radii
;               to be used for the sky annulus
;     BADPIX - Two element vector giving the minimum and maximum value
;               of a good pix (Default [-32765,32767]).    If BADPIX[0] is
;               equal to BADPIX[1] then it is assumed that there are no bad
;               pixels.
;
; OPTIONAL KEYWORD INPUTS:
;     /EXACT -  By default, APER counts subpixels, but uses a polygon 
;             approximation for the intersection of a circular aperture with
;             a square pixel (and normalize the total area of the sum of the
;             pixels to exactly match the circular area).   If the /EXACT 
;             keyword, then the intersection of the circular aperture with a
;             square pixel is computed exactly.    The /EXACT keyword is much
;             slower and is only needed when small (~2 pixels) apertures are
;             used with very undersampled data.    
;     /FLUX - By default, APER uses a magnitude system where a magnitude of
;               25 corresponds to 1 flux unit.   If set, then APER will keep
;              results in flux units instead of magnitudes.
;     PRINT - if set and non-zero then APER will also write its results to
;               a file aper.prt.   One can specify the output file name by
;               setting PRINT = 'filename'.
;     /SILENT -  If supplied and non-zero then no output is displayed to the
;               terminal.
;     SETSKYVAL - Use this keyword to force the sky to a specified value 
;               rather than have APER compute a sky value.    SETSKYVAL 
;               can either be a scalar specifying the sky value to use for 
;               all sources, or a 3 element vector specifying the sky value, 
;               the sigma of the sky value, and the number of elements used 
;               to compute a sky value.   The 3 element form of SETSKYVAL
;               is needed for accurate error budgeting.
;
; OUTPUTS:
;     MAGS   -  NAPER by NSTAR array giving the magnitude for each star in
;               each aperture.  (NAPER is the number of apertures, and NSTAR
;               is the number of stars).   A flux of 1 digital unit is assigned
;               a zero point magnitude of 25.
;     ERRAP  -  NAPER by NSTAR array giving error in magnitude
;               for each star.  If a magnitude could not be deter-
;               mined then ERRAP = 9.99.
;     SKY  -    NSTAR element vector giving sky value for each star
;     SKYERR -  NSTAR element vector giving error in sky values
;
; PROCEDURES USED:
;       DATATYPE(), GETOPT, MMM, PIXWT(), STRN(), STRNUMBER()
; NOTES:
;       Reasons that a valid magnitude cannot be computed include the following:
;      (1) Star position is too close (within 0.5 pixels) to edge of the frame
;      (2) Less than 20 valid pixels available for computing sky
;      (3) Modal value of sky could not be computed by the procedure MMM
;      (4) *Any* pixel within the aperture radius is a "bad" pixel
;
;       APER was modified in June 2000 in two ways: (1) the /EXACT keyword was
;       added (2) the approximation of the intersection of a circular aperture
;       with square pixels was improved (i.e. when /EXACT is not used) 
; REVISON HISTORY:
;       Adapted to IDL from DAOPHOT June, 1989   B. Pfarr, STX
;       Adapted for IDL Version 2,               J. Isensee, July, 1990
;       Code, documentation spiffed up           W. Landsman   August 1991
;       TEXTOUT may be a string                  W. Landsman September 1995
;       FLUX keyword added                       J. E. Hollis, February, 1996
;       SETSKYVAL keyword, increase maxsky       W. Landsman, May 1997
;       Work for more than 32767 stars           W. Landsman, August 1997
;       Converted to IDL V5.0                    W. Landsman   September 1997
;       Don't abort for insufficient sky pixels  W. Landsman  May 2000
;       Added /EXACT keyword                     W. Landsman  June 2000 
;       Allow SETSKYVAL = 0                      W. Landsman  December 2000 
;       Set BADPIX[0] = BADPIX[1] to ignore bad pixels W. L.  January 2001     
;       Fix chk_badpixel problem introduced Jan 01 C. Ishida/W.L. February 2001 
;-
pro x_autophot,image, OUTDIR=outdir

;                                
  if N_params() LT 1 then begin    ;Enough parameters supplied?
      print, $
        'Syntax - X_AUTOPHOT, image, OUTDIR= '
      return
  endif 

  ;; Setup output directory
  if not keyword_set(OUTDIR) then outdir = 'Phot'
  a = findfile(outdir+'/..', count=count)
  if count EQ 0 then file_mkdir, outdir
  a = findfile('Sex/..', count=count)
  if count EQ 0 then file_mkdir, 'Sex'

  ;; Grab image
  img = xmrdfits(image, /fscale, /silent)

  ;; Run sextractor

  return
end
