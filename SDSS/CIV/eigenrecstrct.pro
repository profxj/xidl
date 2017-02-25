;+
; NAME:
;   EIGENRECONSTRUCT
;     Version 2.0
;
; AUTHOR:
;   Melodie M. Kao
;   Massachusetts Institute of Technology
;   Kavli Institute for Astrophysics and Space Research
;   77 Massachusetts Avenue, Building 37-287, Cambridge, MA 02139 
;   melodie.kao@alum.mit.edu
;   mkao@caltech.edu
;
; PURPOSE:
;   Fits predetermined eigenspectra to a spectrum.  Primary function:
;   defines pixel weighting and finds eigen coefficients.
;
;
; CALLING SEQUENCE:
;   eigflux = EIGENRECSTRCT( wave,flux, error, allEigflux, 
;                             [rejectedflux=, eigerror=, weight=,
;                             fail=, _extra=] )
;
;
; DESCRIPTION: 
;
;   EIGENRECSTRCT will shift spectra into a given rest frame and
;   fit given eigenspectra to the spectra.  The eigenspectra must be
;   an ascii table with the first column as z=0.0 wavelengths and the
;   second column as fluxes (example: See Yip C. W. et al., 2004, AJ,
;   128, 2603 for an example of correctly formatted files).  They will
;   be fed into the program as a list of file names that include the
;   complete path to the eigenspectra.
;   The spectra will be returned as 10^(lambda_0 + i*dLoglambda_0).  
;
;
; INPUTS:
;   wave    -- Array of wavelengths associated with fluxes.  Must be same
;              length as the flux array.
;   flux    -- Array of spectrum fluxes associated with each wavelength in
;              the wave array.  The flux and wavelength must be converted 
;              into 10^(lambda_0 + i*dLoglambda_0).
;   error   -- The associated error for each element of the flux array.
;              Must be the same length as the flux and wave arrays.
;   eigfile -- Array of eigenspectra file names to be used in fitting
;
;
; RETURNS:
;   eigflux      -- Double Array of continuum found through
;                   eigenfitting, length is same length as wave/flux
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   rejectedflux=-- Double Array of flux that was not fitted.  Each
;                   fitted pixel has rejectedflux value NaN.  Each
;                   rejected pixel has rejectedflux value equal to
;                   original flux of that pixel.  This is help with
;                   debugging plots in eigqsoconti, so that one can
;                   see exactly which pixels with what flux are not
;                   getting fitted.  If trying to index with pixels
;                   are getting rejected, use: WHERE(FINITE(rejectedflux,/NAN))
;   eigerror= -- errors from eigencontinuum fit ONLY (*not* total
;                errors including flux errors)
;   weight= -- weighting used for fit
;   fail= -- error matrix status
;
; COMMENTS: 
;
; EXAMPLES:
;   See usage in eigqsoconti.pro
;
; PROCEDURES/FUNCTIONS CALLED:
;   eigenrecstrct_mij.pro  
;   interpspec.pro              FIRE/Spextool/pro/
;
; MODIFICATION HISTORY:
;  25 April 2011   Written by M. Kao
;  28 July 2011 Major reorg and overhaul, KLC
;
;-
; Copyright (C) 2011, Melodie Kao
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;------------------------------------------------------------------------------
@eigenrecstrct_mij              ; compile extra functions

function eigenrecstrct, wave, flux, error, allEigflux, rejectedflux=rejectedflux, $
                        EIGERROR=eigerror, WEIGHT=weight, FAIL=fail, _EXTRA=extra

;  if n_params() ne 4 then begin
;     print,'Syntax - eigenrecstrct(wave, flux, error, allEigflux, [rejectedflux=, '
;     print,'                       EIGERROR=, WEIGHT=, FAIL=, _EXTRA=])'
;     return, -1
;  endif 

  sz = size(allEigflux,/dim)
  neig = sz[0]
  npix = sz[1]
  
  ;; ;;;;;;;;;;;;;
  ;; Set weighting  
  ;; weight = 0 if there are gaps in data; if data is too negative
  ;; flux; and if error is 0 or infinite
  weight = DBLARR(npix,/nozero)                      ; weight of each pixel
  rejectedflux = replicate(!VALUES.D_NAN, npix) ; For plotting purposes. NaN = not fitted

  zeroes = WHERE(error LE 0.0, nzeroes)    
  IF nzeroes NE 0 THEN $ 
     rejectedflux[zeroes] = flux[zeroes]
  notzeroes    = WHERE(error NE 0.0,complement=bd) ; weight = 1/error^2  for all else
  weight[notzeroes] = 1.0/error[notzeroes]^2
  if bd[0] ne -1 then weight[bd] = 0.0 ; save time by not instantiating all
  
  nans = WHERE(FINITE(error, /NAN) or FINITE(flux, /NAN), nnans)
  IF nnans NE 0 THEN BEGIN
     weight[nans] = 0.0
     rejectedflux[nans]   = flux[nans]
  ENDIF 

  negs = WHERE(flux LT -3.0*error, nnegs)  ; Too negative flux
  IF nnegs NE 0 THEN BEGIN
     weight[negs] = 0.0   
     rejectedflux[negs]    = flux[negs] ; KLC: I added this; isn't it true?
  ENDIF 
  
  
  ;; ;;;;;;;;;;;;;;;;;; 
  ;; CONSTRUCT SUPERPOSITION OF EIGENSPECTRA 
  m_invert = EIGENRECSTRCT_MIJ(weight, allEigflux, FAIL=fail, _extra=extra) ; /silent
  
  
  ;;; f_j defined in Connolly, A.J. & Szalay, A.S. 1999, AJ, 117, 2052
  ;;; between Eq (4) and Eq (5) 
  tmp = weight*flux
  if nnans ne 0 then tmp[nans] = 0. 
  f_j = alleigflux # tmp       ; total( 1st array's 1st dimen X 2nd array's 2 dimen)
  
  ;;; all_coeff: Connolly, A.J. & Szalay, A.S. 1999, AJ, 117, 2052
  ;;; a_i [Eq (5)]
  f_jarr = replicate(1.d,neig) # f_j ; repeat f_j in the 1st dimension
  all_coeff = m_invert*f_jarr

  ;;; Reconstruction of eigfluxes to determine continuum
  coeff = TOTAL( all_coeff, 2 ) 
  tmp = coeff # allEigFlux     ; KLC: becomes a [1,npix] array
  eigflux = tmp[0,*]
  

  ;;; Find errors due to fitting
  IF KEYWORD_SET(eigerror) THEN BEGIN
     ;;; Connolly, A.J. & Szalay, A.S. 1999, AJ, 117, 2052
     ;;; Eq. (6)
     n = TOTAL(weight, /NAN)
     IF FINITE(n, /NAN) THEN stop,       $
        "eigenrecstrct: Cannot calculate covariance matrix"
     covar = (1./n)*m_invert
     
     ;;; Error propagation equation (3.14 in Bevington and Robinson,
     ;;; Data Reduction and Error Analysis for the Physical Sciences,
     ;;; 2nd edition)
     eigerror  = DBLARR(npix,/nozero)
     for ii=0L,npix-1 do $
        eigerror[ii] = total(covar*(allEigflux[*,ii] # allEigflux[*,ii]))

  ENDIF                         ; /eigerror

;  bd = where(finite(eigflux,/nan),nbd)
;  if nbd eq npix then $
;     stop,'eigenrecstrct() stop: failed to have finite eigflux'
  return, eigflux

end
