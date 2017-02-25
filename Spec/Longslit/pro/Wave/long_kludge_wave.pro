;+
; NAME:
;   long_kludge_wave
;
; PURPOSE:
;   Create dummy files given a wavelength solution (lambda vector) and
;   a template file
;   
;
; CALLING SEQUENCE:
;   long_kludge_wave, lambda, templ_fil, out_root
;
; INPUTS:
;   lambda       - Array of wavelength values (Ang)
;   templ_fil    - Wavelength FITS file used for pixset and fwhmset
;   out_root     - ROOT name for output files
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;  x=findgen(4096L)
;  lambda = 9345.6 + 5.7939e-1 * x + 3.1125e-6 * x^2 - 7.9696e-10 * x^3
;  long_kludge_wave,  lambda,  'wave_aug2010.fits'
;  touch Raw/wave_kludge.fits
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Jan-2011  Written by J.X.P.
;-
;------------------------------------------------------------------------------
pro long_kludge_wave, lambda, templ_fil, out_root, OFIT=ofit

  if not keyword_set(OFIT) then ofit = x_setfitstrct()
  if not keyword_set(OUT_ROOT) then out_root = 'wave-arc_kludge'

  ;; Sets
  pixset = xmrdfits(templ_fil, 1) ;; Has to have the same number of slits as your kludge
  pixset[0].coeff2d[*] = 0.  
  pixset[1].coeff2d[*] = 0.
  fwhmset = xmrdfits(templ_fil, 2)

  ;; Fit
  npix = n_elements(lambda)
  fit = x_fit(findgen(npix), lambda, fitstr=ofit)
  outfil = OUT_ROOT+'.sav'
  print, 'long_kludge_wave: Writing ', outfil
  xfit = replicate(ofit, n_elements(pixset))
  save, xfit, filename=outfil


  ;; Write
  outfil = OUT_ROOT+'.fits'
  print, 'long_kludge_wave: Writing ', outfil
  mwrfits, pixset, outfil, /create
  mwrfits, fwhmset, outfil

  return
end
;------------------------------------------------------------------------------
