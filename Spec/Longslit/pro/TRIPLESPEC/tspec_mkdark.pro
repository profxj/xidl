;+
; NAME:
;        TSPEC_MKDARK
;
;
; PURPOSE:
;        Removes quadrant-to-quadrant crosstalk and capacitive
;        coupling effects in raw Palomar TSPEC spectra.
;
;
; CALLING SEQUENCE:
;        result = tspec_clean(image[,parameters accepted by readfits])
;
;
; INPUTS:
;        image: A 2048x1024 Palomar TSPEC image OR a Palomar TSPEC
;        fits file.  This way you can use it in the place of readfits
;
;
; OPTIONAL INPUTS:
;        Any named parameters accepted by readfits with be passed.
;
;
; KEYWORD PARAMETERS:
;        sigma: optinally NAN all deviant pixels changed by
;        sigma_filter from the Astro IDL library
;
;
; OUTPUTS:
;        result: A 2048x1024 DOUBLE array containing the cleaned
;        image, with crosstalk and capacitive coupling removed.
;
;
; OPTIONAL OUTPUTS:
;        header: String containing the header of the fits file passed,
;        if a fits file is passed for image.
;
;        variance: A 2048x1024 DOUBLE array containing the calculated
;        photon-noise variances of the image (sigma^2).
;
; RESTRICTIONS
;        Requires the IDL Astronomy User's Library: http://idlastro.gsfc.nasa.gov/contents.html
;
; SIDE EFFECTS:
;        Converts image to double.
;
;
; MODIFICATION HISTORY:
;  Written by JXP (April 2011)
;       Version 1.0
;-
;; tspec_mkdark, ['file1', 'file2'], 'dark_300s.fits'
pro tspec_mkdark, files, outfil 


  ;; Input the files (clean along the way)
  nfil = n_elements(files)
  all_img = fltarr(1024L,2048L,nfil)

  ;; 
  for qq=0L,nfil-1 do $
     all_img[*,*,qq] = tspec_clean(files[qq])

  ;; Avg + sig clip
  comb = x_avsigclip(all_img, siglo=2.5, sighi=2.5, maxiter=3)

  bias = djs_median(comb[7:37,*],1)
  savfilt = savgol(60,60,0,4 < (12 - 1))
  bias_fct = convol(bias, savfilt,/EDGE_TRUNCATE)
  bias_img = replicate(1., 1024L) # bias
  comb = comb - bias_img
;  xatv, comb, min=-50., max=500, /bloc


  ;; Transpose and Write
  mwrfits, comb, outfil, /create

  return
end
