;+ 
; NAME:
; fuse_bindat
;   Version 1.1
;
; PURPOSE:
;    Rebin FUSE data from 'Raw' to 'useful'
;
; CALLING SEQUENCE:
;   
;   fuse_bindat, orig_dat, new_dat, nbin
;
; INPUTS:
;   orig_dat -- FITS file of original FUSE data
;   new_dat -- Name of new data file
;   [nbin] -- Number of pixels to smooth over [default: 3L]
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_bindat, 'PKS0405_orig.fits', 'PKS0405_new.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fuse_bindat, orig_dat, new_dat, nbin


;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'fuse_bindat, orig_dat, new_dat, [nbin] (V1.1)'
    return
  endif 

; Optional keywords
  if not keyword_set( NBIN ) then nbin = 3L

; Open data

  orig = xmrdfits(orig_dat, 1, /silent)

  ;; Smooth
  s_fx = smooth(orig.flux, nbin)
  s_wv = smooth(orig.wave, nbin)
  s_err = smooth(orig.error, nbin) / sqrt(nbin)  ;; Not quite right but ok

  ;; Pick out the right pixels
  npix = n_elements(orig.flux)
  strt = (nbin-1)/2
  all_pix = strt + lindgen(npix/nbin - 1)*nbin

  ;; Final
  fin_str = { $
              wave: s_wv[all_pix], $
              flux: s_fx[all_pix], $
              error: s_err[all_pix] $
            }

  ;; 
  print, 'fuse_bindat: Creating ', new_dat
  mwrfits, fin_str, new_dat, /create
  

return
end
