;+ 
; NAME:
; fuse_bindat
;   Version 1.0
;
; PURPOSE:
;    Rebin a single data set to a new wavlength scale
;      Simple adding (no weighting by S/N)
;
; CALLING SEQUENCE:
;   
;   fuse_bindat, gdpix, orig_wv, orig_fx, newwv, newfx
;
; INPUTS:
;   orig_wv
;   orig_fx
;   newwv
;
; RETURNS:
;
; OUTPUTS:
;   newfx
;
; OPTIONAL KEYWORDS:
;  VAR         
;
; OPTIONAL OUTPUTS:
;  NEWVAR      
;
; COMMENTS:
;
; EXAMPLES:
;   x_specrebin, gdpix, orig_wv, orig_fx, newwv, newfx
;
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
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'fuse_bindat, orig_dat, new_dat, nbin [V1.0]'
    return
  endif 

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
