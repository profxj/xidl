;+ 
; NAME:
; fuse_norm
;   Version 1.1
;
; PURPOSE:
;    Normalize a set of FUSE data by a continuum array
;
; CALLING SEQUENCE:
;   
;   fuse_norm, orig_dat, cont_fil, fin_fil 
;
; INPUTS:
;   orig_dat -- FITS file of original FUSE data
;   cont_fil -- Continuum FITS file
;   new_dat -- Name of new data file
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fuse_norm, orig_dat, cont_fil, fin_fil


;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'fuse_norm, orig_dat, cont_fil, fin_fil [V1.1]'
    return
  endif 

; Open data

  orig = xmrdfits(orig_dat, 1, /silent)
  conti = xmrdfits(cont_fil, /silent)
  
  zro = where(conti LE 0., nzro)
  if nzro NE 0 then conti[zro] = 1.

  ;; Normalize
  fin_fx = orig.flux / conti
  fin_err = orig.error / conti

  ;; Final
  fin_str = { $
              wave: orig.wave, $
              flux: fin_fx,$
              error: fin_err $
            }

  ;; 
  mwrfits, fin_str, fin_fil, /create
  print, 'fuse_norm: Created ', fin_fil
  

return
end
