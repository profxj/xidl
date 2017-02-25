;+ 
; NAME:
; lco_guidestr
;   Version 1.1
;
; PURPOSE:
;    Inputting a photometry structure, parse out a set of stars
;  useful for guide stars for WFCCD masks.
;
; CALLING SEQUENCE:
;   lco_guidestr, phot, mmin, mmax
;
; INPUTS:
;   phot -- Photometry structure
;   mmin -- Min mag for guide star
;   mmax -- Max mag for guide star
;
; RETURNS:
;
; OUTPUTS:
;   'Masks/guidestr.dat' -- Guide star output file
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
;   ?? Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lco_guidestr, phot, mmin, mmax

  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'lco_guidestr, phot, magmin, magmax [v1.1]' 
    return
  endif 

 close, /all

 cutstr = where(phot.R LT mmax AND phot.R GT mmin, ncut)
 plot, phot[cutstr].xpix, phot[cutstr].ypix, psym=1
 openw,1,'Masks/guidestr.dat'
  for j=0,ncut-1 do begin
    printf, 1, FORMAT = '(f,f)', phot(cutstr(j)).xpix, phot(cutstr(j)).ypix
  endfor
 close, 1
 return

end
