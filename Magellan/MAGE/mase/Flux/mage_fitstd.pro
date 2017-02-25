; + 
; NAME:
; mage_fitstd
; Version 0.1
;
; PURPOSE:
;  Determines a fit to a standard star and saves it so it can be 
;  used to flux other observations.  
;
; CALLING SEQUENCE:
;
;  mage_fitstd,fluxtable,objstr,outfil 
;
; INPUTS:
;   fluxtable - The file path of an ESO standard star file to
;               calculate the fit from
;   objstr    - The object structure generated from the mage_script
;               extraction routines
;
; RETURNS:
;
; OUTPUTS:
;
;   outfil    - The file path of an IDL save file that contains the
;               fit results
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mage_fitstd,'fgd108.dat',obj_strct,'gd108cal.sav'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;    mage_echfitstd
;
; REVISION HISTORY:
;   16-Jun-2008 CLW


pro mage_fitstd,fluxtable,objstr,outfile

;Read in the calibration table
readcol,fluxtable,swv,sfx,format='D,D'

objstr.ORDER = -objstr.ORDER + 21L
;Normalize the exposure time
for i=0,N_elements(objstr)-1 do objstr[i].fx/=objstr[i].exp
for i=0,N_elements(objstr)-1 do objstr[i].var/=((objstr[i].exp)^2)
for i=0,N_elements(objstr)-1 do objstr[i].sig/=objstr[i].exp
nordr = 15L
;; If the standard does not go red enough, extrapolate it as a power law
maxwv = max(objstr[nordr-1L].WAVE)
IF maxwv GT max(swv) THEN BEGIN ; & $
   ipix = WHERE(swv GT 5000.0 AND sfx GT 0.0) ; & $
   dswv = djs_median((abs(swv - shift(swv, 1)))[ipix]) ; & $
   nadd = ceil((maxwv - max(swv))/dswv) ; & $
   swv_add = max(swv) + (dindgen(nadd) + 1.0d)*dswv ; & $
   ;; fit 
   sfx_coeff = ladfit(alog10(swv[ipix]), alog10(sfx[ipix])) ; & $
   logsfx_add   = poly(alog10(swv_add), sfx_coeff) ; & $
   sfx_add = 10.0d^logsfx_add ; & $
   swv = [swv, swv_add] ; & $
   sfx = [sfx, sfx_add] ; & $
ENDIF

mage_echfitstd,objstr,swv,sfx,outfile, NFITORD = 15, FUNC = 'CHEBY'


end
