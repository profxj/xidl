;+ 
; NAME:
; lco_targets
;   Version 1.1
;
; PURPOSE:
;    Turns a Sextractor output file into an IDL structure
;
; CALLING SEQUENCE:
;   lco_targets, name, galstr
;
; INPUTS:
;   name -- Name of Sextractor file
;
; RETURNS:
;   galstr -- IDL galaxy structure
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
;   ?? Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lco_targets, name, galstr

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'lco_targets, targ_fil, targ_str' 
    return
  endif 

 readcol, name,idg,ragal,decgal,Rgal,xgal,ygal,sgal,kgal,algal, $
   FORMAT='I,F,F,F,F,F,F,F,F,F'

 tmp = { $
         idg: 0L, $
         ra: 0., $
         dec: 0., $
         R: 0., $
         xpix: 0., $
         ypix: 0., $
         sa: 0., $
         ka: 0., $
         al: 0. $
       }

 galstr = replicate(tmp, n_elements(idg))

 galstr.idg = idg
 galstr.ra = ragal
 galstr.dec = decgal
 galstr.R = Rgal
 galstr.xpix = xgal
 galstr.ypix = ygal
 galstr.sa = sgal
 galstr.ka = kgal
 galstr.al = algal

 return

end
