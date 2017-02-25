;+ 
; NAME:
; hamspec_makeflats   
;     Version 1.1
;
; PURPOSE:
;    Combines all flats for a given setup into one final image.
;    ;
; CALLING SEQUENCE:
;   
;  hamspec_makeflats, hamspec, setup, /CLOBBER, /USEBIAS, /SMOOTH, /AVG
;
; INPUTS:
;   hamspec     -  HIRES structure
;   setup    -  Setup identifier 
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup and side
;  (e.g. 'Flats/Flat/_B_01_T.fits.gz')
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite the final fits file
;   /USEBIAS - Use the bias frame in bias subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_mktflat, hamspec, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_subbias
;  xcombine
;  hamspec_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_makeflats, hamspec, setup,CLOBBER=clobber, USEBIAS=usebias, $
                       AVG=avg, SMOOTH=smooth, _EXTRA=extra


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_makeflats, hamspec, setup, /CLOBBER [v1.0]'
      return
  endif 

  hamspec_mkflats,hamspec,setup,'TFLT', CLOBBER=clobber, USEBIAS=usebias, $
                   AVG=avg, _EXTRA=extra
  hamspec_mkflats,hamspec,setup,'TFLT', CLOBBER=clobber, USEBIAS=usebias, $
                   AVG=avg, _EXTRA=extra
  hamspec_mkmflat,hamspec,setup, CLOBBER=clobber, USEBIAS=usebias, $
                   SMOOTH=smooth, _EXTRA=extra

end
