;+ 
; NAME:
; esi_thruput
;    Version 1.1
;
; PURPOSE:
;    Passes back the end-to-end ESI+Keck thruput as a function of
;    wavelength.  Uses a sensitivity function generated from
;    observations of a Standard star through the 6" slit.
;
; CALLING SEQUENCE:
;  thru = esi_thruput(wave)
;
; INPUTS:
;  wave=  -- Wavelength to calculate throughputs at
;
; RETURNS:
; Thruput -- 0-1 value with 1=100% efficiency
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
;       23-Mar-2010 Written by JXP
;       24-May-2011     GDW     Implement default value for XIDL_DIR
;                               if that envar is not defined (required for 
;                               functionality with ION);
;                               add SILENT keyword on XMRDFITS
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;
function esi_thruput, wave

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'thru = esi_thruput(wave) [v1.0]'
      return,-1
  endif 

  ;; define location of output directory...
  dir = getenv('XIDL_DIR')
  if dir eq '' then dir = '/local/home/randyc/idl/xidl/'

  ;; Read sensitivity file
  if not keyword_set(SENS_FIL) then $
    sens_fil = dir + 'Obs/S2N/THRU_PUT_DATA/sens_ESI.fits'
  sens = xmrdfits(sens_fil, 2, /SILENT)
  thru = interpol(sens.eff, sens.wav, wave)
 
  return, thru
end

