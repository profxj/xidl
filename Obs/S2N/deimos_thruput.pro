;+ 
; NAME:
; deimos_thruput
;    Version 1.1
;
; PURPOSE:
;    Passes back the end-to-end ESI+Keck thruput as a function of
;    wavelength.  Uses a sensitivity function generated from
;    observations of a Standard star through the 6" slit.
;
; CALLING SEQUENCE:
;  thru = deimos_thruput(wave)
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
;   23-Mar-2010 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;
function deimos_thruput, wave, str_instr

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'thru = deimos_thruput(wave) [v1.0]'
      return,-1
  endif 

  ;; Read sensitivity file
  case strtrim(str_instr.grating,2) of 
     '600': begin
        case fix(str_instr.cwave) of
           5000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_600_550nm.fits'
           6000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_600_650nm.fits'
           7000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_600_750nm.fits'
           8000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_600_750nm.fits'
           else: stop
        endcase
     end
     '1200': begin
        case fix(str_instr.cwave) of
           5000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_1200_500nm.fits'
           6000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_1200_600nm.fits'
           7000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_1200_700nm.fits'
           8000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_1200_800nm.fits'
           else: stop
        endcase
     end
     '900': begin
        case fix(str_instr.cwave) of
           5000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_900_500nm.fits'
           6000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_900_600nm.fits'
           7000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_900_700nm.fits'
           8000: sens_fil = getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_DEIMOS_900_800nm.fits'
           else: stop
        endcase
     end
     else: stop
  endcase

  ;; Apply
  sens = xmrdfits(sens_fil, 1)
  thru = interpol(sens.eff, sens.wav, wave) > 0.
 
  return, thru
end

