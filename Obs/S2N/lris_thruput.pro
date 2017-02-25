;+ 
; NAME:
; lris_thruput
;    Version 1.1
;
; PURPOSE:
;    Passes back the end-to-end ESI+Keck thruput as a function of
;    wavelength.  Uses a sensitivity function generated from
;    observations of a Standard star through the 6" slit.
;
; CALLING SEQUENCE:
;  thru = lris_thruput(wave)
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
function lris_thruput, wave, str_instr, BIDX=BIDX, RIDX=RIDX

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'thru = lris_thruput(wave, str_instr) [v1.0]'
      return,-1
  endif 

  ;; Initialize
  bidx = -1
  ridx = -1
  thru = fltarr(n_elements(wave))

  ;; Dichroic first
  case str_instr[0].dichroic of
     ;'d46': bidx = where(wave LT 4600., nbwv, complement=ridx, ncomplement=nrwv)
     'D560': bidx = where(wave LT 5600., nbwv, complement=ridx, ncomplement=nrwv)
     else: stop
  endcase

  ;; ;;;;;;;;
  ;; Blue side
  if nbwv GT 0 then begin
     case str_instr[0].grating of
        'B600': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISb_600_4000_D560.fits'
              else: stop
           endcase
        end
        'B300': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISb_300_5000_D560.fits.gz'
              'D680': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISb_300_5000_D680.fits.gz'
              else: sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_Kastb830_3460_d46.fits'
           endcase
        end
        else: stop
     endcase
     sens = xmrdfits(sens_fil, 2)
     thru[bidx] = interpol(sens.eff, sens.wav, wave[bidx])
     ;; No extrapolation
     extrap = where(wave[bidx] LT sens.wav[0], next)
     if next GT 0 then thru[bidx[extrap]] = sens.eff[0]
  endif

  ;; ;;;;;;;;
  ;; Red side
  if nrwv GT 0 then begin
     case str_instr[1].grating of
        '600/7500': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISr_600_7500_D560.fits'
              else: stop 
           endcase
        end
        '400/8500': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISr_400_8500_D560.fits'
              else: stop 
           endcase
        end
        '600/10000': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISr_600_10000_D560.fits'
              else: stop 
           endcase
        end
        '1200/9000': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                    getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISr_1200_9000_D560.fits'
              else: stop 
           endcase
        end
        '831/8200': begin
           case str_instr[0].dichroic of 
              'D560': sens_fil = $
                 getenv('XIDL_DIR')+'/Obs/S2N/THRU_PUT_DATA/sens_LRISr_831_8200_D560.fits'
              else: stop 
           endcase
        end
        else: stop
     endcase
     sens = xmrdfits(sens_fil, 2)
     thru[ridx] = interpol(sens.eff, sens.wav, wave[ridx])
     extrap = where(wave[ridx] LT sens.wav[0], next)
     if next GT 0 then thru[ridx[extrap]] = sens.eff[0]
  endif

  thru = thru > 1e-5  ;; Avoid negative thruputs
 
  return, thru
end

