;+ 
; NAME:
; lris_calcs2n
;    Version 1.1
;
; PURPOSE:
;     This program computes count rates and expected S/N for
;     Keck spectrographs.
;
; CALLING SEQUENCE:
;  keck_calcs2n, iwv, flg
;
; INPUTS:
; iwv -- Input wavelengths
; flg -- Flag for the instrument (1=Old HIRES, 2=New HIRES, 3=ESI)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   imsize - Arcmin of image (default is 5')
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndchrt, 'targets.list'
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;     23/08/11 Written by JXP
;-
;------------------------------------------------------------------------------
pro lris_calcs2n, iwv, flg, INFIL=infil, NOPRINT=noprint, PLOT=plot, $
                   STATE=state, S2N=sn, IORDER=iorder, $
                   B_FSTRCT=b_fstrct, R_FSTRCT=r_fstrct, RED_IDX=ridx, BLUE_IDX=bidx

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'lris_calcs2n, wave, flg, [v1.1]'
      return
  endif 

  str_obs = state.str_obs
  str_instr = state.str_instr
  str_tel = state.str_tel
 
  nwv = n_elements(iwv)
  wave = iwv[sort(iwv)]
  
  ;;  grab through put from within subroutine and communicate
  thru = lris_thruput(wave, str_instr, BIDX=bidx, RIDX=ridx)  ;; This does blue and red 

  ;; S/N calculation
  ;; Blue
  if bidx[0] GE 0 then begin
     b_fstrct =  spec_calcs2n(wave[bidx], thru[bidx], str_tel, str_instr[0], str_obs) 
     sn = b_fstrct.sn
  endif
  ;; Red
  if ridx[0] GE 0 then begin
     r_fstrct =  spec_calcs2n(wave[ridx], thru[ridx], str_tel, str_instr[1], str_obs)
     if bidx[0] GE 0  then sn = [sn, r_fstrct.sn] $
     else sn = r_fstrct.sn
  endif

  return
end
