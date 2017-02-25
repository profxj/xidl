;+ 
; NAME:
; keck_calcs2n
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
;     HISTORY
;     ??/??/?? Written by G. Donald Penrod 
;     11/13/89 Modified by S. Vogt       - made inputs less confusing
;     06/08/92 Modified by M. Keane      - changed from Hamilton to HIRES 
;     02/07/96 Modified by C. Churchill  - structured queries by function
;                                        - set defaults for Decker C1
;                                        - added comments
;     20-Oct-2005 Ported to IDL by JXP
;     23-Mar-2011 Generatlized for "all" Keck spectrometers
;-
;------------------------------------------------------------------------------
pro keck_calcs2n, iwv, flg, INFIL=infil, NOPRINT=noprint, PLOT=plot, $
                  PSFILE=psfile, STATE=state, S2N=sn, IORDER=iorder, $
                  BLAZE=blaze, PIXEL=pixel, FSTRCT=fstrct


  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'keck_calcs2n, wave, flg, [v1.1]'
      return
  endif 

  ;; Parse infil
  if not keyword_set(CSZ) then csz = 1.5
  if not keyword_set(FLG) then flg = 1


  if not keyword_set(STATE) then begin
      ;; Instrument + Telescope
     case flg of
        1: x_inithires, str_instr, flg, STR_TEL=str_tel, $
                        DECKER=hires_decker, INFIL=infil ;; OLD HIRES
        2: x_inithires, str_instr, flg, STR_TEL=str_tel, $
                        DECKER=hires_decker, INFIL=infil ;; OLD HIRES
        3: x_initesi, str_instr, flg, STR_TEL=str_tel, $
                        DECKER=hires_decker, INFIL=infil ;; ESI
        4: x_initdeimos, str_instr, STR_TEL=str_tel, $
                         INFIL=infil ;; DEIMOS
        else: stop
     endcase
      ;; Observing conditions
      str_obs = x_obsinit( infil, /MAUNAKEA )
  endif else begin
      str_obs = state.str_obs
      str_instr = state.str_instr
      str_tel = state.str_tel
  endelse
 
  nwv = n_elements(iwv)
  wave = iwv[sort(iwv)]
  
  ;; Spectral stuff
  if flg LE 2 then begin
     m      = long(str_instr.MLAMBDA/wave)
     center = str_instr.MLAMBDA/m
     fsr = center/m
     
     low = where(center LT 3800.,nlow,complement=high,ncomplement=nhigh)
     sep = fltarr(nwv)
     iorder = lonarr(nwv)
     if nlow NE 0 then begin
        sep[low]    = 2.*str_instr.DELY*(str_instr.MLAMBDA/m[low] - $
                                         str_instr.MLAMBDA/(m[low]+1)) 
        iorder[low] = 2
     endif
     if nhigh NE 0 then begin
        sep[high]  = str_instr.DELY*(str_instr.MLAMBDA/m[high] $
                                     - str_instr.MLAMBDA/(m[high]+1)) 
        iorder[high] = 1
     endif
  endif

  ;;  grab through put from within subroutine and communicate
  case flg of
     1: thru = hires_thruput(wave, center, iorder, fsr, BLAZE=blaze) ; Old HIRES
     2: thru = hires_thruput(wave, center, iorder, fsr, flg=1) ; New HIRES
     3: thru = esi_thruput(wave) ; ESI
     4: thru = deimos_thruput(wave, str_instr)  ; DEIMOS
     else: stop
  endcase


  fstrct =  spec_calcs2n(wave, thru, str_tel, str_instr, str_obs)
  sn = fstrct.sn
  
  return
end
