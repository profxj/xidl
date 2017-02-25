;+ 
; NAME: 
; x_pltarc
;    Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
; x_templarc, spec, lines, guessfit, /FFT, MSK=, MXSHFT=, $
;               MOCK_FFT=, SHFT=, ALL_PK=, PKWDTH=, $
;               /THIN, FORDR=, PKSIG=, SKIPLIN=, $
;               MXOFF=, /LOG
;
; INPUTS:
;   spec       - Input arc image or spectrum
;   lines      - Arc line structure
;  guessfit    - FIT to archived arc
;
; RETURNS:
;
; OUTPUTS:
;   lines     - Sets flg_plt to 1 those lines which are ID'd
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  ALL_PK= --  Pix values of the arc line peaks
;
; COMMENTS:
;
; EXAMPLES:
;   x_templarc, spec, lines, guessfit, SHFT=shft
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_fndpeaks
;
; REVISION HISTORY:
;   26-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_pltarc, tmplfil, ordr

;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_pltarc, tmplfil, ordr [v1.0]'
    	return
  endif 

  ;; Read in FIT file
  restore, tmplfil

  ;; Find order
  gdo = where(guess_ordr EQ ordr, ngd)
  if ngd EQ 0 then return

  ;; Plot
  sz = size(sv_aspec, /dimen)
  wv = x_calcfit(dindgen(sz[0]),fitstr=all_arcfit[gdo])
  if wv[0] LT 10. then wv = 10^wv
  x_splot, wv, sv_aspec[*,gdo], /bloc

      
  return
end
