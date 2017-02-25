;+ 
; NAME:
; hires_mktempl
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  hires_mktempl, fitfil, ordrs, outfil
;
; INPUTS:
;  fitfil - IDL save file output by hires_fitarc
;  ordrs  - Orders to be archived
;  outfil - Archive file
;
; RETURNS:
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
;   22-Aug-2005 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_mktempl, fitfil, ordrs, outfil

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_mktempl, fitfil, ordrs, outfil [v1.1]'
      return
  endif 
  
  ;;  Optional Keywords

  ;; Fitfil
  restore, fitfil

  ;; Good
  gd = where(guess_ordr GE ordrs[0] and $
             guess_ordr LE ordrs[1], ngd)

  ;; Make
  if ngd EQ 0 then return
  guess_ordr = guess_ordr[gd]
  rejstr[0:ngd-1] = rejstr[gd]
  all_arcfit[0:ngd-1] = all_arcfit[gd]
  sv_lines[0:ngd-1] = sv_lines[gd]
  for ss=0L,ngd-1 do begin
      sv_aspec[*,ss] = sv_aspec[*,gd[ss]]
  endfor

  ;; Save
  save, guess_ordr, rejstr, all_arcfit, sv_lines, sv_aspec, $
    filename=outfil

  return

end

