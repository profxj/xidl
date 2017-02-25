;+ 
; NAME:
; x_obsinit   
;    Version 1.1
;
; PURPOSE:
;  Initialize observing structures for an observatory
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /MAUNAKEA -- Use values for MK
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP
;-
;------------------------------------------------------------------------------
function x_obsinit, infil, FLG=flg, MAUNAKEA=maunakea

  str_obs = {obsstruct}

  flg = 0

  if keyword_set(MAUNAKEA) then flg=0

  case flg of
      0: begin  ; Mauna Kea
          str_obs.seeing = 0.7
          str_obs.airmass = 1.1
          str_obs.mphase = 7
          str_obs.exptime = 3600.
          str_obs.mstar = 17.
          str_obs.mtype = 1.
      end
      else: stop
  endcase

  ;; INFIL?
  if keyword_set(INFIL) then begin
      readcol, infil, card, val, FORMAT='A,A' 
      ;; SEEING
      mtch = where(strtrim(card,2) EQ 'SEEING',nmt)
      if nmt NE 0 then str_obs.seeing = float(val[mtch[0]])
      ;; AIRMASS
      mtch = where(strtrim(card,2) EQ 'AIRMASS',nmt)
      if nmt NE 0 then str_obs.airmass = float(val[mtch[0]])
      ;; MOON
      mtch = where(strtrim(card,2) EQ 'MPHASE',nmt)
      if nmt NE 0 then str_obs.mphase = long(val[mtch[0]])
      ;; EXPOSURE TIME
      mtch = where(strtrim(card,2) EQ 'EXPTIME',nmt)
      if nmt NE 0 then str_obs.exptime = float(val[mtch[0]])
      ;; MAGNITUDE
      mtch = where(strtrim(card,2) EQ 'MAGNITUDE',nmt)
      if nmt NE 0 then str_obs.mstar = float(val[mtch[0]])
      ;; MTYPE
      mtch = where(strtrim(card,2) EQ 'MAGTYPE',nmt)
      if nmt NE 0 then str_obs.mtype = float(val[mtch[0]])
  endif

  return, str_obs
end


