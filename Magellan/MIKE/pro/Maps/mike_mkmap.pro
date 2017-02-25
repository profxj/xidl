;+ 
; NAME:
; mike_mkmap   
;     Version 1.0
;
; PURPOSE:
;    Process arc file
;
; CALLING SEQUENCE:
;   
;  mike_mkmap, mike, slit
;
; INPUTS:
;   mike     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mkmap, mike
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_mkmap, ordr_fil, map, CHK=chk

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_mkmap, ordr_fil, [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set(SZ_MAP) then sz_map = [1014L, 2047L]
  if not keyword_set(OUTFIL) then outfil = 'Flats/Map_01.fits'

; Map
  map = dblarr(sz_map[0], sz_map[1]) - 99999.

; Read order structure
  ordr_str = xmrdfits(ordr_fil, 1, /silent) 
  nordr = n_elements(ordr_str)

; Loop on Orders

; Main Loop
  for q=0L,nordr-1 do begin
      ;; Smash row
      yval = ordr_str[q].ycen
      ;; Edges of order
      ltrc = ordr_str[q].lhedg
      rtrc = ordr_str[q].rhedg

      ;; Separation
      lhs = (long(ltrc) - 1) > 0L
      rhs = (long(rtrc) + 2) < (sz_map[0]-1)
      sep = rhs - lhs  ; Check this

      ;; Just the good ones
      a = where(rhs GE 0L AND lhs LE (sz_map[0]-1) AND $
                rhs GE lhs, na)
      ;; Calculate
      for qq=0L,na-1 do begin
          j = a[qq]
          frac = (findgen(rhs[j]-lhs[j]+1) + lhs[j] - ltrc[j]) / sep[j]
          map[lhs[j]:rhs[j],j] = (findgen(rhs[j]-lhs[j]+1) + lhs[j] - $
            (ltrc[yval]+frac*sep[yval]))*(-1.)
      endfor
  endfor

; CHK
  if keyword_set( CHK ) then xatv, map, /block
  
  
; Output
  print, 'mike_echmkmap: Map is in ', outfil
  mwrfits, map, outfil, /create
  ;; compress
  spawn, 'gzip -f ', outfil

  return
end

