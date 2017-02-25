;+ 
; NAME:
; esi_echmkmap   
;     Version 1.0
;
; PURPOSE:
;    Create a Map for rectification of the Arc files
;
; CALLING SEQUENCE:
;   
;  esi_echmkmap, esi, slit
;
; INPUTS:
;   esi     -  ESI structure
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
;   esi_echmkmap, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmkmap, esi

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echmkmap, esi, [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set(HFIT) then hfit = 'Maps/hole_fit.idl'
  if not keyword_set(SZ_MAP) then sz_map = [2048L, 4096L]
  if not keyword_set(OUTFIL) then outfil = 'Maps/ECH_map.fits'
  if not keyword_set(YVAL) then yval = 2048L

; Map
  map = dblarr(sz_map[0], sz_map[1])

; Restore traces
  restore, hfit
  nhtrc = n_elements(fin_fit)

; Main Loop
  for q=0L,nhtrc-2 do begin
      ltrc = x_calcfit(findgen(sz_map[1]), FITSTR=fin_fit[q])
      rtrc = x_calcfit(findgen(sz_map[1]), FITSTR=fin_fit[q+1])
      sep = rtrc - ltrc
      lhs = (long(ltrc) + 1) > 0L
      rhs = long(rtrc) < (sz_map[0]-1)
      ;; Just the good ones
      a = where(rhs GE 0L AND lhs LE (sz_map[0]-1) AND $
                rhs GE lhs, na)
      for qq=0L,na-1 do begin
          j = a[qq]
          frac = (findgen(rhs[j]-lhs[j]+1) + lhs[j] - ltrc[j]) / sep[j]
          map[lhs[j]:rhs[j],j] = (findgen(rhs[j]-lhs[j]+1) + lhs[j] - $
            (ltrc[yval]+frac*sep[yval]))*(-1.)
      endfor
  endfor
  
; Output
  print, 'esi_echmkmap: Map is in ', outfil
  mwrfits, map, outfil, /create
  spawn, 'gzip -f '+outfil

  return
end

