;+ 
; NAME:
; dla_sdssrich   
;   Version 1.0
;
; PURPOSE:
;    Inputs a list of absorption lines, and the quasar emission
;  redshift and then prints all probably DLA (quality > 0.7)
;
; CALLING SEQUENCE:
;   
;   dla_sdssrich, zem, obslin
;
; INPUTS:
;  zem -- Emission redshift of the quasar
;  obslin -- List of observed lines
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; 
; GRAN= -- List of metal-line transitions to key on
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_sdssrich
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   09-Dec-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_sdssrich, zem, obslin, GTRAN=gtran

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'dla_sdssrich, zem, obslin, GTRAN= [v1.1]'
    return
  endif 

; Optional Keywords

  if not keyword_set( GTRAN ) then begin
      gtran = [1302.1685, $
               1304.3702, $
               1334.5323, $
               1526.7066, $
               1608.4511, $
               1670.7874, $
               2344.214, $
               2382.765, $
               2600.1729, $
               2796.352]
      gtran = alog10(gtran)
  endif
  ntran = n_elements(gtran)
  if not keyword_set( TOLER ) then toler = replicate(1.e-4, ntran)

; LOOP ON INPUT

  nlin = n_elements(obslin)
  for i=0L,nlin-1 do begin
      ztemp = 10^(obslin[i]-gtran) - 1.
      gdz = where(ztemp GT 1.6, ngd)
      if ngd EQ 0 then continue
      ;; LOOP ON z
      for j=gdz[0],ntran-1 do begin
          zabs = ztemp[j]
          lrest = obslin - alog10(1.+zabs)
          tst = where(gtran GT alog10( (1.+zem)*1215.6701/(1.+zabs)) $
                      AND gtran LT alog10(9200./(1.+zabs)), ntst)
          if ntst LT 3 then continue
          ;; LOOP on test
          ngd = 0L
          for k=0L,ntst-1 do begin
              mn = min( abs(gtran[tst[k]] - lrest), imn)
              if mn LT toler[tst[k]] then ngd = ngd + 1
          endfor
          if float(ngd)/float(ntst) GT 0.70 then $
            print, zabs, ngd, ntst
      endfor
  endfor

  return
          
          
end
