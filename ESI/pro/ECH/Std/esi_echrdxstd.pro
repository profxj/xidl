;+ 
; NAME:
; esi_echrdxstd   
;     Version 1.0
;
; PURPOSE:
;    Trace a standard star in each ordrer
;
; CALLING SEQUENCE:
;   
;  esi_echrdxstd, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with rdxstdered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echrdxstd, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Dec-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echrdxstd, esi, indx, guide, CHK=chk, CLOBBER=clobber,$
                   SCICLM=sciclm

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'esi_echtrcstd, esi, indx, guide [v1.0]'
      return
  endif 


  for q=0L, n_elements(indx)-1 do begin
      ii = indx[q]

; Set cards
      esi[ii].flat_fil = esi[guide].flat_fil
      esi[ii].arc_fil = esi[guide].arc_fil

; PROC
      esi_echproc, esi, [ii]

; FNDOBJ
      esi_echfndobj, esi, [ii], /std, SCICLM=sciclm, CHK=chk, CLOBBER=clobber
      
; SKYSUB
      esi_echskysub, esi, [ii], /std
      
; EXTRACT
      esi_echextobj, esi, [ii], /std
  endfor

end
              
      
      
