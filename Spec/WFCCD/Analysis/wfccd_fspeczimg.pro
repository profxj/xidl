;+ 
; NAME:
; wfccd_fspeczimg
;    Version 1.0
;
; PURPOSE:
;   Overplot z,id on xatv image
;
; CALLING SEQUENCE:
;   
;   wfccd_fspeczimg, fspec_fil
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_fspeczimg, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_fspeczimg, fspec_fil, PLTALL=pltall, NOERASE=noerase

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_fspeczimg, fspec_fil [v1.0]'
    return
  endif 

;  Optional Keywords

;  Read in

  wfccd_wrfspec, wffspec, fspec_fil, /read

; Get list
  list = wfccd_getobjnm(wffspec, /LST)

;  Parse

  if keyword_set( PLTALL ) then gdwf = lindgen(n_elements(wffspec_fil)) $
  else gdwf = where(wffspec.flg_anly NE 0 AND wffspec.zans.z_err GT 0 AND $
                   strtrim(wffspec.obj_id,2) EQ 'a')

; Erase
  if not keyword_set( NOERASE ) then xatverase
; Plot
  xatvplot, wffspec[gdwf].xyimg[0]-1., wffspec[gdwf].xyimg[1]-1., psym=1

;  Label

  for q=0L,n_elements(gdwf)-1 do xatvxyouts, wffspec[gdwf[q]].xyimg[0]+2., $
    wffspec[gdwf[q]].xyimg[1]+2., string(wffspec[gdwf[q]].zans.z,FORMAT='(f6.4)'),$
    charsize=2.0, color='red'

  if not keyword_set(NONM) then begin
      for q=0L,n_elements(gdwf)-1 do begin
          xatvxyouts, wffspec[gdwf[q]].xyimg[0]+2., $
            wffspec[gdwf[q]].xyimg[1]+12., string(list[q],FORMAT='(a5)'),$
        charsize=2.0, color='blue'
      endfor
  endif

  return
end

