;+ 
; NAME:
; wfccd_fspecxypix
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   indx = wfccd_fspecxypix(wffspec, obj_nm)
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
;   wfccd_fspecxypix, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_fspecxypix, fspec_fil, obj_list

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_fspecxypix, fspec_fil, obj_list [v1.0]'
    return
  endif 

; Optional Keywords

; Open fspec
  wfccd_wrfspec, wffspec, fspec_fil, /read

; Parse obj_list

  readcol, obj_list, id, xpix, ypix, Bmag, Bsig, Rmag, Rsig, $
    FORMAT='(L,F,F,F,F,F,F)'

; Loop

  nobj = n_elements(wffspec)
  for q=0L,nobj-1 do begin
      idx = where(wffspec[q].slit_id EQ id, nidx)
      if nidx NE 1 then message, 'wfccd_fspecxypix: Uh oh'
      wffspec[q].xyimg[0] = xpix[idx[0]]
      wffspec[q].xyimg[1] = ypix[idx[0]]
      wffspec[q].mag = Rmag[idx[0]]
  endfor

; Write fspec
  wfccd_wrfspec, wffspec, fspec_fil

  return
end
