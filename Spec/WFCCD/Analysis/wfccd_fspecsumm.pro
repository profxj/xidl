;+ 
; NAME:
; wfccd_fspecsumm
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   wfccd_fspecsumm, wffspec
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
;   wfccd_fspecsumm, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_fspecsumm, fspec_fil, out_fil

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_fspecsumm, fspec_fil, [out_fil] [v1.0]'
    return
  endif 

; Optional Keywords

; Open fspec
  wfccd_wrfspec, wffspec, fspec_fil, /read

; Outfil

  if not keyword_set( OUTFIL ) then outfil = strmid(fspec_fil, 0, 16)+'_summ.txt'
  close, 56
  openw, 56, outfil

; Output

  nfspec = n_elements(wffspec)

  for q=0L,nfspec-1 do begin

      a = where(wffspec[q].var GT 0.)
      mnwv = min(wffspec[q].wave[a], max=mxwv)
      printf, 56, $
        FORMAT='(a12,a5,1x,i2,1x,3f9.5,1x,2f9.2,2f7.1,1x,f6.3,1x,f6.1)', $
        strtrim(wffspec[q].field,2)+': ',$
        strtrim(wffspec[q].slit_id,2)+strtrim(wffspec[q].obj_id,2), $
        wffspec[q].flg_anly, $
        wffspec[q].zans.z, $
        wffspec[q].zans.z_err, $
        wffspec[q].zans.rchi2, $
        mnwv, mxwv, $
        wffspec[q].xyimg, $
        wffspec[q].mag, $
        total(wffspec[q].texp)
  endfor
  close, 56
  return
end
        
