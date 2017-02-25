;+ 
; NAME:
; sdss_repeats   
;    Version 1.1
;
; PURPOSE:
;  Identify repeats in DLA/MSDLA candidates in the the checked file (default: DR2)
;
; CALLING SEQUENCE:
; sdss_repeats, fil
;
; INPUTS:
;   fil  -- Name of checked file
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
;   Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_repeats, fil

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_repeats, fil [v1.1]'
      return
  endif 

;if not keyword_set(fil) then $
;  fil = '/data/SDSS/DR2_QSO/DR2_chkd.fits'

;read in chkd file
dlastr = xmrdfits(fil, 1, /silent)
tot = n_elements(dlastr)

;loop

flg_NHI = intarr(1)
flg_mtl = intarr(1)
elmt = intarr(1)
svnm = strarr(1)
svz = fltarr(1)

for i=0, tot-1 do begin

    same=where(dlastr[i].qso_name eq dlastr.qso_name and $
               abs(dlastr[i].zabs - dlastr.zabs) lt .01, nsame)  
    
    if nsame gt 1 then begin

        flg_NHI = [flg_NHI, dlastr[i].flg_NHI]
        flg_mtl = [flg_mtl, dlastr[i].flg_mtl]
        elmt = [elmt, i]
        svnm = [svnm, dlastr[i].qso_name]
        svz = [svz, dlastr[i].zabs]
    endif

endfor

srt = sort(svnm)
printcol, elmt[srt], svnm[srt], svz[srt], flg_NHI[srt], flg_mtl[srt] 

;stop
end
