;+ 
; NAME:
; sdss_uniqdla
;    Version 1.1
;
; PURPOSE:
;    Finds those entires unique in data set 2
;
; CALLING SEQUENCE:
;  sdss_uniqdla, data1, data2, uni
;
; INPUTS:
;  data1 -- DLA structure
;  data2 -- DLA structure
;
; RETURNS:
;  uni -- Unique index in data2
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
;   sdss_uniqdla, dr1, dr2, uni
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Nov-2004 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_uniqdla, d1, d2, uni, D2RA=d2ra

  if  N_params() LT 3 then begin 
      print,'Syntax - ' + $
        'sdss_uniqdla, d1, d2, uni, [v1.1]'
      return
  endif 
  
  ;; Get RA, DEC
  nd2 = n_elements(d2)
  if not keyword_set( D2RA ) then begin
      d2ra = dblarr(nd2)
      d2dec = dblarr(nd2)
      for jj=0L,nd2-1 do begin
          x_radec, d2[jj].qso_ra, d2[jj].qso_dec, radd, decdd
          d2ra[jj] = radd
          d2dec[jj] = decdd
      endfor
  endif
  

  ;; Loop on d1 (presumed smaller)
  msk = bytarr(nd2)

  for qq=0L,n_elements(dr1)-1 do begin
      x_radec, d1[qq].qso_ra, d1[qq].qso_dec, radd, decdd
      ;; Search for a match
      mtch = where(abs(d2ra-radd) LT 0.001 AND $
                   abs(d2dec-decdd) LT 0.001 AND $
                   abs(d1[qq].zabs-d2.zabs) LT 0.01, nm)
      if nm NE 0 then msk[mtch] = 1B
  endfor
  uni = where(msk NE 1B)

  return
end
