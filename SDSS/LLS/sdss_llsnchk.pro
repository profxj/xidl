;+ 
; NAME:
; sdss_llsnchk
;
; PURPOSE:
;    Launches x_specplot to fiddle (and set) the redshift of a LLS in
;    SDSS given the LLS structure.
;
; CALLING SEQUENCE:
;   sdss_llszchk
;
; INPUTS:
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   30-May-2007 Modified by JXP
;-
;------------------------------------------------------------------------------
pro sdss_llszchk, lls_fil, new_fil

common x_specplot_lines, $
  flg_lines, $
  lines, $
  zabs

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_llszchk, lls_fil, new_fil [v1.0] '
    return
  endif 

  lls = xmrdfits(lls_fil,1)
  zchk = where(lls.flg_tweak[0] GT 0, nchk) 

  if nchk EQ 0 then begin
      print, 'sdss_llszchk: No LLS to check'
      return
  endif

  for qq=0L,nchk-1 do begin
      ;; Plot
      sdss_objinf, [ lls[zchk[qq]].plate, lls[zchk[qq]].fiber ], $
                   /plot, /lls, zin=lls[zchk[qq]].zlls
      ;; Update
      print, 'sdss_llszchk: Setting zlls = ', zabs
      lls[zchk[qq]].zlls = zabs
  endfor
  
  ;; Write
  mwrfits,lls,new_fil,/create
  spawn, 'gzip -f '+new_fil

end

