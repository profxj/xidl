;+ 
; NAME:
; sdss_llsstrct
;    Version 1.0
;
; PURPOSE:
;    Produces a LLS structure of SDSS LLS 
;
; CALLING SEQUENCE:
;  sdss_llsstrct, lls
;
; INPUTS:
;
; RETURNS:
;  lls -- LLS structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  QSOS=  -- Structure of statistical QSOs in the sample
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
;   16-Feb-2010 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_llsstrct, lls, QSOS=qsos, INIT=init


  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_llsstrct, lls, [v1.1]'
      return
  endif 
  
  sdss_llsinit, init
  if not keyword_set( LLSFIL ) then llsfil = init.llsfil
  if not keyword_set( QSOFIL ) then qsofil = init.qsofil 
  if not keyword_set( MAXDZ ) then maxdz = init.maxoff 
  if not keyword_set( VPROX ) then vprox = init.vprox
  if not keyword_set( XZFIL ) then xzfil = init.xzfil
  if not keyword_set( ZEM_MIN ) then zem_min = init.zem_min
  if not keyword_set( BINS ) then BINS = init.all_bins
  if not keyword_set( ZABS_MAX ) then zabs_max = init.zabs_max


  ;; Find start/end redshift of search for the quasars
  qsos = xmrdfits(qsofil, 1, /silent)
  zmax = x_relvel(qsos.zem,VPROX) < ZABS_MAX

  all_lls = xmrdfits(llsfil, 1)
  idx = sdss_llsstat(all_lls, qsos, VPROX=vprox, PROX=prox, $
                     MAXDZ=maxdz, ZEM_MIN=zem_min)
  lls = all_lls[idx]
  nlls = n_elements(lls)
  x_radec, all_lls.qso_ra, all_lls.qso_dec, llsra, llsdec

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Adjust the starting redshift for tau>2 LLS
  nalls = n_elements(all_lls)
  for jj=0L,nalls-1 do begin
      mtq = where(abs(qsos.ra-llsra[jj]) LT 1e-4 AND $
                  abs(qsos.dec-llsdec[jj]) LT 1e-4, nmt)
      if nmt NE 1 then stop
      ;; Reset lower search to z_LLS
      if qsos[mtq].zt2 GT 0 and all_lls[jj].NHI GT 17.5 then $
        qsos[mtq].zt2 = qsos[mtq].zt2 > all_lls[jj].zabs
  endfor

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Cut on the 'good' quasars
  goodq = where(qsos.zt2 LT ZMAX AND qsos.zt2 GT 0. AND $
                qsos.flg_qso EQ 0 AND qsos.zem GE ZEM_MIN AND $
                qsos.flg_LLS NE 5,ngq) 

  qsos = qsos[goodq]

  return
end
  
