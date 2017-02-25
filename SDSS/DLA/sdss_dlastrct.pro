;+ 
; NAME:
; sdss_dlastrct
;    Version 1.1
;
; PURPOSE:
;    Produces a DLA structure of SDSS DLAs
;
; CALLING SEQUENCE:
;  sdss_dlastrct, dla
;
; INPUTS:
;  GZSTR= -- Structure summarizing g(z) for the quasars
;  GZFIL= -- Filename of the g(z) file for the quasars
;
; RETURNS:
;  dla -- DLA structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /ACAND  -- Return all absorbers (not just DLA)
;  /ALL   -- Return all statistical DLA from all DRs
;  vprox= -- Proximity velocity (to avoid in calculation)
;  BUFF=  -- Require DLA occur BUFF km/s to the red of z1
;  /DR3   -- Return DLA from DR3 (but not DR2 or DR1 unless indicated)
;  /DR5   -- Return DLA from DR4,DR5 (but not early data releases
;            unless indicated)
;  /OLD_QSOZ -- Use 'original' redshifts instead of the updated ones
;               for the quasar [not recommended]
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

pro sdss_dlastrct, dla, DR1=flg_DR1, DR2=flg_DR2, DR3=flg_dr3, PDR5=pdr5, $
                   ALL=all, STAT=stat, ORIG=orig, GZSTR=gzstr, GZFIL=gzfil, $
                   ADR3=adr3, VPROX=vprox, _EXTRA=extra, ACAND=acand, $
                   GZDR3=gzdr3, DR5=flg_dr5, VMIN=vmin, OLD_QSOZ=old_qsoz


  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_dlastrct, dla, /DR1, /DR2, /DR3, /DR5, /ALL, /STAT, /ORIG, /ACAND, /OLD_QSOZ [v1.1]'
      return
  endif 
  
  if not keyword_set( GZFIL ) then begin
      if keyword_set(ORIG) then  $
        gzfil = getenv('SDSSPATH')+'/DR3_QSO/old_dr3_dgz_s2n4.fits' 
      if keyword_set(GZDR3) then $
        gzfil = getenv('SDSSPATH')+'/DR3_QSO/dr3_dlagz.fits' 
      if not keyword_set( GZFIL ) then $
        gzfil = getenv('SDSSPATH')+'/DR5_QSO/dr5_dlagz_s2n4.fits' 
  endif
  print, 'sdss_dlastrct: Using ', gzfil

  if not keyword_set( DLA5LST ) then dla5lst = getenv('DLA')+'Lists/dr5_alldla.lst'
  if not keyword_set( DLA3LST ) then dla3lst = getenv('DLA')+'Lists/dr3_dla.lst'
  if not keyword_set( DLA2LST ) then dla2lst = getenv('DLA')+'Lists/dr2_dla.lst'
  if not keyword_set( DLA1LST ) then dla1lst = getenv('DLA')+'Lists/dr1_dla.lst'
;  if not keyword_set( PDR5LST ) then pdla5lst = getenv('DLA')+'Lists/dr5_proxdla.lst'

  ;; Read in 
  if keyword_set( ALL ) or keyword_set(ACAND) then begin
      flg_dr1 = 1
      flg_dr2 = 1
      flg_dr3 = 1
      flg_pdr5 = 1
      flg_dr5 = 1
  endif

  if keyword_set( ADR3 ) then begin
      flg_dr1 = 1
      flg_dr2 = 1
      flg_dr3 = 1
  endif

  if keyword_set(flg_dr1) then parse_dlalst, dr1, dla1lst, /noelm, ROOT=getenv('SDSSPATH')
  if keyword_set(flg_dr2) then parse_dlalst, dr2, dla2lst, /noelm, ROOT=getenv('SDSSPATH')
  if keyword_set(flg_dr3) then parse_dlalst, dr3, dla3lst, /noelm, ROOT=getenv('SDSSPATH')
;  if keyword_set(flg_pdr5) then parse_dlalst, pdr5, pdla5lst, /noelm, ROOT=getenv('SDSSPATH')
  if keyword_set(flg_dr5) then parse_dlalst, dr5, dla5lst, /noelm, ROOT=getenv('SDSSPATH')

  tmp = {dlanoestruct}
  dla = [tmp]
  if keyword_set(DR1) then dla  = [dla, dr1]
  if keyword_set(DR2) then dla  = [dla, dr2]
  if keyword_set(DR3) then dla  = [dla, dr3]
; if keyword_set(PDR5) then dla  = [dla, pdr5]
  if keyword_set(DR5) then dla  = [dla, dr5]
  ndla = n_elements(dla)
  dla = dla[1:*]

  if keyword_set(STAT) then begin
      if not keyword_set(GZSTR) then gzstr = xmrdfits(gzfil,1,/silent)
      indx = sdss_dlastat(dla,gzstr,vprox, ALL=ACAND, VMIN=vmin, _EXTRA=extra) 
      dla = dla[indx]
   endif

  ;; Update the redshifts
  new_zem = sdss_newzem(dla.sdss_plate, dla.sdss_fibid)
  newi = where(new_zem GT 0., nnew)
  if nnew GT 0 then dla[newi].qso_zem = new_zem[newi]

  return
end
  
