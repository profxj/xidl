;+ 
; NAME:
; sdss_llsdx
;  Version 1.1
;
; PURPOSE:
;    Calculate the cosmological path length dX for a given
;  cosmology for the sample of quasars searched for LLS in the SDSS.
;  This code can also return dz, the redshift path length.
;
; CALLING SEQUENCE:
;   dX = sdss_llsdx(zmin, zmax, GZFIL=, XZFIL=)
;
; INPUTS:
;   zmin -- Min redshift of the bin
;   zmax -- Max redshift of the bin
;  [vprox] -- Proximity velocity (to avoid in calculation)
;  GZSTR= -- Structure summarizing g(z) for the quasars
;  GZFIL= -- Filename of the g(z) file for the quasars
;  XZFIL= -- Filename of the dX calculation as a function of z.
;
; RETURNS:
;   dX  -- The cosmological pathlength 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BUFF=  -- Require DLA occur BUFF km/s to the red of z1
;  ZLLS_IDX= -- Sets the index used for modifying the search path.
;              Default=0
;
; OPTIONAL OUTPUTS:
;  dz=  -- Redshift pathlength instead of dX
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   writecol, 'arrays.dat', array1, array2
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Feb-2009 Written by JXP
;-
;------------------------------------------------------------------------------
function sdss_llsdx, iqsos, b1, b2, VPROX=vprox, XZFIL=xzfil, $
                     VMAX=vmax, BUFF=buff, DZ=dz, PARTIAL=partial, $
                     PROX=prox, MAXDZ=maxdz, ZEM_MIN=zem_min, $
                     ZLLS_IDX=zlls_idx, FLG_ZSRCH=flg_zsrch

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
          'dX = sdss_llsdx(qsos, zmin, zmax, lls, ' + $
          '[vprox], GZFIL=, XZFIL=) [v1.1]'
    return, -1
  endif 

  if not keyword_set(VPROX) then vprox = 3000.  ; km/s
  if not keyword_set(BUFF) then buff = 3000.
  if not keyword_set(MAXDZ) then maxdz = 99.99
  if not keyword_set( ZEM_MIN ) then zem_min = 0.
  if not keyword_set( ZLLS_IDX ) then zlls_idx = 0L
  if not keyword_set( FLG_ZSRCH ) then flg_zsrch = 0L

  ;; Get structure if necessary
;  if not keyword_set( GZFIL ) then $
;    gzfil = getenv('SDSSPATH')+'/DR1_QSO/dr1_dlagz.fits'
  if not keyword_set( XZFIL ) then $
    xzfil = getenv('SDSSPATH')+'/DR7_QSO/xz_val.fits'

  ;; XZ
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)

  ;; Cut the QSOs
  gdqso = where(iqsos.zem GT ZEM_MIN)
  qsos = iqsos[gdqso]

  ;; Search redshift
  case flg_zsrch of
     0: zsrch = qsos.zt2
     1: zsrch = qsos.zt1
     2: zsrch = qsos.zt0
     else: stop
  endcase

  ;; v/c
  zmax = x_relvel(qsos.zem, vprox)

  ;; tau>2 LLS terminates the search
  zmin = (zsrch > (qsos.zem-MAXDZ)) > qsos.zlls[ZLLS_IDX]  

  ;; Proximates?
  if keyword_set(PROX) then begin
      zmin = zsrch > zmax
      zmax = qsos.zem
  endif

  ;; 
  good = where(zmin LT zmax AND qsos.flg_qso EQ 0 $
               and zsrch GT 0. $
               and zsrch LT zmax, ngd)
  if ngd EQ 0 then return, 0.

  ;; Bin
  dX = 0.
  dz = 0.

  for jj=0L,ngd-1 do begin
      ii = good[jj]
      if zmin[ii] LE 0. then continue
      ;; Calculate dX
      mn = min(abs(xzz-(zmin[ii]>b1)),imn1)  ;; Discrete not exact
      mn = min(abs(xzz-(zmax[ii]<b2)),imn2)
      if imn1 GE imn2 then continue
      ;; Add it in
      dX = dX + xzv[imn2]-xzv[imn1]
      dz = dz + xzz[imn2]-xzz[imn1]
  endfor

  return, dX
end
