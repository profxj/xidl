;+ 
; NAME:
; peroux_dx
;  Version 1.1
;
; PURPOSE:
;    Calculate the cosmological path length dX for a given
;  cosmology for the sample of quasars searched for DLA as
;  compiled by Peroux et al. 2003.  Right now, the code defaults to
;  point at all the DR1 values.
;
; CALLING SEQUENCE:
;   
;   dX = peroux_dx(zmin, zmax, /UNIQ, FIL1=, FIL2=, FILRA=, GZFIL=,
;   XZFIL=)
;
; INPUTS:
;   zmin -- Min redshift of the bin
;   zmax -- Max redshift of the bin
;
; RETURNS:
;   dX  -- The cosmological pathlength 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /UNIQ -- Restrict to quasars not found in the SDSS
;  GZFIL -- SDSS gzfil;  Archived g(z) function for the SDSS
;  XZFIL -- SDSS xzfil; contains dX for a given cosmology as a
;           function of redshift
;  FIL1  -- File of DLA listed in Peroux et al. 2003
;  FIL2  -- File of quasars without DLA listed in Peroux et al. 2003
;  FILRA -- File of RA/DEC values for all the quasars compiled by
;           Peroux
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   dx = peroux_dx(2., 2.5)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-March-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function peroux_dx, b1, b2, UNIQ=uniq, FIL1=fil1, FIL2=fil2, FILRA=filra, $
                    GZFIL=gzfil, XZFIL=xzfil, GZSTR=gzstr

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dX = peroux_dx(zmin, zmax, /UNIQ, FIL1=, FIL2=, FILRA=, GZFIL=, XZFIL= [v1.1]'
    return, -1
  endif 
      
  ;; Peroux
  if not keyword_set( FIL1 ) then $
    fil1 = getenv('XIDL_DIR')+'/SDSS/DLA/Peroux/peroux_dlaqso.dat'
  if not keyword_set( FIL2 ) then $
    fil2 = getenv('XIDL_DIR')+'/SDSS/DLA/Peroux/peroux_nodla.dat'
  if not keyword_set( FILRA ) then $
    filra = getenv('XIDL_DIR')+'/SDSS/DLA/Peroux/peroux_radec.dat'
  if not keyword_set( GZFIL ) then $
    gzfil = getenv('SDSSPATH')+'/DR1_QSO/dr1_dlagz.fits'
  if not keyword_set( XZFIL ) then $
    xzfil = getenv('SDSSPATH')+'/DR1_QSO/xz_val.fits'

  ;; Read
  readcol, fil1, a1,a2, zem, z11, z22, format='a,a,f,f,f', /silent
  readcol, fil2, a1,a2, zem, z1, z2, format='a,a,f,f,f', /silent
  z1 = [z11,z1]
  z2 = [z22,z2]
  nqso = n_elements(z1)
  readcol, filra, nam, pra, pdec, format='a,a,a', /silent

  ;; Unique
  if keyword_set( UNIQ ) then begin
      ;; SDSS
      if not keyword_set(GZSTR) then $
        gzstr = xmrdfits(gzfil, 1, /silent)
      ;; Match
      flg_mtch = lonarr(nqso)
      for jj=0L,nqso-1 do begin
          ;; Check ra,dec
          x_radec, strtrim(pra[jj],2), strtrim(pdec[jj],2), rad, decd
          mtch = where(abs(rad-gzstr.ra) LT 0.001 AND $
                       abs(decd-gzstr.dec) LT 0.001, nmtch)
          flg_mtch[jj] = mtch
      endfor
  endif

  ;; XZ
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)

  ;; Bin
  dX = 0.
  gd = where(z2 GE b1 and z1 LE b2, ngd)
  if keyword_set( UNIQ ) then begin
      nmtch = where(flg_mtch[gd] GE 0,nm)
      for kk=0L,nm-1 do begin
          z2[gd[nmtch[kk]]] = $
            gzstr.z1[flg_mtch[gd[nmtch[kk]]]] < z2[gd[nmtch[kk]]] 
      endfor
  endif
;  print, 'nqso = ', ngd
  for ii=0L,ngd-1 do begin
      ;; Calculate dX
      mn = min(abs(xzz-(z1[gd[ii]]>b1)),imn1)
      mn = min(abs(xzz-(z2[gd[ii]]<b2)),imn2)
      if imn2 LT imn1 then continue

      dX = dX + xzv[imn2]-xzv[imn1]
;      printcol, (z1[gd[ii]]>b1), (z2[gd[ii]]<b2), xzv[imn2]-xzv[imn1]
  endfor

;  print, 'total dX = ', dX
  return, dX

end
