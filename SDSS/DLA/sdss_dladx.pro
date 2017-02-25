;+ 
; NAME:
; sdss_dladx
;  Version 1.1
;
; PURPOSE:
;    Calculate the cosmological path length dX for a given
;  cosmology for the sample of quasars searched for DLA in the SDSS.
;  For now, the code defaults to point at the dX vector calcualted
;  when we did the DR1 search (PH04).
; 
;  This code can also return dz, the redshift path length.
;
; CALLING SEQUENCE:
;   dX = sdss_dladx(zmin, zmax, GZFIL=, XZFIL=)
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
;   15-March-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function sdss_dladx, b1, b2, vprox, GZFIL=gzfil, XZFIL=xzfil, GZSTR=gzstr, $
                     VMAX=vmax, BUFF=buff, DZ=dz, VMIN=vmin

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dX = sdss_dladx(zmin, zmax, [vprox], GZFIL=, XZFIL=) [v1.1]'
    return, -1
  endif 

  if not keyword_set(BUFF) then buff = 3000.

  ;; Get structure if necessary
  if not keyword_set( GZFIL ) and not keyword_set(GZSTR) then stop
;    gzfil = getenv('SDSSPATH')+'/DR1_QSO/dr1_dlagz.fits'
  if not keyword_set( XZFIL ) then $
    xzfil = getenv('SDSSPATH')+'/DR7_QSO/xz_val_L74_M26.fits'


  ;; Read gz
  if not keyword_set( GZSTR ) then $
    gzstr = xmrdfits(gzfil, 1, /silent)
  dgzz = gzstr.gzz - shift(gzstr.gzz,1)

  ;; XZ
  print, 'sdss_dladx: Using ', xzfil
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)

  ;; v/c
  if keyword_set(vprox) then begin
      c = x_constants()
      vb = (vprox+buff) / (c.c / 1e5)
      vc = vprox / (c.c / 1e5)
      ovc = 1. - vc
      ovb = 1. - vb

      zbuff = x_relvel(gzstr.zem, (vprox+buff))
      good = gzstr.z1 LE zbuff
      zmin = x_relvel(gzstr.zem, vprox)
      
      if keyword_set(VMAX) then begin
          zmax = x_relvel(gzstr.zem,vmax) 
      endif else zmax = gzstr.zem
  endif else begin
      if keyword_set(VMIN) then zmin = x_relvel(gzstr.z1, vmin) $
        else zmin = gzstr.z1
      zmax = gzstr.z2
      good = gzstr.z2 GT gzstr.z1
  endelse

  ;; Bin
  dX = 0.
  dz = 0.
  gd = where(zmax GE b1 and zmin LE b2 AND gzstr.flg_bal NE 2 $
             and good, ngd)

  for ii=0L,ngd-1 do begin
      ;; Calculate dX
      mn = min(abs(xzz-(zmin[gd[ii]]>b1)),imn1)
      mn = min(abs(xzz-(zmax[gd[ii]]<b2)),imn2)
      ;; Safety check (superfluous now)
      if imn2 LT imn1 then continue
      ;; Add it in
      dX = dX + xzv[imn2]-xzv[imn1]
      dz = dz + xzz[imn2]-xzz[imn1]
  endfor

  return, dX
end
