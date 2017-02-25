;+ 
; NAME:
; slls_dx
;  Version 1.1
;
; PURPOSE:
;    Calculate the cosmological path length dX for a given
;  cosmology for the sample of quasars searched for SLLS.
;
; CALLING SEQUENCE:
;   dX = sdss_dx(zmin, zmax, gzstr, XZFIL=)
;
; INPUTS:
;   zmin -- Min redshift of the bin
;   zmax -- Max redshift of the bin
;   gzstr -- Max redshift of the bin
;
; RETURNS:
;   dX  -- The cosmological pathlength 
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2006 Written by JXP
;-
;------------------------------------------------------------------------------
function slls_dx, b1, b2, gzstr, XZFIL=xzfil

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'dX = slls_dx(zmin, zmax, gzstr, XZFIL=) [v1.1]'
    return, -1
  endif 

  ;; Get structure if necessary
  if not keyword_set( XZFIL ) then $
     xzfil=getenv('SDSSPATH')+'DR7_QSO/xz_val_L74_M26.fits.gz' 

  ;; Read gz
  dgzz = gzstr.gzz - shift(gzstr.gzz,1)

  ;; XZ
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 1, /silent)

  ;; Bin
  dX = 0.
  sz = size(gzstr.z1)
  for qq=0L,sz[2]-1 do begin
      for jj=0L,sz[1]-1 do begin
          if gzstr.z1[jj,qq] LE 0. then continue
          mn = min(abs(xzz-(gzstr.z1[jj,qq]>b1)),imn1)
          mn = min(abs(xzz-(gzstr.z2[jj,qq]<b2)),imn2)
          ;; Safety check (superfluous now)
          if imn2 LT imn1 then continue
          ;; Add it in
          dX = dX + xzv[imn2]-xzv[imn1]
      endfor
  endfor

  return, dX
end
