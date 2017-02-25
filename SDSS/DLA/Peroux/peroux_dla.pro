;+ 
; NAME:
; peroux_dla
;  Version 1.1
;
; PURPOSE:
;    Returns the Peroux compilation of DLA
;
; CALLING SEQUENCE:
;   
;   dla = peroux_dla([zmnx])
;
; INPUTS:
;   [zmnx] -- Min/max redshift of the bin
;
; RETURNS:
;   dla  -- Structure containgin the DLA info
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
;   15-March-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function peroux_dla, zmnx, ALLDR=alldr, FIL1=fil1, FIL2=fil2, FILRA=filra, $
                    GZFIL=gzfil, XZFIL=xzfil, GZSTR=gzstr

;  if (N_params() LT 2) then begin 
;    print,'Syntax - ' + $
;             'dla = peroux_dla([zmnx], ALLDR=) [v1.1]'
;    return, -1
;  endif 
      
  ;; Peroux
  if not keyword_set( FIL1 ) then $
    fil1 = getenv('XIDL_DIR')+'/SDSS/DLA/Peroux/peroux_upddla.dat'
  if not keyword_set( FILRA ) then $
    filra = getenv('XIDL_DIR')+'/SDSS/DLA/Peroux/peroux_radec.dat'

  ;; Read
  readcol, fil1, a1,a2, zem, z11, z22, pzabs, pNHI, $
    format='a,a,f,f,f,f,f', /silent
  readcol, filra, nam, pra, pdec, format='a,a,a', /silent
  gd = where(pNHI GE 20.3,nprx)

  ;; Create the structure
  dlas = replicate({dlanoestruct}, nprx)

  ;; Fill up
  dlas.zabs = pzabs[gd]
  dlas.NHI = pnhi[gd]
  pra = pra[gd]
  pdec = pdec[gd]

  ;; SDSS ??
  if keyword_set(ALLDR) then begin
      ;; Setup the RA/DEC
      ndr = n_elements(alldr)
      drra = dblarr(ndr)
      drdec = dblarr(ndr)
      for jj=0L,ndr-1 do begin
          x_radec, alldr[jj].qso_ra, alldr[jj].qso_dec, radd, decdd
          drra[jj] = radd
          drdec[jj] = decdd
      endfor
      ;; Match?
      pmsk = bytarr(nprx)
      for jj=0L,nprx-1 do begin
          x_radec, strtrim(pra[jj],2), strtrim(pdec[jj],2), rad, decd
          ;; Check ra,dec
          mtch = where(abs(rad-drra) LT 0.001 AND $
                       abs(decd-drdec) LT 0.001 AND $
                       abs(pzabs[jj]-alldr.zabs) LT 0.005, nmtch)
          if nmtch NE 0 then pmsk[jj] = 1B
      endfor
      ;; Parse
      gd = where(not pmsk)
      dlas = dlas[gd]
  endif

  ;; Redshift?
  if keyword_set(zmnx) then begin
      gd = where(dlas.zabs GE zmnx[0] AND $
                 dlas.zabs LT zmnx[1])
      dlas = dlas[gd]
  endif

  return, dlas

end
