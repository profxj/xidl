;+ 
; NAME:
; sdss_llslozx
;
; PURPOSE:
;    Calculate dn/dz or dn/dX for a set of redshift bins
;
; CALLING SEQUENCE:
;   sdss_llslozx
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; NHI_CUT=  See sdss_lssstat
; LLS_CUT=  Use LLS redshifts (in QSOs structure) to cut search path
; /BOOT -- Allows for doubles in the QSO list
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
;   Jul-2009 Written by JXP
;-
pro sdss_llslozx, qsofil, llsfil, strct, VPROX=vprox, PROX=prox, $
                  XZFIL=xzfil, BINS=bins, MOCK=mock, ZEM_MIN=zem_min, $
                  QSOS=qsos, MAXDZ=maxdz, NHI_CUT=nhi_cut, ALL_LLS=all_lls, $
                  ZLLS_IDX=zlls_idx, FLG_ZSRCH=flg_zsrch, LLS_CUT=lls_cut, $
                  DZ_TOLER=DZ_TOLER, BOOT=BOOT
 
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
          'sdss_llslozx, qsofil, llsfil, strct, VPROX=, /PROX, ' + $
          'XZFIL=, MOCK=mock, /BOOT [v1.1]'
    return
  endif 

  if not keyword_set( vprox ) then vprox = 3000.
  if not keyword_set( XZFIL ) then $
    xzfil = getenv('SDSSPATH')+'/DR7_QSO/xz_val.fits'
  ;; Bins
  if not keyword_set(BINS) then begin
      bins = [ [3.28,3.4], [3.4,3.55], [3.55, 3.7], $
               [3.7, 3.85], [3.85, 4.0], [4.0, 4.5], [4.5,5.0] ]
  endif
  sz = size(bins,/dimensions)
  nbin = sz[1]

  ;; QSOs
  if not keyword_set(QSOS) then qsos = xmrdfits(qsofil, 1, /silent)
  if keyword_set(ZEM_MIN) then begin
      badq = where(qsos.zem LT ZEM_MIN, nbad)
      if nbad NE 0 then qsos[badq].zt2 = 99.99
  endif
  nqso = n_elements(qsos)

  ;; Set a maximum dz to a given sightline (as requested)
;  if keyword_set(MAXDZ) then $
;    qsos.zt2 = qsos.zt2 > (qsos.zem - MAXDZ)

  ;; LLS
  if not keyword_set(ALL_LLS) then all_lls = xmrdfits(llsfil, 1)
  idx = sdss_llsstat(all_lls, qsos, VPROX=vprox, PARTIAL=PLLS, PROX=prox, $
                     MOCK=mock, MAXDZ=maxdz, ZEM_MIN=zem_min, NHI_CUT=nhi_cut, $
                     FLG_ZSRCH=flg_zsrch, LLS_CUT=LLS_CUT, DZ_TOLER=dz_toler, $
                     BOOT=boot, IDX_BOOT=idx_boot)
  if keyword_set(BOOT) then lls = all_lls[idx_boot] $
  else lls = all_lls[idx]
  ;print, 'sdss_llslozx: N_LLS = ', n_elements(idx)


  tmp = { $
        medz: 0., $
        meanzabs: 0., $
        meanbin: 0., $
        mlls: 0L, $
        dz: 0., $
        dX: 0., $
        loz: 0., $
        loX: 0., $
        siglz: fltarr(2), $ 
        siglX: fltarr(2) }
  strct = replicate(tmp, nbin)
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   LOOP on bins

  p_dndz = fltarr(nbin)
  p_dndx = fltarr(nbin)
  p_medz = fltarr(nbin)
  p_meanz = fltarr(nbin)
  p_xplt = fltarr(nbin)
  p_errz = fltarr(nbin)
  p_errnz1 = fltarr(nbin)
  p_errnz2 = fltarr(nbin)
  p_errnx1 = fltarr(nbin)
  p_errnx2 = fltarr(nbin)
  p_svnd = intarr(nbin)

  ;; Calculate l(X)
  for ii=0L,nbin - 1 do begin
      ;; Bin
      gd = where(lls.zabs GE bins[0,ii] AND $
                 lls.zabs LT bins[1,ii], nlls_bin)
      if nlls_bin EQ 0 then continue
      p_medz[ii] = median([lls[gd].zabs])
      p_meanz[ii] = mean([lls[gd].zabs])
      p_xplt[ii] = mean(bins[*,ii])
      p_dX = sdss_llsdx( qsos, bins[0,ii], bins[1,ii], $
                         DZ=p_dz, PROX=prox, VPROX=vprox, MAXDZ=maxdz, $
                         ZEM_MIN=zem_min, ZLLS_IDX=zlls_idx, FLG_ZSRCH=flg_zsrch) 
      
      p_dndz[ii] = nlls_bin/p_dz
      p_dndx[ii] = nlls_bin/p_dX
;      print, 'bins: ', bins[*,ii], nlls_bin
      
      ;; Error
      val = x_poisscl(double(nlls_bin), sigma=1., /silent)
      p_errnz1[ii] = val[0]/p_dz - p_dndz[ii]
      p_errnz2[ii] = p_dndz[ii] - val[1]/p_dz
      p_errnx1[ii] = val[0]/p_dX - p_dndx[ii]
      p_errnx2[ii] = p_dndx[ii] - val[1]/p_dX
      
      p_svnd[ii] = nlls_bin
      p_errz[ii] = p_medz[ii]-bins[0,ii]
      
      ;; Save
      strct[ii].dz = p_dz 
      strct[ii].dX = p_dX 
      strct[ii].mlls = nlls_bin
  endfor
      
  ;; Pass back
  strct.medz = p_medz
  strct.meanzabs = p_meanz
  strct.meanbin = p_xplt
  strct.loz = p_dndz 
  strct.lox = p_dndx 
  strct.siglz[0] = p_errnz1
  strct.siglz[1] = p_errnz2
  strct.siglX[0] = p_errnx1
  strct.siglX[1] = p_errnx2
      
  return
end
      
      
