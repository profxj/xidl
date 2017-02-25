;+ 
; NAME:
; lls_stat
;
; PURPOSE:
;    Given a LLS struct and gz list, determine the indices of those
;    LLS satisfying the statistical sample.
;
; CALLING SEQUENCE:
;   sdss_llsstat
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; NHI_CUT =  Criterion for LLS to be 'statistical'.  Default is 17.5
;        to match the tau=2 criterion of the SDSS survey.  Lower to 17.2 for
;        tau=1
; LLS_CUT = Cut search off by presence of an LLS (new with WFC3; not
;          used in original SDSS analysis)
; DZ_TOLER = Tolerance in zabs for pairs of LLS along a sightline (default=0.04)
; /BOOT -- Allows for doubles in the QSO list
;
; OPTIONAL OUTPUTS:
; IDX_BOOT -- Indices for Bootstrap analysis
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Feb-2009 Written by JXP
;-
function sdss_llsstat, llsstr, qsos, PARTIAL=partial, VPROX=vprox, PROX=prox, $
                       OIDX=odix, MOCK=mock, MAXDZ=maxdz, ZEM_MIN=zem_min, $
                       NHI_CUT=nhi_cut, FLG_ZSRCH=flg_zsrch, LLS_CUT=lls_cut, $
                       DZ_TOLER=dz_toler, BOOT=boot, IDX_BOOT=idx_boot

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'indx = lls_stat(llsstr, qsos, /PARTIAL, VPROX=, /PROX, /BOOT), [v1.1]'
    return, -1
  endif 
  if not keyword_set( VPROX ) then vprox = 3000.
  if not keyword_set( MAXDZ ) then maxdz = 99.99
  if not keyword_set( ZEM_MIN ) then zem_min = 0.
  if not keyword_set( NHI_CUT ) then nhi_cut = 17.5
  if not keyword_set( FLG_ZSRCH ) then flg_zsrch = 0
  if not keyword_set( DZ_TOLER ) then dz_toler = 0.04

  if arg_present(idx_boot) then idx_boot = [-1]

  ;; Search redshift
  case flg_zsrch of
     0: zsrch = qsos.zt2
     1: zsrch = qsos.zt1
     2: zsrch = qsos.zt0
     else: stop
  endcase
  
  ;; Modify by LLS along QSO sightline as required
  if size(LLS_CUT,/type) GT 0 then $
     zsrch = zsrch > qsos.zlls[LLS_CUT]

  ;; LLS
  nlls = n_elements(llsstr)
  lls = lindgen(nlls)

  msk_smpl = bytarr(nlls)  ; 0=bad, 1=good

  zmax = x_relvel(qsos.zem,VPROX) 
  if not keyword_set(MOCK) then $
    x_radec, llsstr.qso_ra, llsstr.qso_dec, rad, decd

  for ii=0L,nlls-1 do begin
      qq = lls[ii]

      ;; Two LLS on one sightline?
      if not keyword_set(MOCK) then begin
          mtlls = where(abs(rad-rad[qq]) LT 0.001 AND $
                        abs(decd-decd[qq]) LT 0.001 AND $
                        abs(llsstr.qso_zem-llsstr[qq].qso_zem) LT 0.03 AND $
                        abs(llsstr.zabs-llsstr[qq].zabs) LT DZ_TOLER, nmLLS)
          if nmLLS NE 1 then stop  ;; LLS are probably too close in z
      endif 
      
      ;; Cut on NHI
      if keyword_set(PARTIAL) and llsstr[qq].NHI GT NHI_CUT then continue
      if not keyword_set(PARTIAL) and llsstr[qq].NHI LE NHI_CUT then continue

      ;; Match to QSO RA, DEC
      if not keyword_set(MOCK) then begin
          idx = where(abs(qsos.ra-rad[qq]) LT 0.001 AND $
                      abs(qsos.dec-decd[qq]) LT 0.001 AND $
                      abs(qsos.zem-llsstr[qq].qso_zem) LT 0.03, nidx)
      endif else begin
          idx = where(qsos.plate EQ llsstr[qq].sdss_plate AND $
                      abs(qsos.zem-llsstr[qq].qso_zem) LT 0.03, nidx)
       endelse
      if nidx NE 1 then begin
         if nidx GT 1 and not keyword_set(BOOT) then stop $
         else idx=idx[0]
         if nidx LT 1 then begin
            if keyword_set(BOOT) then continue else stop
         endif
      endif
      
      ;; Query redshift
      if zsrch[idx] GT 0. AND $
        llsstr[qq].zabs GT ((zsrch[idx] > (qsos[idx].zem - MAXDZ)) - 1e-4) AND $
        qsos[idx].zem GT ZEM_MIN then begin 
          if not keyword_set(PROX) and $ ; Intervening
            llsstr[qq].zabs LT zmax[idx] then msk_smpl[qq] = 1B
          if keyword_set(PROX) and $ ; Proximate
            llsstr[qq].zabs GE zmax[idx] then msk_smpl[qq] = 1B
       endif
      ;; Deal with doubles
      if keyword_set(BOOT) and msk_smpl[qq] then begin
         ;; 
         idx_boot = [idx_boot, replicate(lls[qq],nidx)]
      endif
  endfor

  ;; Return
  gd = where(msk_smpl,complement=OIDX,ngd)
  if keyword_set(idx_boot) then idx_boot = idx_boot[1:*]
  if ngd EQ 0 then return, -1 else begin
;      if arg_present(IDXA) then idxa = idxa[gd]
      return, lls[gd]
  endelse

end
