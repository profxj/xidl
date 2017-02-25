;+ 
; NAME: 
; igm_lya_monte
;    Version 1.1
;
; PURPOSE:
;    Generate a Monte-Carlo line list of Lya lines given a redshift range 
;    and a range of N_HI
;
; CALLING SEQUENCE:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Feb-2012 Written by JXP
;-
;------------------------------------------------------------------------------
function igm_lya_monte, zmnx, NHI_min, NHI_max, NEVAL=neval, FN_STRCT=fn_strct, NMONTE=nmonte

  if (N_params() LT 2) then begin 
     print,'Syntax - ' + $
           'igm_lya_monte, zmnx, NHI_min, [NHI_max], FN_STRCT=, NMONTE= [v1.0]'
     return, -1
  endif 

  ;; Initialize
  if not keyword_set(NMONTE) then nmonte = 100L
  if not keyword_set(NEVAL) then neval = 10000L
  if not keyword_set(BVAL_MNX) then BVAL_MNX = [10., 80] ; km/s
  if not keyword_set(NHI_MAX) then begin 
     NHI_MAX = 23.
     infinity=1
  endif

  ;; Use HST f(N) as default
  if not keyword_set(FN_STRCT) then $
     fn_strct = xmrdfits(getenv('XIDL_DIR')+'/IGM/fN_empirical/fn_z2.4_5param.fits', 1,/silen)

  ;; Limited usage thus far
  idx = where(fn_strct.npivot GT 0, na)
  if na GT 1 then stop ;; Not ready for this

  if zmnx[0] LT fn_strct[idx].zmnx[0] $
     OR zmnx[1] GT fn_strct[idx].zmnx[1] then stop ;; Not ready for this either

  if fn_strct[idx].zpivot[0] LT fn_strct[idx].zmnx[0] or $ 
     fn_strct[idx].zpivot[0] GT fn_strct[idx].zmnx[1] then stop ;; Not expecting this

  b = where(abs(fn_strct[idx].zpivot - fn_strct[idx].zpivot[0]) GT 1e-3, nb)
  if nb NE 0 then stop

  ;; Get l(X) at z_pivot
  lox = igm_calc_lox(fn_strct, fn_strct[idx].zpivot[0], NHI_min, NHI_max, $
                     CUMUL=cumul, LGNHI=lgnhi)

  ;; Get average # of lines
  gp1 = fn_strct[idx].gamma[0] + 1
  avg_lin = lox * cosm_dxdz(fn_strct[idx].zpivot[0], /W05MAP, /silent) * $
            ( (1+zmnx[1])^gp1 - (1+zmnx[0])^gp1) / $
            (gp1 * (1+fn_strct[idx].zpivot[0])^fn_strct[idx].gamma[0])

  ;; Monte Carlo
  seed = -1244L
  nlines = round(randomn(seed, NMONTE, POISSON=avg_lin))

  ;; Generate N_HI values
  ntot_lines = total(nlines)
  ranx = randomu(seed, ntot_lines)
  all_NHI = interpol(lgNHI, CUMUL/lox, ranx)

  ;; Generate z values
  ranx = randomu(seed, ntot_lines)
  zdum = zmnx[0] + findgen(10001L)*(zmnx[1]-zmnx[0])/10000.
  cumulz = (1+zdum)^gp1 - (1+zmnx[0])^gp1
  all_z = interpol(zdum, CUMULz/max(cumulz), ranx)

  ;; Generate b values (Hui & Rutledge 1999)
  bsig = 24.
  ranx = randomu(seed, ntot_lines)
  bdum = bval_mnx[0] + findgen(10001L)*(bval_mnx[1]-bval_mnx[0])/10000.
  cumulb = total( bdum^(-5) * exp(-1*bsig^4 / bdum^4),  /cumul )
  all_b = interpol(bdum, CUMULb/max(cumulb), ranx)
  
  ;; Fill
  NHI_array = fltarr(NMONTE, max(nlines))
  b_array = fltarr(NMONTE, max(nlines))
  z_array = fltarr(NMONTE, max(nlines))

  i0 = 0
  for qq=0L,NMONTE-1 do begin
     if nlines[qq] EQ 0 then continue
     b_array[qq,0:nlines[qq]-1] = all_b[i0+lindgen(nlines[qq])]
     z_array[qq,0:nlines[qq]-1] = all_z[i0+lindgen(nlines[qq])]
     NHI_array[qq,0:nlines[qq]-1] = all_NHI[i0+lindgen(nlines[qq])]
     i0 = i0+nlines[qq]
  endfor

  strct = { $
          zabs: z_array, $
          NHI: NHI_array, $
          bval: b_array, $
          nlines: nlines $
          }

  return, strct
end
