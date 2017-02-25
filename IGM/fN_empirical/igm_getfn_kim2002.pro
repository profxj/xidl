;+ 
; NAME: 
; igm_getfn_kim2002
;    Version 1.1
;
; PURPOSE:
;    Create f(N,X) measurements from Kim et al. 2002 line lists
;
; CALLING SEQUENCE:
;
; INPUTS:
;  FLG_COSM=  -- Flag establishing the comoslogy for the calculations
;
; RETURNS:
;
; OUTPUTS:
;  FITS structures containing the empirical constraints
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
;   July-2011 Written by JXP
;-
;------------------------------------------------------------------------------
;; Following Neil Creighton's analysis of Kim et al. 2002
;; Test against Neil's numbers
;;  igm_getfn_kim2002, [2.1, 2.5], [[13.0, 14.], [13.,14.]], f_strct,
;;  flg_cosm=1, /NO_CULL
;; 
pro igm_getfn_kim2002, z_mnx, NHI_mnx, f_strct, FLG_COSM=flg_cosm, $
                       NO_CULL=no_cull

  ;; Set Cosmology
  if not keyword_set( FLG_COSM ) then flg_cosm = 0L

  case flg_cosm of 
     0: begin
        cosmc = 'WMAP5'
        cosm_common, /w05map ;; Dunkley et al. 2009
     end
     1: begin
        cosmc = 'VANILLA'
        cosm_common, /vanilla
     end
     2: begin
        cosmc = 'Neil_Test'
        cosm_common, Omegavac=0.73, OmegaDM=0.27, H0=73.
     end
     else: stop
  endcase

  if not keyword_set(z_mnx) then  z_mnx = [3.5, 3.9]
  if not keyword_set(NHI_MNX) then $
     NHI_mnx = [ [12.5, 13.], [13., 13.5], [13.5, 14.], [14., 14.5], $
                 [14.5, 15.], [15., 15.5]]

  ;; Kim 2002 info
  kim02_qso = ['HE0515-4414','Q1101-264','J2233-606','HE1122-1648', $
              'HE2217-2818','HE1347-2457','Q0302-003','Q0055-269','Q0000-26']
  nqso = n_elements(kim02_qso)
  ;; note 2nd QSO here has a gap 3400-3500 (dla)
  ;; Also Kim 2002 table says 3510-4100, but I'm only using
  ;; 3510-4050, it looks like some lines are missed from 4050-4100.
  kim02_wvlim = [[3080.d,3270], [3230,3778], [3400,3890], [3500,4091], $
                 [3510,4050], [3760,4335], [4808,5150], [4852,5598], [5450,6100]]

  ;; Read in Full linelist
  lines = xmrdfits('Data/z1p5_4.fits',1)
  idx_allkim_HI = where(strtrim(lines.trans,2) EQ 'H I' AND $
                        ( strtrim(lines.ref,2) EQ 'Kim01' OR $
                          strtrim(lines.ref,2) EQ 'Kim02' OR $
                          strtrim(lines.ref,2) EQ 'Lu96'), nall_kimHI)

  allkim02_HI = lines[idx_allkim_HI]
  msk = replicate(1B,nall_kimHI)

  ;; ;;;;;;;;;
  ;; Cull Kim2002
  bad_Q0000 = where(strtrim(allkim02_HI.qso,2) EQ 'Q0000-26' AND $
                    allkim02_HI.wave LT 5450. or $
                    allkim02_HI.wave GT 6100.)
  msk[bad_Q0000] = 0B

  if not keyword_set(NO_CULL) then begin
     bad_bval = where(allkim02_HI.bsig GT (0.25 * allkim02_HI.bval))
     bad_nval = where(allkim02_HI.sigN GT 0.2)

     msk[[bad_Q0000,bad_bval,bad_nval]] = 0B
  endif

  ;; Analysis regions (from wvlim) -- This should be redundant (it is)
;  all_badl = [-1L]
;  for qq=0L,nqso-1 do begin
;     badl = where(strtrim(allkim02_HI.qso,2) EQ kim02_qso[qq] AND $ 
;                  ( allkim02_HI.wave LT kim02_wvlim[0,qq] OR $
;                    allkim02_HI.wave GT kim02_wvlim[1,qq]), nbadl)
;     if nbadl GT 0 then all_badl = [all_badl, badl]
;  endfor
;  bmsk = where(all_badl GT -1L, nbmsk)
;  if nbmsk GT 0 then msk[all_badl[bmsk]] = 0B

  ;; Redshift
  bad_z = where(allkim02_HI.zabs LT z_mnx[0] OR allkim02_HI.zabs GT z_mnx[1])
  msk[bad_z] = 0B

  kim02_HI = allkim02_HI[where(msk)]

  ;; ;;;;;;;
  ;;  Begin f(N,X) analysis
  sz_NHI = size(NHI_mnx, /dimen)

  ;; Calculate DX
  DX = 0.d
  nsight = 0L
  for qq=0L,nqso-1 do begin
     ;; z_lim
     zmin = (kim02_wvlim[0,qq]/1215.6701 -1.) > z_mnx[0]
     zmax = (kim02_wvlim[1,qq]/1215.6701 -1.) < z_mnx[1]
     if zmax LE zmin then continue  ;; No path
     DX_qq = cosm_xz(zmax, /noinit) - cosm_xz(zmin, /noinit)
     DX = DX + DX_qq
     nsight = nsight + 1
     ;; Q1101-264 has a DLA
     if strmatch(strtrim(kim02_qso[qq],2),'Q1101-264') then begin
        z_dla = [3400, 3500]/1215.6701 -1
        if z_dla[0] LT zmax and z_dla[0] GT zmin then begin
           sub_DX = cosm_xz( (zmax < z_dla[1]), /noinit ) - cosm_xz(z_dla[0], /noinit)
           DX = DX - sub_DX
        endif else if z_dla[1] LT zmax and z_dla[1] GT zmin then begin
           sub_DX = cosm_xz( z_dla[1], /noinit ) - cosm_xz(zmin, /noinit)
        endif
     endif
  endfor
           
  print, 'igm_getfn_kim2002: DX = ', DX, nsight

  ;; Loop on NHI bins
  fNX = fltarr(sz_NHI[1])
  sig_fNX = fltarr(sz_NHI[1],2)
  DN = fltarr(sz_NHI[1])
  allz = [0.]

  for qq=0L,sz_NHI[1]-1 do begin
     ;; Just count lines
     gdl = where(kim02_HI.logN GE NHI_mnx[0,qq] AND $
                 kim02_HI.logN LT NHI_mnx[1,qq], nline)
     if nline GT 0 then allz = [allz, kim02_HI.zabs]
     dN[qq] = 10.d^(NHI_mnx[1,qq]) - 10.d^(NHI_mnx[0,qq])
     fNX[qq] = nline / dN[qq] / DX
     sig_fNX[qq,*] = x_poisscl(nline,sigma=1) / dN[qq] / DX
     ;; Log
     fNX[qq] = alog10(fNX[qq])
     sig_fNX[qq,*] = abs( alog10(sig_fNX[qq,*]) - fNX[qq])
     print, nline, NHI_mnx[*,qq], fNX[qq], sig_fNX[qq,*]
  endfor

  print, 'igm_getfn_kim2002: <z> = ', mean(allz[1:*])

  f_strct = { $
            fNX: fNX, $
            sig_fNX: sig_fNX, $
            DX: DX, $
            DN: DN  $
            }
                 
  return
end

  
