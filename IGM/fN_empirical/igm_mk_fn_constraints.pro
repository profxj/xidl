;+ 
; NAME: 
; igm_mk_fn_constraints   
;    Version 1.1
;
; PURPOSE:
;    Generate a FITS files containing the constraints on f(N,X)
;    for the IGM for a set of redshifts
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro igm_mk_fn_constraints_z25, cosmc, LYA_FLG=lya_flg, OUTFIL=outfil

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  z=2.5
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  
  ;; Generate FITS file on z=2.5 constraints for f(N)
  ;;
  ;; EXTEN=1 is the fN constraints (Lya forest, SLLS, DLAs)
  ;;    Lya:  Kim et al. 2001 or Kim et al. 2002
  ;;    SLLS: OMeara et al. 2006, 
  ;;    DLA:  PW09
  ;; EXTEN=2 is the LLS constraints (ACS+WFC3+ ??)  [OPW13]
  ;; EXTEN=3 is the MFP constraints (WFC3)  [OPW13]
  ;; EXTEN=4 is a teff^Lya (D_A) constraint [K05]
  ;; EXTEN=5 is a beta constraint [Kim et al. 2001]

  if n_elements(LYA_FLG) EQ 0 then lya_flg = 1L  ;; Kim et al. 2002

  if not keyword_set(OUTFIL) then begin
     case lya_flg of
        0: outfil = 'fn_constraints_z2.5_beta.fits'
        1: outfil = 'fn_constraints_z2.5.fits'
        2: outfil = 'fn_constraints_z2.5_nolyafN.fits'
        else: stop
     endcase
  endif

  ;; XZ
  case cosmc of 
     'WMAP5': begin
        xzfil=getenv('SDSSPATH')+'DR7_QSO/xz_val_L74_M26.fits.gz' 
        flg_cosm = 0
     end
     'VANILLA': begin
        xzfil=getenv('SDSSPATH')+'DR7_QSO/xz_val_L70_M30.fits.gz' 
        flg_cosm = 1
     end
     else: stop
  endcase

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lya Forest
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lya_fn = {fnconstraint}

  case LYA_FLG of
     0: begin ;; Kim et al. 2001
        lya_fn.zeval = 2.13  
        lya_fn.ref = 'K01'
        lya_fn.type = 'Lya Forest'
        lya_fn.comment = 'Read from Figure 4'
        lya_fn.cosm = cosmc ;; Forced
        lya_fn.DX = 3. ;; Estimated
        
        lya_fn.npt = 1
        
        kim_NHI = [13.0, 13.85]
        kim_lfN  = [-11.4, -12.6]
        for ii=0L,lya_fn.npt-1 do begin

           lya_fn.NHI[ii] = kim_NHI[ii]
           N1 = kim_NHI[ii] - 0.2
           N2 = kim_NHI[ii] + 0.2
           dN = 10.d^N2 - 10.d^N1
           
           lya_fn.dN[ii] = dN
           lya_fn.bins[ii,0] = N1
           lya_fn.bins[ii,1] = N2
           
           kim_fN = 10.d^(kim_lfN[ii])
           dXdz_q0 = (1+lya_fn.zeval)
           dXdz_cosm = cosm_dxdz(lya_fn.zeval, /noinit) ;; This is very close to q=0, actually
           new_kim_fN = kim_fN / dXdz_q0 * dXdz_cosm
           
           lya_fn.fN[ii] = alog10(new_kim_fN)
           lya_fn.sig_fN[ii,*] = 0.1 ;; By-eye estimate
        endfor
     end 
     1: begin ;; Kim et al. 2002
        lya_fn.zeval = 2.34  ;; <z> for the lines in the interval  
        lya_fn.ref = 'K02'
        lya_fn.type = 'Lya Forest'
        lya_fn.comment = 'Recalculated for our Cosmology'
        lya_fn.cosm = cosmc

        NHI_mnx = [ [12.5, 13.], [13., 13.5], [13.5, 14.], [14., 14.5]]
        sz_NHI = size(NHI_mnx, /dimen)

        ;; Calculate on-the-fly
        igm_getfn_kim2002, [2.2, 2.6], NHI_mnx, f_STRCT, FLG_COSM=flg_cosm
        lya_fn.DX = f_strct.DX
        
        lya_fn.npt = sz_NHI[1]
        lya_fn.NHI[0:lya_fn.npt-1] = djs_median(NHI_mnx,1)
        lya_fn.bins[0:lya_fn.npt-1,0] = NHI_mnx[0,*]
        lya_fn.bins[0:lya_fn.npt-1,1] = NHI_mnx[1,*]
        lya_fn.dN[0:lya_fn.npt-1] = f_strct.dN

        lya_fn.fN[0:lya_fn.npt-1] = f_strct.fNX
        lya_fn.sig_fN[0:lya_fn.npt-1,*] = f_strct.sig_fNX
     end
     2: delvarx, lya_fN ;; No Lya f(N) !!
     else: stop
  endcase
        

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SLLS
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; OMeara et al. 2007
  ;; UVES/MIKE for N_HI > 10^19 AND z<3
  ;; Only 24 systems total
  if not keyword_set( MIKELST ) then mikelst = $
     getenv('LLSTREE')+'/Lists/mike_stat.lst'
  if not keyword_set( PRXLST ) then prxlst = $
     getenv('LLSTREE')+'/Lists/prx_stat.lst'
  if not keyword_set( FGZECH )  then $
     fgzech = getenv('XIDL_DIR')+'/LLS/fN/slls_echelle_goz.fits'

  llslst = [mikelst,prxlst]
  nbin = 2
  stp = 0.6
  slls_fnfit_z2, llslst, fgzech, NBIN=nbin, NMIN=19.0, STRCT=slls_strct, $
                 XZFIL=xzfil, ZBIN=[1.7, 3.], STP=stp
  slls_fn = {fnconstraint}
  slls_fn.zeval = slls_strct.zmean
  slls_fn.ref = 'OPB07'
  slls_fn.cosm = cosmc
  slls_fn.type = 'SLLS'
  slls_fn.comment = 'Only 30 systems total'
  slls_fn.DX = slls_strct.dX
  slls_fn.npt = nbin

  slls_fn.dN[0:slls_fn.npt-1] = slls_strct.dN
  slls_fn.NHI[0:slls_fn.npt-1] = slls_strct.logNHI
  slls_fn.bins[0:slls_fn.npt-1,0] = slls_strct.bins
  stp = slls_strct.bins[1]-slls_strct.bins[0]
  slls_fn.bins[0:slls_fn.npt-1,1] = slls_strct.bins+stp+slls_strct.svex ;; Uneven bins

  slls_fn.fN[0:slls_fn.npt-1] = slls_strct.fN
  slls_fn.sig_fN[0:slls_fn.npt-1,0] = slls_strct.sigfn1-slls_strct.fN
  slls_fn.sig_fN[0:slls_fn.npt-1,1] = slls_strct.fN-slls_strct.sigfn2

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; DLA
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  sdss_fndla, GZFIL=GZFIL, STRCT=dr5_fnstr,  STP=0.2, XZFIL=xzfil, $
              /SALL, ALLDR=alldr, CL=0.95, NPLT=10L, ZBIN=[2.3, 2.7]

  dla_fn = {fnconstraint}
  dla_fn.zeval = dr5_fnstr.zmean
  dla_fn.ref = 'PW09'
  dla_fn.cosm = cosmc
  dla_fn.type = 'DLA'
  dla_fn.comment = '$z=[2.3,2.7]$; modest SDSS bias \citep[see][]{np+09}?'
  dla_fn.DX = dr5_fnstr.dX
  dla_fn.npt = n_elements(dr5_fnstr.xval)

  dla_fn.dN[0:dla_fn.npt-1] = dr5_fnstr.dN
  dla_fn.NHI[0:dla_fn.npt-1] = dr5_fnstr.xplt
  dla_fn.bins[0:dla_fn.npt-1,0] = dr5_fnstr.bins
  stp = dr5_fnstr.bins[1]-dr5_fnstr.bins[0]
  dla_fn.bins[0:dla_fn.npt-1,1] = dr5_fnstr.bins+stp

  dla_fn.fN[0:dla_fn.npt-1] = dr5_fnstr.yplt
  dla_fn.sig_fN[0:dla_fn.npt-1,0] = dr5_fnstr.yerr1-dr5_fnstr.yplt
  dla_fn.sig_fN[0:dla_fn.npt-1,1] = dr5_fnstr.yplt-dr5_fnstr.yerr2

  nsys = fix(10.d^dr5_fnstr.yplt * dr5_fnstr.dN * dr5_fNstr.dx)

  ;; Non-detections?
  nod = where(dr5_fnstr.yplt LT -90, nnod)
  if nnod GT 0 then $
     dla_fn.sig_fN[nod] = dr5_fnstr.yerr1[nod]  ;; 2sigma limit
  
  if keyword_set(lya_fn) then $
     all_fn = [lya_fn, slls_fn, dla_fn] $
  else all_fn = [slls_fn, dla_fn]  ;; No Lya!

  ;;
  mwrfits, all_fn, outfil, /create

  print, 'Wrote ', outfil


  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LLS
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
  ;; Should put the following into XIDL once published
  cd, getenv('LLSPAP')+'HST/Science/Figures/pro', curr=curr
  resolve_routine, 'fig_lox', /compile_full_file
  cd, '../'
  fig_lox, all_strct=lls_strct, /nops, XZFIL=xzfil
  cd, curr

  ;;  Parse
  idx=1
  if abs(mean(lls_strct.acs_init.all_bins[*,idx]) - 2.3) GT 0.1 or $
     abs(mean(lls_strct.wfc3_init.all_bins[*,idx]) - 2.3) GT 0.1  then stop
  lls_constr = {llsconstraint}
  lls_constr.z_lls = lls_strct.wfc3_strct[idx].medz ;; Set by HST survey
  lls_constr.ref = 'OPW13'
  lls_constr.type = '\tlox'
  lls_constr.cosm = cosmc
  lls_constr.tau_lim = 2.
  lls_constr.lX = lls_strct.comb_lox[idx]
  lls_constr.sig_lX = lls_strct.sig_lox[idx]

  mwrfits, lls_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; MFP

  ;; HST/WFC3 MFP from IDL file 
  cd, getenv('LLSPAP')+'HST/Science/Analysis/pro/', curr=curr
  hst_sci_initparm, hst_sci
  cd, curr
  ;restore, getenv('LLSPAP')+'/HST/Science/Figures/wfc3_mfp.idl' ;; Turn this off eventually
  ;djs_iterstat, wfc3_mfp, mean=wfc3_mean, sigma=wfc3_sigm, median=wfc3_median

  nreal = hst_sci.nreal ;; Number of realizations (could do 500)
  wfc3_mfp = fltarr(NREAL)
  if not keyword_set( INFIL ) then begin
     case flg_cosm of
        0: infil = getenv('LLSPAP')+'HST/Science/Analysis/boot_mfp_wfc3.fits'
        1: infil = getenv('LLSPAP')+'HST/Science/Analysis/boot_mfp_wfc3_L70_M30_h70.fits'
        else: stop
     endcase
  endif
  lun = fxposit(INFIL+'.gz',1, /compress, /readonly)
  for kk=0L,NREAL-1 do begin
     sub_str = mrdfits(lun, 0)
     wfc3_mfp[kk] = sub_str.best_mfp
  endfor
  djs_iterstat, wfc3_mfp, mean=wfc3_mean, sigma=wfc3_sigm, median=wfc3_median

  mfp_constr = {mfpconstraint}
  mfp_constr.ref = 'OPW13'
  mfp_constr.cosm = cosmc
  mfp_constr.type = '\lmfp'
  mfp_constr.z_mfp = hst_sci.mfp_zeval 
  mfp_constr.mfp = wfc3_mean  ;; This is higher than the median
  mfp_constr.sig_mfp = wfc3_sigm

  mwrfits, mfp_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; tau_eff

  ;; D_A  -- Everything but metals
  teff_constr = {teffconstraint}
  teff_constr.ref = 'K05'  ;; Kirkman et al. 2005
  teff_constr.type = '\tlya'  ;; teff_Lya  (no LLS, no metals, no DLAs) [aka DA8s]
  teff_constr.comment = 'Converted to \tlya\ from $D_A$.  No LLS, no metals.'
  teff_constr.NHI_mnx = [12., 17.]
  teff_constr.z_teff = 2.4
  teff_constr.z_em = 2.5
  DA8s_z = 0.0062d * (1+teff_constr.z_teff)^(2.75)
  sig_DA = 0.006  
  teff_constr.teff = -1 * alog( 1. - DA8s_z)
  teff_constr.sig_teff = (1. / (1-DA8s_z) ) * sig_DA

  ;; Tytler et al. 2004 (good for z<2.2?)
;  teff_constr.ref = 'T04'  ;; Tytler et al. 2004
;  teff_constr.z_teff = 2.3
;  teff_constr.z_em = 2.44
;  DA7_z19 = .128  ; (z=1.9)
;  DA7_z = DA7_z19 * ( (1+teff_constr.z_teff)/(1+1.9))^2.57
;  sig_DA = 0.015  ;; Estimated (10%)
;  teff_constr.NHI_mnx = [12., 23.]
;  teff_constr.teff = -1 * alog( 1. - DA7_z)
;  teff_constr.sig_teff = (1. / (1-DA7_z) ) * sig_DA

  ;; Write
  mwrfits, teff_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; beta (Lya forest)

  if LYA_FLG LT 1 then begin
     beta_str = {betaconstraint}
     beta_str.z_beta = 2.13  
     
     beta_str.ref = 'K02'
     beta_str.type = '$\beta$'
     beta_str.comment = ''
     beta_str.cosm = cosmc ;; Forced
     
     beta_str.NHI = [12.5,14]
     beta_str.beta = -1.42
     beta_str.sig_beta = 0.05 ;; This is larger than Kim et al. 2002 report
     
     ;; Write
     mwrfits, beta_str, outfil
  endif


  print, 'igm_mk_fn_constraints: Wrote ', outfil

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro igm_mk_fn_constraints_z37, cosmc, TEFF_flg=TEFF_flg, LYA_FLG=lya_flg

  ;; Generate FITS file on z=3.7 constraints for f(N)
  ;; EXTEN=1 is the fN constraints (Lya forest, SLLS, DLAs)
  ;;    Lya:  Eyeball one of Kim's data points (dodgy)
  ;;    SLLS: O'Meara et al. 2006, Petitjean group
  ;;    DLA:  PW09
  ;; EXTEN=2 is the LLS constraint  POW10
  ;; EXTEN=3 is the MFP constraint PWO09
  ;; EXTEN=4 is a teff^Lya (D_A) constraint [FG08, Dall'Aglio]
  ;; EXTEN=5 is a beta constraint 

  outfil = 'fn_constraints_z3.7.fits'
  if not keyword_set(TEFF_flg) then teff_flg = 0L

  if not keyword_set(LYA_FLG) then lya_flg = 0L  ;; Kim et al. 2002

  ;; XZ
  case cosmc of 
     'WMAP5': xzfil=getenv('SDSSPATH')+'DR7_QSO/xz_val_L74_M26.fits.gz' 
     'VANILLA': xzfil=getenv('SDSSPATH')+'DR7_QSO/xz_val_L70_M30.fits.gz' 
     else: stop
  endcase

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lya Forest
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lya_fn = {fnconstraint}
  lya_fn.zeval = 3.75  

  case lya_flg of
     0: begin
        lya_fn.ref = 'K01'
        lya_fn.type = 'Lya Forest'
        lya_fn.comment = 'Read from Figure 4'
        lya_fn.cosm = cosmc ;; Forced
        lya_fn.DX = 3. ;; Estimated
        
        lya_fn.npt = 1
        
        kim_NHI = [13.85]
        kim_lfN  = [-12.37]
        for ii=0L,lya_fn.npt-1 do begin
           
           lya_fn.NHI[ii] = kim_NHI[ii]
           N1 = kim_NHI[ii] - 0.2
           N2 = kim_NHI[ii] + 0.2
           dN = 10.d^N2 - 10.d^N1
           
           lya_fn.dN[ii] = dN
           lya_fn.bins[ii,0] = N1
           lya_fn.bins[ii,1] = N2
           
           kim_fN = 10.d^(kim_lfN[ii])
           dXdz_q0 = (1+lya_fn.zeval)
           dXdz_cosm = cosm_dxdz(lya_fn.zeval, /noinit) ;; This is very close, actually
           new_kim_fN = kim_fN / dXdz_q0 * dXdz_cosm
           
           lya_fn.fN[ii] = alog10(new_kim_fN)
           lya_fn.sig_fN[ii,*] = 0.1 ;; By-eye estimate
        endfor
     end 
     1: begin ;; Kim et al. 2002
        lya_fn.zeval = 3.67   ;; <z> for the lines in the interval  
        lya_fn.ref = 'K02'
        lya_fn.type = 'Lya Forest'
        lya_fn.comment = 'Recalculated for our Cosmology'
        lya_fn.cosm = cosmc

        NHI_mnx = [ [12.5, 13.], [13., 13.5], [13.5, 14.], [14., 14.5], [14.5,15]]
        sz_NHI = size(NHI_mnx, /dimen)

        ;; Calculate on-the-fly
        igm_getfn_kim2002, [3.4, 4.0], NHI_mnx, f_STRCT, FLG_COSM=0L
        lya_fn.DX = f_strct.DX
        
        lya_fn.npt = sz_NHI[1]
        lya_fn.NHI[0:lya_fn.npt-1] = djs_median(NHI_mnx,1)
        lya_fn.bins[0:lya_fn.npt-1,0] = NHI_mnx[0,*]
        lya_fn.bins[0:lya_fn.npt-1,1] = NHI_mnx[1,*]
        lya_fn.dN[0:lya_fn.npt-1] = f_strct.dN

        lya_fn.fN[0:lya_fn.npt-1] = f_strct.fNX
        lya_fn.sig_fN[0:lya_fn.npt-1,*] = f_strct.sig_fNX
     end
     else: stop
  endcase
  

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SLLS
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; O'Meara et al. 2007
  ;; UVES/MIKE for N_HI > 10^19 AND z<3
  ;; Only 24 systems total
  if not keyword_set( MIKELST ) then mikelst = $
     getenv('LLSTREE')+'/Lists/mike_stat.lst'
  if not keyword_set( PRXLST ) then prxlst = $
     getenv('LLSTREE')+'/Lists/prx_stat.lst'
  if not keyword_set( FGZECH )  then $
     fgzech = getenv('XIDL_DIR')+'/LLS/fN/slls_echelle_goz.fits'

  llslst = [mikelst,prxlst]
  nbin = 2
  stp = 0.6
  slls_fnfit_z2, llslst, fgzech, NBIN=nbin, NMIN=19.0, STRCT=slls_strct, $
                 XZFIL=xzfil, ZBIN=[3.1,4.5], STP=stp
  slls_fn = {fnconstraint}
  slls_fn.zeval = slls_strct.zmean
  slls_fn.ref = 'OPB07'
  slls_fn.cosm = cosmc
  slls_fn.type = 'SLLS'
  slls_fn.comment = '$z=[3.1,4.5]$'
  slls_fn.DX = slls_strct.dX
  slls_fn.npt = nbin

  slls_fn.dN[0:slls_fn.npt-1] = slls_strct.dN
  slls_fn.NHI[0:slls_fn.npt-1] = slls_strct.logNHI
  slls_fn.bins[0:slls_fn.npt-1,0] = slls_strct.bins
  stp = slls_strct.bins[1]-slls_strct.bins[0]
  slls_fn.bins[0:slls_fn.npt-1,1] = slls_strct.bins+stp+slls_strct.svex ;; Uneven bins

  slls_fn.fN[0:slls_fn.npt-1] = slls_strct.fN
  slls_fn.sig_fN[0:slls_fn.npt-1,0] = slls_strct.sigfn1-slls_strct.fN
  slls_fn.sig_fN[0:slls_fn.npt-1,1] = slls_strct.fN-slls_strct.sigfn2

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; DLA
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if not keyword_set( GZFIL ) then gzfil = '~/SDSS/DR5_QSO/dr5_dlagz_s2n4.fits'
  sdss_fndla, GZFIL=GZFIL, STRCT=dr5_fnstr,  STP=0.2, XZFIL=xzfil, $
              /SALL, ALLDR=alldr, CL=0.95, NPLT=10L, ZBIN=[3.3, 4.2]

  dla_fn = {fnconstraint}
  dla_fn.zeval = dr5_fnstr.zmean
  dla_fn.ref = 'PW09'
  dla_fn.cosm = cosmc
  dla_fn.type = 'DLA'
  dla_fn.comment = '$z=[3.3,4.2]$; modest SDSS bias \citep[see][]{pwo09}'
  dla_fn.DX = dr5_fnstr.dX
  dla_fn.npt = n_elements(dr5_fnstr.xval)

  dla_fn.dN[0:dla_fn.npt-1] = dr5_fnstr.dN
  dla_fn.NHI[0:dla_fn.npt-1] = dr5_fnstr.xplt
  dla_fn.bins[0:dla_fn.npt-1,0] = dr5_fnstr.bins
  stp = dr5_fnstr.bins[1]-dr5_fnstr.bins[0]
  dla_fn.bins[0:dla_fn.npt-1,1] = dr5_fnstr.bins+stp

  dla_fn.fN[0:dla_fn.npt-1] = dr5_fnstr.yplt
  dla_fn.sig_fN[0:dla_fn.npt-1,0] = dr5_fnstr.yerr1-dr5_fnstr.yplt
  dla_fn.sig_fN[0:dla_fn.npt-1,1] = dr5_fnstr.yplt-dr5_fnstr.yerr2

  nsys = fix(10.d^dr5_fnstr.yplt * dr5_fnstr.dN * dr5_fNstr.dx)

  ;; Non-detections?
  nod = where(dr5_fnstr.yplt LT -90, nnod)
  if nnod GT 0 then $
     dla_fn.sig_fN[nod] = dr5_fnstr.yerr1[nod]  ;; 2sigma limit
  
  all_fn = [lya_fn, slls_fn, dla_fn]
  mwrfits, all_fn, outfil, /create

  print, 'Wrote ', outfil


  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LLS (SDSS)
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cd, getenv('LLSPAP')+'SDSS/DR7/Analysis/pro/', curr=curr2
  RESOLVE_ROUTINE, 'max_lozpower'
  lls_dr7_initparm, init
  cd, curr2

  llsfil = init.llsfil
  qsofil = init.qsofil 
  maxdz = init.maxoff 
  vprox = init.vprox
;  xzfil = init.xzfil
  zem_min = init.zem_min
  BINS = init.all_bins

  ;; Calculate from the data
  sdss_llslozx, qsofil, llsfil, lls_strct, BINS=bins, XZFIL=xzfil, VPROX=vprox, $
                ZEM_MIN=zem_min, MAXDZ=maxdz

  lls_constr = {llsconstraint}
  idx = 3L
  lls_constr.z_lls = lls_strct[idx].meanzabs ;; Set by HST survey
  lls_constr.ref = 'POW10'
  lls_constr.type = '\tlox'
  lls_constr.cosm = cosmc
  lls_constr.tau_lim = 2.
  lls_constr.lX = lls_strct[idx].lox
  lls_constr.sig_lX = mean(lls_strct[idx].siglx)

  mwrfits, lls_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; MFP

  ;; ;;;;
  ;; SDSS
  ;sdss_mfpfil = getenv('LLSPAP')+'/taueff/Analysis/mfp_summ.fits'
  sdss_mfpfil = getenv('LLSPAP')+'/HST/Science/Analysis/sdss_mfp_summ.fits'
  sdss_str = xmrdfits(sdss_mfpfil, 1)
  idx = where(sdss_str.zmnx[0] LT 3.7 and sdss_str.zmnx[1] GT 3.7, nidx)
  if nidx NE 1 then stop

  mfp_constr = {mfpconstraint}
  mfp_constr.ref = 'PWO09'
  mfp_constr.cosm = cosmc
  mfp_constr.type = '\lmfp'
  mfp_constr.z_mfp = mean(sdss_str[idx].zmnx)
  mfp_constr.mfp = sdss_str[idx].mfp
  mfp_constr.sig_mfp = sdss_str[idx].sig_mfp

  mwrfits, mfp_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; tau_eff

  ;; D_A  -- Everything but metals
  teff_constr = {teffconstraint}

  case teff_flg of 
     0: begin ;; Faucher-Giguerre
        teff_constr.ref = 'FG08' ;; Faucher-Giguerre et al. 2008
        teff_constr.type = '\tlya' ;; teff_Lya  (no LLS, no metals, no DLAs) [aka DA8s]
        teff_constr.comment = 'Higher than other evaluations.'
        teff_constr.NHI_mnx = [12., 19.]
        teff_constr.z_teff = 3.7
        teff_constr.z_em = 3.8
        teff_constr.teff = 0.795
        teff_constr.sig_teff = 0.062 ;; Includes systematics
     end
     1: begin ;; Kim et al. 2007
        teff_constr.ref = 'K07' ;; 
        teff_constr.type = '\tlya' ;; teff_Lya  (no LLS, no metals, no DLAs) [aka DA8s]
        teff_constr.comment = 'Not much data at this redshift'
        teff_constr.NHI_mnx = [12., 19.]
        teff_constr.z_teff = 3.7
        teff_constr.z_em = 3.8
        teff_constr.teff = 0.653
        teff_constr.sig_teff = 0.06 ;; (My error estimate, not hers)
     end
     else: stop
  endcase

  ;; Write
  mwrfits, teff_constr, outfil

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; beta (Lya forest)

  beta_str = {betaconstraint}
  beta_str.z_beta = 3.7  
  
  beta_str.ref = ' '  ; Done by eye-ball by JXP
  beta_str.type = '$\beta$'
  beta_str.comment = 'Shallower than z=2, but very uncertain.'
  beta_str.cosm = cosmc ;; Forced

  beta_str.NHI = [12.5,14]
  beta_str.beta = -1.3
  beta_str.sig_beta = 0.15

  ;; Write
  mwrfits, beta_str, outfil


  print, 'igm_mk_fn_constraints: Wrote ', outfil

  return
end

;; ;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;
;;  toutfil=  -- Outfil name for the z=2.5 constraints
;; igm_mk_fn_constraints, flg_cosm=1, toutfil='fn_constraints_z2.5_vanilla.fits', z2_lflg=1, /nozt
pro igm_mk_fn_constraints, FLG_COSM=flg_cosm, TEFF_FLG=teff_flg, z2_LFLG=z2_LFLG, $
                           toutfil=toutfil, NOZT=nozt
  compile_opt strictarr

  if not keyword_set( FLG_COSM ) then flg_cosm = 0L

  ;; Set up the cosmology
  case flg_cosm of 
     0: begin
        cosmc = 'WMAP5'
        cosm_common, /w05map ;; Dunkley et al. 2009
     end
     1: begin
        cosmc = 'VANILLA' ;; 0.7, 0.7, 03
        cosm_common, /vanilla ;; Dunkley et al. 2009
     end
     else: stop
  endcase

  ;; z=2.5
  igm_mk_fn_constraints_z25, cosmc, LYA_FLG=z2_LFLG, outfil=toutfil

  ;; z=3.7
  if not keyword_set(NOZT) then $
     igm_mk_fn_constraints_z37, cosmc, teff_flg=teff_flg ;, LYA_FLG=1
       
  return
end

