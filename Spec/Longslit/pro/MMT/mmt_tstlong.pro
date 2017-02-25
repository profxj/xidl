;+ 
; NAME:
; mmt_tstlong
;     Version 1.1
;
; PURPOSE:
;  Tests the Longslit routines on the BCS spectrometer on the MMT
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-May-2007 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mmt_tstlong

;
;  if  N_params() LT 1  then begin 
;      print,'Syntax - ' + $
;        'hires_tstrdx, afil, ffil, ofil, CHIP=, SETUP= [v1.1]'
;      return
;  endif 


  ;; Create plan file
  long_plan, '*.fits.gz', 'Raw/'

  ;; Check the plan file
  spawn, 'diff plan.par save_plan.par > ../diff.log' 

  ;; Create Arc file
  mmt_waveblue, 'Raw/a.0029.fits', 'Raw/a.0026.fits', 'Raw/a.9901.fits'

  ;; Use new plan file
  long_reduce, 'plan_0912.par'
  date = string(bin_date())
  cdate = date[0]+'_'+date[1]+'_'+date[2]
  spawn, 'cp plan.log ../plan'+cdate+'.log'

  ;; Examine the output
  spawn, 'gv wave-a.9901.ps'
  spawn, 'gv Science/sci-a.0015.ps'

  ;; COADD 
  spawn, 'mkdir FSpec'
  files = ['Science/sci-a.0015.fits']
  long_coadd, files, 1, OUTFIL='FSpec/J0912_x.fits'

  ;; Flux
  spawn, 'mkdir Flux'
  sci = 'Science/sci-a.0006.fits'
  sens = long_sensfunc(sci, 'g191b2b_005', 'Flux/g191b2b_bsens.fits')

  long_fluxcal, 'FSpec/J0912_x.fits', 'Flux/g191b2b_bsens.fits', $
                OUTFIL='FSpec/J0912_XF.fits'

  x_specplot, 'FSpec/J0912_XF.fits', inflg=2, /block

  return
end
