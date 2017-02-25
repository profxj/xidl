;+ 
; NAME:
; kast_tstlong  
;     Version 1.1
;
; PURPOSE:
;  Tests the Longslit routines on the Kast spectrometer
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
pro kast_tstlong

;
;  if  N_params() LT 1  then begin 
;      print,'Syntax - ' + $
;        'hires_tstrdx, afil, ffil, ofil, CHIP=, SETUP= [v1.1]'
;      return
;  endif 


  ;; Create plan file
  long_plan, '*.ccd*', 'Raw/'

  ;; Check the plan file
  spawn, 'diff plan.par save_plan.par > ../diff.log' 

  ;; Use new plan file
  long_reduce, 'good_plan.par'
  date = string(bin_date())
  cdate = date[0]+'_'+date[1]+'_'+date[2]
  spawn, 'cp plan.log ../B1200_'+cdate+'.log'

  ;; Examine the output
  spawn, 'gv wave-b223.ccd.ps'
  spawn, 'gv wave-b223.ccd.ps'
  spawn, 'gv Science/sci-b259.ccd.ps'
  spawn, 'gv Science/sci-r260.ccd.ps'

  ;; COADD 
  spawn, 'mkdir FSpec'
  files = findfile('Science/sci-b26*gz')
  long_coadd, files, 1, OUTFIL='FSpec/J0905+1014_xb.fits'
  files = findfile('Science/sci-r26*gz')
  long_coadd, files[1:*], 1, OUTFIL='FSpec/J0905+1014r_x.fits'

  ;; Flux
  spawn, 'mkdir Flux'
  sci = 'Science/sci-b257.ccd.fits'
  bsens = long_sensfunc(sci, 'feige34_005', 'Flux/feige34_bsens.fits')
  scir = 'Science/sci-r258.ccd.fits'
  rsens = long_sensfunc(scir, 'feige34_005', 'Flux/feige34_rsens.fits')
  
  long_fluxcal, 'FSpec/J0905+1014_xb.fits', 'Flux/feige34_bsens.fits', OUTFIL='FSpec/J0905+1014b_XF.fits'
  long_fluxcal, 'FSpec/J0905+1014r_x.fits', 'Flux/feige34_rsens.fits', OUTFIL='FSpec/J0905+1014r_XF.fits'

  x_specplot, 'FSpec/J0905+1014b_XF.fits', inflg=2, /block
  x_specplot, 'FSpec/J0905+1014r_XF.fits', inflg=2, /block

  return
end
