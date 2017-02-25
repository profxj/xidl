;+ 
; NAME:
; lris_tstlong  
;     Version 1.1
;
; PURPOSE:
;  Tests the Longslit routines on LRIS
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
pro lris_tstlong

;
;  if  N_params() LT 1  then begin 
;      print,'Syntax - ' + $
;        'hires_tstrdx, afil, ffil, ofil, CHIP=, SETUP= [v1.1]'
;      return
;  endif 


  ;; Create plan file
  long_plan, '*.fits*', 'Raw/'

  ;; Check the plan file
  spawn, 'diff plan.par save_plan.par > ../diff.log' 

  ;; Use new plan file
  long_reduce, 'plan_B1200.par'
  date = string(bin_date())
  cdate = date[0]+'_'+date[1]+'_'+date[2]
  spawn, 'cp plan.log ../B1200_'+cdate+'.log'

  ;; Examine the output
  spawn, 'gv wave-lblue0136.ps'
  spawn, 'gv Science/sci-lblue1103.ps'

  ;; COADD 
  spawn, 'mkdir FSpec'
  files = ['Science/sci-lblue1103.fits', 'Science/sci-lblue1104.fits'] 
  long_coadd, files, 2, OUTFIL='FSpec/A0235+1624b_x.fits'

  ;; Flux
  spawn, 'mkdir Flux'
  sci = 'Science/sci-lblue1101.fits'
  bsens = long_sensfunc(sci, 'g191b2b_005', 'Flux/g191b2b_bsens.fits')
  
  long_fluxcal, 'FSpec/A0235+1624b_x.fits', 'Flux/g191b2b_bsens.fits', OUTFIL='FSpec/A0235+1624b_XF.fits'

  x_specplot, 'FSpec/A0235+1624b_XF.fits', inflg=2, /block

  return
end
