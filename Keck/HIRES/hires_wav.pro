;+ 
; NAME:
; hires_wav
;
; PURPOSE:
;   Create a line list appropriate for the HIRES emulator 
; provided the redshift.
;
; CALLING SEQUENCE:
;  hires_wav, z, fil, /LLS   
;  
; INPUTS:
;   z  -- Redshift
;  fil -- Output fil
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /LLS  -- Create a line list for LLS
;   /shrt -- shorter line list
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_wav, 2.2, 'q0440-221_wav.dat'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro hires_wav, z, fil, LLS=lls, GRB=grb, PRINT=print, shrt=shrt

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'hires_wav, z, fil, /LLS, [v1.1]'
    return
  endif 

  print, 'hires_wav: Note the wavelengths are converted to air'

  if keyword_set( LLS ) then linfil = getenv('XIDL_DIR')+ $
    '/Spec/Lines/Lists/lls.lst' 
  if keyword_set( GRB ) then linfil = getenv('XIDL_DIR')+ $
    '/Spec/Lines/Lists/qal.lst' 
  if keyword_set( shrt ) then linfil = getenv('XIDL_DIR')+ $
    '/Spec/Lines/Lists/dla_shrt.lst' 
  if not keyword_set(LINFIL) then $
	linfil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/dla.lst' 

  close, /all
  openr, 20, linfil
  openw, 2, fil
  readf, 20, nln, FORMAT='(i3)' 
  wave = 0.d
  mtlnm=' '
  ;; Read
  allwv = dblarr(nln)
  restwv = dblarr(nln)
  allmtl = strarr(nln)
  for i=0,nln-1 do begin
      readf, 20, FORMAT='(f9.4,x,a10)', wave, mtlnm
      allwv[i] = wave
      restwv[i] = wave
      allmtl[i] = mtlnm
  endfor

  ;; Redshift + Air
  allwv = allwv*(1.+z)
  vactoair, allwv

  ;; Print
  for i=0,nln-1 do begin
      printf, 2, allwv[i], 'F 1', allmtl[i], FORMAT='(1x,f8.3,6x,a,2x,a)'
      if keyword_set(PRINT) then print, restwv[i], ' ', allwv[i]
  endfor

  close, /all

  return
end
