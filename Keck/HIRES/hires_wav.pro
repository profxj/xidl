;+ 
; NAME:
; x_radec
;
; PURPOSE:
;    Turns RA, DEC in XX:XX:XX.X -DD:DD:DD.D format to decimal deg
;      or VICE VERSA
;
; CALLING SEQUENCE:
;   
;   x_radec, ra, dec, rad, decd, /ARCS, /FLIP
;
; INPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; RETURNS:
;
; OUTPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; OPTIONAL KEYWORDS:
;   ARCS - Outputs in arcseconds
;   FLIP - Gives RA and DEC from decimal RA,DEC
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_wav, 'q0440-221_wav.dat'
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro hires_wav, z, fil, LLS=lls

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'hires_wav, z, fil, /LLS, (v1.0)'
    return
  endif 

  if keyword_set( LLS ) then linfil = getenv('XIDL_DIR')+ $
    '/Spec/Lines/Lists/lls.lst' $ 
  else linfil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qal.lst' 

  close, /all
  openr, 20, linfil
  openw, 2, fil
  readf, 20, nln, FORMAT='(i3)' 
  wave = 0.d
  mtlnm=' '
  ;; Read
  for i=0,nln-1 do begin
      readf, 20, FORMAT='(f9.4,x,a10)', wave, mtlnm
      obs = (z+1.0)*wave
      ;; Print
      printf, 2, obs, 'F 1', mtlnm, FORMAT='(1x,f8.3,6x,a,2x,a)'
  endfor

  close, /all

  return
end
