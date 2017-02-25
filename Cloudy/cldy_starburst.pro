;+ 
; NAME:
; cldy_starburst
;   Version 1.1
;
; PURPOSE:
;   Creates a Cloudy input file from the Starburst99 template
;
; CALLING SEQUENCE:
;   
;   cldy_cuba, fil, age, outfil, FIXG=fixg
;
; INPUTS:
;   fil  - Filname
;   age  - Age of the starburst
;
; RETURNS:
;
; OUTPUTS:
;   outfil  - Cloudy output file
;
; OPTIONAL INPUTS:
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
;   19-Jun-2006 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_starburst, fil, age, outfil, FIXG=fixg

  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'cldy_cuba, fil, age, outfil, /FIXG [v1.1]'
      return
  endif 

  ;; Parse starburst file
  nlin = 1221L
  close, /all
  openr, 1, fil
  wv = fltarr(nlin)
  flux = dblarr(nlin,36)
  duma = fltarr(36)
  dumf = 0.
  dumc = ''
  aage = [findgen(20)+1,30.,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900]

  ;; Loop
  for ii=0L,2 do readf, 1, dumc
  for ii=0L,nlin-1 do begin
      readf, 1, dumf, duma
      wv[ii] = dumf
      flux[ii,*] = duma
  endfor
  close, 1

  openw, 2, outfil
  printf, 2, 'interpolate (0.00001 -30.0)'

  ;; Interpolate at each wavelength!
  esv = 0.
  for ii=nlin-1,0,-1 do begin

      fx = interpol(flux[ii,*], aage, age, /spline)
      fx = fx[0]
      energy = 912./wv[ii]
;      if fx GT 1E-30 then jnu = alog10(fx) else jnu = -30.
      jnu= fx - 56.
      if keyword_set(FIXG) and wv[ii] LT 50 then jnu = -30.
      if energy NE esv then $  ; Multiple energies in a few spots
        printf, 2, 'continue ('+strtrim(energy,2)+' '+ $
        string(jnu,FORMAT='(f7.3)')+')'
      esv = energy
  endfor
  printf, 2, 'continue (7400000 -30.0)'
  close, 2

  return
  
end
