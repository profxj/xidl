;+ 
; NAME:
; cldy_cuba
;   Version 1.0
;
; PURPOSE:
;    Creates a Uplot for a given NHI, FeH, nH
;
; CALLING SEQUENCE:
;   
; cldy_Uplot, grid, NHI, FeH, nH, ions
;
; INPUTS:
;   grid  - CLOUDY grid
;   NHI - Can be an array of values
;   FeH
;   nH
;   ions  - Array of [Z,ion] vectors
;
; RETURNS:
;   
;
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_cuba, fil, z, outfil
; cldy_cuba, '/u/xavier/Cloudy/Spec/Data/CUBA/Q1G0/bkgthick.out',
; 0.35, '/u/xavier/Cloudy/Spec/Output/q1g0_z035.spec'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_cuba, fil, z, outfil, FIXG=fixg

  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'cldy_cuba, fil, z, outfil, [v1.0]'
      return
  endif 

  ;; Open Madau file
  close, /all
  openr, 1, fil
  zval = fltarr(10)
  wv = dblarr(432L)
  flux = dblarr(432L,10)
  dumf = dblarr(11)

  ;; Loop
  for ii=0L,999 do begin
      readf, 1, zval, FORMAT='(11x,10f11.4)' 
      for jj=0L,431 do begin
          readf, 1, dumf;, FORMAT='(10f11.4)'
          flux[jj,*] = dumf[1:*]
          wv[jj] = dumf[0]
      endfor
      if (z GE zval[0] OR z LE zval[9]) then break
  endfor

  close, 1

  openw, 2, outfil
  printf, 2, 'interpolate (0.00001 -30.0)'

  ;; Spline at each wavelength!
  esv = 0.
  for ii=431L,0,-1 do begin

      fx = interpol(flux[ii,*], zval, z, /spline)
      energy = 912./wv[ii]
      if fx GT 1E-30 then jnu = alog10(fx) else jnu = -30.
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
