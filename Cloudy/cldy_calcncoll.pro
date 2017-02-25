;+ 
; NAME:
; cldy_calcncoll
;  V1.1
;
; PURPOSE:
;    Calculates the column density of a series of ions
;     for CIE given a temperature
;
; CALLING SEQUENCE:
;   
;  cldy_calcncoll, temp, ions, colms, NHI=, MTL=
;
; INPUTS:
;   temp -- Temperature of the gas (K)
;   ions -- Ions to calculate columns for [[Z,i]]
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NHI - N(HI) value  [default = 15]
;  MTL - Metallicity value  [default = 0. (solar)]
;  INFIL - Name of fits file containing Cloudy output of 
;      collsional ionization calculation (e.g. cloudy_collisions.fits)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; cldy_calcncoll, 1e5, [[14,2], [8,6]], colms, NHI=14.8
;
; PROCEDURES CALLED:
;  prs_cldycoll
;  printcol
;
; REVISION HISTORY:
;   04-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_calcncoll, temp, ions, colms, NHI=nhi, MTL=mtl, INFIL=infil

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'cldy_calcncoll, temp, ions, colms, NHI=, MTL=, INFIL= (v1.1)'
    return
  endif 

  if not keyword_set( NHI ) then NHI = 15.
  if not keyword_set( MTL ) then MTL = 0.
  if not keyword_set( INFIL ) then $
    infil = getenv('XIDL_DIR')+'/Cloudy/cloudy_collisions.fits'

  print, 'cldy_calcncoll:  Assuming NHI=', NHI, ' and MTL=', mtl

  ;; Read in structure
  prs_cldycoll, coll, infil

  ;; Loop
  sz = size(ions,/dimensions)
  if n_elements(ions) EQ 2 then nion = 1 else nion = sz[1]
  colms = fltarr(nion) 
  for qq=0L,nion-1 do begin
      
      ;; Grab abundances  
      getabnd, elm, ions[0,qq], abd, flag=1

      ;; Calculate model ratios
      all_colm = coll.X[ions[0,qq],ions[1,qq]] - coll.X[1,1] + NHI + abd - 12. $
        + MTL*(ions[0,qq] NE 1)

      ;; Interpolate
      colms[qq] = interpol(all_colm, coll.T, temp)

      ;; Print
      printcol, ions[0,qq], ions[1,qq], colms[qq]

  endfor
  return
end
