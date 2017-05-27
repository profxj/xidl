function deimos_applyflat, spec, flat2d, slitfn, invvar= invvar, quick=quick
;+
; NAME:
;   Deimos_applyflat
;
; PURPOSE:
;    Apply flat field to a spectrum, 
;     for purposes of removing effects of variations in slit thickness
;     and fringing.
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   flatspec = Deimos_applyflat(spec, flat2d,slitfn, [ invvar= invvar] )
; 
; INPUTS:
;   spec -- 2d rectified spectrum with lambda as fast dimension (transposed
;           DEIMOS data)
;   flat2d -- corresponding rectified 2d flat(MINUS ONE) from DEIMOS_makeflat
;   slitfn -- corresponding slit function from DEIMOS_makeflat; either
;             1d or 2d (flat1d or varslitfn)
;
; KEYWORD PARAMETERS:
;   invvar -- (input/output) inverse variance of input spectrum.
;       rescaled according to the flat field to keep S/N constant.
;
; BUGS:
;
; OUTPUTS:
;   flatspec -- resulting spec has been normalized by flat2d & slitfn while
;              holding saturated lines below 65K
;
; MODIFICATION HISTORY:
;   C. Marinoni, Spring 2001
;   MD Jan-Apr-2002 
;   JAN 3 Jul 2002
;-

  if (NOT keyword_set(invvar)) then invvar=float(spec)*0.+1.
           

  nn = (size(spec))[2] ;spatial range
  mm = (size(spec))[1] ;spectral range

  slitfndim = (size(slitfn))[0]

; Interpolate over bad pixels/columns in flat - due to their
;   flatness in spatial direction, this should be fine                                        
  if quick le 0 then begin
     badpix = (flat2d eq 0.)
     flat2 = djs_maskinterp(flat2d+1.,badpix,iaxis=1,/const)
  endif else flat2 = 1.

; expand slitfn to match the data
  if slitfndim eq 1 then sfall = slitfn ## (fltarr(mm)+1) $
    else sfall = rebin(slitfn, mm, nn)
; Did not use to flat field saturated lines!!!  
;  ii = where(spec gt 65000,  nii) 


; DO NOT WISH TO prevent lines being corrected to a value > 65535
  flatspec = spec/flat2/sfall ;< 65535

; don't make big changes in the inverse variance due to the slit function!
	sfall= (sfall > 0.8) < 1.2

; Note that we have to rescale the inverse variance as well...
  invvar = invvar*(flat2*sfall)^2

;  if (nii gt 0)  then flatspec[ii] = 65535. 

;  not sure if the following is desirable. . .
;  if (nii gt 0)  then invvar[ii] = 0. 

  return, flatspec


end
