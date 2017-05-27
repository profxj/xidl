function deimos_flatfield, spec,  flat_image, flat2d, invvar= invvar, $
     flat1d=flat1d,  twod= twod, mask = mask
;+
; NAME:
;   Deimos_flatfield
;
; PURPOSE:
;    Flat field only in spatial direction only, or in 2d, 
;     for purposes of removing effects of variations in slit thickness
;     and fringing.
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   flatspec = Deimos_flatfield(spec, flat_image, [flat2d, invvar= invvar, $ 
;                     mask=mask, flat1d= flat1d, /twod])
; 
; INPUTS:
;   spec -- 2d spectrum with lambda as fast dimension (transposed
;           DEIMOS data)
;   flat_image -- corresponding flat
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;   invvar -- (input) inverse variance of input flat image
;   mask -- (output) 1 for good data, 0 for bad
;   twod   -- (input) switch to use 2-D flatfielding,  rather than default 1-D
;   flat1d -- (output) derived 1d flat (SLITFN)
;
; BUGS:
;
; OUTPUTS:
;   flatspec -- resulting spec has been normalized, by flat while
;              holding saturated lines below 65K
;   flat2d  -- 2d-flat generated
;
; MODIFICATION HISTORY:
;   C. Marinoni, Spring 2001
;   MD Jan-Apr-2002 
;-

  if (NOT keyword_set(invvar)) then invvar=float(flat_image)*0.+1.
           
; 1] Correction for gradient along dispersion direction
  nn = (size(flat_image))[2] ;spatial range
  mm = (size(flat_image))[1] ;spectral range
  amp =  round(0.1*nn) > 5; width 10%
  inf = amp
  sup = nn-1 - amp
  flat1d = fltarr(nn)

  weight = fltarr(mm) +1.
; Select some good rows (central ones) in order to determine the mean
;  cent = round(nn/2.)
  s1=floor(mm/4)
  s2=ceil(3*mm/4) ;sum over middle half of array for normalization

  weight[0:s1-1] = (findgen(s1)/s1)
  weight[s2:mm-1] = ((mm+s2-findgen(mm-s2))/(mm-s2))

; to make bad columns less relevant, normalize spectrum
   illum = median(flat_image[*, round(nn/2.)-5:round(nn/2.)+5], 11)
   illum = illum[*, 5]
   illum[0:10] = illum[11]
   illum[mm-11:mm-1] = illum[mm-12]
   illum = smooth(illum, 5)
   illumcorr = (fltarr(nn)+1.) ## (1/(illum >  1.))
   
  inv = invvar* ( (fltarr(nn)+1.) ## (weight) )*(flat_image lt 65000) ;weighted matrix


; to deal with bad columns - attempt to normalize in spec direction
;                           initially
 
   flat1d = total(flat_image*illumcorr*inv, 1)/total(inv, 1) 
;   flatnorm = total(flat1d[inf:sup])/(sup-inf+1.)
   flatnorm = median(flat1d[inf:sup], /even)
   flat1d = flat1d/flatnorm

  sfall = flat1d ## (fltarr(mm)+1)
  mask = sfall GT .5 AND sfall LT 1.5

  flatspec = spec*0.
; Do not flat field saturated lines!!!  
  ii = where(spec gt 65000,  nii) 

if NOT keyword_set(twod) then begin ;1-D flatfielding
; to prevent that some lines be corrected to a value > 65535
;  for i=0, nn-1 do flatspec[*, i] = [spec[*, i]/flat1d[i]  <  65535.]
   flatspec = spec/sfall <  65535.

; other are already saturated 
  if (nii gt 0)  then flatspec[ii] = 65535. 

  return, flatspec

endif else begin ;2-D flatfielding

   normflat = flat_image/sfall
   medflat = median(normflat, 5) ;median smooth map
; NEED TO BE CAREFUL AT THE EDGES!
   for i=0, 4 do medflat(*, i) = median(normflat(*, i), 5)
   for i=nn-5, nn-1 do medflat(*, i) = median(normflat(*, i), 5)

  good =  where(medflat gt 100 AND invvar ne 0) ;acceptable regions
  flat2d = medflat*0. +1.
  mask=bytarr(mm,nn)
  mask[good] = 1
  mask[*, 0:inf-1] = 0
  mask[*, sup+1:nn-1] = 0

; the old way, improved
;  rowsum = total(medflat*mask, 2)/total(mask, 2)

;this works a little better still
  summed =  (fltarr(nn)+1.)## illum

  flat2d[good] = medflat[good]/summed[good]
  
  mask=bytarr(mm,nn)
  mask[good]=1B
;not a very sophisticated mask -- should use invvar TBD!

; to prevent that some lines be corrected to a value > 65535
  flatspec = spec/flat2d/sfall < 65535

; Note that we have to rescale the inverse variance as well!
  invvar = invvar*(flat2d*sfall)^2
  if (nii gt 0)  then flatspec[ii] = 65535. 
  if (nii gt 0)  then invvar[ii] = 0. 

;set saturated lines at saturation

  return, flatspec
endelse

end
