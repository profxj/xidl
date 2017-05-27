;+
; NAME:
;   deimos_slit_function
;
; PURPOSE:
;   Derive slit function from rectified flat
; 
; CALLING SEQUENCE:
;   deimos_slit_function, flat, slitfn, sigma, mask
;
; INPUTS:
;   flat(nx,ny) -- rectified flat, with x the lambda dimension
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;   slit (ny) -- slit function along spatial direction
;   sigma (ny) -- 
;   noflat (byte array) describing absence of data
;
; COMMENTS:
;  probably has serious bugs introduced by MD
;
; REVISION HISTORY:
;
;       Fri Feb 22 14:48:01 2002, Douglas Finkbeiner (dfink)
;        MD - transposed to Lambda first
;
;----------------------------------------------------------------------
pro deimos_slit_function, flat, slitfn, slitfnsig, noflat

  ny    = (size(flat, /dimens))[1]
  nx    = 1024
  r     = rebin(flat,nx,ny)
  slitfn     = fltarr(ny)
  slitfnsig  = fltarr(ny)
  noflat = bytarr(nx, ny)
  rt    = total(r,2)
  rflat = r/((rt##(fltarr(ny)+1)/ny) > 1)
  
  for i=0, ny-1 do begin 
     djs_iterstat, rflat[*, i], mean=mn, sigma=sig, mask=mask, sigrej=3.5
     slitfn[i]     = mn
     slitfnsig[i]  = sig
;     noflat[*, i]  = mask eq 0
  endfor 

  if max(slitfn) eq 0 then message, 'slit function is zero ?!?'
  slitfn = slitfn / max(slitfn)

  sfall = slitfn # (fltarr(ny)+1)
  mask = sfall GT .5
  sflat = mask*r/(sfall*mask + (1B-mask))

  mask = mask AND (r GT max(r) / 5)

  noflat = (1B-mask)

; should really make a better mask -- this is TBD!


  return
end
