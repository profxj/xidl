;+
; NAME:
;   findpeaks
; 
; PURPOSE:
;   Find peaks in a spectrum 
;
; CALLING SEQUENCE:
;   peaky = findpeaks(y, mask, nsig=nsig)
;
; INPUTS:
;   y     - dependent variable array containing peaks
;   mask  - mask of pixels to ignore (0=bad, 1=good)
;
; KEYWORDS:
;   nsig  - threshhold for peaks in number of sigma-clipped sigma
;
; OUTPUTS:
;   bitmask: (1=peak, 0=nopeak)
;
; EXAMPLES:
;   sat = where y GT 60000
;   peaky         = findpeaks(y, (1B-sat), nsig=50)
;   peaklocations = where(peaky)
;
; COMMENTS:
;   To be a peak, a pixel must
;     1.  have a larger value than its neighbors
;     2.  have a value nsig sigma above the median
;     3.  not be at either end of the array
;
; REVISION HISTORY:
;
;       Wed Apr 24 22:13:38 2002, Doug Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
function findpeaks, y, mask, nsig=nsig

  if NOT keyword_set(nsig) then nsig = 10
  if NOT keyword_set(mask) then mask = 1B

  djs_iterstat,y,median=y0,maxiter=20
  ysig = djsig(y,maxiter=20)
  ny=n_elements(y)

; downward curvature and 50 sig detection
  yl = shift(y, -1)
  yr = shift(y, +1)

  d2 = (2*y-yr-yl)
  curve = (y GT y0+ysig*nsig) AND (d2 GE 0)
  pk = (y GE yr) AND (y GE yl)

  peaky = curve and pk and mask
  peaky[0] = 0B
  peaky[n_elements(peaky) -1] = 0B
  peakyinit=peaky
  whpeak=where(peaky ne 0,peakct)

; make sure this is a real local maximum (not just a fluctuation in
; vignetted region)
  
  for i=0,peakct-1 do begin
      localmedian=median(y[(whpeak[i]-80)>10:(whpeak[i]+80)<(ny-11)])
      localsig=djsig(y[(whpeak[i]-80)>10:(whpeak[i]+80)<(ny-11)])
      peaky[whpeak[i]]=y[whpeak[i]] gt localmedian+localsig*nsig
  endfor

  if total(peaky) lt 10 then peaky=peakyinit

  return, peaky
end
