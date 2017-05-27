function peaklines, spec, yin, width=width, nsmooth=nsmooth
;+
; NAME:
;   peaklines
;
; PURPOSE:
;   to tune up the initial estimates of line centers from optical
;   model so that they better match the data!  Undetected lines flagged
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   yout = peaklines( spectrum, yin, width=width)
;
; INPUTS:
;   spectrum -- input spectrum for which we seek peaks
;   yin      -- list of centers of lines from optical model?
;
; KEYWORD PARAMETERS:
;    width = width -- search radius allowed to find peak (default =10 pixels)
;    nsmooth = nsmooth -- number of pixels to smooth function by before search
;
; OUTPUTS:
;   yout    --  improved line centers based on nearest high peak to incident.
;               Bad lines are flagged as yout=-1
;
; PROCEDURE:
;   spectrum is smoothed and then highest point within +- width is selected
;
; EXAMPLE:
;   arcline_x = peaklines(spec, arcline_x)
;
; MODIFICATION HISTORY:
;   md 12feb02
;   md 27may02 flag bad lines
;   dpf 17oct02 added nsmooth keyword
;-


  if NOT keyword_set(width) then width = 10
  if n_elements(nsmooth) EQ 0 then nsmooth = 5

  width=width > round(nsmooth*2./3.)
  yinf = round(yin)
  yout = yin
  if nsmooth GE 2 then temp = smooth(spec, nsmooth) else temp = spec
  nnn = n_elements(spec)-1
  
  for i=0,  n_elements(yin)-1 do begin ;find local peak
     tempj = temp[yinf[i]-width > 0:yinf[i]+width < nnn]
     nw = n_elements(tempj)
 ;    rms = stddev(tempj)
 ;    mmean = mean(tempj)
     tempjj = [1.0e6, tempj, 1.0e6] ;pad buffer with large values
     tl=shift(tempjj,1)
     tr=shift(tempjj,-1)
     tl=tl[1:nw] ;cut buffer back to original
     tr=tr[1:nw]
     ipeak = where(tempj ge tl and tempj ge tr, npeak) 
     if npeak eq 0 then begin
       yout[i] = -1        
;       print, 'deleting weak line: ', i
     endif else begin
; local peaks, but not end points, then select maximum of these
       if npeak eq 1 then ii = ipeak else begin
         xx = max(tempj[ipeak],ip) ;if multiple peaks, select maximum
         ii = ipeak[ip]
       endelse 
       yout[i] = yin[i] + ii-width
       if ((tempj[ii])[0] lt 1.5*min(tempj)) then begin ;weird that this construction is needed!
          yout[i] = -1.
;          print, 'deleting weak line: ', i
;        stop
       endif
     endelse
  endfor

;  window, 1
;  plothist, yin-yout, bin=.25

  return, yout
end








