;+
; NAME:
;   deimos_arcwave_qa
;
; PURPOSE:
;   Look at arcwave_qa files
;
; CALLING SEQUENCE:
;   deimos_arcwave_qa, chip
; 
; INPUTS:
;   chip - chip
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   2002-Oct-25   D. Finkbeiner
;
;----------------------------------------------------------------------
pro deimos_arcwave_qa, chip

  flist = findfile('arcwave*.fits*')

  if keyword_set(chip) then i = chip-1 else $
    message, 'deimos_arcwave_qa, chip'

  print, 'Reading ', flist[i]
  dline = mrdfits(flist[i], 0, /silent)
  dlinesig = mrdfits(flist[i], 1, /silent)
  lamps = mrdfits(flist[i], 2, /silent)
  
  nslit = (size(dline, /dim))[1]
  rms = fltarr(nslit)
  for j=0, nslit-1 do begin 
     w = where(dlinesig[*, j] NE 0)     
     dlamp = dline[w, j]
     dlampsig = dlinesig[w, j]
     lambda = lamps[w].lambda
     foo = max(abs(dline[*, j]), indmax)
     
     if mean(lambda) lt 7800 then xrange = [6200, 8200] else $
       xrange = [7400, 9400]

     if mean(lambda) lt 6200 then xrange=[4000,8000]

     plot, lambda, dlamp, ps=1, yr=[-.1,.1]*2, $
       title='Arc wavelength residual', xtit='Ang', ytit='Delta [Ang]', $
       xr=xrange, /xst
     errplot, lambda, dlamp-dlampsig, dlamp+dlampsig
     oplot, [0, 1e4], [1, 1]*.01, line=1
     oplot, [0, 1e4], -[1, 1]*.01, line=1
;     oplot, chipedge[0]*[1, 1], [-1, 1], line=2
;     oplot, chipedge[1]*[1, 1], [-1, 1], line=2

     rms[j] = stdev(dlamp)*1000
     print, j, rms[j], max(abs(dlamp))*1000
     if max(abs(dlamp))*1000 GT 100 then begin 
        wait, 5
        print, 'outlier: ', lamps[indmax].lambda, lamps[indmax].element
     endif
     wait, .2
  endfor 

  print, 'Chip ', chip, '    Median RMS [mA]:', median(rms), '   Worst [mA]: ', max(rms)

  return
end
