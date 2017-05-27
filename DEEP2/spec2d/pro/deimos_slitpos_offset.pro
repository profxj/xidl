
; mask  - mask definition structure
; X direction is the spatial direction in this routine

; COMMENTS:
;   We do NOT need to know the chip number (or numbering scheme!)
;    in order to find the correct offset.  Just try everything...
;   We should probably throw out bad "found" slits very carefully
;    before calling this routine to avoid getting fooled. 
; D. Finkbeiner 23 Aug 2001
function deimos_slitpos_offset, xfound, xpred

  nfound = n_elements(xfound)
  ind = lonarr(nfound)

; -------- Generate set of offsets to probe
  step = 5
  minx = min(xpred, max=maxx)
  maxf = max(xfound)
  xoffs = findgen((maxx-minx+maxf)/step)*step+minx-maxf

; -------- Store results in sdev
  sdev = fltarr(n_elements(xoffs))

; -------- Loop over chips to find best match

  for j=0, n_elements(xoffs)-1 do begin 
;    for each FOUND slit, get the nearest predicted slit position

     for i=0, nfound-1 do begin
        junk = min(abs((xfound[i]+xoffs[j])-xpred), indi)
        ind[i] = indi
     endfor
     sdev[j] = stdev(xfound-xpred[ind]) ; big if match is bad
  endfor

  jmean = mean(where(sdev EQ min(sdev)))

  return, xoffs[jmean]
end

