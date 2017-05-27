function gfit_line, xcen, rad_in, spec, specivar

  rad = round(rad_in) ; integer
  g1 = {amp: 0., cen: 0., sig: 0., base: 0., mask: 0B}
  if size(spec, /n_dim) eq 1 then begin 
     nrow = 1
     nline = n_elements(xcen) 
  endif else begin 
     nrow = (size(spec, /dim))[1]
     ncol = (size(spec, /dim))[0]
     nline = (size(xcen, /dim))[1]
  endelse 

  xx = findgen(rad*2+1)-rad

  g = make_array(nrow, nline, value=g1)
  for i=0, nrow-1 do begin 
     for j=0, nline-1 do begin 
        xmid = round(xcen[i, j])
        minx = (xmid-rad) > 0
        maxx = (xmid+rad) <  (ncol-1)
        if total(specivar[minx:maxx, i] eq 0) eq 0 then begin 
           yfit = gaussfit(xx, spec[minx:maxx, i], a, $
                           estimates=[10000, 0, 1.5, 0], nterms=4)
           g[i, j].amp = a[0]
           g[i, j].cen = a[1]+xmid
           g[i, j].sig = a[2]
           g[i, j].base = a[3]
           
        endif 
     endfor 
  endfor 
  g.mask = (g.base LT 1000) AND (abs(g.sig) GT .75) AND (abs(g.sig) LT 4) AND g.amp GE 0

  return, g
end
