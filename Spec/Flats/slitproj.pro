
function slitproj, ystart, ywave, arc_slope, ordr_str

      nrow = n_elements(ordr_str.lhedg)
      ycol = dindgen(nrow)
      oo = 0.5*(ordr_str.lhedg + ordr_str.rhedg)
    

      xy2traceset, ycol[1:*]-0.5, double(oo[1:*] - oo), oset, ncoeff=4, /silent
      traceset2xy, oset, ywave, order_slope  

      length = ordr_str.rhedg - ordr_str.lhedg
      xy2traceset, ycol, length, lset, ncoeff=4, /silent
      traceset2xy, lset, ywave, slit_length

      trouble = where(arc_slope* order_slope GE 0.5, nt)
      if nt GT 0 then begin
         print, 'WARNING: Slopes are very large, check order traces and and arc line slopes'
      endif

      slit_proj = slit_length / (1.0 - arc_slope*order_slope)

      return, slit_proj
end


