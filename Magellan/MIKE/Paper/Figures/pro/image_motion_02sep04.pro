
mike = mike_ar('data/mikestrct_02sep04.fits')
bstr = mrdfits('data/OStr_B_02sep04.fits',1)
rstr = mrdfits('data/OStr_R_02sep04.fits',1)

bblue = max(bstr.order)
bred  = min(bstr.order)

rblue = max(rstr.order)
rred  = min(rstr.order)

bb = where(mike.side EQ 1 AND total(mike.arc_xyoff,1) NE 0)
bbuniq = uniq(mike[bb].arc_xyoff[0])
bs = bb[bbuniq]
nb = n_elements(bs)
utb = lonarr(nb)-1
bblue_shift = fltarr(nb)
bred_shift = fltarr(nb)
for i=0,nb -1 do begin
  if mike[bs[i]].type EQ 'ARC' then utb[i] = mike[bs[i]].ut $
  else begin
    af = mike[bs[i]].arc_fil
    arc_frame = long(strmid(af, strpos(af, 'mb')+2,4))
    utb[i] = mike[where(mike.frame EQ arc_frame AND mike.side EQ 1)].ut
  endelse
  bblue_shift[i] = poly([bblue],mike[bs[i]].arc_xyoff)
  bred_shift[i] = poly([bred],mike[bs[i]].arc_xyoff)
endfor

rr = where(mike.side EQ 2 AND total(mike.arc_xyoff,1) NE 0)
rruniq = uniq(mike[rr].arc_xyoff[0])
rs = rr[rruniq]
nr = n_elements(rs)
utr = lonarr(nr)-1
rblue_shift = fltarr(nr)
rred_shift = fltarr(nr)
for i=0,nr -1 do begin
  if mike[rs[i]].type EQ 'ARC' then utr[i] = mike[rs[i]].ut $
  else begin
    af = mike[rs[i]].arc_fil
    arc_frame = long(strmid(af, strpos(af, 'mr')+2,4))
    utr[i] = mike[where(mike.frame EQ arc_frame AND mike.side EQ 2)].ut
  endelse
  rblue_shift[i] = poly([rblue],mike[rs[i]].arc_xyoff)
  rred_shift[i] = poly([rred],mike[rs[i]].arc_xyoff)
endfor

tflat_ut = 38350.

x_psopen, 'image_motion.ps', /color, /square, /portrait
djs_plot, (utb - 86400.*(utb GT 40000) -tflat_ut)/3600., bblue_shift*2, $
       yr=[-2,6], thick=8, /nodata, $
       xtitle='Time difference from trace flat calibration (hr)', $
       ytitle='Order shift in native CCD pixels', chars=1.5
djs_oplot, (utb - 86400.*(utb GT 40000)-tflat_ut)/3600., bblue_shift[0:10]*2, thick=12, color='blue'
djs_oplot, (utb - 86400.*(utb GT 40000)-tflat_ut)/3600., bred_shift[0:10]*2, color='blue'
djs_oplot, (utr - 86400.*(utr GT 40000)-tflat_ut)/3600., rblue_shift[0:10]*2, thick=12., color='red', lines=2
djs_oplot, (utr - 86400.*(utr GT 40000)-tflat_ut)/3600., rred_shift[0:10]*2, color='red', lines=2
x_psclose

end
