;;pro oh_lineplot

  data = mrdfits("FIRE_modelsky_R6000.fits", 0, hdr)

  naxis1 = sxpar(hdr,"NAXIS1")
  crval1 = sxpar(hdr,"CRVAL1")
  cdelt1 = sxpar(hdr,"CDELT1")

  wv = 10^((crval1+4) + cdelt1*findgen(naxis1))

  lines_per_micron = 54.5 / 1000.
  blaze_angle = 46.0 * 3.14159 / 180.

  set_plot, "ps"
  device, file='arc_key.ps', /landscape

  openw, 10, "FIRE_OH_R6000.lst"

  !p.multi=[0,1,3]

  colors=getcolor(/load)

  for order=31, 11, -1 do begin

     cw  = (2 / lines_per_micron * sin(blaze_angle)) / float(order) * 10000.
     fsr = cw / float(order)

     in = where(wv GT cw-0.8*fsr AND wv LT cw + 0.8*fsr)
     plot, wv, data, xrange=[cw-0.8*fsr, cw+0.8*fsr], /xsty, title='Order '+strtrim(order,2), yrange=[0,1.4*max(data[in])]

     x_fndpeaks, data[in], pks, peak=ypeak, pkwdth=4.0D, /thin, nsig=1.0

     npk = n_elements(pks)
     for ipk=0, npk-1 do begin
        plots, [wv[in[pks[ipk]]],wv[in[pks[ipk]]]], [data[in[pks[ipk]]]*1.05,data[in[pks[ipk]]]*1.30], color=colors.red
        xyouts, wv[in[pks[ipk]]], data[in[pks[ipk]]]*1.15, strtrim(wv[in[pks[ipk]]],2), orientation=90, charsize=0.5

        printf, 10, wv[in[pks[ipk]]], data[in[pks[ipk]]], "   OH"

     endfor

  endfor
  
  close, 10

  device, /close
  set_plot, "x"

end
