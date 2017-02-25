pro order_width

x_psopen, 'slit_length.ps', /square, /color, /portrait
blue = mrdfits('data/OStr_B_02sep04.fits',1)
diff = blue.rhedg - blue.lhedg
hypotenuse = sqrt(1.0 + blue.arc_m^2) * diff
full = hypotenuse[1024,*]
mid = diff[1024,*]
up = diff[2047,*]
dn = diff[0,*]
djs_plot, blue[*].order, mid,yr=[18,22], xr=[38,110], /xs, $
      chars=1.5, xtitle='Echelle Order', $
                 ytitle='Slit Lengths (binned 2 pixels)'
errplot, blue[*].order, dn[*], up[*], color=djs_icolor('blue')
oplot, blue[*].order, full, ps=-1, color=djs_icolor('blue')

red = mrdfits('data/OStr_R_02sep04.fits',1)
red = red[2:32]
diff = red.rhedg - red.lhedg
hypotenuse = sqrt(1.0 + red.arc_m^2) * diff
full = hypotenuse[1024,*]
mid = diff[1024,*]
up = diff[2047,*]
dn = diff[0,*]
oplot, red[*].order, mid
errplot, red[*].order, dn[*], up[*], color=djs_icolor('red')
oplot, red[*].order, full, ps=-1, color=djs_icolor('red')

x_psclose
end

