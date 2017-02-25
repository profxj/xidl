pro fwhm_and_atmos

restore, 'data/fracstuff.dat'
objred = mrdfits('data/Obj_mr0006.fits.gz',1)
objblue = mrdfits('data/Obj_mb0006.fits.gz',1)

x_psopen, 'fwhm_frac.ps', /color, /square, /portrait
djs_plot, ordr_b.order, frac_blue, xr=[35, 110], yr=[0,1], ps=4, $
      xtitle='Echelle Order', chars=1.5, /xs, $
      ytitle='Object centroid in Slit & FWHM (")'
djs_oplot, ordr_b.order, frac_blue_fit, color='blue' 
djs_oplot, ordr_r.order, frac_red, ps=4
djs_oplot, ordr_r.order, frac_red_fit, color='red' 

djs_oplot, objred.order, objred.spatial_fwhm, ps=1
djs_oplot, objblue[0:10].order, objblue[0:10].spatial_fwhm, ps=1
djs_oplot, objblue[11:*].order, objblue[11:*].spatial_fwhm

x_psclose
end
