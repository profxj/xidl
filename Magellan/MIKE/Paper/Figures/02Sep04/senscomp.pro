pro senscomp


eg274_1=mrdfits('Extract/Obj_mb0003.fits.gz',1)
sens_274_1=mrdfits('Extract/Sens_b0003.fits',1)
good_274_1 = where(eg274_1.wave GT 2500.0)

n = n_elements(eg274_1[0].wave)
ord = good_274_1 / n

x274_1 = 2.0d*((eg274_1.wave)[good_274_1] - (sens_274_1.wmin)[ord]) / $
	((sens_274_1.wmax)[ord] - (sens_274_1.wmin)[ord]) - 1.0
ncoeff_274_1 = n_elements(sens_274_1[0].lcoeff)

s_274_1 = exp(total(flegendre(x274_1, ncoeff_274_1) * $
                 transpose((sens_274_1.lcoeff)[*, ord]),2))
stop
eg274_2=mrdfits('Extract/Obj_mb0004.fits.gz',1)
eg274_3=mrdfits('Extract/Obj_mb0009.fits.gz',1)






return
end
