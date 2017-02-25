;pro senscomp2,objfil,fluxfil,objfil2,fluxfil2,order
pro senscomp2,order

	objfil='Extract/Obj_mb0013.fits.gz'
	objfil2='Extract/Obj_mb0039.fits.gz'
	fluxfil='Extract/Sens_b0013.fits'
	fluxfil2='Extract/Sens_b0039.fits'
	flat=mrdfits('Flats/Tflt_SPEC_B_01.fits')

      objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)

      good = where(objstr[order].wave GT 2500.0)
      n = n_elements(objstr[order].wave)
      ord = good / n
      svset = xmrdfits(fluxfil, 1, /silent)
      x = 2.0d*((objstr[order].wave)[good] - (svset[order].wmin)[ord]) / $
                ((svset[order].wmax)[ord] - (svset[order].wmin)[ord]) - 1.0

  
      ncoeff = n_elements(svset[order].lcoeff)
      sens = exp(total(flegendre(x, ncoeff) * $
                 transpose((svset[order].lcoeff)[*, ord]),2))

;;;;;;FILE 2
      objstr2 = xmrdfits(objfil2, 1, STRUCTYP='mikeobjstrct', /silent)

      good = where(objstr2[order].wave GT 2500.0)
      n = n_elements(objstr2[order].wave)
      ord = good / n
      svset2 = xmrdfits(fluxfil2, 1, /silent)
      x = 2.0d*((objstr2[order].wave)[good] - (svset2[order].wmin)[ord]) / $
                ((svset2[order].wmax)[ord] - (svset2[order].wmin)[ord]) - 1.0

  
      ncoeff = n_elements(svset2[order].lcoeff)
      sens2 = exp(total(flegendre(x, ncoeff) * $
                 transpose((svset2[order].lcoeff)[*, ord]),2))

	point=sens[floor(n_elements(sens)/2)]
	point2=sens2[floor(n_elements(sens2)/2)]
	scale=point/point2
	sens2=sens2*scale


	tflat=flat[*,order]
	point3=tflat[floor(n_elements(tflat)/2)]
	scale3=point/point3
	tflat=tflat*scale

;;;;;;ANALYSIS
	plot,objstr[order].wave[good],alog10(sens),xr=[svset[order].wmin, $
						svset[order].wmax]
	oplot,objstr2[order].wave[good],alog10(sens2),color=255
	xyouts,(svset[order].wmax-svset[order].wmin)/2 + svset[order].wmin, $
        (max(alog10(sens)) - min(alog10(sens)))/2 + min(alog10(sens)), $
	'sens2 scaled by '+ string(scale)
	oplot,objstr[order].box_wv,alog10(tflat),color=-255

stop


return
end
