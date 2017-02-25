pro bluecomb,bluefile,zlls=zlls,TITLE=title,OUTFIL=outfil

	spec=mrdfits(bluefile,1)
	norders=n_elements(spec.order)
	npix = norders * 2048

	velpix= 1.50d * 2
	cdelt = alog10(1.0d + velpix / 299792.458d)
	crval1=3000.0d
	wave0 = alog10(3000.0d)
	ioffset = floor((alog10(crval1)-alog10(3000.0d))/(cdelt))
	tot_wave = 10^(wave0 + (dindgen(npix)+ioffset)*cdelt)

	tot_flux = fltarr(npix)
	weight = fltarr(npix)
	sig=replicate(-1.,npix)

	for qq=(norders-1),0,-1 do begin
		indx=where(abs((tot_wave - spec[qq].wave[0])/spec[qq].wave[0]*299792.458) LT 0.1d, nindx)
		;non-zero flux points (good points)
		a = where(spec[qq].flux[*] NE 0.00000, ngd)

		;inverse wieghts for the good points
		wtmp = (spec[qq].sig[a])^2

;		tmp = tot_flux[indx+a[0]:indx+a[ngd-1]] + spec[qq].flux[a]/wtmp
		tmp=fltarr(ngd)
		for j=0L,ngd-1 do begin
			tmp[j] = tot_flux[indx+a[j]] + $
				spec[qq].flux[a[j]]/(spec[qq].sig[a[j]])^2
		endfor

		for i=0L,ngd-1 do begin
			tot_flux[indx+a[i]] = tmp[i]
			weight[indx+a[i]] = weight[indx+a[i]] + 1./wtmp[i]
		endfor
		;stop
	endfor

	a = where(weight NE 0.)
	tot_flux[a] = tot_flux[a] / weight[a]
	sig[a] = sqrt(1./weight[a])
	;good=where(tot_wave LT 4650.0,gcount)
	;plot,tot_wave[good],tot_flux[good],xrange=[3800.0,8500.0],yr=[0,20], $
	;     xstyle=1

	max10=floor(npix/15)
	
	tmpflux=fltarr(max10+1)
	tmpwav=fltarr(max10+1)
	for qq=0L,max10 do begin
		tmpflux[qq]=tot_flux[15*qq]
		tmpwav[qq]=tot_wave[15*qq]
	endfor
	good=where(tmpwav LT 4780.0,gcount)
	;clean=where(tmpflux GT 13,ccount)
	;tmpflux[clean]=0.0
	;plot,tmpwav[good],medsmooth(tmpflux[good],5),xrange=[3800.0,8900.0],yr=[0,20], $
	     ;xstyle=1
	;x_psopen,'specplot.ps'
	;set_plot,'ps'
	;;;;;plot,tmpwav[good],tmpflux[good],xrange=[3500.0,8900.0],yr=[0,1.3*max(tmpflux[good])], $
	;	xstyle=1,xtitle='Wavelength',ytitle='Flux'

;;LABELS
	if keyword_set(title) then begin
		xyouts,7400,0.95*max(tmpflux[good]),title,charsize=1.5
	endif

	;xyouts,7400,15,'V=17.3',charsize=1.5
	;xyouts,7400,14,'2400 seconds',charsize=1.5
	;xyouts,7400,13,'z(LLS)=3.5505',charsize=1.5


	;zlls=3.55046
	if keyword_set(zlls) then begin
		posvec=(1.0 + zlls)*[914.,1025.7,1215.67,1393.76,1548.20]
		namevec=['Ly-Limit','Ly-!7b!3 (1026)','Ly-!7a!3 (1216)', $
			'Si IV (1394)','C IV (1548)']
		yline=findgen(12)
		xline=fltarr(12)
		for qq=0,4 do begin
			xyouts,posvec[qq],0.6*max(tmpflux[good]), $
			namevec[qq],charsize=1.0, $
				orientation=90.
			xline[*]=posvec[qq]
			oplot,xline,yline,linestyle=1
		endfor
	endif
	
	;x_psclose
;	device,/close
;	set_plot,'x'

spec1d=temporary(tot_flux)
;; HEADER

if not keyword_set(OUTFIL) then outfil = 'out.fits'


;head = xheadfits(outfil)
mkhdr,head,spec1d
sxaddpar, head, 'CRVAL1', alog10(tot_wave[0])
sxaddpar, head, 'CDELT1', cdelt
sxaddpar, head, 'CRPIX1', 1
sxaddpar, head, 'CTYPE1', 'LINEAR'
sxaddpar, head, 'DC-FLAG', 1
sxaddpar, head, 'BITPIX', -32
sxaddpar, head, 'NAXIS', 1
sxaddpar, head, 'NAXIS1', n_elements(spec1d)
sxdelpar, head, 'NAXIS2'

mwrfits,spec1d,outfil,head,/create,/silent
;stop	
		
return
end
