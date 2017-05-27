pro make_superdark, files

	directory=getenv('DEIMOS_DATA')

        if n_elements(files) eq 0 then files=['/2002aug14/d0814_0084.fits','/2002aug14/d0814_0085.fits','/2002aug14/d0814_0086.fits']

	if strpos(files[0],'/') gt 0 then directory=directory+'/'
        
        searchstring = str_sep(files[0], '/')
        hasdate = strpos(searchstring, '200') eq 0
        whdate = where(hasdate, datect)
	if datect gt 0 then outfile='superdark.'+searchstring[whdate[0]]+'.fits' $
            else outfile = 'superdark.fits'


        datename=searchstring[whdate[0]]
        print,'processing: ',datename

        nfiles=n_elements(files)

	array=fltarr(2048,4096,nfiles)
	invvar=array+1

        header=headfits(directory+files[0])
	mwrfits, junk, outfile, header,/create ;dummy primary




	for i=1,8 do begin
		badmask=mrdfits(getenv('CALIB_DATA')+'/*bad*.Z',i)
		pixflat=mrdfits(getenv('CALIB_DATA')+'/proc*.gz',i)
		wh=where(pixflat eq -1,whct)
		if whct gt 0 then pixflat[wh]=0.
		mask=badmask*0
		if i eq 5 AND (datename eq  '2002aug14' $
                               OR datename eq '2002sep15') $
                               then mask[1320:1475,*]=1


		for j=0,nfiles-1 do begin
			tempim=deimos_read_chip(directory+files[j],i) > (-5)
                        
			tempim=tempim*(pixflat+1.)
			badmask=((badmask AND 1) eq 1) OR pixflat gt 0.25 OR pixflat lt -0.3
; fix the different BG level in frame 84
			if j eq 0 AND files[0] eq '/2002aug14/d0814_0084.fits'then tempim=tempim-0.75
			medtemp=djs_median(tempim,width=9,boundary='reflect')
			crpix = (tempim gt 70. and medtemp lt 2*median(medtemp)) $
				OR tempim gt 200. 
			tempim=djs_maskinterp(tempim, $
				crpix OR mask OR badmask OR $
				dilatemask(badmask OR crpix,2) OR $
				dilate(badmask,fltarr(9)+1), $
				iaxis=0)
			array[*,*,j]=tempim
			invvar[*,*,j]=1/((tempim>5)+2.3^2)*(1.-crpix)
	
		endfor
		if j eq 0 AND files[0] eq '/2002aug14/d0814_0084.fits'then invvar[*,*,0]=0.3*invvar[*,*,0]
		clipdark=avsigclip(array,invvar,1,5)
		medflux=djs_median(clipdark.flux,width=5,boundary='reflect')
		highmask=medflux gt 20 AND ((badmask AND 1) eq 0)
		whhigh=where(highmask AND NOT mask,highct)
		output=djs_median(medflux,width=9,boundary='reflect')
;		output=smooth(output,35,/edge_truncate)
                output=boxcar2d(output,[32,32],/ALL)
                medflux=clipdark.flux
		if highct gt 50 then output[whhigh]=(medflux)[whhigh]

                mkhdr,  header,  output
                mosaic_keywords,  header,  i
		output=float(output)
		mwrfits,floatcompress(output,ndig=8),outfile, header

	endfor

	spawn,'gzip -f '+outfile
return
end









