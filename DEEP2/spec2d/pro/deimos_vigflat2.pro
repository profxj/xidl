;+
; NAME:
;   deimos_vigflat2
;
; PURPOSE:
;   Generate vignetting correction flat - WORKING VERSION!
;
; CALLING SEQUENCE:
;   deimos_vigflat2, path
; 
; INPUTS:
;   path   - path for output file
;
; OUTPUTS:
;     file vignetting_flat.fits.gz in directory path.
;     This file should be put into the CALIB_DATA directory
;      (or simply downloaded from aquila.berkeley.edu)
;     Divide by this file to perform a vignetting correction for science 
;      data (NOT INTERNAL CALIBRATIONS!!!!!)
;
; MODIFICATION HISTORY:
;    21-Oct-2002 JAN - first version
;-

pro fitvignetting,chip,poly_level_to_coeff,fitlevels=fitlevels,file=file

; routine to fit contours of constant vignetting in one chip on a dome flat.
;  the code is slightly messy, but seems to work well.



; IF YOU WANT TO USE A DIFFERENT DOME FLAT, CHANGE THE FOLLOWING!!!!!!!
dir='/cosmo3/marc/DEIMOS_data/2002aug14/'
if n_elements(file) eq 0 then file='d0814_0047.fits'


lowvig=(chip LE 4)

fitlevels=findgen(21)*.025+0.45

if chip eq 4 then fitlevels=fitlevels(where(fitlevels gt 0.53 and fitlevels lt 0.87))
if chip eq 5 then fitlevels=fitlevels(where(fitlevels gt 0.45 and fitlevels le 0.925))
if chip eq 8 then fitlevels=fitlevels(where(fitlevels le 0.925))
nlevels=n_elements(fitlevels)


flat=deimos_read_chip(dir+file,chip)

; subtract off scattered light

			 badchip=mrdfits(getenv('CALIB_DATA') + $
				'/deimos_badmask.fits.Z',chip)

                         deimos_findmins,flat,badchip NE 0,planepars

                         flat=deimos_submins(flat,planepars)
                         
                         if min(flat) lt 0 then flat=flat-min(flat)+1.


if NOT lowvig then flat=rotate(flat,7)

flat=flat/median(flat)
nx=674
deltax=3
ncoeff=3-(chip eq 4); OR chip eq 5)
xcoords=findgen(nx)*deltax+5*deltax
transguess=xcoords*0-1.
allcuts=fltarr(nx,nlevels)
fits=dblarr(nlevels,3)-1

yindex=findgen(4096)

for i=0,nx-1 do begin
	column=flat[xcoords[i],*]
	column=reform(column)
	nderiv=40
	derivcols=findgen(nderiv)*100+100
	smder=deriv( (smooth(column,75))[derivcols])
	goodcol=max(smder) gt 1.5*abs(min(smder))
	if goodcol then begin
		neartransition=max(smder[0:15],whmax)
		realtransition=min(where(findgen(nderiv) gt whmax and smder lt 0.4*neartransition))*100
		transguess[i]=realtransition
; avoid things that cross slit boundaries
		column=median(column,5)
		if realtransition gt 300 then minvig = 0.6*realtransition else minvig=20
		maxvig=0.85*realtransition
		mincon=1.2*realtransition > 400.
		maxcon=1.6*realtransition > 700.			
		vigpars=linfit(yindex[minvig:maxvig],column[minvig:maxvig])
		vigfit=vigpars[0]+vigpars[1]*yindex
		continpars=linfit(yindex[mincon:maxcon],column[mincon:maxcon])
		continfit=continpars[0]+continpars[1]*yindex
		ratio=vigfit/continfit
		whinf=where(finite(ratio) eq 0,infct)
		if infct gt 0 then ratio[whinf]=0
	
		for j=0,nlevels-1 do begin
			allcuts(i,j)=min(where(ratio ge fitlevels[j] AND yindex lt (realtransition*1.+(fitlevels[j]-0.75)*.5)))
	
		endfor
	endif
	
endfor


;	splot,xcoords,lowvig*transguess+(1-lowvig)*(4095-transguess),psym=3,/nodata

	for j=0,nlevels-1 do begin
		cut=reform(allcuts(*,j))
		cutgood=(abs(cut-shift(cut,1)) lt 50) and abs(cut-shift(cut,-1)) lt 50 and cut gt 75. and xcoords gt 75 and xcoords lt 4020
		whok=where(cutgood,okct)
		if okct gt 40 then begin

			if NOT lowvig then cut=4095-cut
			poly_iter,xcoords[whok],cut[whok],2,2.,curve,coeff=coeff			
			fits[j,*]=coeff
;			soplot,xcoords[whok],curve
		endif else print,'bad cut',fitlevels[j]
	
	endfor
	smfits=fits

	goodfit=where(fits[*,0] ne -1,ngood)
	
	fitlevels=fitlevels[goodfit]
	fits=fits[goodfit,*]
	smfits=smfits[goodfit,*]

	poly_level_to_coeff=dblarr(3,3);ncoeff)
	
	poly_iter,fitlevels,fits[*,0],ncoeff-1,3.,curve,coeff=coeff

	poly_level_to_coeff[0,0:ncoeff-1]=coeff

	smfits[*,0]=curve
	poly_iter,fitlevels,fits[*,1],ncoeff-1,3.,curve,coeff=coeff
	smfits[*,1]=curve
	poly_level_to_coeff[1,0:ncoeff-1]=coeff
	poly_iter,fitlevels,fits[*,2],ncoeff-1,3.,curve,coeff=coeff
	smfits[*,2]=curve
	poly_level_to_coeff[2,0:ncoeff-1]=coeff

;	for j=0,ngood-1 do soplot,xcoords,poly(xcoords,smfits[j,*]),color=1
;	splot,fitlevels,fits[*,2],psym=4
;	soplot,fitlevels,curve	

	;splot,xcoords,transguess

end


pro deimos_vigflat2,path

xind=lindgen(2048,4096)
yind=xind / 2048
xind = xind MOD 2048

	if n_elements(path) eq 0 then path = ''

  nfiles = n_elements(files)

  if NOT keyword_set(path) then message, 'please call with a path!', /informational
  
  outfilename = path+'vignettingflat.fits'
 

  mwrfits, junk, outfilename, /create ;dummy primary
	for chip=1,8 do begin
		needcorr=(chip MOD 4 eq 1 or chip MOD 4 eq 0)

		if needcorr then begin

			 badchip=mrdfits(getenv('CALIB_DATA') + $
				'/deimos_badmask.fits.Z',chip)
			  nearvig=(badchip AND 4b) eq 4

			lowvig=(chip LE 4)
			leftvig=(chip MOD 4) eq 1
			fitvignetting,chip,coeff
			a=(coeff[2,2]*xind^2+coeff[1,2]*xind+coeff[0,2])
			b=(coeff[2,1]*xind^2+coeff[1,1]*xind+coeff[0,1])
			c=(coeff[2,0]*xind^2+coeff[1,0]*xind+coeff[0,0])-yind

			fofxy_1=(-b+sqrt(b^2-4*a*c))/(2*a)
			fofxy_2=(-b-sqrt(b^2-4*a*c))/(2*a)

			if min(a) eq 0. and max(a) eq 0. then begin
				fofxy_1=-c/b
				whinf=where(finite(fofxy_1) eq 0,infct)	
				
				vigflat=fofxy_1
				if infct gt 0 then vigflat[whinf]=1.
				vigflat=(vigflat < 1.) > 1.e-2
			endif else begin
				vigflat=( (fofxy_1*(lowvig eq 1)+fofxy_2*(lowvig eq 0)) < 1.) > 1.e-2
		
				if NOT lowvig AND leftvig then begin
					vigrange=fltarr(2048)
					for i=0,2047 do begin
						vigrange[i]=max(where(vigflat[i,*] eq 1.))
						if vigrange[i] gt -1 then vigflat[i,0:vigrange[i]-1]=1.
					endfor
				endif 

				if NOT leftvig then begin
					vigrange=fltarr(4096)
					for i=0,4095 do begin
						vigrange[i]=max(where(vigflat[*,i] eq 1.))
						if vigrange[i] gt -1 then vigflat[0:vigrange[i]-1,i]=1.
					endfor
				endif 

				badcalc=finite(vigflat) eq 0 OR (b^2-4*a*c) lt 0
				whinfh=where((badcalc OR vigflat lt 0) and NOT nearvig,infct)
				if infct gt 0 then vigflat[whinfh]=1.
				whinfl=where(badcalc and nearvig,infct)
				if infct gt 0 then vigflat[whinfl]=1.e-2

		endelse

		endif else vigflat=fltarr(2048,4096)+1
		 mkhdr, header, vigflat
	         mosaic_keywords, header, chip
	         mwrfits, floatcompress(float(vigflat)), outfilename, header; chip8
		
	endfor
	spawn,'gzip -f '+outfilename
return
end

