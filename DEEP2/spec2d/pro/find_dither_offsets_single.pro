function find_dither_offsets_single,slitfile,Nexp,poff=poff,foff=foff,mfwhm=mfwhm

fin_fluxoff=fltarr(Nexp)
fin_offsets=fltarr(Nexp)
fin_fwhm = fltarr(Nexp)
fin_weight = fltarr(Nexp)

  filename = slitfile
  print,'Using Slit File =', filename

  pixpos=fltarr(Nexp)
  fluxoff=fltarr(Nexp)
  offsets=fltarr(Nexp,2)
  fwhm = fltarr(Nexp)
  weight = fltarr(Nexp)

  FOR I=0,(Nexp-1) DO BEGIN

   ;check if file exists
   goodfile= file_test(filename)
   if goodfile LT 1 THEN good=0
   if goodfile LT 1 THEN BREAK

   ext=I*2+1
   slit = mrdfits(filename,ext,/silent)

   mask=where(slit.mask GT 0)
   slit.flux(mask)=0
   profile = find_object(slit)
   Npix = size(profile,/DIMENSIONS)
  
   cprofile=profile(5:Npix-6) ; clip the slit edges
   Npix = size(cprofile,/DIMENSIONS)
   index = indgen(Npix)
   fitprofile = gaussfit(index,cprofile,fitparam)
   plot,cprofile
   oplot,fitprofile
;   print,'flux = ',fitparam(0)
   weight(I) = sigma(cprofile-fitprofile)

   fluxoff(I)=fitparam(0)
   pixpos(I)=fitparam(1)
   fwhm(I) = fitparam(2)*sqrt(2*alog(2))*0.1185
ENDFOR

FOR I=0,(Nexp-1) DO BEGIN
    offsets(I,0)=-1*(pixpos(I)-pixpos(0))
    offsets(I,1)=fluxoff(0)/fluxoff(I)
    
    fin_fluxoff(I)=offsets(I,1)
    fin_offsets(I)=offsets(I,0)
    fin_fwhm(I)=fwhm(I)
ENDFOR


 
maxflux = max(fin_fluxoff)

print,'Position Offsets =',fin_offsets
print,'Flux Scaling =', fin_fluxoff/maxflux
print,'FWHM (asec) =', fin_fwhm


poff=fin_offsets
foff=fin_fluxoff/maxflux
mfwhm=fin_fwhm

offsets=fltarr(Nexp,3)
offsets(*,0)=poff
offsets(*,1)=foff
offsets(*,2)=mfwhm

RETURN, offsets

end
