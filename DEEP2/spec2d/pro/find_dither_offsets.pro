function find_dither_offsets,filename,Nexp,offset

zeropos=offset/0.1185;

pixpos=fltarr(Nexp)
fluxoff=fltarr(Nexp)
offsets=fltarr(Nexp,2)

FOR I=0,(Nexp-1) DO BEGIN
  ext=I*2+1
  slit = mrdfits(filename,ext)
  mask=where(slit.mask GT 0)
  slit.flux(mask)=0
  profile = find_object(slit)
  Npix = size(profile,/DIMENSIONS)
  
  cprofile=profile(5:Npix-6) ; clip the slit edges
  Npix = size(cprofile,/DIMENSIONS)

;  print,Npix,size(cprofile,/DIMENSIONS)

  index = indgen(Npix)
  fitprofile = gaussfit(index,cprofile,fitparam)
  plot,cprofile
  oplot,fitprofile
  print,'flux = ',fitparam(0)
  fluxoff(I)=fitparam(0)
  pixpos(I)=fitparam(1)
ENDFOR

;print,pixpos

bestexp = max(fluxoff)

FOR I=0,(Nexp-1) DO BEGIN
  offsets(I,0)=-1*(pixpos(I)-pixpos(0)-zeropos)
  offsets(I,1)=bestexp/fluxoff(I)
ENDFOR

print,offsets

RETURN, offsets

end
