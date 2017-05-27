function find_dither_offsets_auto,bintab,maskname,Nexp,poff=poff,foff=foff,mfwhm=mfwhm

color = ['B','R']



srctable = mrdfits(bintab,1)

brtobjects = where(srctable.mag lt 24 and srctable.mag gt 0 and srctable.objclass eq 'Program_Target      ')

Nbrtobjects = n_elements(brtobjects)-1

if Nbrtobjects  LT 1 THEN brtobjects = where(srctable.objclass eq 'Program_Target      ')

Nbrtobjects = n_elements(brtobjects)-1

fin_fluxoff=fltarr(Nexp)
fin_offsets=fltarr(Nexp)
fin_fwhm = fltarr(Nexp)
fin_weight = fltarr(Nexp)

FOR ifile=0, Nbrtobjects DO BEGIN
 FOR ic=0, 1 DO BEGIN


  slitno = brtobjects[ifile];
  good=1 ; assume good, flag if bad

  if slitno LT 1000 THEN slitname = '.' + strtrim(slitno,2)
  if slitno LT 100 THEN slitname = '.0' + strtrim(slitno,2)
  if slitno LT 10  THEN slitname = '.00' + strtrim(slitno,2)

  filename = 'spSlit.' + maskname + slitname + color(ic) + '.fits'
  print,filename

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
;   plot,cprofile
;   oplot,fitprofile
;   print,'flux = ',fitparam(0)
   weight(I) = sigma(cprofile-fitprofile)

   fluxoff(I)=fitparam(0)
   pixpos(I)=fitparam(1)
   fwhm(I) = fitparam(2)*sqrt(2*alog(2))*0.1185

   ;do some flagging
   ;first check the baseline
   if fitparam(3) GT 2 THEN good=0
   if fitparam(3) LT -2 THEN good=0
   ;check for positive flux
   if fitparam(0) LT 0 THEN good=0
   ;check for good FWHM
   if fwhm(I) GT 1.5 THEN good=0
   if fwhm(I) LT 0.2 THEN good=0
   ;check if centriod off of slit
   if pixpos(I) LT 3 THEN good=0
   if pixpos(I) GT NPIX - 3 THEN good=0

   if good EQ 0 THEN BREAK
 
  ENDFOR


  IF GOOD GT 0 THEN BEGIN
   FOR I=0,(Nexp-1) DO BEGIN
    offsets(I,0)=-1*(pixpos(I)-pixpos(0))
    offsets(I,1)=fluxoff(0)/fluxoff(I)

    w = 1.0/(weight(I)*weight(I))

    fin_fluxoff(I)+=offsets(I,1)*w
    fin_offsets(I)+=offsets(I,0)*w
    fin_fwhm(I)+= fwhm(I)*w
    fin_weight(I)+=w
   ENDFOR
 
   print,weight
   print,offsets(*,0)
   print,offsets(*,1)
   print,fwhm
  ENDIF

 ENDFOR
ENDFOR


maxflux = max(fin_fluxoff/fin_weight)

print,fin_offsets/fin_weight
print,(fin_fluxoff/fin_weight)/maxflux
print,fin_fwhm/fin_weight


poff=fin_offsets/fin_weight
foff=(fin_fluxoff/fin_weight)/maxflux
mfwhm=fin_fwhm/fin_weight

offsets=fltarr(Nexp,3)
offsets(*,0)=poff
offsets(*,1)=foff
offsets(*,2)=mfwhm

RETURN, offsets

end
