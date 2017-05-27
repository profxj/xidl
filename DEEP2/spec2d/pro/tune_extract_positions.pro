function tune_extract_positions,bintab,maskname


color = ['B','R']

srctable = mrdfits(bintab,1)
brtobjects = where(srctable.mag lt 24.5 and srctable.mag gt 0 and srctable.objclass eq 'Program_Target      ')

Nbrtobjects = n_elements(brtobjects)-1

pos = fltarr(Nbrtobjects*2+2)
predpos = fltarr(Nbrtobjects*2+2)
weight = fltarr(Nbrtobjects*2+2)

I=0
FOR ifile=0, Nbrtobjects DO BEGIN
 FOR ic=0, 1 DO BEGIN


  slitno = brtobjects[ifile];
  good=1 ; assume good, flag if bad

  if slitno LT 1000 THEN slitname = '.' + strtrim(slitno,2)
  if slitno LT 100 THEN slitname = '.0' + strtrim(slitno,2)
  if slitno LT 10  THEN slitname = '.00' + strtrim(slitno,2)

  filename = 'slit.' + maskname + slitname + color(ic) + '.fits'
  print,filename

  goodfile= file_test(filename)
  if goodfile GT 0 THEN BEGIN

   slit = mrdfits(filename,1,/silent)

   mask=where(slit.mask GT 0)
   slit.flux(mask)=0

   p = total(slit.flux*slit.ivar,1)
   w = total(slit.ivar,1)
   profile = p/w
   mask = where( w LE 0)
   profile(mask) = 0

   Npix = size(profile,/DIMENSIONS)

   ;get slit dimensions
   slitsize=slit.rawsize
   dithersize=slit.dithersize/2.0

   index = indgen(Npix)

   fitprofile = gaussfit(index,profile,fitparam)
   noise = sigma(profile-fitprofile)

   plot,profile
   oplot,fitprofile
   print,noise   



   fwhm=fitparam(2)*sqrt(2*alog(2))*0.1185

   ;do some flagging
   ;first check the baseline
   if fitparam(3) GT 2 THEN good=0
   if fitparam(3) LT -2 THEN good=0
   ;check for positive flux
   if fitparam(0) LT 0 THEN good=0  
   ;check for good FWHM
   if fwhm GT 1.5 THEN good=0
   if fwhm LT 0.2 THEN good=0
   ;check if centriod off of slit
   if fitparam(1) LT dithersize THEN good=0
   if fitparam(1) GT NPIX - dithersize THEN good=0
  ENDIF

  IF goodfile LT 1 THEN good=0

  IF GOOD GT 0 THEN BEGIN  
   weight(I) = noise
   pos(I)=fitparam(1)
   predpos(I)=find_objpos(slitno,slitsize)+dithersize
   print,pos(I),predpos(I)
  ENDIF

  IF GOOD LT 1 THEN BEGIN
    weight(I) = -1
  ENDIF
 
  I+=1

 ENDFOR
ENDFOR

plot,pos,predpos,psym=3

RETURN,pos-predpos

end
