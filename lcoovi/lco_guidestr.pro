pro lco_guidestr, phot, mmin, mmax

  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'lco_guidestr, phot, magmin, magmax' 
    return
  endif 

 close, /all

 cutstr = where(phot.R LT mmax AND phot.R GT mmin, ncut)
 plot, phot[cutstr].xpix, phot[cutstr].ypix, psym=1
 openw,1,'Masks/guidestr.dat'
  for j=0,ncut-1 do begin
    printf, 1, FORMAT = '(f,f)', phot(cutstr(j)).xpix, phot(cutstr(j)).ypix
  endfor
 close, 1
 return

end
