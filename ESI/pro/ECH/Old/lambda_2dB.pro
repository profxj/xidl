;+
; NAME:
;       lambda_2d
;
; PURPOSE:
;       calculate the two dimensional wavelength function for an arclamp 
;       order of ESI echelle data. Requires the 
;	data to be rectified first. requires a linelist file with the 
;        following columns:
;	   wavelength (in angstroms)	intensity	good or not
;	program only uses 'good' lines
;
; CALLING SEQUENCE:
;      lambda_2d,arcim,model,wave,coeffarr,filename=filename,order=order(,navg=navg,residuals=residuals, peaks=peaks,fit=fit,ncfit=ncfit,print=print,bin=bin,smooth=smooth,nfit=nfit,fit=fit,ncfit=ncfit,print=print,reject=reject,range=range)
;
; INPUTS:
;	arcim 	- 2d arclamp image , one order, 4096 x 170
;	model 	- array to hold the wavelength vs pixel number array, 4096 x 170
;	wave 	- initial model,  wavelength vs pixel number array, 4096 x 10, I got the coefficents off the 
;		keck ESI website and calculated...
;	coeffarr- the array to hold the legendre coefficents of the new fitted wavelength scale, numcoeff x4096
;	navg	- number of orws to average the arclamp image over when calculating the fit. 
;
;
; OPTIONAL INPUTS:
;
;      
;
; REQUIRED KEYWORDS:
;      
;	filename- name for the linelist file
;	order - echellette order number
;
; OPTIONAL KEYWORDS: 
;	navg - number of rows to average
;	residuals - if residuals ne 0 then lammbda_2d will plot residuals as it goes..
;	peaks	-  if peaks ne 0 then lammbda_2d will plot peaks as it goes.
;	bin - search window = =/- bin in lambawpeak
;	smooth - smooths the lambda model in the lambda direction only
;	nfit - order for the legendra fit
;       fit - does a linear fit to  each COEFFICIENT as a function of slit height (relative to the ensemble of coeffs)
;	ncfit - order for the wavelength model coefficeint fit ncfit=1 = linear
;	print - print out diagnostics
;	reject - reject fits greter than 1 sigma from mean      
;       range - the number of rows to calculate a lambda model over   
;
; OUTPUTS:
;   
;
; OPTIONAL OUTPUTS:
;   
;   
; COMMENTS:
;
;
; EXAMPLES:
;
;
; BUGS:
;
;
; PROCEDURES CALLED:
;
;	lambdawpeak
;	lamp
;	svdfit	   
;
; REVISION HISTORY:
;    written, 30-Aug-01, Andrew Sheinis, Lick observatory, UCSC
;   revised,19 mar 2002 , AIS. calculate coeffs first then fit them to a low order poly as a function of slit height
;   revised 29 march 2002 , AIS, add reject feature to  coeffs fit...
;-
;--------------------------------------------------------------------------- 



pro lambda_2dB,arcim,model,wave,coeffarr,navg=navg,residuals=residuals,peaks=peaks,filename=filename,stop=stop,order=order,bin=bin,smooth=smooth,nfit=nfit,fit=fit,ncfit=ncfit,print=print,reject=reject,range=range

;initailize********************

   if NOT keyword_set(navg) then navg = 1
   if NOT keyword_set(residuals) then residuals = 0
   if NOT keyword_set(peaks) then peaks = 0
   if NOT keyword_set(filename) then print,'error, no lamp filename entered'
   if NOT keyword_set(stop) then stop = 0
   if NOT keyword_set(order) then print,'no order selected'
   if NOT keyword_set(bin) then bin = 7
   if NOT keyword_set(smooth) then smooth = 0
   if NOT keyword_set(nfit) then nfit = 6
   if NOT keyword_set(ncfit) then ncfit = 1
   if NOT keyword_set(fit) then fit = 0
   if NOT keyword_set(print) then print = 0
   if NOT keyword_set(range) then range = 100

smsz=15			;smooth size


if print ne 0 then print,'nfit =  ',nfit

coeff=fltarr(7)
model=fltarr(4096,range)
coeffarr=fltarr(nfit,range)
arctemp=fltarr(4096)

;create wavelength structure*****

	;order number = 15-i

i=15-order


xl=-1+2*findgen(4096)/4096.
yl=wave(*,i)

lcoeff=svdfit(xl,yl,4,/legendre)
wset={func:'legendre',xmin:0,xmax:4096,coeff:lcoeff}

;read in arc and linelist*************

;.run lamp_esi
;.run lambdawpeak_esi

pixel=xl



lampsold=lamp(filename,  quality='GOOD')

for j=0,range-1 do begin ; run through each central row, call lambdawpeak and fit a model to the wavelength scale

if (print ne 0) then   print,'j =  ',j


  ;couple of differant averaging schemes. for the arc rows********************************************


  if (j ge navg)  then arc=rebin(arcim(*,i*180+(j-navg):i*180+j),4096,1) else arc=rebin(arcim(*,i*180+j:i*180+(j+navg)),4096,1) & if (print ne 0) then  print,'averaging',navg,'   rows'


  ;find peaks******************************************
  lambdawpeak, arc, wset,lampsold, 0,xstart,lamps, parab=1,bin=bin
  if peaks ne 0 then   plot,arc,yrange=[0,200] & if peaks ne 0 then oplot,xstart,replicate(50,n_elements(xstart)),psym=5

  ;make a model***************************************


  ;or...........

  ;yl=model from keck website
  ;yy=actual values found from arc image using lambdawpeak
  ;yyy=model I calculate from yy using svdfit
  ;yt=my model values at yy locations i.e. yt-yy = residuals
  ;yyt=double check that my model values at yy locations

  
  xx=-1+2*xstart/4096
  yy=lamps.lambda
  coef=svdfit(xx,yy,nfit,/legendre,/double,yfit=yt)
  xxx=-1+2*findgen(4096)/4096.
  xxx=double(xxx)
  coeff(0:nfit-1)=coef
  


;  yyy=coeff(0)+xxx*coeff(1)+xxx^2*coeff(2)+xxx^3*coeff(3)+xxx^4*coeff(4)+xxx^5*coeff(5) ;polynomials
;  yyt=coeff(0)+xx*coeff(1)+xx^2*coeff(2)+xx^3*coeff(3)+xx^4*coeff(4)+xx^5*coeff(5)

;  yyy=coeff(0)*(nfit ge 1)+coeff(1)*(nfit ge 2)*xxx+coeff(2)*(nfit ge 3)*(-.5+1.5*xxx^2)+coeff(3)*(nfit ge 4)*(-1.5*xxx+2.5*xxx^3)+coeff(4)*(nfit ge 5)/8*(3-30*xxx^2+35*xxx^4)+coeff(5)*(nfit ge 6)/8*(15*xxx-70*xxx^3+63*xxx^5)+coeff(6)*(nfit ge 7)/16*(-5+105*xxx^2-315*xxx^4+231*xxx^6); legendre polynomials

  yytt=coeff(0)*(nfit ge 1)+coeff(1)*(nfit ge 2)*xx+coeff(2)*(nfit ge 3)*(-.5+1.5*xx^2)+coeff(3)*(nfit ge 4)*(-1.5*xx+2.5*xx^3)+coeff(4)*(nfit ge 5)/8*(3-30*xx^2+35*xx^4)+coeff(5)*(nfit ge 6)/8*(15*xx-70*xx^3+63*xx^5)+coeff(6)*(nfit ge 7)/16*(-5+105*xx^2-315*xx^4+231*xx^6); legendre polynomials

ord=nfit-1
yyy = coeff(0)*(ord ge 0)+$
coeff(1)*legendre(xxx, 1)*(ord ge 1)+ $
coeff(2)*legendre(xxx, 2)*(ord ge 2)+ $ 
coeff(3)*legendre(xxx, 3)*(ord ge 3)+ $
coeff(4)*legendre(xxx, 4)*(ord ge 4)+ $
coeff(5)*legendre(xxx, 5)*(ord ge 5)

 
  xxnorm=findgen(4096)



   if residuals ne 0 then plot,xstart,.153*(yt-yy),psym=5,title='residuals',xtitle='pixels',ytitle='angstroms',subtitle=j  & if residuals ne 0 then xyouts,xstart,.153*(yt-yy),alignment=0.5,string(lamps.lambda),charsize=1.5  ;& if residuals ne 0 then wait,.2


  model(*,j)=yyy
  coeffarr(*,j)=coef

if (print ne 0) then   print,'coeff ='
if (print ne 0) then   print, coeff


endfor

;we can take the ensemble of coeffs and fit them to a function, then recalcualet the lambdascale based on that function

if fit ne 0 then begin
	
	for i = 0,nfit-1 do begin
	  xx=findgen(range)
	  yy=coeffarr(i,*)
	  momt=moment(yy)
	  mean=momt(0)
	  stdev=(momt(1))^0.5
	  if reject ne 0 then begin

             xx=xx(where(yy le mean+stdev and yy ge mean-stdev))
	     yy=yy(where(yy le mean+stdev and yy ge mean-stdev))
	  endif
	  cfit=poly_fit(xx,yy,ncfit,/double,yfit=yfit)
	  xx=findgen(range)
	  dfit=fltarr(4)
	  dfit(0:ncfit)=cfit
	  coeffarr(i,*)=dfit(0)*(ncfit ge 0)+$
	  dfit(1)*xx*(ncfit ge 1)+ $
	  dfit(2)*xx^2*(ncfit ge 2)+ $
	  dfit(3)*xx^3*(ncfit ge 3)
;	  coeffarr(i,*)=yfit
	endfor
   	

	for j=0,range-1 do begin
	  ord=nfit-1
	  doeffarr=fltarr(6,range)
	  doeffarr(0:ord,*)=coeffarr	
	  model(*,j) = doeffarr(0,j)*(ord ge 0)+$
	  doeffarr(1,j)*legendre(xxx, 1)*(ord ge 1)+ $
	  doeffarr(2,j)*legendre(xxx, 2)*(ord ge 2)+ $ 
	  doeffarr(3,j)*legendre(xxx, 3)*(ord ge 3)+ $
	  doeffarr(4,j)*legendre(xxx, 4)*(ord ge 4)+ $	  
	  doeffarr(5,j)*legendre(xxx, 5)*(ord ge 5)
	endfor

endif




return

end








