pro findcoeff,image,coeff

;this routine finds the maximum of a series of 5 star spectra per order and fits
;those maxima to a polynomial. It is used to input coeff to a_maporders.pro to
;rectify ESI images....note we've rotaed the esi images to look like lris and diemos...
;x and y are reversed

;modified 4/10/01..to use more thqan three spectra pointas per fit and therefor higher
;order than quadratic.

 
;plot
window,1,xsize=1024,ysize=512
tv,rebin(image,1024,512)/100

;tv,image/100

lns=10; 1/2 the number of columns to average
;arrays holding the positions of the maxima at each of 5 positions, 10 orders

center = fltarr(10,5,20);10 orders, 5 stars per order, 20 points per star


coeff=fltarr(4,10,5)  ;4 polynomialfit coefficients for 10 orders, 20 y positiosn, 5 x positions within an order

ypos=findgen(4096) ; y position variable

;yfitpos=findgen(20)*3095/20+1000 ;20 points, from 1000 to 4095
yfitpos=(findgen(20)+1)*4095/21

;I went through by hand and found the top and bottoms of each order....

ylow=[1750,1750,1650,1400,1100,925,600,200,10,10];blue edge
ymid=replicate(2048,10) ;in the middle
yhigh=[3400, 3700, 3500, 4000, 4000, 4000, 4000, 4000, 4000, 2500] ;y boundary at red edge

 ;at blue edge, beginning of order
;xlow is the bottom of the slit
;xhigh is the top of the slit

xlowbeg=[110,400,648,866,1060,1240,1394,1526,1664,1822]
xhighbeg=[294,566,806,1018,1208,1384,1536,1664,1800,1958]


; x position of exiting edge of order, at y=2048
xlowmid=[128,402,644,860,1054,1228,1388,1536,1682,1832]
xhighmid=[302,568,802,1012,1202,1374,1528,1676,1816,1960]

;at the red end of the order
xlowend=[10,195,463,535,709,873,1017,1165,1305,1781]
xhighend= [184,363,622,691,861,1013,1165,1305,1443,1907]




for i=1,9 do begin ;find centers for all three y positions over all orders
  ;at 5 positions within an order.

     ;do a quadratic fit to roughly find the orders
	
    rufxhigh=[xhighbeg(i),xhighmid(i),xhighend(i)]
    rufxlow=[xlowbeg(i),xlowmid(i),xlowend(i)]
    rufypos=[ylow(i),ymid(i),yhigh(i)]

    yfithigh=poly_fit(rufypos,rufxhigh,2)  ;quadratic fit coefficients
    yfitlow=poly_fit(rufypos,rufxlow,2) 

 



  for j=0,4 do begin ;go star by star...

  xhigh= yfithigh(0)+ yfithigh(1)*yfitpos +yfithigh(2)*yfitpos*yfitpos
  xlow= yfitlow(0)+ yfitlow(1)*yfitpos +yfitlow(2)*yfitpos*yfitpos
xyouts,yfitpos/4,xlow/4,'*',/device,alignment=0.5
xyouts,yfitpos/4,xhigh/4,'o',/device,alignment=0.5


  
       for k=0,19 do begin; 20 points along the spectra , onestar, one order
 
   	;cut out a small piece of the order in spatial direction
	;qq=the offset within each slitlet...

	qq=0
	
    	dx=round((xhigh(k)-xlow(k)+2*qq)/5.) 

    	ystrip=image((yfitpos(k)-lns):(yfitpos(k)+lns),(xlow(k)+j*dx +qq):(xlow(k)+(j+1)*dx-1+qq));pick y
   	ystrip=rebin(ystrip,1,n_elements(ystrip(0,*)));rebin to 1 x strip
	ystrip=smooth(ystrip,3)




    	xstrip=findgen(dx)+xlow(k)+(j*dx);pick x


	;my usual trick, blow up the image by 1000, find the max of the expanded imagetv,rebin(image,1024,512)/100
    	ytemp=rebin(ystrip,1,(1000*dx))
    	xtemp=rebin(xstrip,1000*dx,1)
	temp=(max(ytemp))
	temp=where(ytemp ge temp)
	center(i,j,k)=congrid(xtemp(temp),1)

	centery=congrid(ytemp(where(max(ytemp))),1)


;	diagnostic plots	
;	plot,findgen(5*dx)+xlow(k),image(yfitpos(k),xlow(k):xlow(k)+5*dx-1)
;	oplot,xstrip,0.9*ystrip


;	xyouts,center(i,j,k),centery,'+',alignment=0.5
;	wait,0.5



      endfor

    	;take the 20 center points for each star in each order...

    	yy=yfitpos
    	xx=center(i,j,*)

    	;fit to a ploynomial
    	a=poly_fit(yy,xx,3)  ;quadratic fit coefficients
    	coeff[*,i,j]=a
    	;plot out the fits


xyouts,yy/4,xx/4,'.',/device,alignment=0.5



 


  endfor
endfor


return

end
