;+
;pro scatter_model2,image,scattermod
;
;;
; NAME:
;       scatter_model2
;
; PURPOSE:
;      
;this routine uses the inter-order intensties to fit a scattered light model. The fit uses a 2-dimensional 
;b_spline, which calculates a B-spline in the least-squares sense 
;     based on two variables: x which is sorted and spans a large range
;                                 where bkpts are required
;               and           x2 which can be described with a low order
;                                 polynomial    
;
;
;then subtracts that model from the data.
;it expects, trimmed, bias subtratcted data, rotated 90 degrees and NOT rectified.
;
; CATEGORY:
;       ESI reduction
;
; CALLING SEQUENCE:
;       scatter_model,image,imageout,scattermod
; INPUTS:
;	image = the image to model scattered light from
; KEYWORD PARAMETERS:
;      none
; OUTPUTS:
;       imageout =image - scatter model
;	scattermod = the inra-order scattered light model
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;       don't mix with alcohol or barbituates
; RESTRICTIONS:
;       None.
; MODIFICATION HISTORY:
;
;	written by A. Sheinis lick Observatory 7/26/01
; PROCEDURES CALLED:
;	bspline_iterfit
;	bspline_valu
;	
; PROCEDURE:
;
pro scatter_model2,image,scattermod
;
;
;- 

;I went through by hand and found the top and bottoms of each order....

xblue=[1750L,1350,1280,1150,880,490,100,10,10,10];blue edge
xgreen=replicate(2048L,10) ;in the middle
xred=[3100L,3300,3500,3650,4086,4086,4086,4086,4086,2940] ; boundary at red edge

;Y position at the blue  end of the order

ybluelow=[130L,387,638,864,1055,1216,1352,1508,1668,1826]
ybluehigh=[285L,535,782,1004,1190,1348,1484,1637,1796,1949]

; y position of center of order, at x=2048
ygreenlow=[135L,406,644,860,1054,1230,1390,1537,1686,1836]
ygreenhigh=[291L,561,794,1005,1196,1364,1517,1667,1807,1951]


;y position  at the red end of the order

yredlow=[64L,284,466,630,686,849,998,1142,1286,1786]
yredhigh= [230L,442,618,775,832,988,1134,1273,1413,1904]


;make sure we only look in between the orders
dn=6

yredlow=yredlow-dn
ygreenlow=ygreenlow-dn
ybluelow=ybluelow-dn

yredhigh=yredhigh+dn
ygreenhigh=ygreenhigh+dn
ybluehigh=ybluehigh+dn





yy=findgen(2048L)
xx=findgen(4096L)
yylow=fltarr(10,4096L)
yyhigh=fltarr(10,4096L)

scatterdat=0*image+1

for i=0,9 do begin ;fit the orders


     ;do a quadratic fit to roughly find the orders
	
    rufyhigh=[ybluehigh(i),ygreenhigh(i),yredhigh(i)]
    rufylow=[ybluelow(i),ygreenlow(i),yredlow(i)]
    rufxpos=[xblue(i),xgreen(i),xred(i)]

 ;   xyouts,rufxpos/4,rufyhigh/4,'+',/device,alignment=0.5
 ;   xyouts,rufxpos/4,rufylow/4,'o',/device,alignment=0.5


    cf_hi=poly_fit(rufxpos,rufyhigh,2)  ;quadratic fit coefficients
    cf_low=poly_fit(rufxpos,rufylow,2)

    yyhigh(i,*)=cf_hi(0)+cf_hi(1)*xx+cf_hi(2)*xx^2
    yylow(i,*)=cf_low(0)+cf_low(1)*xx+cf_low(2)*xx^2 
 
	
	
;    xyouts,xx/4,yyhigh(i,*)/4,'-',/device,alignment=0.5
;    xyouts,xx/4,yylow(i,*)/4,'-',/device,alignment=0.5

    
yylow=yylow * (yylow ge 0 and yylow le 2047)
yyhigh=yyhigh * (yyhigh ge 0 and yyhigh le 2047)

endfor

for i=0,4095 do begin
  for j=0,9 do begin
   scatterdat(i,yylow(j,i):yyhigh(j,i))=0.0
  endfor
endfor



scatterdat=scatterdat* image

npix=1024		;processor can't handle the whole image so well, so rebin to a smaller image


scatterdat=rebin(scatterdat,npix,npix/2.)

x=findgen(npix)
x=rebin(x,npix,npix/2.,/sample)

x2=findgen(npix/2)
x2=transpose(x2)
x2=rebin(x2,npix,npix/2.,/sample)



x=x(where(scatterdat gt 0))
x2=x2(where(scatterdat gt 0))
scatterdat=scatterdat(where(scatterdat gt 0))
temp=fltarr(npix,npix/2.)
scattermod=fltarr(npix,npix/2.)

;everyn= something like number of x elements of scatterdat /(10 or 20).

;everyn=fix(npix/10.)
  
sset = bspline_iterfit(x2, scatterdat,x2=x, invar=ivar, everyn=25000, yfit=yfit,npoly=4)

for i=0L,n_elements(x)-1  do begin
;	scattermod(x(i),x2(i))=yfit(i)
	temp(x(i),x2(i))=scatterdat(i)	

endfor

;x=findgen(4096)/16.
;x=rebin(x,4096,2048)

;x2=findgen(2048)/16.
;x2=transpose(x2)
;x2=rebin(x2,4096,2048)

x=findgen(npix)
x=rebin(x,npix,npix/2.,/sample)

x2=findgen(npix/2)
x2=transpose(x2)
x2=rebin(x2,npix,npix/2.,/sample)


scattermod  = bspline_valu( x2, sset, x2=x)



scattermod=rebin(scattermod,4096,2048)

scatterdat=rebin(temp,4096,2048,/sample)

imageout=image-scattermod



return

end

