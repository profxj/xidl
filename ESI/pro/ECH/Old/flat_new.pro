; NAME:
;       flat_new
;
; PURPOSE:
;      
;
;	This routine flattens the flat. It extracts the full order 
;       based on a pre-determined order
;	location (maybe something more flexible can be placed here, 
;       like unsharp masking) and 
;	then averages (should median?) the order down to a single slit
;       location. It then smooths these 
;	10 order averages and subtracts from the original flat field. 
;
;	NOTE: works remarkably well, excpet for one problem:
;	tends to amplify  ERRORS in DEBIASED areas, mostly the red end
;	of orders 10 and 11.
;
; CATEGORY:
;       ESI reduction
;
; CALLING SEQUENCE:
;	flat_new,flat,flatout

; INPUTS:
;	flat - flat field image

; KEYWORD PARAMETERS:
;      none
; OUTPUTS:
;       flatout =flattened flat
;
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;       don't mix with alcohol or barbituates
; RESTRICTIONS:
;       None.
; MODIFICATION HISTORY:
;
;	written by A. Sheinis lick Observatory 3/20/2002
; PROCEDURES CALLED:
;	smooth
;	poly_fit
;	
; PROCEDURE:

pro flat_new,flat,flatout	


flatin=flat
flatout=flatin*0.0
flattmp=flatin*0.0

;I went through by hand and found the top and bottoms of each order....

xblue=[1750,1350,1280,1150,880,490,100,10,10,10];blue edge
xgreen=replicate(2048,10) ;in the middle
xred=[3100,3300,3500,3650,4086,4086,4086,4086,4086,2940] ;y boundary at red edge

;Y position at the blue  end of the order

ybluelow=[120,380,638,864,1052,1216,1352,1508,1668,1826]
ybluehigh=[290,540,788,1004,1194,1348,1484,1642,1796,1952]

; x position of center of order, at x=2048
ygreenlow=[132,406,648,862,1056,1230,1390,1542,1686,1832]
ygreenhigh=[300,564,798,1010,1194,1366,1524,1670,1810,1958]


;y position  at the red end of the order

yredlow=[64,284,466,624,686,844,998,1142,1284,1714]
yredhigh= [230,442,620,780,832,988,1134,1276,1418,1836]


;make sure we get the entire order
dn=-8

yredlow=yredlow-dn
ygreenlow=ygreenlow-dn
ybluelow=ybluelow-dn

yredhigh=yredhigh+dn
ygreenhigh=ygreenhigh+dn
ybluehigh=ybluehigh+dn





yy=findgen(2048)
xx=findgen(4096)
yylow=fltarr(10,4096)
yyhigh=fltarr(10,4096)


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

flat_dat=fltarr(4096,10)


for j=0,9 do begin
  for i=0L,4095 do begin
;   flat_dat(i,yylow(j,i):yyhigh(j,i))=0.0
    flat_dat(i,j)=median(flatin(i,yylow(j,i):yyhigh(j,i)))
  endfor
  flat_dat(*,j)=smooth(flat_dat(*,j),311)
  for i=0L,4095 do begin
;   flat_dat(i,yylow(j,i):yyhigh(j,i))=0.0
    flattmp(i,yylow(j,i):yyhigh(j,i))=flat_dat(i,j)
  endfor

endfor

;JXP
flattmp=abs(flattmp) 
a = where(flattmp NE 0)

flatout[a]=flat[a]/flattmp[a]

;flatout=flatout * (flatout gt 0 and flatout lt 10)

; Release memory
delvarx, flattmp, flatin

return
end
