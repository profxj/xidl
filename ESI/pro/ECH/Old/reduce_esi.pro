;
; NAME:
;       reduce_esi
;
; PURPOSE:
;       reduce and rectify esi data
;
; CALLING SEQUENCE:
;
;
;reduce_esi, data, data_nosky, datafile=datafile,
;      dataoutfile=dataoutfile, lambdafile=lambdafile,
;      flatfile=flatfile, darkfile=darkfile,
;      (reffile=reffile,refoutfile=refoutfile,ref_datafile=ref_datafile,
;      mansky=mansky)
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;      
;
; REQUIRED KEYWORDS:
;      
;	datafile - 	filename for input data
;	dataoutfile - 	filename for output data
;	lambdafile - 	filename for input wavelength model
;	flatfile - 	filename for input reference flats
;	
;
; OPTIONAL KEYWORDS: 
;	   
;	reffile	-	filename for input standard star
;	refoutfile - 	filename for output standard star
;	ref_datafile -	filename for input standard star photometry
;	
;
; OUTPUTS:
;	data -		rectified,flattened,debised data
;   	data_nosky - 	data with sky subtracted
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
;	rotate
;	congrid
;   	qzap
;	scatter_model2	
;	flat_new
;	maporders
;	skysub_esi
;		
;	
;	
;
; REVISION HISTORY:
;  
;-
;--------------------------------------------------------------------------- 
pro reduce_esi, data, data_nosky, datafile=datafile, dataoutfile=dataoutfile, $
                lambdafile=lambdafile, flatfile=flatfile, darkfile=darkfile, $
                reffile=reffile, refoutfile=refoutfile, $
                ref_datafile=ref_datafile,mansky=mansky,nframe=nframe,zap=zap


  if NOT keyword_set(reffile) then reffile = 0
  if NOT keyword_set(mansky) then mansky = 0
  if NOT keyword_set(nframe) then nframe = 1
  if NOT keyword_set(zap) then zap = 0

;read in data ***********************************************************
print,'reading datafile'
data=readfits(datafile,datahd, /no_unsigned)
print,'reading flatfile'
flat=readfits(flatfile,flat_hd, /no_unsigned)
print,'reading darkfile'
dark=readfits(darkfile,dark_hd, /no_unsigned)
if reffile ne 0 then ref=readfits(reffile,ref_hd)


;trim overscan**********************************************************
print,'trimmimg overscan'
   a = systime(1)
szd=size(data)			;make sure not already trimmed
if szd(1) gt 2048 then data=data(25:2072,*)
szf=size(flat)
if szf(1) gt 2048 then flat=flat(25:2072,*)
szdd=size(dark)
if szdd(1) gt 2048 then dark= dark(25:2072,*)
if reffile ne 0 then ref=ref(25:2072,*)
   print,  'calculation time', systime(1)-a

;or fix the size************************************************************
   a = systime(1)
if szd(1) lt 2048 then data=congrid(data,2048,4096)
if szdd(1) lt 2048 then dark=congrid(dark,2048,4096)
if szf(1) lt 2048 then flat=congrid(flat,2048,4096)
   print,  'calculation time', systime(1)-a
;rotate***********************************************************

print,'rotate data sets'
   a = systime(1)
data=rotate(data,4)
flat=rotate(flat,4)
dark=rotate(dark,4)
if reffile ne 0 then ref=rotate(ref,4)
   print,  'calculation time', systime(1)-a
;subtract bias***********************************************************

print,'subtract bias for ',nframe,'frames'
if szd(1) gt 2048 then data=data-nframe*dark	;usually if trimmed then biased also
flat=flat-dark
if reffile ne 0 then ref=ref-dark


;create hot pixel mask and multiply*************************************
print,'hot pixel mask'
   a = systime(1)
;mask=dark*0+(dark lt 800)
;mask = (mask gt 0.9)

mask=data*0+1
;mask(2646:4095,421:441)=0	;hot column
mask(2646:4095,410:441)=0
;mask(3875:3914,58:110)=0	;hot pixel
mask(3782:4003,0:217)=0
   print,  'calculation time', systime(1)-a

data=data*mask
flat=flat*mask
if reffile ne 0 then ref=ref*mask

;zap if neccesary *****************************************
if zap ne 0 then begin
	  a = systime(1)
	print,'qzap'
	qzap,data,data
	if reffile ne 0  then  qzap,ref,ref

	   print,  'calculation time', systime(1)-a
endif
; model and remove intra-order scattered light*************************

print,'scatter subtract'
;.run scatter_model2
   a = systime(1)
scatter_model2, flat, scattermod
flat=flat-scattermod
   print,  'calculation time', systime(1)-a
;   a = systime(1)
;scatter_model2,data,scattermod
;data=data-scattermod
;   print,  'calculation time', systime(1)-a
;if reffile ne 0 then scatter_model2,ref,scattermod
;if reffile ne 0 then  ref=ref-scattermod

;data=data*mask
;flat=flat*mask
;if reffile ne 0 then  ref=ref*mask

; flatten the flat and divide out ***********************************************

print,'flatten data'

;.run flat_new
   a = systime(1)
flat_new,nframe*flat,flatout
;JXP
flatout=flatout+ 1e10 * (flatout lt 0.01)
data=data/flatout*mask
;data[where(data GT 60000. OR data LT -1000.)] = 0.
data=data *( data le 60000 and data ge -1000)
if reffile ne 0 then ref=ref/flatout*mask
if reffile ne 0 then  ref=ref* (ref le 60000 and ref ge -1000)
   print,  'calculation time', systime(1)-a



;rectify***********************************************************

print,'rectify data'

;this is where you would include a better coefficient fit......

coeff=fltarr(4, 10, 5)

openr,1,'rectify_coeff.txt'
readf,1,coeff
close,1

;.run maporders
   a = systime(1)
maporders,data,temp,coeff
data_rect=temp
   print,  'calculation time', systime(1)-a
if reffile ne 0 then maporders,ref,temp,coeff
if reffile ne 0 then ref=temp
;  a = systime(1)
;maporders,mask,temp,coeff
;mask=temp

 print,  'calculation time', systime(1)-a
;read in lambda model :***********************************************************
print,'read in lambda model'

lambda_model=readfits(lambdafile,lambda_modelhd)
;calculate inverse variance************************************************************************************
print,'calculate inverse variance'

invvar=1/data *(data ne 0)
invvar=invvar*mask
;hard clip some bad areas
i=9
invvar(2880:4095,i*180:i*180+169)=0  
i=3
invvar(0:980,i*180:i*180+169)=0
i=1
invvar(0:1000,i*180:i*180+169)=0
invvar(2601:4095,i*180:i*180+169)=0


;subtract the sky****************************************************************************************

print,'subtract the sky'
data_nosky=fltarr(4096,2048)
;.run skysub_esi.pro

if mansky ne 0 then begin
  i=7
  wset,0
  temp=rebin(data(*,(i*180):(i*180+169)),1,170)
  plot,temp

  for i=0,9 do oplot,rebin(data(*,(i*180):(i*180+169)),1,170)

  print,"select lower position for skyrows1"
  cursor,x1,y1
  wait,.3
  print,"select upper position for skyrows1"
  cursor,x2,y2
  wait,.3

  skyrows1=findgen(x2-x1)+x1

  print,"select lower position for skyrows2"
  cursor,x1,y1
  wait,.3
  print,"select upper position for skyrows2"
  cursor,x2,y2
  wait,.3
  skyrows2=findgen(x2-x1)+x1

  oplot,[skyrows1],temp,color=200,thick=5
  oplot,[skyrows2],temp,color=200,thick=5

endif


if mansky eq 0 then skyrows1=findgen(20)+10
if mansky eq 0 then skyrows2=findgen(20)+130


;hard clip some bad areas
i=9
data(2880:4095,i*180:i*180+169)=0  
i=3
data(0:980,i*180:i*180+169)=0
i=1
data(0:1000,i*180:i*180+169)=0
data(2601:4095,i*180:i*180+169)=0

for i=2,9 do begin
  print,'order number ',i
  galspec=data(*,i*180:i*180+169)
  model=lambda_model(*,i*180:i*180+169)
  skysub_esi,galspec,galspec, model, model, skyrows=[skyrows1,skyrows2],objnosky=objnosky,invvar=invvar
data_nosky(*,i*180:i*180+169)=objnosky

endfor

writefits,dataoutfile,data_nosky,datahd

return

end
