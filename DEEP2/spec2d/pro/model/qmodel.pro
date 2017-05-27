pro qmodel,sys,xmm,ymm,lambda,xics,yics,ccdpix,xpix,ypix, $
   alphafile,betafile,goodlocation=goodlocation,cubic=cubic

;+
; NAME:
;        QMODEL
;
;
; PURPOSE:
;        Calculate pixel locations on the DEIMOS mosaic for a given
;        slit position on a mask and wavelength
;
; CATEGORY:
;        DEIMOS optical model
;
;
; CALLING SEQUENCE:
;        qmodel,sysstructure,xmm,ymm,lambda,xics,yics,ccdpix,xpix,ypix,[alphafile,betafile,goodlocation=goodlocation,cubic=cubic]
; 
; INPUTS:
;        sysstructure - structure containing DEIMOS system parameters,
;                       provided by SYSINIT
;        xmm - array of x positions on the mask for which output is
;              desired (mm)
;        ymm - array of y positions on the mask for which output is
;              desired (mm)
;        lambda - array of wavelengths for which output is desired
;                 (Angstroms)
;
; OPTIONAL INPUTS:
;        alphafile - Name of an IDL binary file containing the amap
;            (pre-grating transformation) interpolation table.  
;            Default: amap.sav
;        betafile  - Name of an IDL binary file containing the bmap 
;            (post-grating transformation) interpolation table.  
;            Default: bmap.sav
;	
; KEYWORD PARAMETERS:
;        goodlocation - an array containing 1 for locations with good
;                       results, 0 for all others
;                       bad results are also flagged with XICS=YICS=-1E4
;        CUBIC - Sets value of the CUBIC keyword to use in INTERPOLATE
;
; OUTPUTS:
;        xics - Array of X locations corresponding to xmm,ymm, and
;               lambda, in the ICS coordinate system
;        yics - Array of X locations corresponding to xmm,ymm, and
;               lambda, in the ICS coordinate system
;        ccdpix - Array containing the CCD numbers on which the points
;                 defined by xmm,ymm, and lambda fall
;        xpix,ypix - Array containing the X/Y pixels on which the points
;                 defined by xmm,ymm, and lambda fall, on the ccd
;                 given by ccdpix
;
; OPTIONAL OUTPUTS:
;        None
; COMMON BLOCKS:
;        None
; SIDE EFFECTS:
;        None
; RESTRICTIONS:
;        SYSINIT must be run before to initialize system variables.
;
;        QTRMAP must have been run at some point (for the same values
;        of sys, other than grating parameters) to generate the
;           transformation tables.
;
;        xmm,ymm, and lambda must all have the same dimensions.
;
; EXAMPLE:
;        sys=sysinit(10.,gmm=1200)
;        qtrmap, sys
;        qmodel,sys,xmm,ymm,lambda,xics,yics,ccdnum,xpix,ypix,/cubic
;
;
; MODIFICATION HISTORY:
;        Based on a routine by Drew Phillips; see QMODEL in qmodel.x
;        Finished testing 12/5/01 JAN
;
;-







;real	xmm, ymm, wave			# X,Y in slitmask, wave
;double	a3[3,3]				# grating transform
;real	xics, yics			# pixel values in ICS
;real	xpix, ypix			# pixel values on CCD

initdeimos

xmm=double(xmm)
ymm=double(ymm)
wave=lambda*1.D-4	; in microns
npt=n_elements(xmm)
ccdpix=intarr(npt)-1
xpix=fltarr(npt)
ypix=fltarr(npt)


if n_elements(cubic) eq 0 then cubic=-0.5
if cubic gt 0 then cubic=-0.5

if n_elements(alphafile) eq 0 then alphafile = 'amap.sav'
if n_elements(betafile) eq 0 then betafile = 'bmap.sav'
amap = alphafile
bmap = betafile

; read amap
restore,amap

; set up the grating transform
a3=gsetup(sys)

; get mapping and convert to r
s = size(tanx)
sx = s(1)
sy = s(2)

xindx=(xmm-xmin)/xstep
yindx=(ymm-ymin)/ystep
tanxx=interpolate(tanx,xindx,yindx,cubic=cubic, missing=-1E10)
tanyy=interpolate(tany,xindx,yindx,cubic=cubic, missing=-1E10)
whbad = where(xindx lt 4 or xindx gt (sx-4) or yindx lt 4 or yindx gt (sy-4), badct)
if badct gt 0 then tanxx(whbad) = -1E10
if badct gt 0 then tanyy(whbad) = -1E10


; NOTE: THE CODE HAS BEEN VECTORIZED.  ALL OPERATIONS ARE VECTOR OPERATIONS
; SAVE THE XFORMS.

r=dblarr(npt,3)
r[*,2]=-1.d0 / sqrt (1. + tanxx*tanxx + tanyy*tanyy)
r[*,0]=r[*,2]*tanxx
r[*,1]=r[*,2]*tanyy

; xform into grating system

for i=0L,npt-1 do begin 
	rg=reform(r[i,*],3)    ; the reform forces rg to be a vector, not 1x3
				; matrix
	gen_xfm,rg,a3,/forward 
	r[i,*]=rg 
endfor

; convert to alpha,gamma

alpha=-atan(-r[*,1],-r[*,2])
gamma=atan(r[*,0],sqrt(r[*,2]*r[*,2]+r[*,1]*r[*,1]))

; Apply the grating equation

beta = asin ((sys.ORDER*sys.GRLINES*wave / cos (gamma)) - sin (alpha))

; convert beta, gamma into x,y,z (cf Schroeder p259); note sign reversal of beta


; added 3/5/03 - match drew's code for negative-lambda case
wavesign= 1-2*(wave lt 0)

r[*,0] = sin (gamma)
r[*,1] = sin (-beta*wavesign) * cos (gamma)
r[*,2] = cos (-beta*wavesign) * cos (gamma)

; xform out of grating system.  Must do separately for each r

for i=0L,npt-1 do begin 
	rg=reform(r[i,*],3)    				
	gen_xfm,rg,a3
	r[i,*]=rg 
endfor

; convert to tanx,tany

tanxi=r(*,0)/r(*,2)
tanyi=r(*,1)/r(*,2)

; get mapping into ICS pixels

; but first read bmap
restore,bmap

; use grids to determine values of xics,yics from tanxi,tanyi
xindx=(tanxi-txmin)/txstep
yindx=(tanyi-tymin)/tystep

s = size(gridx)
sx = s(1)
sy = s(2)

xics=interpolate(gridx,xindx,yindx,cubic=cubic, missing=-1E10)
yics=interpolate(gridy,xindx,yindx,cubic=cubic, missing=-1E10)
whbad = where(xindx lt 4 or xindx gt (sx-4) or yindx lt 4 or yindx gt (sy-4), badct)
if badct gt 0 then xics(whbad) = -1E10
if badct gt 0 then yics(whbad) = -1E10

wh = where(abs(xics) gt 1E9 or abs(yics) gt 1E9, badct)
if badct gt 0 then xics(wh) = -1E4
if badct gt 0 then yics(wh) = -1E4
goodloc = xics*0b+1b
if badct gt 0 then goodloc(wh) = 0

goodlocation=goodloc

; setup ccd geometry
ccd=ccd_geom(sys)

; determine locations on ccds
for n=0,7 do begin
	ics_to_ccd,xics,yics,ccd,n,xn,yn
	
	wh=where(xn gt 0 and xn le 2048 and yn gt 0 and yn le 4096)
;                AND goodloc)
	if total(wh) gt -1 then begin
		ccdpix(wh)=n	
		xpix(wh)=xn(wh)
		ypix(wh)=yn(wh)
	endif
endfor

; convert arrays to 0-indexed
xpix=xpix-1
ypix=ypix-1

return
end
