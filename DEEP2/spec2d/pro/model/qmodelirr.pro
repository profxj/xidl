pro qmodelirr,amap,bmap,sys,xmm,ymm,lambda,xics,yics,ccdpix,xpix,ypix,text=text,cubic=cubic

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

if n_elements(text) eq 0 then text=0
if n_elements(cubic) eq 0 then cubic=0
if cubic gt 0 then cubic=-0.5


dat=1-text

; read amap
if dat then restore,amap

; set up the grating transform
a3=gsetup(sys)

; get mapping and convert to r

xindx=(xmm-xmin)/xstep
yindx=(ymm-ymin)/ystep
tanxx=interpolate(tanx,xindx,yindx,cubic=cubic)
tanyy=interpolate(tany,xindx,yindx,cubic=cubic)


; NOTE: THE CODE HAS BEEN VECTORIZED.  ALL OPERATIONS ARE VECTOR OPERATIONS
; SAVE THE XFORMS.

r=dblarr(npt,3)
r[*,2]=-1.d0 / sqrt (1. + tanxx*tanxx + tanyy*tanyy)
r[*,0]=r[*,2]*tanxx
r[*,1]=r[*,2]*tanyy

; xform into grating system

for i=0,npt-1 do begin 
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

r[*,0] = sin (gamma)
r[*,1] = sin (-beta) * cos (gamma)
r[*,2] = cos (-beta) * cos (gamma)

; xform out of grating system.  Must do separately for each r

for i=0,npt-1 do begin 
	rg=reform(r[i,*],3)    				
	gen_xfm,rg,a3
	r[i,*]=rg 
endfor

; convert to tanx,tany

tanxi=r(*,0)/r(*,2)
tanyi=r(*,1)/r(*,2)

; get mapping into ICS pixels

; but first read bmap
if dat then restore,bmap

; set up arrays of x,y at each point
xmat=xarr # (fltarr(n_elements(yarr))+1)
ymat=(fltarr(n_elements(xarr))+1) # yarr

; determine range in tanx,tany needed for interpolations
nsamp=50.
mmx=minmax(tanxi)
xdiff=(mmx(1)-mmx(0)) > 1.E-2
mmx(0)=mmx(0)-.05*xdiff
mmx(1)=mmx(1)+.05*xdiff

mmy=minmax(tanyi)
ydiff=(mmy(1)-mmy(0)) > 1.E-2
mmy(0)=mmy(0)-.05*ydiff
mmy(1)=mmy(1)+.05*ydiff

; if things go crazy in the middle, probably need to alter this
if xdiff gt 0.1 and ydiff gt 0.1 then nsamp=200.

wh=where(tanx gt mmx(0)-xdiff and tanx lt mmx(1)+xdiff and tany gt mmy(0)-ydiff and tany lt mmy(1)+ydiff)

; perform delaunay triangulation to be used by trigrid
triangulate,tanx(wh),tany(wh),triang;,b

;interpolate to regular grid in tanx & tany
interpx=trigrid(tanx(wh),tany(wh),xmat(wh),triang,[(mmx(1)-mmx(0))/(nsamp-1),(mmy(1)-mmy(0))/(nsamp-1)],[mmx(0),mmy(0),mmx(1),mmy(1)]);,extrap=b)
interpy=trigrid(tanx(wh),tany(wh),ymat(wh),triang,[(mmx(1)-mmx(0))/(nsamp-1),(mmy(1)-mmy(0))/(nsamp-1)],[mmx(0),mmy(0),mmx(1),mmy(1)]);,extrap=b)

; use grids to determine values of xics,yics from tanxi,tanyi
xindx=(tanxi-mmx(0))/(1.1*xdiff)*(nsamp-1.)
yindx=(tanyi-mmy(0))/(1.1*ydiff)*(nsamp-1.)

xics=interpolate(interpx,xindx,yindx,cubic=cubic)
yics=interpolate(interpy,xindx,yindx,cubic=cubic)

; setup ccd geometry
ccd=ccd_geom(sys)

; determine locations on ccds
for n=0,7 do begin
	ics_to_ccd,xics,yics,ccd,n,xn,yn
	
	wh=where(xn gt 0 and xn le 2048 and yn gt 0 and yn le 4096)
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
