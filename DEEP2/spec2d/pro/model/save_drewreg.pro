pro save_drewreg,infile,outfile,oversamp

; ingests drew's table of x,y,tanx,tany and interpolates to a regular
; grid; use for bmap only!

if n_elements(oversamp) eq 0 then oversamp=1.
readcol,infile,x1,y1,tanx1,tany1,format='d,d,d,d'

nrow=n_elements(x1)
xmin=min(x1)
xmax=max(x1)
ymin=min(y1)
ymax=max(y1)

xvals=[xmin]
yvals=[ymin]

for i=0,nrow-1 do begin
	whx=where(xvals eq x1(i))
	why=where(yvals eq y1(i))
	if total(whx) eq -1 then xvals=[xvals,x1(i)]
	if total(why) eq -1 then yvals=[yvals,y1(i)]
endfor

nx=n_elements(xvals)
ny=n_elements(yvals)


; determine range in tanx,tany needed for interpolations
;nsamp=50.
mmx=minmax(tanx1)
xdiff=(mmx(1)-mmx(0)) 
mmx(0)=mmx(0)+1.E-4*xdiff
mmx(1)=mmx(1)-1.E-4*xdiff

mmy=minmax(tany1)
ydiff=(mmy(1)-mmy(0)) 
mmy(0)=mmy(0)+1.E-4*ydiff
mmy(1)=mmy(1)-1.E-4*ydiff

nx=float(nx)*oversamp
ny=float(ny)*oversamp

txstep=(mmx(1)-mmx(0))/(nx-1)
tystep=(mmy(1)-mmy(0))/(ny-1)

; perform delaunay triangulation to be used by trigrid
triangulate,tanx1,tany1,triang,b

;interpolate to regular grid in tanx & tany
gridx=trigrid(tanx1,tany1,x1,triang,[(mmx(1)-mmx(0))/(nx-1),(mmy(1)-mmy(0))/(ny-1)],[mmx(0),mmy(0),mmx(1),mmy(1)],missing=-1.E10,/quintic);,extrap=b)
gridy=trigrid(tanx1,tany1,y1,triang,[(mmx(1)-mmx(0))/(nx-1),(mmy(1)-mmy(0))/(ny-1)],[mmx(0),mmy(0),mmx(1),mmy(1)],missing=-1.E10,/quintic);,extrap=b)

txmin=mmx(0)
txmax=mmx(1)
tymin=mmy(0)
tymax=mmy(1)

print,(mmx(1)-mmx(0))/txstep
print,(mmy(1)-mmy(0))/tystep

;gridx2=min_curve_surf(x1,tanx1,tany1, $ 
                      ;xout=mmx(0)+findgen(nx-1)*txstep, $
                   ;yout=mmy(0)+findgen(ny-1)*tystep)
;                      gs=[txstep,tystep],$
;                      bounds=[mmx(0),mmy(0),mmx(1),mmy(1)])

;gridy2=min_curve_surf(y1,tanx1,tany1,gs=[txstep,tystep],$
;                      bounds=[mmx(0),mmy(0),mmx(1),mmy(1)])


save,gridx,gridy,txmin,txmax,tymin,tymax,txstep,tystep,gridx2,gridy2,f=outfile

;save,xarr,yarr,tanx,tany,xmin,xmax,xstep,ymin,ymax,ystep,f=outfile

return
end
