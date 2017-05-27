pro save_drewtxt,infile,outfile

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
xstep=(xmax-xmin)/(nx-1.)
ystep=(ymax-ymin)/(ny-1.)

xarr=xmin+dindgen(nx)*xstep
yarr=ymin+dindgen(ny)*ystep

tanx=dblarr(nx,ny)
tany=dblarr(nx,ny)

for i=0,nx-1 do begin
	for j=0,ny-1 do begin
		wh=where(abs(x1-xarr(i)) lt 1.D-4 and abs(y1-yarr(j)) lt 1.D-4)
		if total(wh) gt -1 then tanx(i,j)=tanx1(wh)
		if total(wh) gt -1 then tany(i,j)=tany1(wh)

	endfor
endfor

save,xarr,yarr,tanx,tany,xmin,xmax,xstep,ymin,ymax,ystep,f=outfile

return
end
