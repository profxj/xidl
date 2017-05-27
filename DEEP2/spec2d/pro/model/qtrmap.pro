pro qtrmap,sys,alphafile,betafile,nxstep,nystep,nextrastep,text=text, oversamp=oversamp
;+
; NAME:
;        QTRMAP
;
;
; PURPOSE:
;        Calculate "Quick-Trace-Maps" for the DEIMOS optical path to be
;        used by QMODEL
;
;
; CATEGORY:
;        DEIMOS optical model
;
;
; CALLING SEQUENCE:
;        qtrmap,sysstructure,[alphafile,betafile,nxstep,nystep,nextrastep,/TEXT,oversamp=]
; 
; INPUTS:
;        sysstructure - structure containing DEIMOS system parameters,
;                       provided by SYSINIT
;
;
; OPTIONAL INPUTS:
;        alphafile - Name of an IDL binary file to contain the amap
;            (pre-grating transformation) interpolation table.  
;            Default: amap.sav
;        betafile  - Name of an IDL binary file to contain the bmap 
;            (post-grating transformation) interpolation table.  
;            Default: bmap.sav
;        nxstep - roughly, half the number of intervals to use in
;                 the X grid.  Default: 50
;        nystep - roughly, half the number of intervals to use in
;                 the Y grid for the bmap; roughly the whole number
;                 for the amap.  Default: 50
;        nextrastep - a number of extra intervals to add to the grid
;                 in each dimension on either end, going outside
;                 the nominal range.  Provides a buffer for
;                 interpolation, etc.  Default: 5
;	
; KEYWORD PARAMETERS:
;        /TEXT - write output as a text file rather than an IDL binary file.
;        OVERSAMP= Factor by which to over/under-sample the data grid
;        in producing the bmap interpolation table.
;
; OUTPUTS:
;        amap, bmap interpolation tables saved in IDL binary .sav
;        files (unless /TEXT is set, in which case the input data for
;        the tables is recorded)
;
; OPTIONAL OUTPUTS:
;        amap/bmap text files if /TEXT is set.
;
; COMMON BLOCKS:
;        None
; SIDE EFFECTS:
;        None
; RESTRICTIONS:
;        SYSINIT must be run before to initialize system variables.
;
; EXAMPLE:
;        sys=sysinit(10.,gmm=1200)
;        qtrmap, sys
;        qmodel,sys,xmm,ymm,lambda,xics,yics,ccdnum,xpix,ypix,/cubic
;
;
; MODIFICATION HISTORY:
;        Based on a routine by Drew Phillips; see QTRMAP in qtrmap.x
;        Finished testing 12/5/01 JAN
;
;-



; see qtrmap.x

; STATUS: maps OK -- now need extension to individual CCDs and a wrapper script.

; QTRMAP: "Quick Trace Mapping" -- generate the data and/or mappings needed
; for the "quick-trace".  In this case, all pre- and post-grating elements
; are used to produce a single map each; only the grating x-form and grating
; equation are needed with these maps to produce all info.
;
; based largely on TRACE; is TRACE is changed we need to review here.
; Also, this is a good opportunity to start to "clean up" these codes.  Eg,
;build "spec_optics.x", break "pt_xfm" into the three parts required here, etc

if n_elements(text) eq 0 then text=0
if n_elements(alphafile) eq 0 and text then alphafile='amap.txt'
if n_elements(betafile) eq 0 and text then betafile='bmap.txt'
if n_elements(alphafile) eq 0 then alphafile='amap.sav'
if n_elements(betafile) eq 0 then betafile='bmap.sav'
if n_elements(nxstep) eq 0 then nxstep=50
if n_elements(nystep) eq 0 then nystep=50
if n_elements(nextrastep) eq 0 then nextrastep=5
if n_elements(nskip) eq 0 then nskip=1
if n_elements(oversamp) eq 0 then oversamp=1.
amap=alphafile
bmap=betafile
nx=nxstep
ny=nystep
nextra=nextrastep
nex=nextra

; set up the transforms

initdeimos

setup,e1,a2,a3,a4,ccd,sys

; write the amap

xstep=double(364. / (nx-1))
ystep=double(220. / (ny-1))

indexarr=lindgen(2*nx+2*nex+1,ny+2*nex+1)
xarr=(indexarr mod (2*nx+2*nex+1)) -nx-nex
xarr=xarr*xstep
yarr=(indexarr / (2*nx+2*nex+1)) -nex
yarr=yarr*ystep
;tanx=indexarr*0d0
;tany=indexarr*0d0
tmp=minmax(xarr)
xmin=tmp(0)
xmax=tmp(1)
tmp=minmax(yarr)
ymin=tmp(0)
ymax=tmp(1)



; calculate the amap 
;for i=0,n_elements(indexarr)-1,nskip do begin
	pre_grating,xarr,yarr,e1,a2,sys,r

	tanx = reform(r[*,0] / r[*,2],2*nx+2*nex+1,ny+2*nex+1)
	tany = reform(r[*,1] / r[*,2],2*nx+2*nex+1,ny+2*nex+1)
;endfor




; write the amap to file
if text then begin
	openw,2,amap
	for i=0,n_elements(indexarr)-1 do printf,2,xarr(i),yarr(i),tanx(i),tany(i),format='(%"%10.5f %10.5f %10.7f %10.7f")'      
	close,2
endif else save,xarr,yarr,tanx,tany,xmin,xmax,xstep,ymin,ymax,ystep,f=amap


; now calculate the bmap

ny = nx
xstep = double(4596. / (nx-1))
ystep = double(4596. / (ny-1))

indexarr=lindgen(2*nx+2*nex+1,2*ny+2*nex+1)
xarr=(indexarr mod (2*nx+2*nex+1)) -nx-nex
xarr=xarr*xstep
yarr=(indexarr / (2*nx+2*nex+1)) -ny -nex
yarr=yarr*ystep
;rad=(xarr*xarr)+(yarr*yarr)
;useable= rad le 5080.^2
tmp=minmax(xarr)
xmin=tmp(0)
xmax=tmp(1)
tmp=minmax(yarr)
ymin=tmp(0)
ymax=tmp(1)

;tanx=dblarr(2*nx+2*nex+1,2*ny+2*nex+1)
;tany=tanx

; calculate the bmap
;for i=0,n_elements(indexarr)-1 do begin
	;if useable(i) then begin
		ics_post_grating,xarr,yarr,a4,sys,r
		
		tanx = (r[*,0] / r[*,2])
		tany = (r[*,1] / r[*,2])
	;endif
;endfor

; write the bmap to file
if text then begin
	openw,2,bmap
	for i=0,n_elements(indexarr)-1 do if useable(i) then printf,2,xarr(i),yarr(i),tanx(i),tany(i),format='(%"%10.4f %10.4f %10.7f %10.7f")'      
	close,2
endif else begin

	; determine range in tanx,tany needed for interpolations

	mmx=minmax(tanx)
	xdiff=(mmx(1)-mmx(0)) 
	mmx(0)=mmx(0)+1.E-4*xdiff
	mmx(1)=mmx(1)-1.E-4*xdiff

	mmy=minmax(tany)
	ydiff=(mmy(1)-mmy(0)) 
	mmy(0)=mmy(0)+1.E-4*ydiff
	mmy(1)=mmy(1)-1.E-4*ydiff
	
	nx=2*float(nx)*oversamp
	ny=2*float(ny)*oversamp

	txstep=(mmx(1)-mmx(0))/(nx-1)
	tystep=(mmy(1)-mmy(0))/(ny-1)

	; perform delaunay triangulation to be used by trigrid
	triangulate,tanx,tany,triang,b

	; interpolate to regular grid in tanx & tany
	gridx=trigrid(tanx,tany,xarr,triang,[(mmx(1)-mmx(0))/(nx-1),(mmy(1)-mmy(0))/(ny-1)],[mmx(0),mmy(0),mmx(1),mmy(1)]);,extrap=b)
	gridy=trigrid(tanx,tany,yarr,triang,[(mmx(1)-mmx(0))/(nx-1),(mmy(1)-mmy(0))/(ny-1)],[mmx(0),mmy(0),mmx(1),mmy(1)]);,extrap=b)

	txmin=mmx(0)
	txmax=mmx(1)
	tymin=mmy(0)
	tymax=mmy(1)
	save,gridx,gridy,txmin,txmax,tymin,tymax,txstep,tystep,f=bmap

endelse

return
end
