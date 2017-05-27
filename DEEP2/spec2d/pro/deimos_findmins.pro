;+
; NAME:
;   deimos_findmins
;
; PURPOSE:
;    Determine minima levels in slit gaps to allow us to subtract them off
;
; CALLING SEQUENCE:
;   deimos_findmins,image,mask,planeparams[,minarr,xarr,yarr,nx=nx,ny=ny]
; 
; INPUTS:
;   image  - flat/data frame to find minima in
;   mask   - bad pixel mask for the image frame (0 on good pixels)
;
;
; OUTPUTS:
;   planeparams - parameters defining a planar fit to the set of minima.  
;   This is a 3-element x 2 array of the form 
;    [constant, x coefficient, y coefficient (where e.g. x=0...2047)], 
;   where the 0th row is for the bottom half of the chip, and the 1st row is fit for the top half.
;
;   minarr- array containing the minima found
;   xarr  - mean x index at each point of minarr
;   yarr  - mean y index at each point of minarr
;
; KEYWORDS:
;   nx - number of regions to use across the x direction
;   ny - number of regions to use across the y direction
;        onlymeds -  if set, only the region-by-region median is calculated 
;        (useful to run on a superdark)
; MODIFICATION HISTORY:
;    
;-
pro deimos_findmins,image,mask,planeparams,minarr,xarr,yarr,nx=nx,ny=ny,planefit=planefit,onlymeds=onlymeds
 

nrows=4096
ncols=2048

medmiddle = median(image[1000:1100,2000:2200])

isdata=(medmiddle lt 300)

nbin=4+28*isdata

if n_elements(nx) eq 0 then nx=8
if n_elements(ny) eq 0 then ny=8
if n_elements(mask) eq 0 then mask=image*0
if n_elements(onlymeds) eq 0 then onlymeds=0
	

	newimage=rebin(image,ncols,nrows/nbin)
	newmask=rebin(mask,ncols,nrows/nbin)

	minarr=fltarr(nx-1,ny)
	minima=-1
	ctarr=minarr
	xarr=minarr
	yarr=minarr
	sm2=shift(newimage,-2,0)
	sm1=shift(newimage,-1,0)
	sp2=shift(newimage,2,0)
	sp1=shift(newimage,1,0)

	locmin=(newimage lt sm2 and newimage lt sm1 and $
		newimage lt sp2 and newimage lt sp1 and newmask eq 0)
;	fracderiv=(sp1-sm1)/(newimage > 1E-10)
;	fracderiv2=(sp2-sm2)/(newimage > 1E-10)
;	derivlim=1E10;0.45+0.3*isdata

	fracderiv=(sm1-sm2)/(newimage > 1E-10)
	fracderiv2=(sp1-sp2)/(newimage > 1E-10)
	derivlim=0.10+0.3*isdata

	flatarea=(((fracderiv) lt derivlim) and ((fracderiv2) lt derivlim)) 
		;$OR image lt 6

	delvarx,sp1,sp2,sm1,sm2

	for j=0,ny-1 do begin
		for i=0,nx-2 do begin
			minrow=fix(nrows/nbin/ny*j)
			maxrow=fix(nrows/nbin/ny*(j+1)-1)<4095
			nrow=(maxrow-minrow+1)
		
			mincol=fix(ncols/nx*i)

;			maxcol=fix(ncols/nx*(i+1)-1)<2047
			maxcol=(mincol+2*ncols/nx)<2047

			region=newimage[mincol:maxcol,minrow:maxrow]
			regmin=locmin[mincol:maxcol,minrow:maxrow]
			regflat=flatarea[mincol:maxcol,minrow:maxrow]
			regmask=newmask[mincol:maxcol,minrow:maxrow]
			medlevel=median(region)

                        if  onlymeds eq 0 then begin
                            cutoff = 0.15+0.4*isdata
                            if NOT isdata then begin
                                if medmiddle lt 1.2E4 then cutoff = 5E3/medlevel >  cutoff
                            endif
                            minimapix=(regmin and region lt cutoff*medlevel)


; 		this lets us avoid narrow gaps  
                            broadminima=(minimapix and regflat)
			
                            broadct=total(broadminima,2)
                            broadarr=(smooth(broadct,3) gt 0.4/3*nrow/nbin) # (fltarr(nrow)+1.)

                            whminima=where(broadminima and broadarr,minct)
;		atv,broadminima and broadarr
                            if minct ge 1 then minima=region[whminima] 	
;			if minct ge 1 then plothist,minima,xrange=[0,30.*(100.-99.*isdata)],xa,ya,bin=(10.-9*isdata)
		
                            s=sort(minima)

                            if minct ge 1 then $
	minarr[i,j]=minima[s[(n_elements(s)/50. > 10) <(n_elements(s)-1)/5]] $
                            else minarr[i,j]=-1
                            ctarr[i,j]=minct
                            
                    endif else begin
                        minarr[i,j]=medlevel
                        ctarr[i,j]=1

                    endelse

                    xarr[i,j]=[mincol+maxcol]/2.
                    yarr[i,j]=[minrow+maxrow]/2.*nbin

              endfor
        endfor

        if onlymeds eq 0 then begin
            err=1/sqrt(ctarr > 1E-10)

; deweight highest pixel in each row/column
            for j=0,ny-1 do begin
                maxpix=max(minarr[*,j],whmax)
                err[whmax,j]=err[whmax,j]*10.
            endfor

            for i=0,nx-2 do begin
                maxpix=max(minarr[i,*],whmax)
		err[i,whmax]=err[i,whmax]*10.
		minpix=min(minarr[i,*],whmax)
		err[i,whmax]=err[i,whmax]*3.
		badcol=total(minarr[i,*] eq -1)
		if badcol gt 3 then err[i,*]=err[i,*]*5.
            endfor
	
            top=where(yarr ge nrows/2)
            bot=where(yarr le nrows/2)

            planefit=xarr*0.
            planeparams=fltarr(3,2)

            fitx=xarr[bot]
            fity=yarr[bot]
            fitxy=transpose([[fitx],[fity]])

            fitmins=minarr[bot] ;reform(minarr,n_elements(minarr))
            fiterr=err[bot]
            wh=where(fitmins eq -1,whct)
            if whct gt 0 then fitmins[wh]=median(fitmins)
            if whct gt 0 then fiterr[wh]=1.E20

            botparams=regress(fitxy, $
                  fitmins,measure=fiterr,yfit=planefit,sigma=sigparams,const=botconst)

            planefit[bot]=botconst+botparams[0]*xarr[bot]+botparams[1]*yarr[bot]

            planeparams[*,0]=[botconst,botparams[0],botparams[1]]

            fitx=xarr[top]
            fity=yarr[top]
            fitxy=transpose([[fitx],[fity]])
            fitmins=minarr[top] ;reform(minarr,n_elements(minarr))
            fiterr=err[top]
            wh=where(fitmins eq -1,whct)
            if whct gt 0 then fitmins[wh]=median(fitmins)
            if whct gt 0 then fiterr[wh]=1.E20

            topparams=regress(fitxy, $
		fitmins,measure=fiterr,yfit=planefit,sigma=sigparams,const=topconst)

            planefit[top]=topconst+topparams[0]*xarr[top]+topparams[1]*yarr[top]

            planeparams[*,1]=[topconst,topparams[0],topparams[1]]
     endif else planeparams=0


  return 
end








