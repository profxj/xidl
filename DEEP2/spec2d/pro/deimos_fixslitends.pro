;+
; NAME:
;
;  DEIMOS_FIXSLITENDS
;
; PURPOSE:
;
;  Remove derivative tails in flux at the spatial ends of slits
;
; CATEGORY:
;
;   spslit processing
;
; CALLING SEQUENCE:
;
;   fixedspslit=deimos_fixslitends(spslit,slitfn2d,leftshift,rightshift)
;
; INPUTS:
;     spslit -- an spslit structure
;     slitfn2d -- a 2d slitfunction array
;
; OUTPUTS:
;     fixedspslit - a fixed version of spslit
; OPTIONAL OUTPUTS:
;     leftshift -- shift applied to low-X end of slit
;     rightshift -- shift applied to high-X end of slit
; MODIFICATION HISTORY:
;   jan 2003 mar 26
;-

function deimos_fixslitends, spslit,slitfn2d,shiftl, shiftr

	
	s=size(spslit.flux,/DIM)
	nx=s[0]
	ny=s[1]
	nslitfn=n_elements(slitfn2d)/ny

    derivfn=slitfn2d*0.
    dderivfn=derivfn
    for i=0,nslitfn-1 do derivfn[i,*]=deriv(slitfn2d[i,*])
    for i=0,nslitfn-1 do dderivfn[i,*]=deriv(deriv(slitfn2d[i,*]))


    signature=derivfn/slitfn2d
    signature2=(derivfn/slitfn2d)^2-0.5*dderivfn/slitfn2d

	nsegments=8
    
    profile=fltarr(nsegments-1,ny)

    xarray=lindgen(nx,ny) MOD nx

    signature=rebin(signature,nsegments,ny)
    signature2=rebin(signature2,nsegments,ny)
    fix=fltarr(nsegments-1,2)
	sigfix=fix
    npixi=fltarr(nsegments-1)
    npixsky=npixi
        edge=7
    junk=1.

            tempslit=spslit
; check to see if on bright object
        fullprof=find_object(tempslit,prof=pivar,npix=junk,mode=mode)
        peakinfo,fullprof,pkcol,fwhm,s2n_fwhm=sig2noise, $
          prof=pivar,npix=junk,pk_q=pk_q
	  if n_elements(pk_q) eq 0 then pk_q = -1

	if total(pk_q) lt 0  then model=fltarr(ny) $
	   else begin
	        model=exp(-(findgen(ny)-pk_q[0])^2 / $
			(2.*fwhm[0]^2/2.35^2*1.1^2))
	        model=model*max(fullprof)
	   endelse

        leftmodel=0.
        rightmodel=0.
        if max(fullprof) gt 0.75*mode then begin
            leftmodel=max(model[0:edge-1])/mode
            rightmodel=max(model[ny-edge-2:ny-1])/mode
        endif

    for i=0,nsegments-2 do begin
        profile[i,*]=find_object(tempslit,/NOSUB,prof=pivar, $
                                 pixrange=[(nx/nsegments)*i,(nx/nsegments)*(i+2)-1],npix=npix)
        tempslit=spslit
        npixsky[i]=median(npix)
        
; deweight first and last pixel
        pivar[0]=pivar[0]/10
        pivar[ny-1]=pivar[ny-1]/10

        pivarl=pivar
        pivarr=pivar
; and only use left or right side at one time
        pivarl[edge:ny-1]=1.E-20
        pivarr[0:ny-1-edge]=1.E-20

; fit over 2 bins at a time
        sig=total(signature[i:i+1,*],1)/2.
        sig2=total(signature2[i:i+1,*],1)/2.
        
        prof=reform(profile[i,*])
        resultl=svdfit(sig,prof,2,measure=1/sqrt(pivarl),yfit=yfitl,sigma=sigmal)
;        print,'left: ',resultl[1]/resultl[0]
        fix[i,0]=resultl[1]/resultl[0]*(leftmodel lt 0.05)
	sigfix[i,0]=sqrt(sigmal[1]^2/resultl[1]^2+ $
                         sigmal[0]^2/resultl[0]^2)*abs(fix[i,0])
        resultr=svdfit(sig,prof,2,measure=1/sqrt(pivarr),yfit=yfitr,sigma=sigmar)
;        print,'right: ',resultr[1]/resultr[0]
        fix[i,1]=resultr[1]/resultr[0]*(rightmodel lt 0.05)
        npixi[i]=median(npix)
	sigfix[i,1]=sqrt(sigmar[1]^2/resultr[1]^2+ $
           sigmar[0]^2/resultr[0]^2)*abs(fix[i,1])

    endfor


; combine the fit shifts from the various bins

    wt=(npixi)
    lowsections=where(npixi lt 100,lowct)
    highsections=where(npixi ge 100,highct)

    if lowct gt 0 then wt[lowsections]=1.E-4

    if highct gt 3 then shiftl=median(fix[highsections,0],/even) else shiftl=total(fix[*,0]*wt)/total(wt)
    if highct gt 3 then shiftr=median(fix[highsections,1],/even) else shiftr=total(fix[*,1]*wt)/total(wt)
	
    if highct gt 3 then begin
	sigl=median(sigfix[highsections,0])
	sigr=median(sigfix[highsections,1])
	sigl2=djsig(fix[highsections,0])
	sigr2=djsig(fix[highsections,1])
    endif else begin
	sigl=median(sigfix[*,0])
	sigr=median(sigfix[*,1])
	sigl2=djsig(fix[*,0])
	sigr2=djsig(fix[*,1])
    endelse

	badshifts=(shiftl gt 6 OR shiftr gt 6 OR (sigl > sigr) gt 0.5 OR (sigl2 > sigr2) gt 1.5)

	print,'Left shift :', shiftl, ' +/- ',sigl,' /',sigl2, ' regions: ',highct,form='(A,f8.3,A,f8.4,A,f8.4,A,i2)'

	print,'Right shift:', shiftr, ' +/- ',sigr,' /',sigr2,form='(A,f8.3,A,f8.4,A,f8.4)'



; do not make large shifts
    if badshifts  then shiftl=0.
    if badshifts then shiftr=0.

	if badshifts then print,'Applied no shifts!'

    fixl=rebin( (slitfn2d+derivfn*shiftl)/slitfn2d,nx,ny)

    fixr=rebin( (slitfn2d+derivfn*shiftr)/slitfn2d,nx,ny)

; smoothly transition between left/right/center
    wtl=fixl*0
    wtl[0:edge]=1.
    wtl= ( (1+(edge-findgen(ny))*0.2) < 1.) > 0.
    wtl=(fltarr(nx)+1.) # wtl
    wtr= ( (1+(findgen(ny)-ny+edge+1)*0.2) < 1.) > 0.
    wtr=(fltarr(nx)+1.) # wtr


    fixfull=fixr*wtr+fixl*wtl+(1.-wtl-wtr)

    tempslit.flux=spslit.flux/fixfull
	tempslit.ivar=tempslit.ivar*fixfull^2
	
return,tempslit
end








