;+
; NAME:
;
;  RESUBTRACT
;
; PURPOSE:
;
;  Refine 2d sky subtraction for a slit file 
;
; CATEGORY:
;
;  Kinematics
;
; CALLING SEQUENCE:
;
;  resubtract,filename[,manual=manual,multiplier=multiplier,nrows=nrows,plot=plot]
;
; INPUTS:
;
;  filename - name of a slit file.  Right now, need to do B & R files separately.
;
; KEYWORD PARAMETERS:
;

; OPTIONAL INPUTS:
;
;   MANUAL -- set to a 2d array containing the minimum & maximum row
;             numbers to use for the sky region.  This is a last
;             resort.  Note that the first row is row 0 in IDL, and
;             that the first 5 and last 5 rows of each slit are bad
;             for sky regions (due to PSF effects). If MANUAL has 4
;             elements instead of 2, it will take them to be
;             [min1,max1,min2,max2] for two sky regions, min1:max1 and
;             min2:max2.  MANUAL values
;             will override a MULTIPLIER or NROWS setting
;
;   MULTIPLIER -- set to a value which the routine will multiply the
;                 FWHMs by for the full width of an exclusion region
;                 around each object.  ~2.25 is reasonable,
;                 based on tests.  MULTIPLIER values will override an
;                 NROWS setting.
;   NROWS -- set to the number of sky rows to use; RESUBTRACT will
;            then find the best (in terms of furthest from objects)
;            NROWS rows to use as the sky region.  This seems to work
;            well, and is the default behavior (with NROWS = 10).
;
;   PLOT -- set to produce two diagnostic plots: a 1d plot showing the
;           collapsed light profile along the slit, with original sky
;           rows indicated as diamonds and the new sky region as red
;           stars, and an ATV window showing, from top to bottom, the
;           original 2d spectrum, the new spectrum, and the residual
;           sky model that has been taken out.  Can set to 'B' or 'R'
;           to only plot for B or R.
;
; OUTPUTS:
;
;    for input FITS file XXX.fits[.gz] containing a slit structure, a
;    gzipped FITS file XXX.modified.fits.gz will be produced
;
; MODIFICATION HISTORY:
; written JAN 6/21/05
;-



function resubtract_work,slit,hdrsl,color,manual=manual,multiplier=multiplier,nrows=nrows,plot=plot


if n_elements(manual) eq 0 then manual=0
if n_elements(nrows) eq 0 then nrows=10

if n_elements(multiplier) eq 0 then multi=2.25 else multi=multiplier

if n_elements(plot) eq 0 then plot=0

objinfo=mrdfits('obj_info*',1,/silent)

slitno=sxpar(hdrsl,'SLITNO')
isdata=1

sizey = (size(slit.flux, /dimens))[1]
ncol = (size(slit.flux, /dimens))[0]
wave=lambda_eval(slit,/double)


wh=where(objinfo.slitno eq slitno and $
	strpos(objinfo.objno,'s') ne 0 $
	AND objinfo.color eq color,ct)

nobj=ct

yobj_cat=fltarr(ct)
cat_fwhm=fltarr(ct)
yobj=yobj_cat
corr_fwhm=cat_fwhm
meas_fwhm=cat_fwhm

for i=0,ct -1 do yobj_cat[i]=objinfo[wh[i]].cat_objpos
for i=0,ct -1 do cat_fwhm[i]=objinfo[wh[i]].cat_fwhm
for i=0,ct -1 do yobj[i]=objinfo[wh[i]].objpos
for i=0,ct -1 do corr_fwhm[i]=objinfo[wh[i]].corr_fwhm
for i=0,ct -1 do meas_fwhm[i]=objinfo[wh[i]].fwhm

yobj=yobj+(yobj eq 0)*yobj_cat
meas_fwhm=meas_fwhm+(meas_fwhm eq 0)*corr_fwhm

; original region
offslit = yobj_cat LT 0 OR yobj_cat GT (sizey-1)

; refined region
offslit_new = yobj_cat LT 0 OR yobj GT (sizey-1)

        skyok = bytarr(sizey)
	skyok_new=skyok

	edge=5

badness=skyok*0.
idx=findgen(n_elements(badness))

skyok[edge:sizey-edge-1] = 1B
skyok_new[edge:sizey-edge-1] = 1B
skyok_manual=skyok_new

badness=badness+100*(skyok eq 0)
badness[edge-1]=.2
badness[sizey-edge]=.2

        ; remove rows near objects (within wid pixels)
        wid = 14

; exclude regions near objects                
     for k=0, nobj-1 do begin
         
          IF NOT(offslit[k]) THEN skyok[(yobj_cat[k]-wid/2) >0: $
                                        (yobj_cat[k]+wid/2) < (sizey-1)] = 0B
                     
          IF NOT(offslit_new[k]) THEN badness=badness+exp(-(yobj[k]-idx)^2/meas_fwhm[k]^2)


          IF NOT(offslit_new[k]) THEN skyok_new[(yobj[k]-multi*meas_fwhm[k]/2) >0: $
                               (yobj[k]+multi*meas_fwhm[k]/2) < (sizey-1)] = 0B
     

endfor

; exclude bright regions
        medprof=skyok*0.
        for i=0,sizey-1 do medprof[i]=median(slit.flux[*,i])
        if total(isdata) ge 1 then $
          skyhigh=medprof  gt 15. $
        else skyhigh = 0

        skyind = where(skyok AND (skyhigh eq 0), ct)




; make sure at least 7 sky rows
        while (ct lt 7 AND wid gt 2) DO BEGIN 
	      wid=wid-1 
	      skyok[floor(1+wid/4.)>(edge-1) : $
                    ceil(sizey-1-wid/4.)<(sizey-edge)] = 1B 
	      for k=0, nobj-1 do skyok[(yobj[k]-wid/2) > 0: $
                        (yobj[k]+wid/2) < (sizey-1)] = 0B &$
	      skyind = where(skyok AND isdata, ct) 
	      if wid lt 10 then vprint, 2, $
  	        'Warning: < 10 pix around object excluded on slit ',slitno 
	ENDWHILE   


        skyind_new = where(skyok_new AND (skyhigh eq 0), ct_new)

; and do the same for the new region
        while (ct_new lt 7 AND wid gt 2) DO BEGIN 
	      wid=wid-1 
	      skyok_new[floor(1+wid/4.)>(edge-1) : $
                    ceil(sizey-1-wid/4.)<(sizey-edge)] = 1B 
	      for k=0, nobj-1 do skyok_new[(yobj[k]-wid/14.*multi*meas_fwhm[k]/2) > 0: $
                        (yobj[k]+wid/14.*multi*meas_fwhm[k]/2) < (sizey-1)] = 0B &$
	      skyind_new = where(skyok_new AND isdata, ct_new) 
	      if wid lt 10 then print, $
  	        'Warning: < 10 pix around object excluded on slit ',slitno 
	ENDWHILE   

   specmask=1b-((slit.mask AND 4b) eq 4b)
   vig_mask=((slit.mask AND 8b) eq 8b)


; check for real bad stuff in cases with few good sky rows
; old sky region
	if ct lt 5 then begin

            isdata=specmask[100,*] eq 1 OR specmask[ncol-101,*] eq 1
            isdata=isdata AND (calib.slitfn gt 0.85) AND (vig_mask[ncol/2,*] eq 0)

	     skyok[edge:sizey-edge-1] = 1B
	     skyind = where(skyok AND isdata, ct)
	     vprint,2, 'no good sky region-using whole slit ',slitno
	endif

        if ct lt 5 then begin
             skyok[(edge-1)>0:(sizey-edge)<(sizey-1)] = 1B
             skyind = where(skyok AND isdata, ct)
             vprint,2, 'no good sky region-using even more of slit ',slitno
        endif

	if ct lt 5 then begin
	 	 vprint,2, 'no sky region on slit ',slitno
		 vprint,2, 'bspline will be very poor'
		skyind=where(skyok,ct)
            endif

	if n_elements(skyind) NE total (skyok) $
          AND total(skyind) ge 0 then begin
            skyok = skyok*0B
            skyok[skyind] = 1B
        endif

; new sky region
        if ct_new lt 5 then begin
             skyok_new[(edge-1)>0:(sizey-edge)<(sizey-1)] = 1B
             skyind_new = where(skyok_new AND isdata, ct_new)
             print, 'no good sky region-using even more of slit ',slitno
         endif

	if ct_new lt 5 then begin
	 	 print,2, 'no sky region on slit ',slitno
		 print,2, 'bspline will be very poor'
		skyind_new=where(skyok_new,ct_new)
	endif


	if n_elements(skyind_new) NE total (skyok_new) $
          AND total(skyind_new) ge 0 then begin
            skyok_new = skyok_new*0B
            skyok_new[skyind_new] = 1B
        endif

        flagivar=slit.ivar
     flagivar=flagivar/(1.+99.*vig_mask)


; favor manual over multiplier, nrows default
if total(manual) gt 0 then begin
    if n_elements(manual) lt 4 then $
      skyind=manual[0]+indgen(manual[1]-manual[0]+1) $
    else begin
        skyind=manual[0]+indgen(manual[1]-manual[0]+1) 
        skyind=[skyind,manual[2]+indgen(manual[3]-manual[2]+1)]
    endelse

    skyind=skyind[where(skyind gt edge-2 and skyind lt sizey-edge+1)]

endif else if n_elements(multiplier) gt 0 then begin 
    skyind=skyind_new
endif else begin
    srt=sort(badness)
    skyind=srt[0:nrows-1]

endelse

print,n_elements(skyind),' rows being used!'

           skywave = wave[*, skyind]
           skyflux = slit.flux[*,skyind]

; trying to keep bspline from going crazy at the ends by giving these
; masked pixels nonzero weight
           skyivar = (flagivar*((slit.mask AND 22b) eq 0))[*, skyind]
           skyx=lindgen(ncol,n_elements(skyind)) MOD ncol
           badsky=(skyivar eq 0.)
           whok=where(badsky eq 0b,okct)
           if okct gt 0 then begin
              meanivar=mean(skyivar[whok])
              whbad=where((badsky AND (skyx lt EDGE+1 OR skyx gt (ncol-edge-2))) ,badct)
              if badct gt 0 then skyivar[whbad]=meanivar/250.
           endif

; original bspline
           everyn = 2.*n_elements(skyind)/3.
           bkspace = (max(skywave)-min(skywave))/(1.5*ncol)
           testrange=wave[ncol/2,*]
           testrange=testrange[sort(testrange)]
 ; determine range of wavelengths in a given column
           deltalambda=abs(testrange[fix(0.9*sizey)]-testrange[fix(0.1*sizey)])

           nbkpts = 1.5*ncol      
           outmask1=1
           outmask2=1

; do this because bspline_iterfit is BROKEN and requires sorting
           ind = sort(skywave)
;           ind=ind[where(randomu(s,n_elements(ind)) gt 0.6)]

; only do this if we are very poorly sampled

           if deltalambda lt 0.5*bkspace then begin
              vprint,2, 'doing fixed breakpoints for slit ',slitno
              vprint,2, 'deltalambda: ',deltalambda,'  bkspace: ',bkspace
              nsky=n_elements(skyind)

              if max(offslit) eq 0 then row = floor(mean(yobj)) $
                 else row = sizey/2
              
              if row lt min(skyind) OR row gt max(skyind) $
                then row = median(skyind)

              vprint, 2, 'setting breakpoints to match row ', row

              breakpoints = (wave[0:ncol-1, row:row+1])
              breakpoints=total(breakpoints,2)/2.
	
;          print, 'repeating bspline - slit was too vertical!'
              sset = bspline_iterfit(skywave[ind], skyflux[ind], $
                   invvar=skyivar[ind], upper=20, lower=20, $ 
                   maxiter=3, fullbkpt = breakpoints, $
                   outmask=outmask2, /silent)
              outmask1 = outmask2

           endif else begin
; the old way:
             sset = bspline_iterfit(skywave[ind], skyflux[ind], $
               invvar=skyivar[ind], upper=20, lower=20, $ 
               maxiter=3, everyn=everyn, $
               outmask=outmask1, /silent)
           endelse

     skymodel = bspline_valu(wave, sset)

     if size(plot,/type) eq 7 then doplot=string(plot) eq color $
       else doplot = (plot eq 1)

     if doplot then atv,[[skymodel],[(slit.flux-skymodel)],[slit.flux]],min=-10,max=50



if doplot then begin
    curve=find_object(slit,/nosubtract)
    splot,curve,xtit='Row',ytit='Counts'
    if total(skyok) gt 0 then soplot,where(skyok),curve[where(skyok)],psym=4
    soplot,skyind,curve[skyind],psym=2,color=1
endif

out=slit
out.flux=out.flux-skymodel

return,out

end



pro resubtract,file,manual=manual,multiplier=multiplier,nrows=nrows,plot=plot

filename=findfile(file)
filename=filename[0]

posfits=strpos(filename,'B.fits') > strpos(filename,'R.fits') 
base=strmid(filename,0,posfits)


; do 'B' first
color='B'
fileb=base+'B.fits*'
slit=mrdfits(fileb,1,hdrsl,/silent)  
outfile=strcompress(base+'B.modified.fits',/remove)

fixed=resubtract_work(slit,hdrsl,color,manual=manual,multiplier=multiplier,nrows=nrows,plot=plot)

mwrfits,fixed,outfile,hdrsl,/create
spawn,'gzip -f '+outfile

; now do R
color='R'
filer=base+'R.fits*'
slit=mrdfits(filer,1,hdrsl,/silent)  
outfile=strcompress(base+'R.modified.fits',/remove)

fixed=resubtract_work(slit,hdrsl,color,manual=manual,multiplier=multiplier,nrows=nrows,plot=plot)

mwrfits,fixed,outfile,hdrsl,/create
spawn,'gzip -f '+outfile

return

end
