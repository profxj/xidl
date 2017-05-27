;+
; NAME:
;   deimos_vigflat
;
; PURPOSE:
;   Generate vignetting correction flat
;
; CALLING SEQUENCE:
;   deimos_vigflat, files,path
; 
; INPUTS:
;   files  - list of files to use in making flat
;   path   - path for output file
;
; OUTPUTS:

;
; MODIFICATION HISTORY:
;    14-Oct-2002 JAN - first version
;-
pro deimos_vigflat, files, path


  if n_elements(path) eq 0 then path = ''
  if n_elements(nofix) eq 0 then nofix=1

  nfiles = n_elements(files)


  if NOT keyword_set(path) then message, 'please call with a path!', /informational

  
  outfilename = path+'vignetting_image.fits'

  calib_data = getenv('CALIB_DATA')+'/'
  badmaskname = calib_data+ 'deimos_badmask.fits.Z'
  pixflatname = calib_data+'processed_pix_mult_flat.fits.gz'
  outpixflat = path+'processed_pix_mult_flat.fits'

  mwrfits, junk, outfilename, /create ;dummy primary
  mwrfits, junk, outpixflat, /create ;dummy primary
  

; TO DEBUG CODE: DO ONLY 1 CHIP
  for i=1,8 do begin

     deimos_badchip = mrdfits(badmaskname, i, /silent) 
     deimos_pixflat = mrdfits(pixflatname, i, /silent) 

     mkhdr, header, deimos_pixflat
     mosaic_keywords, header, i
     mwrfits, deimos_pixflat, outpixflat, header; chip8
     whneg=where(deimos_pixflat eq -1.,negct)
     if negct gt 0 then deimos_pixflat[whneg]=0.

 ; set up vignetting mask.  this is done now to
 ; keep djs_maskinterp from doing
 ; anything in the vignetted regions.
   
        mask = (deimos_badchip AND 2B eq 2b)

	nkernel=25
	vigkernel=intarr(nkernel)+1


	vigarr=(deimos_badchip AND 2b) eq 2b

	badlyvig=erode(vigarr,vigkernel) AND erode(vigarr,transpose(vigkernel))


	; check to see what, if any, corners and sides are vignetted
	llvig=i eq 1
	ulvig=i eq 5
	lrvig=i eq 4
	urvig=i eq 8

	cornerdel=nkernel*sqrt(2)-5

	leftvig=(i eq 1 OR i eq 5)
	rightvig=(i eq 4 or i eq 8)

	botrange=[0*llvig+(2047-1023+cornerdel)*lrvig, $
		(1023-cornerdel)*llvig +2047*lrvig]

	toprange=[0*ulvig+(2047-1023+cornerdel)*urvig, $
		(1023-cornerdel)*ulvig +2047*urvig]
	
	if leftvig then leftrange=[0,4095] else leftrange = $
		[0*llvig+(4095-1023+cornerdel)*ulvig, $
		(1023-cornerdel)*llvig +4095*ulvig]
	if rightvig then rightrange=[0,4095] else rightrange = $
		[0*lrvig+(4095-1023+cornerdel)*urvig, $
		(1023-cornerdel)*lrvig +4095*urvig]
		

	badlyvig[botrange[0]:botrange[1],0:nkernel/2]= $
		vigarr[botrange[0]:botrange[1],0:nkernel/2]	
	badlyvig[toprange[0]:toprange[1],4095-nkernel/2:4095]= $
		vigarr[toprange[0]:toprange[1],4095-nkernel/2:4095]	
	badlyvig[0:nkernel/2,leftrange[0]:leftrange[1]]= $
		vigarr[0:nkernel/2,leftrange[0]:leftrange[1]]	
	badlyvig[2047-nkernel/2:2047,rightrange[0]:rightrange[1]]= $
		vigarr[2047-nkernel/2:2047,rightrange[0]:rightrange[1]]	

; done setting up vignetting mask for the pixflat code


     data = fltarr(2048, 4096, nfiles)
;     mdata = data

     ; DIVIDE COLUMNS INTO N SECTIONS FOR THE "FUNNYCORR": THE CORRECTION FOR 
     ; VIGNETTING OF OBJECTS THAT ONLY AFFECT THE PIXFLAT
	; Make sure N divides 4096!
     nsection=16

     funnycorr=fltarr(2048,nsection)


     for j=0, n_elements(files)-1 do begin
           ; read in frame
        tempimage = deimos_read_chip(files[j], i)
	
           ; interpolate over bad, but NOT vignetted, regions
        data[*, *, j] = djs_maskinterp(tempimage, $
               (deimos_badchip AND 1B ) eq 1, iaxis=0)
	delvarx,tempimage
	; funny things happen at the edge in 1/8 and 4/5.  But the 
	; funny business on 3/7 is close enough we don't want to affect it.
	; this seems like the solution for now.

	fixpix = 5; 11+50*(i eq 1 or i eq 4 or i eq 5 or i eq 8)
	data[*,*,j]=deimos_fixglow(data[*,*,j],nfix=fixpix)
	nfix=fixpix ; me and my memory...

; OK, this was a hard choice.  But it looks like the vertical features
;I'm seeing in chip 3 are real.  They are present in the raw frames,
;that is sure.  What I'm going to do now is take it out from the
;median, to make the next processing steps work, but keep it in the
;final pixflat.  Hopefully that will prove to be good.  The 'features'
;don't seem to show up in real flats.  According to Sandy, it's due to
;the data having been taken in slitless mode - so we do have to
;correct them.
	
; Unfortunately, the effect seems to depend slightly on row #, so we can't
; just use one fix for the whole column.  The below code does the job, though.

	for section=0,nsection - 1 do begin

		funnystart=fltarr(2048)-1

		for k=0,2047 do begin
                   colmask = deimos_badchip[k, *]
; watch out for vignetted regions
                   usepix = findgen(4096) gt section*(4096/nsection) $
                     AND findgen(4096) lt (section+1)*(4096/nsection)-1 $
                     AND colmask eq 0
                   if total(usepix) ge 1 then funnystart[k]= $
			median(data[k,where(usepix),j]) 
                endfor

                nodata = funnystart eq (-1)
                whnodata= where(nodata, nodatact)
                funnystart = djs_maskinterp(funnystart, nodata)
                if nodatact eq 4096 then funnystart=fltarr(2048)+1

	; need to make sure we don't throw off the median-smooth normalization
	; from the funny region

		normcorr=funnystart / djs_median(funnystart,width=251,$
		 boundary='reflect')

                normcorr = djs_maskinterp(normcorr, nodata)

		xindx = findgen(2048)
		;nbuf=100

		weirdpix=(normcorr lt .995 OR normcorr gt 1.01)

	; make sure we aren't influenced by the sharp-gradient (vignetted?)
	; regions
		;bufcols=((i eq 1 OR i eq 5) AND xindx lt nbuf) OR $
		;   ((i eq 8 OR i eq 4) AND xindx gt 2048-nbuf)
		weirdpix=weirdpix OR  nodata ; OR bufcols


	; want to make sure the odd regions don't throw off the median levels
	; expand by +/- 3 pix
		weirdplus=dilate(weirdpix,fltarr(7)+1) OR weirdpix
		weirdplus[0:nfix]=0
		weirdplus[2048-nfix-1:2048-1]=0
		tempcorr=djs_maskinterp(funnystart, weirdplus)

; the width in the next line doesn't matter much - 51 or 251 is nearly
; indistinguishable.  We just need to take out the continuum level.
; these few lines of code would have saved me endless hours back when
; I was working for alex (of course, he'd have to have used IDL in those days
; instead of LOLITA, too.  sigh...)

		funnycorr[*,section]=funnystart / $
			djs_median(tempcorr,width=201,$                  
	                boundary='reflect')
                if nodatact gt 0 then funnycorr[whnodata, section] = 1.
	endfor



	; ONLY apply a correction to the 'weird' pixels

	discrepcol=(abs(funnycorr-1.) gt 0.005)

        ; presume that isolated deviations are spurious

	discrepcol=(discrepcol AND erode(discrepcol,[1,1,1]))
	
	; dilate does some funny things, have to deal with it
	discrepcol[0:nfix,*]=1
	discrepcol[2047-nfix:2047,*]=1
	discrepcol=dilate(discrepcol,fltarr(9)+1) OR discrepcol

;	bufarr=bufcols # (fltarr(nsection)+1.)
;	discrepcol=discrepcol OR bufarr

	print,total(discrepcol)/float(nsection),' columns corrected'
	wh=where(NOT discrepcol,whct)
	if whct ge 1 then funnycorr[wh]=1.

	funnyfix=fltarr(2048,4096)
	sectionrows=4096/2./nsection+4096/nsection*findgen(nsection)

	columns=lindgen(2048,4096) MOD 2048
	rows=lindgen(2048,4096)/ 2048
	
	funnyfix=interpolate(funnycorr, columns,rows*float(nsection)/4096.-0.5)
	
; /missing didn't seem to work right - force it to behave

	for k=0,4096./2./nsection do funnyfix[*,k]=reform(funnyfix[*,4096./2./nsection+1])
	for k=4096-4096./2./nsection-1,4095 do funnyfix[*,k]=$
		reform(funnyfix[*,4096-4096./2./nsection-2])

;	savefunny=funnyfix


	data[*,*,j]=data[*,*,j]/funnyfix
;	savedata=data[*,*,j]
;	save,savefunny,savedata,discrepcol,f='tempfunny.sav'

;tempimage=djs_maskinterp(data[*,*,j], $
;              vigarr, iaxis=0,/const)

           ; make median-smoothed version for rescaling
;       tempmed = djs_median(tempimage, width=31, $
;                boundary='reflect')
;delvarx,tempimage

;discreppix=(data[*,*,j]/tempmed lt 0.8 AND NOT vigarr)

; mask out regions around 'dust' for good measure
;dilation=2.
;discreppix=dilatemask(discreppix,dilation)
	
;meddata=djs_maskinterp(data[*,*,j], $
;              discreppix OR badlyvig, iaxis=0)	
;mdata[*,*,j]=djs_median(meddata, width=41, $
;                boundary='reflect')	
     endfor

	delvarx,discreppix,badlyvig,vigarr,meddata

     ndata = data*(deimos_pixflat+1)
     delvarx,data   ; try to conserve memory!
     invdata = 1./ndata   ; == mdata^2/data


 ;    whinf = where(finite(ndata) eq 0, infct)
 ;    if infct ne 0 then ndata[whinf] = 0.
 ;    if infct ne 0 then invdata[whinf] = 0.

;     whbad = where(finite(invdata) eq 0,  badct)
;     if badct ne 0 then invdata[whbad] = 0.

     ; do the avsigclip to build the pixflat
     clipdata = avsigclip(ndata, invdata, 20., 1.5)
     flat = (clipdata.flux)

     delvarx, clipdata 
;     flat = 1./flat
;     whinf = where(finite(flat) eq 0, infct)
;     if infct ne 0 then flat(whinf) = 1.

;     whbad=where((deimos_badchip ne 0) ,badct)
;     if badct ne 0 then flat(whbad)=1.

;     whzero = where(flat eq 0, zeroct)
;     if zeroct ne 0 then flat(whzero) = 1.

	; can't quite trust the outer edges, it looks like.  Err on the side of not correcting them
	; rather than miscorrecting them.
		
	mwrfits, flat, outfilename, header

; USING COMPRESSED FILE - REPLACE WITH THIS TO REVERSE
        ;    mwrfits, (flat*mask-1.), maskfilename 
            print, 'wrote chip:', i
  endfor  

; USING COMPRESSED FILE - COMMENT TO REVERSE
  spawn, 'gzip -f '+outfilename

  return 
end







