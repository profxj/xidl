;+
; NAME:
;    quickcalc
;
; PURPOSE:
;    Perform measurements from quicklook frames
; CALLING SEQUENCE:
;    quickcalc, inputfile
;
; INPUTS:
;    inputfile -- name of science frame previously analyzed with quicklook.pro
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;    outdir -- directory where data were output to (current if not specified)
;    
; EXAMPLES:
;
; COMMENTS:
;    generates several plots (should be only 1)
;    generate calibSlit files first, then call this routine! 
;    Input data comes from current directory.
;    Output data goes into current directory.
;
; REVISION HISTORY:
;   MD, 6Aug2002
;   extensive hacking aferward by JAN, AC, MD
;
;----------------------------------------------------------------------
pro quickcalc,   inputfile,  outdir=outdir

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  calib_data = getenv('CALIB_DATA')+'/'

;  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
   if deimos_data eq '/' then deimos_data = ''

  if calib_data eq '/' then message, 'You need to set $CALIB_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'

  if n_elements(outdir) eq 0 then outdir = ''
  if outdir ne '' then begin
     if strmid(outdir, strlen(outdir)-1, 1) ne '/' then outdir = outdir+'/'
  endif

  outdatadir = outdir

  objectcat = deimos_tables(inputfile, bluslits=slitcoords, $
                 slitobjmap=slitobjmap) 
  slitnums = slitcoords.slitno

  nobj = n_elements(objectcat)
  nslit = n_elements(slitcoords)
 
  isastar = slitcoords.slitwidth gt 2*median(slitcoords.slitwidth)
  xmms = slitcoords.xmm
  ymms = slitcoords.ymm
  
  slitidxforobjects = fltarr(nobj)
  objidxforslits = fltarr(nslit)
  mags = fltarr(nslit)+99.9999


  for i=0, nobj-1 do begin & $
     slitid = slitobjmap[i].dslitid & $
     slitidxforobjects[i] = where(slitid eq slitcoords.dslitid) & $
     if objectcat[i].mag lt mags[slitidxforobjects[i]] then begin & $
        mags[slitidxforobjects[i]] = objectcat[i].mag & $
        objidxforslits[slitidxforobjects[i]] = i & $
     endif & $
  endfor

  objectnums = float(objectcat[objidxforslits].object)
  issky = (objectnums MOD 1E6) ge 5E5
  lengths = slitobjmap[objidxforslits].botdist+slitobjmap[objidxforslits].topdist


  inputstem = strmid(inputfile, strpos(inputfile,'.fits')-10, 10)

;7590-7660
  quickfiles = findfile(outdir+'quick*.'+inputstem+'*')

 
;open PS file
  psname = 'ql.' +inputstem+ '.ps'
;-----
; PostScript

    !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
    set_plot, "PS"
    xsize= 6.0 & ysize= 8.0
    device, file=psname,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
    !P.THICK= 1.0
    !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
    !P.CHARSIZE= 1.0
;    !P.TITLE= 'Aband test results'
    !X.STYLE= 1
    !X.TITLE= ''
    !X.MARGIN= [9,0]
;    !X.MARGIN = [0, 0]
    !X.OMARGIN= [0, 0]
    !X.CHARSIZE= 1.0; 0.01
;    !X.TICKLEN= 0.04
    !Y.STYLE= 1
    !Y.TITLE= ''
    !Y.MARGIN= [0,0]
    !Y.OMARGIN= [8,4]
    !Y.CHARSIZE= 1.0
    !P.MULTI= [0,1,4]
    

; do alignment stars first!
  slitnumpos = strpos(quickfiles, 'B.fits') > strpos(quickfiles, 'R.fits') 
  nquick = n_elements(quickfiles)
  slitnums = fltarr(nquick)
  for i=0, nquick-1 do slitnums[i] =fix( strmid(quickfiles[i], $
                      slitnumpos[i]-3, 3))

  quickstars = isastar(slitnums)
  quickfiles = quickfiles(sort(10-1.*quickstars))

  starcount = 0
  for i=0, n_elements(quickfiles)-1 do begin

     spslit = mrdfits(quickfiles[i], 1, header, /silent)
     
     slitn = sxpar(header, 'SLITNO')
     print, '------------------------------------'
     print, 'processing slit: ', slitn
     print, 'processing file: ', quickfiles[i]

     slitwid = sxpar(header, 'SLITWID')
     thisisastar = isastar(slitn)
     if thisisastar then print, 'Alignment star!'
     mymag = mags[slitn]
     print
     print, 'R mag.: ', mymag
     mylength = lengths[slitn]
;     print, 'slit length: ', mylength
     myx = xmms[slitn]
     myy = ymms[slitn]
;     print, 'xmm/ymm: ', myx, myy
;     print
;     print,  "To view 2d rectified spectrum type: spec=mrdfits( '"+quickfiles[i]+ "' ,1)"
;      print,  'and then: atv,spec.flux'
  
; look for object, while rejecting cosmic rays, using the badmask,
; using the tilts of the lines, and subtracting out the sky with a mode
 
; alternatively, we could do the following stuff at the end, after
; files are written out.  Would probably make a fair amount of sense


    if thisisastar then begin
       lightprofile = find_object(spslit, profivar=ivarout, npix=npixout, $
                        /CR, /BPM, /NOSUBTRACT, /USETILT)
       nrows = n_elements(lightprofile)
       peakinfo, lightprofile, pkcol, fwhm, pk_cent=pkcent, $
          profivar=ivarout, window=7, s2n=sig2noise
       pixel_shift = fix(pkcol-nrows/2) 
;shift of peak pixel from center of slitlet, location of star
        sig2noiseperpix = sig2noise[0]/sqrt(npixout[pkcol]) ; typical signal-to-noise in an extracted spectrum, off sky lines, per (wavelength) pixel.  

 
      lightprofile = lightprofile-min(lightprofile(5:nrows-6))
       wavelengths = lambda_eval(spslit)
;       if max(wavelengths) gt 7590 or min(wavelengths) lt 7660 then $
; make it cover the ENTIRE A-band!
       if max(wavelengths) gt 7630 AND min(wavelengths) lt 7590 then begin
         starcount = starcount + 1
         starexam, quickfiles[i], mymag, spec=spec
         stats = {mag:mymag, lshift:spec.shift, xshift:pixel_shift[0], $
                   xpeak:spec.xpeak, fwhm:fwhm[0], s2n:sig2noiseperpix[0], $
                   transparency:spec.transparency } 
;results of positional analysis
         star_stats = (starcount eq 1) ? stats : [star_stats, stats] ;cumulant 

;         print
;         print, 'To do the A band test again on just this file, type: '
;         print, "starexam, '"+quickfiles[i]+"'"
;plot results for aband stars only--

         plot, spec.lambda, spec.spec, ytit='flux', $
            title='Aband region for file: '+ $
            inputstem + '   Slit: ' + string(slitn, format='(i3.3)')
         oplot, spec.lambda, spec.template, line=2
         xpos = spec.lambda[0] +10.
         dypos = (max(spec.spec) - min(spec.spec))/10.
         ypos = max(spec.spec) - dypos
         xyouts, xpos, ypos, 'R mag of star: '+string(stats.mag, $
            format='(f6.2)'), charthick = 2, charsize = 1.0
         xyouts, xpos, ypos-dypos, 'pixel shift, peak X-corr: ' $
            +strmid(strcompress(spec.shift,$
            /remove_all), 0, 5)+'  '+strmid(strcompress(spec.xpeak, $
              /remove_all), 0, 4),  charsize = 1.0, charthick = 2
         xyouts, xpos, ypos-2.*dypos, 'Spatial shift (pixels): ' + $
              string(stats.xshift, format='(f5.1)'), charthick=2
         xyouts, xpos+50., ypos, ' SNR: ' + string(stats.s2n, $
           format= '(f6.1)') +  ' FWHM (pixels): ' + $
              string(stats.fwhm, format='(f5.1)'), charthick=2
         xyouts, xpos, ypos-8.5*dypos, 'Transparency: ' + $
            string(stats.transparency, format='(f5.2)'), charthick=2

       endif else print, 'Spectrum does not include the (entire) A-band'
     endif else begin
       lightprofile = find_object(spslit, profivar=ivarout, npix=npixout, $
                    /CR, /BPM, /MMM, /USETILT)
      nrows = n_elements(lightprofile)
     endelse   

        print, 'Spatial profile information:'
        peakinfo, lightprofile, pkcol, fwhm, pk_cent=pkcent, $
          profivar=ivarout, window=7, s2n=sig2noise
       
                                ; take only the most significant peak.
                                ; The quantities we have measured are:

        pkcol = pkcol[0]  ; location of maximum in the profile along the
                             ;slit (spatial direction), in pixels
        fwhm = fwhm[0]    ; peak fwhm, in pixels
        pkcent = pkcent[0] ; centroid of the peak profile
        sig2noise = sig2noise[0] ; S/N for detecting the object, given 
                                ;an ivar-weighted integration over the
                                ;7-spatial-pixel window
                                ;defined above
        sig2noiseperpix = sig2noise[0]/sqrt(npixout[pkcol]) ; typical signal-to-noise in an extracted spectrum, off sky lines, per (wavelength) pixel.  


; may want to dump the lightprofile to a PS file, along with the
; predicted location of the peak (yobj, above)
;        print, 'predicted position: ', spslit.yobj
;        print, 'actual location of max & centroid: ', pkcol, pkcent
        offset = spslit.yobj-pkcent
        offset_arcsec = offset*0.119
; note redundancy here, but we needed to get this info for the
; alignment stars,and they are done first!
;        print, '**** Need to SHIFT SPATIALLY in X by '+ $
;            strcompress(offset_arcsec)+' arcsec'
;        print, 'Signal to noise of detection: ', sig2noise
        print, 'Typical S/N/pixel in extracted spectrum: ', $
          sig2noiseperpix
        print, 'Object FWHM along slit (pix): ', fwhm

; we might want to save arrays of these numbers in future, etc. 

; make sure don't go off the slit, and have some room for sky
;        exfwhm = fwhm < (pkcol-6)
;        exfwhm = exfwhm < (nrows-pkcol-7)
;        exfwhm = exfwhm >  1
;
;      extract = extract_1d(spslit, pkcol-exfwhm, pkcol+exfwhm, 4, nrows-5)    
   ;         pkcol-2*fwhm > 3, pkcol+2.*fwhm < nrows-4 )
      
;      mwrfits, extract, quickfiles[i]
;      print
;      print, 'Added simply-extracted spectrum to '+quickfiles[i]+' as HDU 2'
;      print, "To plot extracted spectrum type: extract=mrdfits( '"+quickfiles[i]+ "' ,2)"
;      print, 'and then: splot, extract.lambda, extract.flux'
;      print, ' ( splot,extract.sky would plot sky spectrum)'
;      names=strsplit(quickfiles[i],'.',/extract)
;      ps_open, 'ps/extracted_spectrum_'+names[3]
;      fname = 'ps/extracted_spectrum_'+names[3]+'.ps'
;      set_plot,'PS'
;     device,file=fname,/portrait, /color, xsize=6, ysize=4, /inch, $
;       SET_FONT='Courier', /TT_FONT,SET_CHARACTER_SIZE=[140,180]

;      plot, extract.lambda,  extract.spec, yr = [min(extract.spec) <  (-20),max(extract.spec) < 70]
;      ps_close
;      device,/close_file ; close off the ps file
;      set_plot,'X'
;      mwrfits, lightprofile, quickfiles[i]
;      print
 ;     print, 'Added light profile  along slit to '+quickfiles[i]+' as HDU 3'
 ;     print, "To plot light profile type: splot, mrdfits( '"+quickfiles[i]+"' ,3)"
;      print
;      print
;      print
  endfor

; generate average statistics for alignment stars
     nalign = n_elements(star_stats)
     print
     print,  '****CUMULATIVE STATISTICS FOR ALIGNMENT STARS***
     print, 'Mean spatial shift (pixels): ', mean(star_stats.xshift)
     print, 'Mean spectral shift (pixels): ', mean(star_stats.lshift)
     print, 'Mean transparency: ', mean(star_stats.transparency)
     print, 'Mean FWHM: ', mean(star_stats.fwhm)

  device, /close ;close PS file
  set_plot, 'X'
  !p.multi = [0, 0, 0]

  return
end











