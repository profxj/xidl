;+
;
; NAME
;       quality_test.pro
;
; PURPOSE
;       Create plots detailing the quality of the data for a given
;       mask. 
;
;
; SYNTAX
;       quality_test, result=result, /noplot, /terse, /silent
;
; INPUTS
;       None. But note that the routine should be run in the directory
;       which contains the spec1d.xxxx.nnn.oooooooo.fits files and the
;       obj_info.xxxx.fits file.
;
; KEYWORDS
;       /noplot = set this keyword and the post-script file is not
;                 created. 
;       /terse = set this keyword and only a subset of the plots are
;                made.
;       /silent = set this keyword to stop the routine from spitting
;                 out unwanted warning messages.
;       /speedy = if this keyword is set, then the signal-to-noise
;                 values for each slit will not be determined
;                 empirically (by determining the mean and std dev of
;                 the continuum counts) 
;
; OUTPUTS
;       result = a variable which will be returned as an array of 2
;                elements. The first element in the array will be the
;                s/n value (determined from integrating over the fwhm
;                in each object's spatial profile) at an R magntiude
;                of 23.5. This value is determined by doing a
;                median-median fit to the s/n versus Rmag
;                distribution. The second element is exactly the same
;                except that the s/n in this case is the
;                signal-to-noise as derived empircally from the 1d
;                spectra (by calculating a mean and standard deviation
;                in a given lambda window).
;
;       Also, a single post-script file containing several plots
;       detailing the signal-to-noise (quality) of the data for the
;       given mask. This file will be written to the ps/ subdirectory,
;       thus that subdirectory needs to exist!
;
; PROCEDURES CALLED 
;       fill_gap.pro
;
; COMMENTS
;
;
; EXAMPLES
;       None.
;
; HISTORY
;       Created October 1, 2002 by mcc. 
;
;-


pro quality_test, result=result, noplot=noplot, silent=silent, $
                  terse=terse, speedy=speedy, deep=deep
; find the object info file.
  objfile = findfile('obj_info*.fits*', count=nfiles)
  if nfiles gt 0 then objfile = objfile[0] $
    else message, 'ERROR: No obj_info.xxxx.fits file found!'
; read in the object info file.
  objinfo = mrdfits(objfile, 1, objhdr, /silent)

; find the bin table file.
  binfile = findfile('*.bintabs.*fits*', count=nfiles)
  if nfiles eq 0 then $
    message, 'ERROR: No xxxx.bintabs.fits file found!' $
  else begin
; check if the bin table file is compressed.
      binfile = binfile[0]
      len = strlen(binfile)
      IF strmid(binfile, len-3, 3) eq '.gz' then begin
          spawn, 'gunzip -v ' + binfile
          binfile = strmid(binfile, 0, len-3)
      endif
  endelse
; check if this is a DEEP2-1HS mask.
  fits_open, binfile, fcb
  maskdex = where(fcb.extname eq 'MaskDesign', maskcnt)
  if maskcnt eq 0 then $
    message, '(quality_test.pro) ERROR: No MaskDesign bin table found!'
  masktab = mrdfits(binfile, maskdex[0], /silent)
  project = strcompress(masktab.PROJNAME, /rem)
  if project eq 'DEEP2-1HS' then isdeep = 1 else isdeep = 0
; read in the ObjectCat table from the bin table file.
  catdex = where(fcb.extname eq 'ObjectCat', catcnt)
  if catcnt eq 0 then $
    message, '(quality_test.pro) ERROR: No ObjectCat bin table found!'
  objcat = mrdfits(binfile, catdex[0], /silent)

; the obj_info files contains the s/n numbers for each object as
; derived from the spatial profile (recall that the obj_info file is
; created by the routine do_extract.pro). let's also find the s/n for
; each object in a more empirical manner by determining the mean and
; standard deviation of the 1-d spectra in a given lambda range. only
; this part of the code is SLOW, so only proceed if the /speedy
; keyword is not set.
  if not(keyword_set(speedy)) then begin
      maskname = sxpar(objhdr, 'SLMSKNAM')
; if this is a DEEP2 mask, then trim the maskname.
      if isdeep and strlen(maskname) gt 4 then $
        maskname = strmid(maskname, 0, 4)
; construct a list of all of the spec1d file names using the
; information stored in the objinfo file.
      files = 'spec1d.' + maskname + '.' + $
        string(objinfo.slitno, format='(I3.3)') + '.' + $
        strcompress(objinfo.objno, /rem) + '.fits'
      nfiles = n_elements(files)
; set a threshold value determining the number of pixels needed in the
; wavelength region (given by blam) in order to do statistical
; test. if there isn't enough pixels, then the s/n for that slitlet
; will simply be taken to be 0.
      thresh = 500. ;pixels
; define the wavelength regions in which we will measure s/n. these
; regions are chosen to avoid sky lines but still provide a large
; enough number of pixels (6600 < lambda < 6830 and 6960 < lambda <
; 7250).
      blam = [6600., 6830., 6960., 7250.] ;angstroms 
; define the output array.
      s2n_emp = fltarr(nfiles)
      for i=0,nfiles-1 do begin
          ss1d = fill_gap(files[i], /optimal)
          if size(ss1d, /tname) eq 'STRUCT' then begin
              region = where((ss1d.lambda gt blam[0] and ss1d.lambda lt blam[1]) or $
                             (ss1d.lambda gt blam[2] and ss1d.lambda lt blam[3]), cnt)
              if cnt gt thresh then begin
                  djs_iterstat, ss1d.spec[region], mean=signal, sigma=noise, $
                    sigrej=3.0
                  s2n_emp[i] = signal / noise
              endif
          endif else $ 
              if not(keyword_set(silent)) then $
                print, 'ERROR: possible wavelength overlap in file ' + files[i]
      endfor
  endif

; now sort the ObjectCat according to object number (going object by
; object), and thereby extract the magnitude for each object.
; here we will construct 4 arrays: magB, magR, dexB, and dexR. recall
; that the objinfo file likely has 2 entries per object: 1 for each
; the Blue slitfile and the Red slitfile. 
  rflag = 0
  bflag = 0
  nobj = n_elements(objinfo.objno)
  for i=0,nobj-1 do begin
; for the ith object in the objinfo file, find the corresponding entry
; in the objcat structure. note that serendips will NOT be found in
; the objcat, so exclude them. also exclude alignment stars.
      dex = where(strcompress(objcat.object, /rem) eq $
                  strcompress(objinfo[i].objno, /rem) and $
                  strmid(objcat.object, 2, 1) ne '5' and $
                  objinfo[i].objtype ne 'A', cnt)
      if cnt gt 0 then begin
          if cnt gt 1 then $
            if not(keyword_set(silent)) then $
            print, '(quality_test.pro) WARNING: Multiple objects ' + $
            'with same object number!'
; if the ith entry in the objinfo file corresponds to a Blue slitfile,
; then add the magnitude value to the magB array and vice versa for
; the Red files.
          if objinfo[i].color eq 'B' then begin
              if bflag eq 0 then begin
                  magB = objcat[dex[0]].mag
                  dexB = [i]
                  bflag = 1
              endif else begin
                  magB = [magB, objcat[dex[0]].mag]
                  dexB = [dexB, i]
              endelse
          endif else begin
              if rflag eq 0 then begin
                  magR = objcat[dex[0]].mag
                  dexR = [i]
                  rflag = 1
              endif else begin
                  magR = [magR, objcat[dex[0]].mag]
                  dexR = [dexR, i]
              endelse
          endelse
      endif else if not(keyword_set(silent)) then $
        print, '(quality_test.pro) WARNING: No magnitude found ' + $
        'for object ' + strcompress(objinfo[i].objno, /rem) + '!!!'
  endfor

; use these indexes (dexB and dexR) that we have built up to extract
; the s/n info for all of the objects for which a R magntiude was
; found. that is, to exclude the serendips (as well as alignment stars
; and sky-only slits)!
  s2nB = objinfo[dexB].s2n_fwhm
  s2nR = objinfo[dexR].s2n_fwhm
  s2nB_det = objinfo[dexB].s2n_window
  s2nR_det = objinfo[dexR].s2n_window
; trim the empirical s/n array like we trimmed the other arrays
; (excluding serendips; objects for which we don't have a magnitude in
; the pcat). 
  if not(keyword_set(speedy)) then s2n_emp = s2n_emp[dexB]


;----------------------------------
; plot the results in a single post-script file.
  if not(keyword_set(noplot)) then begin
      set_plot, 'ps'
      device, file='ps/quality_plots.ps', /landscape, /color
      tname = 'mask ' + strcompress( sxpar(objhdr, 'SLMSKNAM'), /remove_all)
      if not(keyword_set(terse)) then begin
; plot histogram of s/n (detection) distribution: s/n over 5 pixel
; window...for both blue and red portions of slit.
          xname = 'detection signal-to-noise'
          subname = 'signal-to-noise from spatial profile'
          binsz = 5.
          plothist, s2nB_det, xhistB, yhistB, bin=binsz, /noplot
          plothist, s2nR_det, xhistR, yhistR, bin=binsz, /noplot
          ymax = MAX([yhistB, yhistR])
          plothist, s2nB_det, bin=binsz, /xsty, xr=[0, 120], xtitle=xname, $
            thick=2, xthick=1.5, ythick=1.5, /ysty, yr=[0, 1.1*ymax], $
            title=tname, subtitle=subname
          plothist, s2nR_det, bin=binsz, line=2, thick=2, /overplot
          legend, ['Solid Line = Blue side of slit', $
                   'Dashed Line = Red side of slit'], /right_legend
      endif
; plot histogram of s/n (fwhm) distribution...for both the blue and
; red portions of slit.
      xname = 'signal-to-noise (per pixel)'
      subname = 'signal-to-noise from spatial profile'
      binsz = 0.25
; catch for places where the signal-to-noise might be a non-finite
; value...just in case.
      nonfin = where(finite(s2nB) eq 0 or finite(s2nR) eq 0, cnt)
      if cnt gt 0 then s2nB[nonfin] = 0.
      if cnt gt 0 then s2nR[nonfin] = 0.
      plothist, s2nB, xhistB, yhistB, bin=binsz, /noplot
      plothist, s2nR, xhistR, yhistR, bin=binsz, /noplot
      ymax = max([yhistB, yhistR])
      plothist, s2nB, bin=binsz, /xsty, xr=[-0.5, 10], xtitle=xname, $
        thick=2, xthick=1.5, ythick=1.5, /ysty, yr=[0, 1.1*ymax], $
        title=tname, subtitle=subname
      plothist, s2nR, bin=binsz, line=2, thick=2, /overplot
      legend, ['Solid Line = Blue side of slit', $
               'Dashed Line = Red side of slit'], /right_legend
; plot histogram of s/n (empirical) distribution...calculated on the
; blue side of slit (6600A - 7300A).
      if not(keyword_set(speedy)) then begin
          xname = 'signal-to-noise (per pixel)'
          subname = 'signal-to-noise from 1-d spectra'
          binsz = 0.5
          plothist, s2n_emp, bin=binsz, xtitle=xname, thick=2, $
            xthick=1.5, ythick=1.5, title=tname, subtitle=subname, $
            /xsty, xr=[-0.5, 10]
; plot s/n (emp) versus s/n (fwhm).
          plot, s2n_emp, s2nB, psym=2, xthick=1.5, $
            ythick=1.5, xtitle='s/n from 1-d spectra (per 1d pixel)', $
            ytitle='s/n from spatial profile (per 1d pixel)', $
            /xsty, /ysty, xr=[-0.5, 2.5], yr=[0.0, 3.0], title=tname
          oplot, findgen(10)-5, findgen(10)-5, line=2, thick=2
      endif

; plot s/n (detection) versus R magnitude.
      if not(keyword_set(terse)) then begin
          xname = strcompress(objcat[0].pband, /rem) + ' magnitude'
          yname = 'signal-to-noise (of detection)'
          plot, magb, s2nB_det, PSYM=6, xtitle=xname, ytitle=yname, $
            /xsty, /ysty, xr=[20, 24.6], yr=[0, 120], xthick=1.5, $
            ythick=1.5, title=tname
          oplot, magr, s2nR_det, PSYM=4
          legend, ['Square = Blue side of slit', $
                   'Diamond = Red side of slit'], /LEFT_LEGEND
      endif
; plot s/n (fwhm) versus R magnitude.
      xname = strcompress(objcat[0].pband, /REMOVE_ALL) + ' magnitude'
      yname = 'signal-to-noise (per pixel)'
      plot, magb, s2nb, PSYM=6, xtitle=xname, ytitle=yname, $
        /xsty, /ysty, xr=[20, 24.6], yr=[-0.25, 6], xthick=1.5, $
        ythick=1.5, title=tname
      oplot, magr, s2nr, PSYM=4
      legend, ['Square = Blue side of slit', $
               'Diamond = Red side of slit'], /LEFT_LEGEND
; plot s/n (empirical) verus R magitude.
      if not(keyword_set(speedy)) then begin
          yname = 'signal-to-noise per pixel (empirical)'
          plot, magb, s2n_emp, psym=5, xtitle=xname, ytitle=yname, $
            /xsty, /ysty, xr=[20, 24.6], yr=[-0.25, 5], xthick=1.5, $
            ythick=1.5, title=tname
      endif
  endif

;---------------------------
; now fit the data s/n versus magnitude data with a median-median
; fit. use the s/n data from the blue side of the slit since it is
; closer to the R band magnitudes in the pcat.
  xdata = magB
  ydata = s2nB
; bin data into 2 bins.
  dex1 = where(magB gt 23.1 and magB le 23.6, cnt1)
  dex2 = where(magB gt 23.6 and magB le 24.1, cnt2)
  if cnt1 eq 0 or cnt2 eq 0 then begin
      if not(keyword_set(noplot)) then begin
          device, /close
          set_plot, 'x'
      endif
      result = [-999, -999]
      return
  endif
; take median values in both X and Y directions.
  x1 = median(xdata[dex1])
  x2 = median(xdata[dex2])
  y1 = median(ydata[dex1])
  y2 = median(ydata[dex2])
; fit a line to the median-median points.
  mcc_polyfit, [x1, x2], [y1, y2], [0, 1], a=a
  xfit = xdata[sort(xdata)]
  yfit = a[0] + xfit*a[1]
; evaluate the linear fit at R=23.5.
  result = [a[0] + (23.5)*a[1], -999]
; now determine distance of point from the fit line (at that
; magnitude). that is, find distance (vertically) in s/n.
  if not(keyword_set(noplot)) then begin
      subset = where(magB gt 23., subcnt)
      if subcnt lt 5 then $
        message, 'ERROR not enough points to continue!'
      fitvals = a[0] + xdata[subset]*a[1]
      dist = ydata[subset] - fitvals
; bin the data.
      binsz = 0.25
      plothist, dist, xhist, yhist, bin=binsz, /noplot
      nbins = n_elements(xhist)
      
; plot more results. 
; plot 1: the plot of s/n vs. R mag for all objects showing the
; median-median fit.
      yname = 'signal-to-noise (per pixel)'
      plot, xdata, ydata, psym=6, xtitle=xname, ytitle=yname, $
        /xsty, /ysty, xr=[20., 24.6], yr=[-0.25, 6], xthick=1.5, $
        ythick=1.5, title=tname
      oplot, xfit, yfit, thick=2
      legend, ['Square = Blue side of slit', $
               's/n (R=23.5) = ' + string(result[0], format='(F6.2)')], $
        /left_legend
; plot 2: plot of s/n vs. R mag for all object with the size of the
; symbol scaled according to distance from the fit line.
      plot, xfit, yfit, thick=2, xtitle=xname, ytitle=yname, $
        /xsty, /ysty, xr=[22.8, 24.2], yr=[-0.25, 4.5], xthick=1.5, $
        ythick=1.5, title=tname
      maxsz = 5.
; make a user-defined plotting symbol: an open circle.
      res = 32
      foo = findgen(res) * (!pi * 2 / (res-1))
      usersym, cos(foo), sin(foo)
      for i=0,nbins-1 do begin
          dex = where(dist lt xhist[i] and $
                      dist ge xhist[i]-binsz, cnt)
          if cnt gt 0 then $
            oplot, (xdata[subset])[dex], (ydata[subset])[dex], psym=8, $
            symsize=(i*maxsz/nbins+0.5)
      endfor
; plot3: plot of deimos x position vs. deimos y position for all the
; objects....use the same symbol size coding.
; read-in the bintab file to get the x and y positions.
      bludex = where(fcb.extname eq 'BluSlits', blucnt)
      desidex = where(fcb.extname eq 'DesiSlits', desicnt)
      if desicnt eq 0 or blucnt eq 0 then $
        message, '(quality_test.pro) ERROR: trouble in indexing ' + $
        'the bin table file!'
      bluslits = mrdfits(binfile, bludex[0], /silent)
      desislits = mrdfits(binfile, desidex[0], /silent)
      xpos = (bluslits.slitx1 + bluslits.slitx2 + $
              bluslits.slitx3 + bluslits.slitx4)/4.
      ypos = (bluslits.slity1 + bluslits.slity2 + $
              bluslits.slity3 + bluslits.slity4)/4
; match the bluslits structure to the desislits via the dslitid tag.
      xpos = xpos[sort(bluslits.dslitid)]
      ypos = ypos[sort(bluslits.dslitid)]
      slitname = desislits[sort(desislits.dslitid)].slitname
; we can now match these xpos and ypos values to our s/n numbers via
; the slitno.
      slitno = objinfo[dexB].slitno
      for i=0,n_elements(slitno)-1 do begin 
          mdex = where(string(slitname, format='(I3.3)') eq $
                       string(slitno[i], format='(I3.3)'), mcnt)
          if mcnt gt 0 then begin
              if keyword_set(srtord) then srtord = [srtord, mdex[0]] $
              else srtord = mdex[0]
          endif
      endfor
      xpos = xpos[srtord]
      ypos = ypos[srtord]
; as an aside, let's store this x and y positions in another
; variables. 
      xx = xpos
      yy = ypos
; finally make the freakin' plot!
      xname = 'Deimos x position'
      yname = 'Deimos y position'
      maxsz = 5.
      pflag = 0
      for i=0,nbins-1 do begin
          dex = where(dist lt xhist[i] and $
                      dist ge xhist[i]-binsz, cnt)
          if cnt gt 0 then begin
              if i ge nbins/4. and i lt nbins/3.*2. then hue = getcolor('red')
              if i ge nbins/3.*2. then hue = getcolor('green')
              if i lt nbins/4. or pflag eq 0 then hue = getcolor('black')
              if pflag eq 0 then begin
                  plot, (xpos[subset])[dex], (ypos[subset])[dex], $
                    xtitle=xname, ytitle=yname, $
                    /xsty, /ysty, xr=[-400, 400], yr=[0, 200], $
                    xthick=1.5, ythick=1.5, title=tname, $
                    symsize=(i*maxsz/nbins+1), psym=8, color=hue
                  pflag = 1
              endif else $
                oplot, (xpos[subset])[dex], (ypos[subset])[dex], psym=8, $
                symsize=(i*maxsz/nbins+0.5)  , color=hue 
          endif
      endfor
  endif

;----------------------------
; now do this all again: (1) do median-median fit and (2) plot
; results....but now for the empircal s/n numbers.
  if not(keyword_set(speedy)) then begin
      xdata = magB
      ydata = s2n_emp
; bin data into 2 bins.
      dex1 = where(magB gt 23.1 and magB le 23.6, cnt1)
      dex2 = where(magB gt 23.6 and magB le 24.1, cnt2)
      if cnt1 eq 0 or cnt2 eq 0 then begin
          if not(keyword_set(noplot)) then begin
              device, /close
              set_plot, 'x'
          endif
          result = [result, 0]
          return
      endif
      x1 = median(xdata[dex1])
      x2 = median(xdata[dex2])
      y1 = median(ydata[dex1])
      y2 = median(ydata[dex2])
      mcc_polyfit, [x1, x2], [y1, y2], [0, 1], a=a
      xfit = xdata[sort(xdata)]
      yfit = a[0] + xfit*a[1]
      result[1] = a[0] + (23.5)*a[1]
; now determine distance of point from the fit line (at that
; magnitude). that is, find distance (vertically) in s/n.
      if not(keyword_set(noplot)) then begin
          subset = where(magB gt 23., subcnt)
          if subcnt lt 5 then $
            message, 'ERROR not enough points to continue!'
          fitvals = a[0] + xdata[subset]*a[1]
          dist = ydata[subset] - fitvals
; bin the data.
          binsz = 0.25
          plothist, dist, xhist, yhist, bin=binsz, /noplot
          nbins = n_elements(xhist)

; plot more results. 
; plot 1: the plot of s/n vs. R mag for all objects showing the
; median-median fit.
          yname = 'signal-to-noise per pixel (empirical)'
          plot, xdata, ydata, psym=6, xtitle=xname, ytitle=yname, $
            /xsty, /ysty, xr=[20., 24.6], yr=[-0.25, 5], xthick=1.5, $
            ythick=1.5, title=tname
          oplot, xfit, yfit, thick=2
          legend, ['Square = Blue side of slit', $
                   's/n (R=23.5) = ' + string(result[1], format='(F6.2)')], $
            /left_legend
; plot 2: plot of s/n vs. R mag for all object with the size of the
; symbol scaled according to distance from the fit line.
          plot, xfit, yfit, thick=2, xtitle=xname, ytitle=yname, $
            /xsty, /ysty, xr=[22.8, 24.2], yr=[-0.25, 3.0], xthick=1.5, $
            ythick=1.5, title=tname
          maxsz = 5.
; make a user-defined plotting symbol: an open circle.
          res = 32
          foo = findgen(res) * (!pi * 2 / (res-1))
          usersym, cos(foo), sin(foo)
          for i=0,nbins-1 do begin
              dex = where(dist lt xhist[i] and $
                          dist ge xhist[i]-binsz, cnt)
              if cnt gt 0 then $
            oplot, (xdata[subset])[dex], (ydata[subset])[dex], psym=8, $
                symsize=(i*maxsz/nbins+0.5)
          endfor

; lastly, make a plot of s2n_emp versus x,y position on the detector.
;          binsz = 0.5
;          plothist, s2n_emp, xhist, yhist, bin=binsz, /noplot
          nbins = n_elements(xhist)
          pflag = 0
          for i=0,nbins-1 do begin
              dex = where(dist lt xhist[i] and $
                          dist ge xhist[i]-binsz, cnt)
              if cnt gt 0 then begin
                  if i ge nbins/4. and i lt nbins/3.*2. then hue = getcolor('red')
                  if i ge nbins/3.*2. then hue = getcolor('green')
                  if i lt nbins/4. or pflag eq 0 then hue = getcolor('black')
                  if xhist[i] lt 0 then symbol = 8 else symbol = 6
                  if pflag eq 0 then begin
                      plot, xx[subset[dex]], yy[subset[dex]], $
                        xtitle=xname, ytitle=yname, $
                        /xsty, /ysty, xr=[-400, 400], yr=[0, 200], $
                        xthick=1.5, ythick=1.5, title=tname, $
                        symsize=(i*maxsz/nbins+1), psym=symbol, color=hue
                      pflag = 1
                  endif else $
                    oplot, xx[subset[dex]], yy[subset[dex]], psym=symbol, $
                    symsize=(i*maxsz/nbins+0.5), color=hue
              endif
          endfor
      endif
  endif

;------------------------
; make plots of z-completeness versus chip position.
  if not(keyword_set(noplot)) and keyword_set(deep) then begin
; define the location of the zcat fits file.
      zfile = '/g/dsm/Redshifts/zcat.fits'
; read in the zcat from this file.
      zz = mrdfits(zfile, 1, /silent)
      spawn, 'pwd', fulldir
      fulldir = strsplit(fulldir, '/', /extract)
      mask = fulldir[n_elements(fulldir)-2]
      dex = where(strcompress(zz.maskname, /rem) eq mask, cnt)
      if cnt gt 0 then begin
          slitno = objinfo[dexB].slitno
          objno = objinfo[dexB].objno
          for i=0,n_elements(slitno)-1 do begin 
              mdex = where(string(zz[dex].slitname, format='(I3.3)') eq $
                           string(slitno[i], format='(I3.3)') and $
                           strcompress(zz[dex].objname, /rem) eq $
                           strcompress(objno[i], /rem), mcnt)
              if mcnt gt 0 then begin
                  if keyword_set(srtord2) then $
                    srtord2 = [srtord2, dex[mdex[0]]] $
                  else srtord2 = dex[mdex[0]]
              endif
          endfor
          zz = zz[srtord2]
      endif else begin
          print, 'No redshift info found for mask!'
          goto, finish
      endelse
; plot Q value as a function of s/n for each object.
      if not(keyword_set(speedy)) then begin
          xname = 'z quality'
          yname = 'signal-to-noise per pixel (empirical)'
          plot, zz.zquality, s2n_emp, psym=4, /xsty, /ysty, xr=[0.5, 4.5], $
            yr=[-0.25, 4], xthick=1.5, ythick=1.5, title=tname, $
            xtitle=xname, ytitle=yname
      endif
      xname = 'z quality'
      yname = 'signal-to-noise per pixel (fwhm)'
      plot, zz.zquality, s2nB, psym=4, /xsty, /ysty, xr=[0.5, 4.5], $
        yr=[-0.25, 5], xthick=1.5, ythick=1.5, title=tname, $
        xtitle=xname, ytitle=yname
      foo = where(zz.zquality eq 1 and s2nB ge 2., cnt)
      if cnt gt 0 then begin
          print, 'Object found with s/n > 2 and Q=1...'
          print, zz[foo].objname
      endif else print, 'No objects with s/n > 2 and Q=1...'
; plot Q value versus position on chip. 
      xname = 'Deimos x position'
      yname = 'Deimos y position'
      plot, [-1], [-1], psym=3, /xsty, /ysty, xr=[-420, 400], $
        yr=[0, 200], xtitle=xname, ytitle=yname, title=tname
      for i=0,n_elements(zz)-1 do begin
          if zz[i].zquality eq 255 then $
            oplot, [xpos[i]], [ypos[i]], psym=2, color=getcolor('black')
          if zz[i].zquality eq 1 then $
            oplot, [xpos[i]], [ypos[i]], psym=4, color=getcolor('red')
          if zz[i].zquality eq 2 then $
            oplot, [xpos[i]], [ypos[i]], psym=5, color=getcolor('blue')
          if zz[i].zquality eq 3 then $
            oplot, [xpos[i]], [ypos[i]], psym=6, color=getcolor('green')
          if zz[i].zquality eq 4 then $
            oplot, [xpos[i]], [ypos[i]], psym=1, color=getcolor('green')
          oplot, [-410], [180], psym=2, color=getcolor('black')
          oplot, [-410], [170], psym=4, color=getcolor('red')
          oplot, [-410], [160], psym=5, color=getcolor('blue')
          oplot, [-410], [150], psym=6, color=getcolor('green')
          oplot, [-410], [140], psym=1, color=getcolor('green')
          xyouts, fltarr(5)-405, [178, 168, 158, 148, 138], $
            ['Stars', 'Junk', 'Q=2', 'Q=3', 'Q=4']
      endfor

; also make plots comparing s/n versus position in slit.
      slitpos = objinfo[dexB].nrows-objinfo[dexB].cat_objpos < $
        objinfo[dexB].cat_objpos
      nonz = where(objinfo[dexB].objpos gt 0., cnt)
      xname = 'Distance from slit edge (pixels)'
      yname = 'signal-to-noise (fwhm)'
      x1 = min(slitpos) - 5
      x2 = max(slitpos) + 5
      plot, slitpos, s2nB, psym=2, xtitle=xname, ytitle=yname, $
        title=pname, /ysty, yr=[-0.5, 10], /xsty, xr=[x1, x2]
      if cnt gt 0 then begin
          slitpos2 = objinfo[dexB].nrows-objinfo[dexB].objpos < $
        objinfo[dexB].objpos
          oplot, slitpos2[nonz], s2nB[nonz], psym=6, color=getcolor('red')
      endif
      if not(keyword_set(speedy)) then begin
          yname = 'signal-to-noise (1-d spectra)'
          plot, slitpos, s2n_emp, psym=2, xtitle=xname, ytitle=yname, $
            title=pname, /ysty, yr=[-2, 8], /xsty, xr=[x1, x2]
          if cnt gt 0 then $
            oplot, slitpos2[nonz], s2n_emp[nonz], psym=6, color=getcolor('red')
      endif
; plot z-completeness versus position in slit.
      yname = 'z quality'
      plot, slitpos, zz.zquality, psym=2, xtitle=xname, $
        ytitle=yname, title=tname, /ysty, yr=[0, 4.5]
      foo = where(zz.zquality gt 200, cnt)
      if cnt gt 0 then $
        oplot, slitpos[foo], fltarr(cnt) + 0.5, psym=2


  endif

; close the ps file.
finish:
  if not(keyword_set(noplot)) then begin
      device, /close
      set_plot, 'x'
  endif
; finally close the open file!
  fits_close, fcb



end




