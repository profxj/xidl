;+
;
; NAME
;     quickslit.pro
;
; PURPOSE
;     quickslit provides the muscle behind the DEIMOS QuickSlit
;     reduction GUI. the routine takes the science frame(s) specified
;     by the user and reduces a single user-defined slitlet. the
;     routine launches "atv" and "splot" widgets to display the
;     reduced data.
;
;
; SYNTAX
;     quickslit, file, slitno, errval=errval, mess=mess
;
; INPUTS
;     file = a string specifying the DEIMOS science (spectral) image
;            file on which to run the QuickSlit package. file can be
;            an array of strings, in which case the spectra from the
;            multiple images will be co-added. the files are assumed
;            to be uncompressed, but the user may use any of the
;            following syntax "44", "0044", "d0522_0044", or
;            "d522_0044.fits" to specify file "d0522_0044.fits". if
;            the user selected a non-standard root filename (e.g. "not
;            d0522_"), then the full image name must be specified
;            (e.g. "d0522_0044" or d0522_0044.fits").
;     slitno = an integer specifying the number of the slit which is
;              to be reduced. this slit number is the number assigned
;              to the slit in the mask design information.
;
; KEYWORDS
;     errval = a variable passed as the errval will be set equal
;              either the integer value 0 when quickslit runs without
;              error and set to unity when quickslit catches an error.
;     mess = a variable passed as the mess will be set equal to a
;            string specifying the exit status of the routine. in
;            cases where an error is caught by the routine, the string
;            gives the error message. when the routine runs
;            successfully, the string specifies the path for the
;            reduced data files (see OUTPUTS).
;
; OUTPUTS
;     see the errval and mess keywords.
;
;     also, the routine generates the reduction for the slit
;     specified. output files are written to a directory on a
;     "scratch" disk. the location of the files are printed to stdout
;     upon completion of the reduction and returned via the "mess"
;     optional parameter. the reduced 2-d spectrum is automatically
;     display using "atv.pro" and the 1-d spectrum is plotted using
;     "splot.pro". 
;
; PROCEDURES CALLED 
;     much of the UC-Berkeley DEEP2 DEIMOS data pipeline code and
;          the IDLUTILS package.
;     atv.pro (in IDLUTILS)
;     splot.pro (in IDLUTILS)
;
; EXAMPLES
;     quickslit, 'd0522_0014', 22, errval=isdone, mess=outstr
;
; COMMENTS
;     the input DEIMOS image files are assumed to reside in the
;     current DEIMOS output directory, as specified by the "outdir" 
;     command (/home/kics/instr/bin/outdir). quickslit checks the
;     outdir and looks for the user-defined image file within the
;     directory specified by the "outdir" command. for this reason,
;     quickslit should not be initiated while a MIRA focus is being
;     run. the MIRA script temporarily points the outdata directory to
;     a separate DEIMOS engineering directory.
;
; HISTORY
;     Created May 19, 2004 by mcc.
;
;-



PRO zspec_atv_overplot_lines, lambdabottom, lambdatop, linewaves, linenames, nbins, height, color
	; iteratively go through the above values
	; and plop them down on the ATV

	right = 1
  	FOR i=0,n_elements(linenames)-1 DO BEGIN
		l = linewaves[i]

      		; convert wavelength into an index:
      		nl0 = where(lambdabottom ge l-.5 and lambdabottom le l+.5, count)
      		if count lt 1 then continue
      		nl0 = nl0[(n_elements(nl0)-1)/2]
      		nlbottom = nl0[0]

      		nl0 = where(lambdatop ge l-.5 and lambdatop le l+.5, count)
      		if count lt 1 then continue
      		nl0 = nl0[(n_elements(nl0)-1)/2]
      		nltop = nl0[0]

      		; figure out the x & y positions
      		x = fix((float(nlbottom)/nbins-nlbottom/nbins)*nbins)
      		y = height*(7-nlbottom/nbins)

		if right then $
      			atvxyouts, x+10, y, linenames[i], color=color, charsize=2 $
		else $
      			atvxyouts, x-10-10*strlen(linenames[i]), y, linenames[i], color=color, charsize=2
		right = not right

      		atvplot, fltarr(7)+x, findgen(7)+y, color=color, thick=2
      		x = fix((float(nltop)/nbins-nl0/nbins)*nbins)
      		y = height*(7-nltop/nbins)
      		atvplot, fltarr(7)+x, height-4-findgen(7)+y, color=color,thick=2

      		if linenames[i] eq 'OII' then i = i+2
		if linenames[i] eq 'CA' then i = i+2
  	ENDFOR
END

function get_froot, outdir
  dirs = strsplit(outdir, '/', /extract)
  date = strcompress(dirs[n_elements(dirs)-1], /rem)
  mon = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', $
         'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
  mn = string((where(strmid(date, 4, 3) eq mon))[0] + 1, format='(I2.2)')
  dy = strmid(date, 7, 2)
  froot =  'd' + mn + dy + '_'
  return,  froot
end

pro quickslit, file, slitno, errval=errval, mess=mess
; define the error value and message.
  errval = 1
  mess =  'ERROR!'

; determine the rawdata directory.
  spawn, 'outdir', rawdatadir, sperr
; determine the standard root filename.
  froot = get_froot(rawdatadir)

; compress the file string.
  file = strcompress(file,  /rem)
  fs =  file[sort(file)]
  file =  fs[uniq(fs)]  
; determine the number of science frames.
  sparts = strsplit(file,  ',',  /extract)
  nframes = n_elements(sparts)
  filename = strarr(nframes)
  for ii=0,nframes-1 do begin
; check if the filename has the root specified.
     if strpos(sparts[ii], froot) eq -1 then $
       sparts[ii] = froot + string(sparts[ii], format='(I4.4)')

; check if the filename has the .fits on the end.
     flen = strlen(sparts[ii])
     if strmid(sparts[ii], flen-5, flen) ne '.fits' then $
       file_ii = strcompress(sparts[ii],  /rem) + '.fits' $
     else file_ii = sparts[ii]
; check that the raw file exists.
     frame = concat_dir(rawdatadir, file_ii)
     filename[ii] = findfile(frame, count=nfile)
     if nfile eq 0 then $
       message, 'Could not find file: ' + frame, /info
  endfor
; remove all unfound files.
  fnd = where(filename ne '', nframes)
  if nframes eq 0 then begin
     mess = 'No frames found!'
     return
  endif
  filename = filename[fnd]

; set the $DEIMOS_DATA environment variable.
  setenv, 'DEIMOS_DATA=/'

; check the clock.
  time0 = systime(/seconds)

; read the header.
  hdr = headfits(filename[0], ext=0, /silent)
  mask = sxpar(hdr, 'SLMSKNAM')
  lamps = sxpar(hdr, 'LAMPS')
  exptime = sxpar(hdr, 'EXPTIME')
  obstype = sxpar(hdr, 'OBSTYPE')
  grating = sxpar(hdr, 'GRATENAM')
  object = sxpar(hdr, 'OBJECT')
  gratepos = sxpar(hdr, 'GRATEPOS')
  hatchpos = sxpar(hdr, 'HATCHPOS')
  if gratepos eq 3 then wave = sxpar(hdr, 'G3TLTWAV') $
  else wave = sxpar(hdr, 'G4TLTWAV') 
  if gratepos eq 3 then tval = sxpar(hdr, 'G3TLTVAL') $
  else tval = sxpar(hdr, 'G4TLTVAL') 

; check that each frame is a science frame and not a direct image,
; arc, etc.
  is_science = strpos(lamps, 'Off') ge 0 and $
    exptime gt 0.0 and strpos(obstype, 'Object') ge 0 and $
    tval gt -18.3 and strpos(hatchpos, 'open') ge 0
  if is_science eq 0 then begin
     mess = 'This is not a science frame: ' + filename[0]
     return
  endif
; make a directory in which the output files will be written.
  spawn, 'best_scratch_disk -h',  sdisk,  sperr
  user = strcompress(getenv('USER'),  /rem)
  if user eq '' then user =  'deimosXX'
  udir =  concat_dir(sdisk, user)
  spawn,  'mkdir -p ' + udir
  outdir = concat_dir(udir,  'QSLIT_' + strcompress(mask, /rem))
  spawn, 'mkdir -p ' + outdir, spres, sperr
  cd, outdir

; create the planfile.
  planfile = strcompress(mask, /rem) + '.plan'
  message, 'Writing planfile: ' + planfile,  /info
  openw, wlun, planfile, /get_lun
  
  printf, wlun, '# Plan file auto-generated by quickslit.pro ', $
    systime()
  printf, wlun, '# Grating: ', grating, '       Wavelength: ', wave
  printf, wlun, mask, format='("MASK: ",A)'
  printf, wlun, rawdatadir, format='("RAWDATADIR: ",A)'
  printf, wlun, 'polyflag   - use polyflag for fitting lambda'
; find the corresponding arcs and flats and add them to the planfile.
  files = findfile(concat_dir(rawdatadir,'*.fits'), count=nfiles)
  for ii=0,nfiles-1 do begin
      h = headfits(files[ii], ext=0, /silent)
      mask_ii = sxpar(h, 'SLMSKNAM')
      gratepos_ii = sxpar(h, 'GRATEPOS')
      if gratepos_ii eq 3 then wave_ii = sxpar(h, 'G3TLTWAV') $
      else wave_ii = sxpar(h, 'G4TLTWAV') 
; if the mask, grating, and grating tilt agree with input frame, then
; check to see if the file is an arc or flat and include in the
; planfile.
      if (strcompress(mask_ii, /rem) eq strcompress(mask, /rem)) and $
        (gratepos_ii eq gratepos) and (wave_ii eq wave) then begin
          lamps_ii = sxpar(h, 'LAMPS')
          hatchpos_ii = sxpar(h, 'HATCHPOS')
          obstype_ii = sxpar(h, 'OBSTYPE')
          is_flat = strpos(lamps_ii, 'Qz') ge 0
          is_arc = strpos(lamps_ii, 'Ar') ge 0 or $
            strpos(lamps_ii, 'Kr') ge 0 or strpos(lamps_ii, 'Ne') ge 0 or $
            strpos(lamps_ii, 'Xe') ge 0 or strpos(lamps_ii, 'Cd') ge 0 or $
            strpos(lamps_ii, 'Hg') ge 0
          fparts = strsplit(files[ii], '/', /extract)
          fname = fparts[n_elements(fparts)-1]
          if is_arc then $
            printf, wlun, fname, format='("ARCNAME: ",A)'
          if is_flat then $
            printf, wlun, fname, format='("FLATNAME: ",A)'
      endif
  endfor
; write the science image file.
  for ii=0,nframes-1 do begin
     if ii gt 0 then begin
        h = headfits(files[ii], ext=0, /silent)
        mask_ii = sxpar(h, 'SLMSKNAM')
        gratepos_ii = sxpar(h, 'GRATEPOS')
        if gratepos_ii eq 3 then wave_ii = sxpar(h, 'G3TLTWAV') $
        else wave_ii = sxpar(h, 'G4TLTWAV') 
        if (strcompress(mask_ii, /rem) eq strcompress(mask, /rem)) and $
          (gratepos_ii eq gratepos) and (wave_ii eq wave) then begin
           fparts = strsplit(filename[ii], '/', /extract)
           fname = fparts[n_elements(fparts)-1]
           printf, wlun, fname, format='("SCIENCENAME: ",A)'
        endif else begin
           message, 'Skipping bad frame: ' + fname,  /info
           if n_elements(badframes) eq 0 then badframes = fname $
           else badframes = [badframes, fname]
        endelse
     endif
     fparts = strsplit(filename[ii], '/', /extract)
     fname = fparts[n_elements(fparts)-1]
     printf, wlun, fname, format='("SCIENCENAME: ",A)'
  endfor
; close the planfile.
  close, wlun
  free_lun, wlun
  nbad = n_elements(badframes)
  nframes = nframes - nbad


; get the object names for the object(s) in slit number slitno.
  fits_open, filename[0], fcb
  extnames = strcompress(fcb.extname, /rem)
  fits_close, fcb
  desidex = where(extnames eq 'DesiSlits', desicnt)
  if desicnt eq 0 then begin
     mess = 'No DesiSlits table found!'
     return
  endif
  objdex = where(extnames eq 'ObjectCat', objcnt)
  if objcnt eq 0 then begin
     mess = 'No ObjectCat table found!'
     return
  endif
  mapdex = where(extnames eq 'SlitObjMap', mapcnt)
  if mapcnt eq 0 then begin
     mess = 'No SlitObjMap table found!'
     return
  endif
  desi = mrdfits(filename[0], desidex[0], /silent)
  objcat = mrdfits(filename[0], objdex[0], /silent)
  smap = mrdfits(filename[0], mapdex[0], /silent)
  dex0 = where(long(desi.slitname) eq slitno, cnt0)
  if cnt0 eq 0 then begin
     mess = 'Slit not found in DesiSlit table!'
     return
  endif
  dex1 = where(desi[dex0[0]].dslitid eq smap.dslitid, cnt1)
  if cnt1 eq 0 then begin
     mess = 'Slit not found SlitObjMap table!'
     return
  endif
  dex2 = where(smap[dex1[0]].objectid eq objcat.objectid, nobj)
  if nobj eq 0 then begin
     mess = 'Slit not found ObjectCat table!'
     return
  endif
  objnames = strcompress(objcat[dex2].object,  /rem)
  objstr = ''
  for ii=0, nobj-1 do objstr = objstr + objnames + ' '


; make the bintab file.
  make_bintab_file, planfile, /quick

; find the chips on which the slit is contained.
  ss6 = mrdfits('*bintab*', 6, /silent)
  ss3 = mrdfits('*bintab*', 3, /silent)
  dex3 = where(long(ss3.slitname) eq slitno, cnt3)
  if cnt3 eq 0 then begin
     mess = 'slit not on mask!'
     return
  endif
  dex6 = where(ss3[dex3].dslitid eq ss6.dslitid, cnt6)
  xmm = (ss6[dex6].slitx1 + ss6[dex6].slitx2 + $
         ss6[dex6].slitx3 + ss6[dex6].slitx4) / 4.0
  scale_asec_per_pix = 0.117371
  scale_mm_per_asec = 0.73
  f1 = 8192 * scale_asec_per_pix * scale_mm_per_asec / 4.0
  gap1 = -1.0 * f1
  gap2 = 0
  gap3 = f1
  if xmm le gap1 then chiplist = [1,5]
  if xmm gt gap1 and xmm le gap2 then chiplist = [2,6]
  if xmm gt gap2 and xmm le gap3 then chiplist = [3,7]
  if xmm gt gap3 then chiplist = [4,8]

; run the pipeline on the slit.
  deimos_mask_calibrate, /noplot, /quick, slitlist=slitno, chiplist=chiplist

  deimos_2dreduce, file=fname, /quick

  qlist = findfile('quickSlit.*' + $
                   string(slitno, format='(I3.3)') + '*.fits', count=nfiles)
  if nfiles eq 0 then begin
     mess = 'No quickSlit files found!'
     return
  endif
  splist = qlist
  for ii=0,nfiles-1 do begin
      spos = strpos(qlist[ii],'Slit')
      splist[ii] = 'sp' + strmid(qlist[ii],spos)
      spawn, 'mv -f ' + qlist[ii] + ' ' + splist[ii]
  endfor
  spslit_combine, splist, nlsky=0
  
  slitfiles = findfile('slit.*' + string(slitno, format='(I3.3)') + $
                       '*.fits',  count=nslit)
  if nslit eq 0 then begin
     mess = 'No slitfiles found!'
     return
  endif
  do_extract, files=slitfiles, nsigma_optimal=1.75, nsigma_boxcar=1.1 

; now display the data.
  spec1d = findfile('spec1d.*.' + string(slitno, format='(I3.3)') + $
                    '.*.fits', count=nfile)
  if nfile gt 0 then begin
      ss1d = fill_gap(spec1d[0], /horne)
      splot, ss1d.lambda, ivarsmooth(ss1d.spec, ss1d.ivar, 11), $
        xtitle='!6wavelength'
      medval =  median(ss1d.spec)
      skyspec0 = 25.0
      skyoff =  -500.0
      skyspec = ss1d.skyspec / median(ss1d.skyspec) * skyspec0 + skyoff
      soplot,  ss1d.lambda, skyspec, color=3
  endif


  blue = mrdfits(slitfiles[0], 1, /silent)
  wv1 = lambda_eval(blue)
  if nslit gt 1 then red = mrdfits(slitfiles[1], 1, /silent)
  wv2 = lambda_eval(red)

  nbins = (size(wv1, /dimensions))[0]
  nbins = nbins/4
  height = n_elements(blue.flux[0,*])+2
  blambda0 = wv1[*,0]
  rlambda0 = wv2[*,0]

  flux1 = blue.flux[0:nbins-1,*]
  flux2 = blue.flux[nbins:2*nbins-1,*]
  flux3 = blue.flux[2*nbins:3*nbins-1,*]
  flux4 = blue.flux[3*nbins:4*nbins-1,*]
  lam1  = blambda0[0:nbins-1]
  lam2  = blambda0[nbins:2*nbins-1]
  lam3  = blambda0[2*nbins:3*nbins-1]
  lam4  = blambda0[3*nbins:4*nbins-1]
  flux_blue = [[flux4],[lam4],[lam4],[flux3],[lam3],[lam3],[flux2],[lam2],[lam2],[flux1],[lam1],[lam1]]

  flux1 = red.flux[0:nbins-1,*]
  flux2 = red.flux[nbins:2*nbins-1,*]
  flux3 = red.flux[2*nbins:3*nbins-1,*]
  flux4 = red.flux[3*nbins:4*nbins-1,*]
  lam1  = rlambda0[0:nbins-1]
  lam2  = rlambda0[nbins:2*nbins-1]
  lam3  = rlambda0[2*nbins:3*nbins-1]
  lam4  = rlambda0[3*nbins:4*nbins-1]
  flux_red = [[flux4],[lam4],[lam4],[flux3],[lam3],[lam3],[flux2],[lam2],[lam2],[flux1],[lam1],[lam1]]


  atv, [[flux_red],[flux_blue]], min=-5, max=25

  hlambda = (size(wv1))[2] < (size(wv2))[2]
  lambdabottom = [blambda0, rlambda0]
  lambdatop = [wv1[*,hlambda-1], wv2[*,hlambda-1]]
  offset = (height - (size(wv1, /dimensions))[1])/2>0  
  telluricnames = ['B-Band', 'A-Band']
  telluricwaves = [6860, 7600]
  zspec_atv_overplot_lines, lambdabottom, lambdatop, $
    telluricwaves, telluricnames, nbins, height, 5

  extraction_offset = ss1d.r1
  extraction_r1 =  ss1d.r1
  extraction_r2 =  ss1d.r2

  for i=0, 7 do begin
     base = fltarr(10) + height*i
     atvplot, findgen(10), base + extraction_r1, color=1, thick=2
     atvplot, findgen(10), base + extraction_r2, color=1, thick=2
     atvplot, 1024-findgen(10), base + extraction_r1, color=1, thick=2
     atvplot, 1024-findgen(10), base + extraction_r2, color=1, thick=2
  endfor

  help,  ss1d,  /str

  time1 = systime(/seconds)
  time = (time1 - time0) / 60.0
  print
  print
  print, '-------------------------------------------------------'
  print, 'QUICKSLIT Reduction Complete'
  print, '  Michael Cooper (cooper@astron.berkeley.edu)'
  print, 'Elapsed time: ' + strcompress(time, /rem) + ' minutes.'
  print, 'Outdata Directory = ' + outdir
  print, 'Mask: ' + strcompress(mask,  /rem)
  print, 'Slit: ' + string(slitno,  format='(I3.3)')
  print, 'Object(s): ' + objstr
  print, strcompress(nframes) + ' science frames reduced;'
  if nbad gt 0 then begin
     print, 'Following frames were ignored: '
     for ii=0,nbad-1 do print, badframes[ii]
  endif
  print, '-------------------------------------------------------'
  print


  errval = 0
  mess = outdir

end
