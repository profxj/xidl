
; routine to tabulate all of the signal-to-noise (quality
; assurance) data for the obererved masks. quality_wrapper calls the
; routine quality_test recursively to generate the data.

pro quality_wrapper, maskno=maskno, _extra=extra

; find the directory in which the results are held.
  result_dir = getenv('D2_RESULTS')
  if result_dir eq '' then $
    message, '(quality_wrapper) Specify D2_RESULTS directory in ' + $
    '.idlenv file!'
  cd, result_dir
; find all directories (masks) for which the code successfully
; completed (made quality plots). recall that the directory structure
; is ~D2_RESULTS/xxxx/yyyymmdd/.
  dirlist = findfile('*/*/ps/quality_plots.ps', count=ndir)
  if keyword_set(maskno) then begin
      length = strlen(maskno)
      dex = where(strmid(dirlist,0,length) eq maskno, ndir)
      if ndir gt 0 then dirlist = dirlist[dex]
  endif
  if ndir eq 0 then $
    message, '(quality_wrapper) ERROR: no quality output found!'
  
  for i=0,ndir-1 do begin
      jump1:
; switch to the ith mask directory.
      dirend = strpos(dirlist[i], '/ps')
      this_dir = strmid(dirlist[i], 0, dirend[0])
      cd, this_dir
; grab the mask name from the directory string.
      subdir = strsplit(dirlist[i], '/', /extract)
      if n_elements(subdir) lt 2 then $
        message, '(quality_wrapper) ERROR: Improper directory structure!'
      this_mask = strmid(dirlist[i], 0, 4)
      this_date = strmid(subdir[1],0,9) ;date of form yyyymmmdd
      print, 'Doing QA analysis of mask ' + this_mask + '.....'
; check that the mask is a DEEP2 mask.
      binfile = findfile('*bintab*.fits', count=bincnt)
      if bincnt eq 0 then begin
          print, 'No bin tab file found!'
          cd, result_dir
          if i eq ndir-1 then goto, jump2
          i = i + 1
          goto, jump1
      endif
      fits_open, binfile[0], fcb
      dex = where(fcb.extname eq 'MaskDesign', maskcnt)
      if maskcnt eq 0 then begin
          print, 'No MaskDesign bin table found!'
          fits_close, fcb
          cd, result_dir
          if i eq ndir-1 then goto, jump2
          i = i + 1
          goto, jump1
      endif
      table = mrdfits(binfile[0], dex[0], /silent)
      project = strcompress(table.projname, /remove_all)
      if project ne 'DEEP2-1HS' then begin
          print, 'Skipping non-DEEP2 mask!'
          fits_close, fcb
          cd, result_dir
          if i eq ndir-1 then goto, jump2 
          i = i + 1
          goto, jump1
      endif
      fits_close, fcb
; get the seeing from the seeing data file (seeing values from the
; observing logs).
      calib_dir = getenv('CALIB_DATA')
      if calib_dir eq '' then $
        message, '(quality_wrapper) Specify CALIB_DATA directory in ' + $
        '.idlenv file!' $
      else seeingfile = filepath('seeing.dat', root_dir=calib_dir)
      readcol, seeingfile, mask, date, seeing, format='A,A,F', $
        comment='#', /silent
      mask = strcompress(mask, /rem)
      date = strcompress(date, /rem)
      dex = where(mask eq this_mask and date eq this_date, cnt)
      if cnt eq 0 then begin
          print, '(quality_wrapper) ERROR: mask ' + this_mask + $
            ' not found in seeing data file.'
          seeval = 99.
      endif else seeval = seeing[dex[0]]
; also get the seeing from the alignment stars...
      objfile = findfile('*obj_info*', count=nobjfile)
      if nobjfile eq 0 then begin
          print, 'No obj_info file found!'
          cd, result_dir
          if i eq ndir-1 then goto, jump2
          i = i + 1
          goto, jump1 
      endif
      objinfo = mrdfits(objfile[0], 1, objhdr, /silent)
      astar = where(objinfo.objtype eq 'A', acnt)
      if acnt eq 0 then begin
          print, '(quality_wrapper) No alignment stars found for mask ' + $
            this_mask + '!!!'
          astar_seeing = 0.
      endif else begin
          pixscl = 0.117371 ;arcsec/deimos pixel
          astar_seeing = median(objinfo[astar].fwhm) * pixscl
      endelse
; also get the seeing from comparison of the fwhm distribution to the
; pcat seeing value. take the seeing difference on the blue side since
; it most closely matches to the photometric R band.
      cfht_pixscl = 0.207 ;arcsec/cfht pixel
      psee = sxpar(objhdr, 'PCAT_SEE') * cfht_pixscl ;in arcseconds
      seediff = sxpar(objhdr, 'SEEDIFFB') ;in arcseconds
      if seediff lt 0 then pcat_seeing = sqrt(psee^2 - seediff^2) $
      else pcat_seeing = sqrt(psee^2 + seediff^2)

; run the quality_test routine without producing any plots.
      quality_test, result=result, /noplot, /silent, _extra=extra
; get the airmass.
      slitfiles = findfile('slit.*', count=nslit)
      if nslit eq 0 then $
        message, '(quality_wrapper) ERROR: no slit files found ' + $
        'for mask ' + this_mask + '!!' $
      else begin
          hdr = headfits(slitfiles[0], ext=1, /silent)
          airmass = sxpar(hdr, 'AIRMASS')
      endelse
; make an output structure.
      if i eq 0 then $
        str = {mask:this_mask, date:this_date, $
               sn1:float(result[0]), sn2:float(result[1]), $
               pcat_seeing:float(pcat_seeing), $
               astar_seeing:float(astar_seeing), $
               log_seeing:float(seeval), airmass:float(airmass)} $
      else begin
          tmp = {mask:this_mask, date:this_date, $
                 sn1:float(result[0]), sn2:float(result[1]), $
                 pcat_seeing:float(pcat_seeing), $
                 astar_seeing:float(astar_seeing), $
                 log_seeing:float(seeval), airmass:float(airmass)}
          str = [str,tmp]
      endelse
; go back to the result directory.
      cd, result_dir
  endfor  

  jump2:
; write the output structure to a file in the $D2_RESULTS directory.
  if keyword_set(maskno) then outfile = 'qa_' + maskno + '_summary.fits' $
  else outfile = 'qa_summary.fits'
  mwrfits, str, outfile, /create, /silent

; determine the number of masks for which quality results were found.
  nres = n_elements(str)

; remove the NaN's...set them equal to zeros.
  nonfinite = where(finite(str.sn1) eq 0, cnt)
  if cnt gt 0 then str[nonfinite].sn1 = 0.
  nonfinite = where(finite(str.sn2) eq 0, cnt)
  if cnt gt 0 then str[nonfinite].sn2 = 0.

; now plot up the results.
  set_plot, 'ps'
  device, file='qa_summary.ps', /landscape
  ptitle='Total Number of Masks = ' + string(nres, format='(I3.3)')
; plot histogram of s/n (fwhm)at R=23.5.
  plothist, str.sn1, bin=0.1, thick=2, xthick=1.5, $
    ythick=1.5, /xsty, xr=[0.3, 2.5], title=ptitle, $
    xtitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5'
; plot histogram of s/n (1-d spectra) at R=23.5.
  if n_elements(uniq(str.sn2)) gt 1 then sflag = 1 else sflag = 0
  if sflag then plothist, str.sn2, bin=0.1, thick=2, xthick=1.5, $
    ythick=1.5, /xsty, xr=[0.3, 2.5], title=ptitle, $
    xtitle='signal-to-noise from 1-d spectra (per 1d pixel) at R=23.5'
; plot s/n (fwhm) vs. s/n (1-d spectra).
  plot, str.sn2, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='s/n from 1-d spectra (per 1d pixel) at R=23.5', $
    ytitle='s/n from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.0, 1.5], yr=[0.5, 2.0], title=ptitle
; plot s/n (fwhm) vs. pcat_seeing.
  plot, str.pcat_seeing, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='pcat seeing (arcsec)', $
    ytitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.5, 1.6], yr=[0.5, 2.0], title=ptitle
; plot s/n (1-d spectra) vs. pcat_seeing.
  if sflag then plot, str.pcat_seeing, str.sn2, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='pcat seeing (arcsec)', $
    ytitle='signal-to-noise from 1-d spectra (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.5, 1.6], yr=[0.3, 1.0], title=ptitle
; plot s/n (fwhm) vs. astar_seeing.
  plot, str.astar_seeing, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='astar seeing (arcsec)', $
    ytitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.6, 2.0], yr=[0.5, 2.0], title=ptitle
; plot s/n (1-d spectra) vs. astar_seeing.
  plot, str.astar_seeing, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='astar seeing (arcsec)', $
    ytitle='signal-to-noise from1-d spectra (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.6, 2.0], yr=[0.5, 2.0], title=ptitle
; plot s/n vs. log_seeing.
  plot, str.log_seeing, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='log seeing (arcsec)', $
    ytitle='signal-to-noise (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.5, 1.5], yr=[0.5, 2.0], title=ptitle
; plot s/n vs. airmass
  plot, str.airmass, str.sn1, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='airmass', $
    ytitle='signal-to-noise (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.9, 2.0], yr=[0.5, 2.0], title=ptitle

; close ps file.
  device, /close
  set_plot, 'x'


stop
end






