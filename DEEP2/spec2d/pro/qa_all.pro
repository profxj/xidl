
pro make_ascii_file, ss, file=file
; determine how many masks we have.
  nmasks = n_elements(ss)
; open our output file for writing.
  openw, unit, file, /get_lun

; print the median stats.
  printf, unit, 'Median sn_sprofB = ' + $
    strcompress(median(ss.sn_sprofB), /rem)
  printf, unit, 'Median sn_spec1d = ' + $
    strcompress(median(ss.sn_spec1d), /rem)
; print the header for the file.
  printf, unit, 'Mask     Date-Obs   Airmass   astar_seeing  ' + $
    'sn_sprofB   sn_spec1d  zcompleteness'
  printf, unit, '--------------------------------------------------' + $
    '----------------------------------------------'

; loop thru the masks and print the data to the file.
  sp = ' '
  frmt = '(A5, A3, A10, A3, F5.3, A3, F5.3, A3, ' + $
    'F5.3, A3, F5.3, A3, F6.2)'
  for ii=0,nmasks-1 do printf, unit, ss[ii].maskname, sp, ss[ii].date, $
    sp, ss[ii].airmass, sp, ss[ii].astar_seeing, sp, ss[ii].sn_sprofB, $
    sp, ss[ii].sn_spec1d, sp, ss[ii].zcom*100.0, format=frmt
; close the output file.
  close, unit
end

pro qa_all

; find the directory in which the spec2d results are located.
  d2dir = getenv('D2_RESULTS')
  if d2dir eq '' then $
    message, '(qa_all.pro) Specify D2_RESULTS directory in ' + $
    '.idlenv file!'
  cd, concat_dir(d2dir,'qa')
; collect the output files from the qa_check routine.
  files = findfile('qa.*.fits*', count=nfiles)
  if nfiles eq 0 then $
    message, '(qa_all.pro) ERROR: no quality output found!'


; read in the latest zcat.
  zfile = findfile('$DEEP2PRODUCTS/zcat*.fits*', count=zcnt)
  if zcnt eq 0 then message, 'No zcat found!!!'
  zfile = reverse(zfile[sort(zfile)])
  zz = mrdfits(zfile[0], 1, /silent)

; loop through the qa output files and tabulate the results.
  tmp = {maskname:'', date:'', airmass:0.0, $
         astar_seeing:0.0, pcat_seeing:0.0, $
         sn_sprofB:0.0, sn_spec1d:0.0, zcom:0.0}
  all = replicate(tmp, nfiles)
  for i=0,nfiles-1 do begin
      ss = mrdfits(files[i], 1, hdr, /silent)
      all[i].maskname = sxpar(hdr, 'maskname')
      all[i].date = sxpar(hdr, 'date-obs')
      all[i].airmass = sxpar(hdr, 'airmass')
      all[i].astar_seeing = sxpar(hdr, 'asee')
      all[i].pcat_seeing = sxpar(hdr, 'psee')
      all[i].sn_sprofB = sxpar(hdr, 'sn_val1')
      all[i].sn_spec1d = sxpar(hdr, 'sn_val2')
      all[i].zcom = sxpar(hdr, 'zge3')
      if all[i].zcom eq 0.0 then begin
          totwh = where(strcompress(zz.maskname, /rem) eq $
                        strcompress(all[i].maskname, /rem), totnum)
          if totnum gt 0 then begin
              q3wh = where(zz[totwh].zquality ge 3, q3num)
              all[i].zcom = float(q3num) / float(totnum)
          endif 
      endif
      if keyword_set(tab) then tab = [tab, ss] $
      else tab = ss
  endfor


; write the output structure to a file in the $D2_RESULTS directory.
  outfile = 'summary.qa.fits'
  mwrfits, all, outfile, /create, /silent

; make an ascii file containing the s/n values for each mask.
  make_ascii_file, all, file='summary.sn.dat'

; determine the number of masks for which quality results were found.
  nres = n_elements(all)

; remove the NaN's...set them equal to zeros.
;  nonfinite = where(finite(str.sn1) eq 0, cnt)
 ; if cnt gt 0 then str[nonfinite].sn1 = 0.
 ; nonfinite = where(finite(str.sn2) eq 0, cnt)
 ; if cnt gt 0 then str[nonfinite].sn2 = 0.

; now plot up the results.
  set_plot, 'ps'
  device, file='ps/summary.qa.ps', /landscape
  ptitle='Total Number of Masks = ' + string(nres, format='(I3.3)')
; plot histogram of s/n (fwhm)at R=23.5.
  plothist, all.sn_sprofB, bin=0.1, thick=2, xthick=1.5, $
    ythick=1.5, /xsty, xr=[0.3, 2.5], title=ptitle, $
    xtitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5'
; plot histogram of s/n (1-d spectra) at R=23.5.
  plothist, all.sn_spec1d, bin=0.1, thick=2, xthick=1.5, $
    ythick=1.5, /xsty, xr=[0.3, 2.5], title=ptitle, $
    xtitle='signal-to-noise from 1-d spectra (per 1d pixel) at R=23.5'
; plot s/n (fwhm) vs. s/n (1-d spectra).
  plot, all.sn_spec1d, all.sn_sprofB, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='s/n from 1-d spectra (per 1d pixel) at R=23.5', $
    ytitle='s/n from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.0, 1.5], yr=[0.5, 2.0], title=ptitle
; plot s/n (fwhm) vs. pcat_seeing.
  plot, all.pcat_seeing, all.sn_sprofB, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='pcat seeing (arcsec)', $
    ytitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.5, 1.6], yr=[0.5, 2.0], title=ptitle
; plot s/n (1-d spectra) vs. pcat_seeing.
  plot, all.pcat_seeing, all.sn_spec1d, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='pcat seeing (arcsec)', $
    ytitle='signal-to-noise from 1-d spectra (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.5, 1.6], yr=[0.3, 1.0], title=ptitle
; plot s/n (fwhm) vs. astar_seeing.
  plot, all.astar_seeing, all.sn_sprofB, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='astar seeing (arcsec)', $
    ytitle='signal-to-noise from spatial profile (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.6, 2.0], yr=[0.5, 2.0], title=ptitle
; plot s/n (1-d spectra) vs. astar_seeing.
  plot, all.astar_seeing, all.sn_spec1d, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='astar seeing (arcsec)', $
    ytitle='signal-to-noise from1-d spectra (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.6, 2.0], yr=[0.5, 2.0], title=ptitle
; plot s/n vs. airmass
  plot, all.airmass, all.sn_sprofB, psym=2, xthick=1.5, $
    ythick=1.5, xtitle='airmass', $
    ytitle='signal-to-noise (per 1d pixel) at R=23.5', $
    /xsty, /ysty, xr=[0.9, 2.0], yr=[0.5, 2.0], title=ptitle

; plot z quality versus distance from edge of slit for all masks.
  xmax = max(tab.edge_dist)
  plot, tab.edge_dist, tab.zquality, /xsty, /ysty, $
    xr=[0, xmax], yr=[0, 4], ytitle='z quality', $
    xtitle='distance from edge of slit (pixels)', $
    xthick=1.5, ythick=1.5, psym=2


; close ps file.
  device, /close
  set_plot, 'x'












end


