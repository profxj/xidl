pro quick_sciqa, quicklevel=quicklevel, file=file, silent=silent

bad = 0
s2nlimit = 0.4

      fileexists = file_test('../science_qa.dat')
      if fileexists then readcol, '../science_qa.dat', filename, masks, align, pm1, alignsig, seeing, pm2, seeingsigma, s2narr, pm3, s2nsigma, $
        skipline=1, format='A,A,F,A,F,F,A,F,F,A,F', /SILENT

      if fileexists then bestseeing = min(seeing) else bestseeing = 999.

if n_elements(quicklevel) eq 0 then quicklevel = 2


if n_elements(silent) eq 0 then silent = 0


qlfiles = findfile('quick*.fits',  count=count)
fits_info, qlfiles[0], /SILENT, n_ext=next
if quicklevel lt 2 then bspline = 1 else bspline = 0

;hdu = next/(1+bspline)
hdu = next - bspline

; Check for buckling
chiparr=fltarr(count)
dxarr=fltarr(count)

for i=0,count-1 do begin &$
  hdr=headfits(qlfiles[i],ext=hdu) &$
  chiparr[i]=sxpar(hdr,'CHIPNO') &$
  dxarr[i]=sxpar(hdr,'DX') &$
endfor

srtchips=sort(chiparr)
uchips=chiparr[uniq(chiparr,srtchips)]
udx=dxarr[uniq(chiparr,srtchips)]
nchips=n_elements(uchips)
nspec=intarr(nchips)

for i=0,nchips-1 do nspec[i]=total(chiparr eq uchips[i])

meddxs = string(median(dxarr,/even), format='(f6.3)')
sigdxs = string(djsig(dxarr), format='(f6.3)')

  bfile = '*bintab*.fits'
  bfile = bfile[0] 

simple_tables, bfile, slitnames=slitnames, mags=magall, slitwid=slitwid, $
  slitlen=slitlen, xmm=xmm, ymm=ymm, objnames=objnames


; find and process alignment stars
whalign = where(slitwid gt ( 1.2*median(slitwid) > 3.), nalign)


for i=0, nalign-1 do begin
   wh = where(strpos(qlfiles, $
          strcompress('.'+slitnames[whalign[i]], /REMOVE)) ge 0, whct)
          if whct gt 0 then if n_elements(alignfiles) eq 0 then $
            alignfiles = qlfiles[wh] $
          else alignfiles = [alignfiles, qlfiles[wh]]
endfor

issci = fltarr(n_elements(qlfiles))

for i=0, n_elements(qlfiles) - 1 do begin
   issci[i] = total(alignfiles eq qlfiles[i]) eq 0.
endfor

scifiles = qlfiles[where(issci)]

nsci = n_elements(scifiles)
nalign = n_elements(alignfiles)
centers = fltarr(nalign)
fwhms = centers
blue = strpos(alignfiles, 'B.fits') ge 0
peaks = centers

for i=0, nalign-1 do begin
; could do a flux vs. mag. here
   slitcenter = 1
   quickstruct = MRDFITS(alignfiles[i], hdu, /SILENT)
   profile = find_object(quickstruct,prof=pivar,align=slitcenter, /quick)

   peakinfo,profile,pkcol,fwhm,prof=pivar,pk_q=pk_q,pk_c=pk_c
   centers[i] = slitcenter
   peaks[i] = pk_q[0]
   fwhms[i] = fwhm[0]

endfor

whb = where(blue gt 0)
whr = where(blue eq 0)

peaks = peaks*0.11
centers = centers*0.11
fwhms = fwhms*0.11

djs_iterstat, peaks[whb]-centers[whb], median=medshiftb, sigma=sigshiftb
djs_iterstat, peaks[whr]-centers[whr], median=medshiftr, sigma=sigshiftr
djs_iterstat, fwhms[whb], median=medfwhmb, sigma=sigfwhmb
djs_iterstat, fwhms[whr], median=medfwhmr, sigma=sigfwhmr
djs_iterstat, peaks-centers, median=medshift, sigma=sigshift
djs_iterstat, fwhms, median=medfwhm, sigma=sigfwhm

medshifts = string(medshift, format='(f6.3)')
sigshifts = string(sigshift/sqrt(2./3.*(nalign-1)), format='(f6.3)')
medfwhms = string(medfwhm, format='(f6.3)')
sigfwhms = string(sigfwhm/sqrt(2./3.*(nalign-1)), format='(f6.3)')

print, 'spatial alignment offset: ', medshifts, ' +/- ', sigshifts, ' arcsec'
print, 'seeing fwhm             : ', medfwhms, ' +/- ', sigfwhms, ' arcsec'
print, 'Median flexure/buckle offset from flats: ', meddxs, '   RMS: ', sigdxs, ' pixels'
bcorr = sqrt(2./3.*(n_elements(whb)-1))
rcorr = sqrt(2./3.*(n_elements(whr)-1))

nstar = n_elements(whb) < n_elements(whr)

   cd,current=cwd
     stringsep = strsplit(strcompress(cwd, /REMOVE), '/', /extract)
     mask = stringsep[n_elements(stringsep)-1]


; if mask might be buckled, alert!
   if median(dxarr,/even) gt 2.5 OR min(dxarr) gt 1.5 OR djsig(dxarr) gt 1 then begin
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/doh.au &'
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/doh.au &'
  
            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2
            printf, 2, 'Warning: Mask may be off FCS position or buckled!'
            printf, 2, 'Mask: '+mask
            if n_elements(file) ne 0 then printf, 2, 'File: ', file
            printf, 2
            printf,2, 'Median offset from flats: ', meddxs, '   RMS: ', sigdxs, ' pixels'
            printf,2
            printf,2,'Deviations of 1-2 pix might be off position; 5+ pix or '
            printf,2,'strongly variable chip-to-chip may indicate buckling.'
            printf,2
            printf, 2, 'Chips:     dx_B      dx_R     # of spectra'
            for i=0, nchips/2-1 > 0 do printf, 2, uchips[i],' /',$
              uchips[i+nchips/2], udx[i], $
              udx[i+nchips/2],nspec[i], form='(I2,A,I2,2f10.3,I10)'
            printf, 2
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning &'
      if silent lt 2 then spawn, spawnstring
      bad = 1
      endif



; if seeing is bad, alert!
   if medfwhm gt 1.1 then begin
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/doh.au &'
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/doh.au &'
  
            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2
            printf, 2, 'Warning: Bad seeing!'
            printf, 2, 'Mask: '+mask
            if n_elements(file) ne 0 then printf, 2, 'File: ', file
            printf, 2
            printf,2, 'Median seeing fwhm: ', medfwhms, ' +/- ', sigfwhms, ' arcsec'
            printf, 2
            printf,2, 'Seeing fwhm (arcsec):'
            printf, 2
            printf, 2, 'Star:    B       R'
            for i=0, nstar-1 > 0 do printf, 2, slitnames[whalign[i]], fwhms[whb[i]], $
              fwhms[whr[i]], form='(A4,2f8.3)'
            printf, 2
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning &'
      if silent lt 2 then spawn, spawnstring
      bad = 1
      endif



; if seeing has degraded, alert!
   if medfwhm gt (bestseeing+0.3) then begin
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/doh2.au &'
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/doh2.au &'
 
            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2
            printf, 2, 'Warning: Seeing has degraded!'
            printf, 2, 'Mask: '+mask
            if n_elements(file) ne 0 then printf, 2, 'File: ', file
            printf, 2
            printf,2, 'Median seeing fwhm: ', medfwhms, ' +/- ', sigfwhms, ' arcsec'
            printf, 2
            printf,2, 'Seeing fwhm (arcsec):'
            printf, 2
            printf, 2, 'Star:    B       R'
            for i=0, nstar-1 > 0 do printf, 2, slitnames[whalign[i]], fwhms[whb[i]], $
              fwhms[whr[i]], form='(A4,2f8.3)'
            printf, 2
            printf, 2, 'History for the night:'
            printf, 2, 'File:', 'Mask', 'Seeing FWHM:', format='(3A20)'
            for i=0, n_elements(masks)-1 do printf, 2, filename[i], masks[i], '   ', seeing[i], ' +/- ',seeingsigma[i], format='(2A20,A,f6.3,A,f6.3)'
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning &'
      if silent lt 2 then spawn, spawnstring
      bad = 1
      endif

; if alignment is bad, alert!
   if medshift gt 0.25 then begin
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/STTNG-redalert.au &'
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/STTNG-redalert.au &'

            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2
            printf, 2, 'Warning: Bad alignment offset!'
            printf, 2, 'Mask: '+mask
            if n_elements(file) ne 0 then printf, 2, 'File: ', file
            printf, 2
            printf,2, 'Median alignment offset: ', medshifts, ' +/- ', sigshifts, ' arcsec'
            printf, 2
            printf,2, 'Alignment offset (arcsec):'
            printf, 2
            printf, 2, 'Star:    B       R'
            for i=0, nstar-1 > 0 do printf, 2, slitnames[whalign[i]], peaks[whb[i]]-centers[whb[i]], $
              peaks[whr[i]]-centers[whr[i]], form='(A4,2f8.3)'
            printf, 2
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning &'
      if silent lt 2 then spawn, spawnstring
      bad = 1
      endif


s2ns = fltarr(nsci)
mags = s2ns
badarr = s2ns
blue = s2ns

; process science objects
for i=0, nsci-1 do begin

   slitpos = strpos(scifiles[i], 'B.fits') > strpos(scifiles[i], 'R.fits')
   whichslit = where(strcompress(slitnames, /REMOVE) $
                     eq strmid(scifiles[i], slitpos-3, 3))

   mag = magall[whichslit]
   quickstruct = MRDFITS(scifiles[i], hdu, /SILENT)
   profile = find_object(quickstruct,prof=pivar,npix=npix, /quick)

   peakinfo,profile,pkcol,fwhm,prof=pivar,pk_q=pk_q,pk_c=pk_c, s2n_fwhm=s2n, npix=npix, s2n_window=s2n_window
;   s2ns[i] = s2n[0] ;> s2n_window[0]
   s2ns[i] = s2n[0] > s2n_window[0]/sqrt(median(npix))

   mags[i] = mag[0]
   badarr[i] = abs(pk_q[0]-pk_c[0]) gt 2 OR abs(pkcol[0]-pk_c[0]) gt 2 OR abs(pkcol[0]-pk_q[0]) gt 2 OR s2ns[i] lt 0
   blue[i] = strpos(scifiles[i], 'B.fits') gt 0
endfor

whgood = where(badarr eq 0)
s2ns = s2ns[whgood]
mags = mags[whgood]
blue = blue[whgood]


    whb = where(blue, comp=whr)

if silent lt 3 then begin

    window, 0, xs=512, ys=512
    loadct, 12, /SILENT
    plot, mags, s2ns, psym=4, xtit='R mag.',  ytit='S/N per pixel (extracted spectrum)',  yrange=[0, max(s2ns)]+0.1
    oplot, mags[whb], s2ns[whb], color=16*6, psym=4
    oplot,mags[whr], s2ns[whr], color=12*16, psym=4

endif

    dex1 = where(mags lt median(mags, /even), comp=dex2)
      x1 = median(mags[dex1], /even)
      y1 = median(s2ns[dex1], /even)
      x2 = median(mags[dex2], /even)
      y2 = median(s2ns[dex2], /even)
      mcc_polyfit, [x1, x2], [y1, y2], [0, 1], a=a
      xfit = [0, 100]
      yfit = a[0] + xfit*a[1]
     if silent lt 3 then  oplot, xfit, yfit

      noms2n = a[0]+23.25*a[1]
      s2nnoise = stdev(s2ns-a[0]-a[1]*mags)/3.

      print, 'Signal-to-noise at mag. 23.25: ',noms2n , ' +/- ', s2nnoise, ' Nominal: ', s2nlimit, format='(A,f7.3,A,f7.3,A,f7.3)'

print


; if signal-to-noise is bad, alert!
  if  (noms2n+s2nnoise) lt s2nlimit then begin
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/Homer_Scream.au &'
      if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/Homer_Scream.au &'
     cd,current=cwd
     stringsep = strsplit(strcompress(cwd, /REMOVE), '/', /extract)
     mask = stringsep[n_elements(stringsep)-1]

            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2
            printf, 2, 'Warning: Low signal-to-noise!'
            printf, 2, 'Mask: '+mask
            if n_elements(file) ne 0 then printf, 2, 'File: ', file
            printf, 2
            printf,2, 'Signal-to-noise at R=23.25: ', noms2n, ' +/- ', s2nnoise, ' arcsec', form='(A,f7.3,A,f7.3,A)'
            printf, 2
            printf, 2, 'It is advisable to take extra frames!'
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning &'
      if silent lt 2 then spawn, spawnstring
      bad = 1
      endif


      if n_elements(file) gt 0 then begin
          openw, 2, '../science_qa.dat', append=fileexists
          if fileexists eq 0 then printf, 2, 'Filename', 'Mask', $
            'Alignment offset', 'Seeing', 'Signal-to-noise', 'FCS offset',format='(6A20)'
          printf, 2, file, mask, medshifts+' +/-'+sigshifts, $
            medfwhms+' +/- '+sigfwhms,  '   ', noms2n, ' +/- ', $
            s2nnoise, meddxs+' +/- '+sigdxs,format='(4A20,A,f6.3,A5,f6.3,A20)'
          close, 2
       endif

      if bad eq 0 then if silent eq 0 then spawn, '/home/kics/instr/bin/play_sound -v 49 -host pohue -account deepteam /home/deepteam/sounds/Homer-Wohoo.au &'

return
end
