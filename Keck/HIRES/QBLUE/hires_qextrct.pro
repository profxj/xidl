pro hires_qextrct, finfil, outfil, sigfil, APER=aper, SKYREG=skyreg, $
                 ARCFIL=arcfil, ARCCALIB=arccalib, LINLST=linlst, $
                 RA=ra, DEC=dec, DISPLAY=display, CHKTRC=chktrc

 if not keyword_set( LINLST ) then linlst = 'Arcs/hires_blue.lst'
 if not keyword_set(CRVAL1) then crval1 = alog10(2700.d)
 if not keyword_set(CDELT) then cdelt = 2.0267028d-6
 if not keyword_set(NPIX) then npix = 83000L
 if not keyword_set(SIGREJ) then sigrej = 2.5

 ;; Read data
 dat = xmrdfits(finfil, /silent)
 var = xmrdfits(finfil, 1, /silent)
 
 ;; Run x_apall
 spec = x_apall(transpose(dat), ERROR=error, /CRUDE, VAR=transpose(var), $
                APER=aper, STRCT=strct,  GAIN=1., RN=2.5, /NOOV, $
                skyreg=skyreg,  CHKTRC=CHKTRC, DISPLAY=display) 

 ;; Extract Arc
 adat = xmrdfits(arcfil,/silen)
 arc = x_apall(adat, CLINE=strct.cline, APER=strct.aper, $
                TRACE=*strct.trace, /NOOV, /NOSKY, /ROT)


 ;; Calib arc
 if not keyword_set( ARCCALIB ) then stop
 restore, arccalib

 ;; Cross-correlate
 step = lindgen(100L) - 50L 
 corr = c_correlate((arc<10000.), $
                    (arcspec<10000.), step, /double)
 mx = max(corr, imx)
 imx = step[imx]
 shft = -imx 
 print, ', Shift', shft, format='(a,i4)'

 ;; Read in line list
 x_arclist, linlst, lines
 x_templarc, arc, lines, calib, $
   SHFT=shft, ALL_PK=all_pk, PKWDTH=3L, FORDR=9, LOG=FLOG, $
   PKSIG=35.

 ;; Fit arc
 gdfit = where(lines.flg_plt EQ 1, ngd)
 fin_fit = x_setfitstrct(NITER=3, HSIG=sigrej, LSIG=sigrej, $
                         NORD=3, /flgrej, MINPT=5)

 wfit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, $
                FITSTR=fin_fit, REJPT=rejpt)
 wav = x_calcfit(dindgen(n_elements(spec)), fitstr=fin_fit)

 dl = wav - shift(wav,1)
 print, 'WAVE FIT: RMS = ', fin_fit.rms/median(dl), ' pix'

 ;; Vacuum
 airtovac, wav

 ;; Heliocentric
 head = xheadfits(finfil)
 mjd = sxpar(head,'MJD')
 x_radec, ra, dec, radeg, decdeg
 helio = x_keckhelio(radeg, decdeg, 2000., jd=mjd)
 hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
 wav = wav * hel_corr

 ;; Rebin to global solution
 tot_wave = 10^(CRVAL1 + dindgen(npix)*cdelt)

 x_specrebin, wav, spec, tot_wave, nwfx, VAR=error, NWVAR=nwvar

 ;; Divide by exposure time
 exptim = sxpar(head,'EXPTIME')
 nwfx = nwfx / exptim
 nwvar = nwvar / exptim^2

 ;; Output
 sxaddpar, head, 'NAXIS', 1
 sxaddpar, head, 'NAXIS1', npix
 sxdelpar, head, 'NAXIS2'
 sxaddpar, head, 'CRVAL1', crval1
 sxaddpar, head, 'CDELT1', cdelt
 sxaddpar, head, 'CRPIX1', 1
 sxaddpar, head, 'CTYPE1', 'LINEAR'
 sxaddpar, head, 'DC-FLAG', 1
 sxaddpar, head, 'BITPIX', -32
 sxdelpar, head, 'BZERO'
 mwrfits, nwfx, outfil, head, /create
 mwrfits, sqrt(nwvar), sigfil, head, /create

 return
end
