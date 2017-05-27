;+
; NAME:
;    abandtest
;
; PURPOSE:
;    determines spectral shift of aband absorption relative to
;    template.  Intended for quick-look analysis
;
; CALLING SEQUENCE:
;    shift = abandtest(flux,lambda, sky,[spec=spec],plot=plot)
;
; INPUTS:
;    flux -- spectrum of object (1d)
;    lambda -- wavelength of each pixel (1d)
;    sky   -- sky spectrum (1d) with same lambda scale
;
; KEYWORD:
;    plot  -- set if output plot desired
;
; OUTPUTS:
;   shift -- (pixels) shift relative to stored template
;
; OPTIONAL OUTPUTS:
;    spec  -- structure giving spectrum, lambda, template
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   does not require the full A-band to be present in spectrum; will
;   use whatever is available.  Present scheme is not properly
;   treating the edges.  Apodizing is applied.  Need to check lambda
;   solution against sky spectrum in order to shift if necessary
;   
; REVISION HISTORY:
;   MD   1jul02
;   MD   9sep02 structure for output, optional plots, apodize spectra
;   MD   11sep02 use sky spectrum to help lambda solution
;----------------------------------------------------------------------


pro getwave,h,wave
;+
; NAME:
; 	GETWAVE
;
; PURPOSE:
;	Create a wavelength vector from a FITS header of a spectrum
;	The FITS header must contain a CRVAL1 keyword and a CRPIX1 or
;	CD1_1 keyword.
;
; CALLING SEQUENCE:
;	GETWAVE, hdr, wave
;
; INPUTS:
;	hdr - FITS header for a 1-d spectrum, string array
;
; OUTPUTS:
;	wave - Wavelength vector.     The number of elements in WAVE is
;		given by the value of NAXIS1 in the FITS header.   The
;		element is given by CRVAL1 and the wavelength interval is
;		taken from CRPIX1 or CD1_1
;
; PROCEDURES CALLED:
;	Procedure: ZPARCHECK
;	Function: SXPAR
; REVISION HISTORY:
;	Written, W. Landsman   (HSTX)            October, 1993
;-
On_Error,2

if N_params() LT 2 then begin
	print,'Syntax - GETWAVE, h, wave'
	return
endif

 zparcheck, 'GETWAVE', h, 1, 7, 1, 'FITS header'

 xsize = sxpar(h,'NAXIS1')
 if xsize LE 0 then message, $
	'ERROR - FITS header does not contain positive value of NAXIS1'
 w0 = sxpar(h,'CRVAL1')
 if !ERR EQ -1 then message,'FITS header does not contain CRVAL1 value'
 wdelt = sxpar(h,'CDELT1')
 if !ERR EQ -1 then begin
	wdelt = sxpar(h,'CD1_1')
        if !ERR EQ -1 then message, $
             'FITS header does not contain CDELT1 or CD1_1 value'
 endif
 
 wave = w0 + wdelt*findgen(xsize)

 return
 end

pro read_aband,aa
  aband_dir = getenv('IDLSPEC2D_DIR')+'/etc/'
  aband=readfits(aband_dir +'aband.fits',heada, /silent)
  getwave,heada,ll
  aa = {aband:aband, ll:ll}

;normal and trim input spectrum
  aband=float(aband)/((total(aband[1400:1499])+total(aband[2500:2599]))/200.)

  aa = {ll:ll[1500:2000], aband:aband[1500:2000] }
return
end




function abandtest, flux, lambda, sky, filename, rmag, plot=plot, spec=spec

;   if (max(lambda) lt 7590. or min(lambda) gt 7660.) then begin
; make it cover the ENTIRE A-band!
   if (max(lambda) lt 7630. AND min(lambda) gt 7590.) then begin
      print, 'spectral coverage does not include A-band region!'
      return, 0
   endif 

;get aband template and wavelength
   read_aband, aa

; original data from file is sampled at .23 Ang/pixel, while DEIMOS
; with 1200 l/mm grating will be sampled at .33 Ang/pixel, so rebin to 
; match
   npix = n_elements(aa.aband)
   inscale = (aa.ll[npix-1]-aa.ll[0])/(npix-1)

   llow = where(lambda lt aa.ll[0], nlow) ;set range to cover aband
   if nlow eq 0 then lmin = 0 else $
     lmin=max(llow)
   llmax = where(lambda gt aa.ll[npix-1], nmax)
   if (nmax eq 0) then lmax = n_elements(flux)-1 else $
       lmax = min(llmax)

   aflux = flux(lmin:lmax)
   alambda = lambda(lmin:lmax)
   scale = alambda[1]-alambda[0] ;Angs/pixel in this spectral region
   npix_data = n_elements(aflux)

;rebin template to match data spectral dispersion
   rr = fix(npix*inscale/scale)
   aaband=congrid(aa.aband,rr, cubic=-.5)
   ll = float(congrid(aa.ll, rr, cubic=-.5))
   aacoeff = svdfit(ll, aaband, 3, yfit=aband_fit) ;quadratic polynomial fit
   aaband = aaband - aband_fit ;remove poly trend
   aaband = aaband/sqrt(variance(aaband)) ;unit normalization

;apodize the ends!
   npix = n_elements(aaband) ;reset number of template pixels
   for i=0, npix/10-1 do begin 
       aaband[i] = aaband[i]*sin(i*!pi/2/(npix/10))
       aaband[npix-1-i] = aaband[npix-1-i]*sin(i*!pi/2/(npix/10))
   endfor


;next check to see if data covers full range of template.  If not, extend.
   if alambda[0] gt ll[0] then begin
      nadd = fix((alambda[0]-ll[0])/scale) +1
      aflux = [(fltarr(nadd)+aflux[0]), aflux] ;prepend with constant value
      alambda = [fltarr(nadd), alambda] ;extend lambda array
      npix_data = n_elements(aflux)
   endif
   if alambda[npix_data-1] lt ll[npix-1] then begin
      nadd = fix((ll[npix-1]-alambda[npix_data-1])/scale) +1
      aflux = [aflux, (fltarr(nadd)+aflux[npix_data-1])] ;append, constant.
      alambda = [alambda, fltarr(nadd) ] ;extend lambda array
      npix_data = n_elements(aflux)
   endif
   
   aflux = smooth(aflux, 3) ;smooth out H-F noise
   aflux = aflux/median(aflux) ;normalize to unity amplitude

   pixs = npix <  npix_data ;take lessor of the two array lengths

   ll    = ll[0:pixs-1]
   aflux = aflux[0:pixs-1]
   aaband = aaband[0:pixs-1] ;force arrays to have same length

   ffcoeff = svdfit(ll, aflux, 3, yfit=aflux_fit) ;get cubic poly trend
   aflux = aflux - aflux_fit ;remove trend in star flux
   aflux = aflux/sqrt(variance(aflux)) ;unit normalization

;apodize the ends!
   for i=0, pixs/10-1 do begin 
       aflux[i] = aflux[i]*sin(i*!pi/2/(pixs/10))
       aflux[pixs-1-i] = aflux[pixs-1-i]*sin(i*!pi/2/(pixs/10))
   endfor


   xnorm = sqrt(variance(aaband)*variance(aflux))  ;use for normalization

   ncor = 21 
   noff = fix(ncor/2)
   xcor=fltarr(ncor)
   xx=findgen(ncor)-noff

   for j=0,ncor-1 do xcor[j]= total(aaband*shift(aflux,j-noff))/pixs/xnorm
;ncor element cross-correlation
   xpeak = max(xcor, ipeak) ;peak of x-corr.
   xxs = xx[ipeak-2 > 0:ipeak+2 < ncor-1]
   xcors = xcor[ipeak-2 > 0:ipeak+2 < ncor-1] ;select subset around peak
   fit=poly_fit(xxs,xcors,2)
   dist=-.5*fit[1]/fit[2] ; get central value of fit
;   print, 'X-correlation: ', xcor
   print, 'pixel shift, peak X-corr: ', dist, xpeak

    if keyword_set(plot) then begin
;   ps_open, 'ps/Aband_'+filename
      fname = 'ps/Aband_'+filename +'.ps'
      set_plot,'PS'
       device,file=fname,/portrait, /color, xsize=6, ysize=6, /inch, $
         SET_FONT='Courier', /TT_FONT,SET_CHARACTER_SIZE=[140,180]

      plot, ll, aflux/sqrt(variance(aflux))
      oplot, ll, aaband, line=2
      xyouts, 150, 150, 'pixel shift, peak X-corr: '+strmid(strcompress(dist,$
          /remove_all), 0, 5)+'  '+strmid(strcompress(xpeak, /remove_all), $
            0, 5), /device, charsize = 1.3, charthick = 2
       xyouts, 450, 700, 'R mag of star: '+strcompress(rmag), /device, $
        charthick = 2, charsize = 1.3
;   ps_close
       device,/close_file ; close off the ps file
      set_plot,'X'
    endif
  
; structure returned to be used by calling routine
 spec = {lambda:ll, spec:aflux/sqrt(variance(aflux)), template:aaband,$
         shift:dist, xpeak:xpeak,  transparency:1.}

 return, dist
end








