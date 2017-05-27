;+
; NAME:
;   starexam
;
; PURPOSE:
;   to plot sky subtracted stellar spectrum -- testing for A-band
;
; CALLING SEQUENCE:
;   starexam, filein,[star=star, lambda=lambda,sky=sky ]
; 
; INPUTS:
;   filein  -- spSlit* type file of star
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   star -- keyword returns array of sky subtracted stellar spectra,
;           one for each exposure
;   lambda -- wavelength solution for spectrum
;   sky   -- sky spectrum returned
;   spec  -- structure containing single spectrum, template,lambda,
;            and shift,rmag,xcor peak
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;    this is intended to make it easier to examine stars, looking for
;    A-band features. calls abandtest
;
; REVISION HISTORY:
;   MD 29Jun02
;   md 7sep02 extensively hacked for quicklook, returns structure
;----------------------------------------------------------------------
pro starexam, filein,  rmag, star=star, lambda=lambda, sky=sky,  spec=spec


  result = 0
  status = 0
  i = 1
; read 
;  while (status GE 0) do begin
    data = mrdfits(filein, 2*i-1 , /silent, status=status)
;    if (size(data_in, /n_dimen) eq 0) then status = -2 else begin 
;concatentate frames
       nrows =  (size(data.flux, /dimen))[1]
;       sset_in =  mrdfits(file, 2*i, /silent) ;get bspline structure
;       sky_in = bspline_valu(data_in.lambda+  $
;             data_in.lambda0#(fltarr(nrows) +1.), sset_in)     
;       if (i eq 1) then data = data_in else $
;          data = [data, data_in]
;       if (i eq 1) then sky = sky_in else $
;          sky = [sky, sky_in]
;    endelse
;    i = i+1 
;  endwhile

;  nframes = (size(data.flux, /dimen))[2]
;  nframes = i-2
          nframes = 1
  star = fltarr(4096, nframes)
  lambda = star

;should probably account for pixel shifts before generating sky-subtraction
  for i=0, nframes-1 do begin
;    summed = total(data[i].flux[1000:3000, *], 1)
;    peak = max(summed, j)
    
;    star[*, i]= total(data[i].flux[*, j-3 > 0:j+3 < nrows-1], 2)
;    sky =  total(data[i].flux[*, j-8:j-4], 2) + $
;               total(data[i].flux[*, j+4:j+8], 2)
;    star[*, i] = star[*, i]- sky*7./10. ;scale to same number of rows
       flag_cr, data, newivar ;flag bad Cosmic Rays
       data.ivar = newivar 
       
       lightprofile = find_object(data, profivar=ivarout, npix=npixout, $
                        /CR, /BPM, /NOSUBTRACT, /USETILT)
       nel = n_elements(lightprofile)
       lightprofile = lightprofile-min(lightprofile(4:nel-5))

       nrows = (size(data.flux, /dimens))[1]

        peakinfo, lightprofile, pkcol, fwhm, pk_cent=pkcent, $
          profivar=ivarout, window=9, s2n=sig2noise
     pkcol = pkcol[0]
     fwhm = fwhm[0]
;       maxlight = max(lightprofile, pkcol)
;       pkcol = (pkcol < (nrows-6)) >  5 ;limit possible solutions
       
       fwhm = 5 ;force this to be small, since alignment box is only 4" long
       extract = extract_1d(data, pkcol-fwhm, pkcol+fwhm, $
            pkcol-2*fwhm > 4, pkcol+2.*fwhm < nrows-5 )
       star[*, i] = extract.spec
       lambda[*, i] = extract.lambda
       sky = extract.sky
 ;fix to get sky at same time!! TBD      
; calculate system throughput using magnitude information 
; Dave Koo calculates 1.0 detected electron/sec/Ang in DEIMOS at
; R=20.0 (assumes 25% efficiency of telescope and detector)
; AB magnitudes
; USE this number to estimate throughput!
; (effective Lambda=6400AA, but our Rmags have larger ll_eff)
; (use 1.0 e-/sec/Ang as effective interpolant at 6900 Ang )
;    
       star_smooth = smooth(star[*, i], 100) ;get smooth spectrum
       lesser = lambda[*, i] lt 6900. ;get last pixel at this lambda
       jll = where(lesser gt 0, nct) ;does range go below our effective ll_R? 
       fluxR =  (nct gt 0) ? star_smooth[max(jll)] : star_smooth[100]
       dispersion = (nct gt 0) ? lambda[max(jll), i]-lambda[max(jll)-1, i] : $
              lambda[1000, i] - lambda[999, i]  ;Ang/pixel

; get integration time TBD!!
       int_time = 1200.
       calibration = 1.0 ;e-/sec/Ang for 20th mag obj at 6700 Ang

       transparency =  fluxR/int_time/dispersion/calibration/$
          10^(-.4*(rmag-20.))   

;    print, 'frame, peak row: ', i, pkcol
   
  endfor

;  splot, lambda, star[*, 0], yr=[-100, 1500]
;  soplot, lambda, sky/10., color=red ;plot sky 
;  for i=1,  nframes-1 do soplot, lambda, star[*, i], lines=i

   filename = strsplit(filein,'.', /extract)
   for i=0, nframes-1 do shift = abandtest(star[*, i], lambda[*, i], sky, $
      filename[3], rmag,  spec=spec)
   spec.transparency = transparency

; compute effective efficiency, relative to measured throughput of DEIMOS



return
end









