;+
; NAME:
;   long_sensfunc
;
; PURPOSE:
;   Use one or more standard star spectra to determine the
;   spectroscopic response function. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;   scifile         - file containing object structure which has spectrum 
;                     for standard star
;
;   standard_name   - name of the standard star
;   
;   sensfuncfile    - File to write the sensitivity function out to
;
;
; OPTIONAL INPUTS:
;   OBJID           - object id of standar star in the object structure. Default
;                     is to the first object. 
;   nresln          - Break point spacing in resolution elemnets
;                     (default=20)
;   /MSK_BALM       - Mask Balmer lines (recommended but not default)
;   /NOGREY         - do not grey-shift (applicable when multiple
;     standards are passed)
;   WRANGE      - Vector with [WAVE_MIN, WAVE_MAX] to be
;                     considered in the fit
;
; OUTPUTS:
;   mag_set        - structure containing b-spline info for sensitivity function
;
; OPTIONAL OUTPUTS:
;  sensfunc         - sensitivity function evaluated
;
; COMMENTS:
;                   See README file in /apps2/iraf211/iraf/noao/lib/onedstds/
;                   for list of standard stars and the names of the
;                   associated files
;
; EXAMPLES:
;
; BUGS:
;                   Does not take into account atmospheric extinction!!!
;                   Leaves out first and last wavelength bins of
;                   sensitivity function
;
; PROCEDURES CALLED:
;   traceset2xy (idlutils)
;   xy2traceset (idlutils)
;   splog       (idlutils)
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   01-Oct-2005  Written by J. Hennawi UC Berkeley
;   08-Jun-2011  J. Moustakas, UCSD - optionally deal with a list of
;     standards, and generate a QAplot
;   08-Feb-2017  K-G Lee, LBL - Added WRANGE
;------------------------------------------------------------------------------
FUNCTION long_sensfunc, scifile, sensfuncfile $
                        , SCIHDR = SCIHDR1, STANDARD_STRUCT = STANDARD_STRUCT1 $
                        , std_name = std_name $
                        , sensfunc = sensfunc, sensfit = sensfit $
                        , wave = wave, flux=flux, nresln = nresln $
                        , chk = chk, MSK_BALM = msk_balm $
                        , LINE_BALM = LINE_BALM $
                        , BALM_MASK_WID = BALM_MASK_WID $
                        , sciind = sciind1, std_flux = flux_std_int $
                        , STDFILE = STDFILE,  INMASK = INMASK $
                        , FINALMASK = FINALMASK, WRANGE=wave_range $
                        , IVAR=ivar, resln=resln, nofit=nofit, nogrey=nogrey

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'sfunc = long_sensfunc(stdfil, outfil, STD_NAME=) [v1.0]'
      return, -1
  endif 

; check whether SCIFILE and STD_NAME vectors
  nstd = n_elements(scifile)
  if (nstd gt 1) then begin
     for ii = 0, nstd-1 do begin
        sens1 = long_sensfunc(scifile[ii],sensfuncfile,inmask=inmask1,$
          std_name=std_name[ii],wave=wave1,flux=flux1,ivar=ivar1,$
          nresln=nresln,MSK_BALM=msk_balm,BALM_MASK_WID=BALM_MASK_WID,$
                              std_flux=std_flux1,STDFILE=STDFILE,/nofit,resln=resln1, $
                              wrange=wave_range)
        if (ii eq 0) then begin
           wave = wave1
           flux = flux1
           ivar = ivar1
           std_flux = std_flux1
           resln = resln1
           if (n_elements(inmask1) ne 0) then begin
              inmask = inmask1
           endif
           npix = n_elements(wave1)
        endif else begin
           wave = [[wave],[wave1]]
           flux = [[flux],[flux1]]
           ivar = [[ivar],[ivar1]]
           std_flux = [[std_flux],[std_flux1]]
           resln = [resln,resln1]
        endelse
     endfor
     meanwave = fix(djs_mean(wave))
     dwave = (max(wave)-min(wave))*0.01
     grey = fltarr(nstd)
; calculate the grey shifts, but dont necessarily apply them
     for ii = 0, nstd-1 do begin
        good = where((wave[*,ii] gt meanwave-dwave) and $
          (wave[*,ii] lt meanwave+dwave) and (std_flux[*,ii] gt 0) and $
          (flux[*,ii] gt 0) and (ivar[*,ii] gt 0),ngood)
        if (ngood eq 0) then message, 'Fix me!'
        grey[ii] = djs_median(flux[good,ii]/std_flux[good,ii])
     endfor
     grey = 2.5*alog10(max(grey)/grey)
     if (keyword_set(nogrey) eq 0) then for ii = 0, nstd-1 do $
       flux[*,ii] = flux[*,ii]*10^(0.4*grey[ii])
; now do the fit
     if (n_elements(inmask) ne 0) then $
       thismask = rebin(reform(inmask,npix,1),npix,nstd) else $
         thismask = wave*0+1
     wsrt = sort(wave)
     wave = reform(wave[wsrt],npix*nstd)
     flux = reform(flux[wsrt],npix*nstd)
     ivar = reform(ivar[wsrt]*thismask[wsrt],npix*nstd)
     std_flux = reform(std_flux[wsrt],npix*nstd)
     
     bkspace = max(resln)*nresln*2

     mag_set = bspline_magfit(wave,flux,ivar,std_flux,/nointerp,$
       bkspace=bkspace,maxiter=10,maxrej=5,upper=2.0,$
       lower=2.0,sensfit=sensfit,sensfunc=sensfunc,$
       wave_min=wave_min,wave_max=wave_max,outmask=outmask)

     if keyword_set(CHK) then x_splot, wave, sensfit, /blo 
     IF KEYWORD_SET(sensfuncfile) THEN mwrfits, mag_set, sensfuncfile, /create

; make a qaplot     
     IF KEYWORD_SET(sensfuncfile) THEN begin
        psfile = repstr(sensfuncfile,'.fits','.ps')
        dfpsplot, psfile, /color, /square
     endif
     sensfunc_plot = 2.5*alog10(sensfunc>0.1)+16 ; standards are in units of 10^-16 erg/s/cm2/A
     sensfit_plot = 2.5*alog10(sensfit>0.1)+16
     djs_iterstat, sensfunc_plot-sensfit_plot, sigrej=2.5, invvar=ivar gt 0, mask=msk

     yrange = minmax(sensfunc_plot[where(msk eq 1)])
     djs_plot, wave, sensfunc_plot, ps=6, sym=0.2, xsty=3, ysty=3, $
       xtitle='Wavelength (\AA)', $ ; yrange=[max(sensfit),min(sensfit)], $
       ytitle='+2.5 log_{10} SF [(e^{-} s^{-1}) / (erg s^{-1} cm^{-2} \AA^{-1})]', $
       yrange=reverse(yrange), charthick=2, color='grey'
     djs_oplot, wave, sensfit_plot, thick=6, color='red';, psym=6, symsize=0.2
     plots, meanwave, interpol(sensfit_plot,wave,meanwave), $
       psym=6, symsize=2, thick=6
     srt = sort(grey)
     legend, repstr(repstr(file_basename(scifile[srt]),'.fits',''),'.gz','')+' '+$
       textoidl('\Delta'+'m=')+strtrim(string(grey[srt],format='(F12.3)'),2), $
       /right, /bottom, box=0, charthick=2, charsize=1.4
     IF KEYWORD_SET(sensfuncfile) THEN dfpsclose
     return, mag_set
  endif

; stay five pixels away from edge
IF NOT KEYWORD_SET(NRESLN) THEN NRESLN = 20  ;; JXP -- Avoid 'narrow' features
IF NOT KEYWORD_SET(LO_BUFFER) THEN LO_BUFFER = 8L
IF NOT KEYWORD_SET(HI_BUFFER) THEN HI_BUFFER = 5L
IF NOT KEYWORD_SET(BALM_MASK_WID) THEN BALM_MASK_WID = 5.0D

IF KEYWORD_SET(scihdr1) AND KEYWORD_SET(STANDARD_STRUCT1) THEN BEGIN
   scihdr = scihdr1
   standard_struct = standard_struct1
ENDIF ELSE BEGIN
   scihdr = xheadfits(scifile)
   temp = xmrdfits(scifile, 5)
   IF NOT KEYWORD_SET(SCIIND1) THEN BEGIN maxf = max(temp.PEAKFLUX, sciind)
   ENDIF ELSE sciind = sciind1 - 1L 
   standard_struct = temp[SCIIND]
ENDELSE

IF KEYWORD_SET(STDFILE) THEN BEGIN
   readcol, stdfile, wave_std, flux_std1, format = 'D,D', /silent
   flux_std = 10.0d*flux_std1 
;; this factor of 10 converts 1.0e-16 erg/cm^2/s/A to units of 1.0e-17
;; erg/cm^2/s/A. Since the xidl stds have that normalization
;; Not all do..  
ENDIF ELSE BEGIN
   longslit_dir = getenv('LONGSLIT_DIR')
   std_file = longslit_dir+ '/calib/standards/calspec/' + std_name + '.fits.gz'
; read in standar star spectrum
   std = xmrdfits(std_file, 1)
   wave_std = std.WAVELENGTH
   flux_std = 1.0d17*std.FLUX   ; fluxes are in units of 1.0e-17 erg/cm^2/s/A
ENDELSE
   
; uncalibrated observed spectrum
wave = standard_struct.wave_box
flux = standard_struct.flux_box
ivar = standard_struct.ivar_box
IF KEYWORD_SET(INMASK) THEN BEGIN
   IF n_elements(INMASK) NE n_elements(ivar) THEN $
      message, 'Your mask must be aligned with the std spectrum'
   badpix = WHERE(INMASK EQ 0, nbad)
   IF nbad GT 0 THEN ivar[badpix] = 0.0
ENDIF
; parse headers and read in extinctino file
ext = long_extinct(wave, scihdr, AIRMASS = AIRMASS, EXPTIME = EXPTIME)
; extinction correct data and divide by exposure time
flux = flux*ext
ivar = ivar/ext^2

; find the min and max of the calibration spectrum
wave_min_std = min(wave_std)
wave_max_std = max(wave_std)
; find the min and max of the standard
if keyword_set(wave_range) then begin
   wave_min_obs = min(wave_range)
   wave_max_obs = max(wave_range)
endif else begin
   ind_sort = sort(wave)
   wave_min_obs = min(wave)     ;wave[ind_sort[LO_BUFFER-1]]
   wave_max_obs = max(wave)     ;wave[ind_sort[nwave-1L-HI_BUFFER]]
endelse
   
wave_min = wave_min_std > wave_min_obs
wave_max = wave_max_std < wave_max_obs 

calib_inds = WHERE(wave GE wave_min AND wave LE wave_max)
wave = wave[calib_inds]
flux = flux[calib_inds]
ivar = ivar[calib_inds]
sort_ind = sort(wave)
wave = wave[sort_ind]
flux = flux[sort_ind]
ivar = ivar[sort_ind]
nwave = n_elements(wave)
;; Dont use the edges in the fit
;IF KEYWORD_SET(LO_BUFFER) THEN ivar[0:LO_BUFFER] = 0
;IF KEYWORD_SET(HI_BUFFER) THEN ivar[nwave-1L-HI_BUFFER:nwave-1L] = 0

;interpolate calbiration spectrum onto observed wavelengths
flux_std_int = interpol(flux_std, wave_std, wave)

; Compute an effective resolution for the standard. This could be improved
; to setup an array of breakpoints based on the resolution. At the 
; moment we are using only one number
std_res = 2.0*djs_median(abs(wave_std - shift(wave_std, 1)))
IF TAG_EXIST(standard_struct, 'PIX_RES') THEN BEGIN
    pix_res = djs_median(standard_struct.PIX_RES)
    disp = djs_median(abs(wave-shift(wave, 1)))
    resln = pix_res*disp >  std_res
ENDIF ELSE resln = std_res


if keyword_set(MSK_BALM) then begin
   IF NOT KEYWORD_SET(LINE_BALM) THEN $
      line_balm = [3836.4 $
                   , 3969.6 $
                   , 3890.1 $
                   , 4102.8 $                        ;; (H9,CaII H, Hf,H-delta)
                   , 4102.8 $
                   , 4341.6 $
;;           , 4687.3 $
                   , 4862.7 $ ;; (H-gamma,H-beta,H-alpha,Fe)
                   , 5407.0 $
                   , 6564.6 $
                   , 8224.8 $
                   , 8239.2] 
   nbalm = n_elements(line_balm)
   for qq = 0L, nbalm-1 do begin
      mskwv = where(abs(wave-line_balm[qq]) LE BALM_MASK_WID*resln, nbad)
      if nbad NE 0 then ivar[mskwv] = 0.
   endfor
endif
;; Mask telluric absorption
tell = (wave GE 7580.0D AND wave LE 7750.0D) OR $
       (wave GE 7160.0D AND wave LE 7340.0D) OR $
       (wave GE 6860.0  AND wave LE 6930.0D)
tell_ind = where(tell, ntell)
IF ntell GT 0 THEN ivar[tell_ind] = 0.0
finalmask = ivar GT 0

if (keyword_set(nofit) eq 0) then begin ; jm09dec19ucsd
   mag_set = bspline_magfit(wave, flux, ivar, flux_std_int $
     , bkspace = resln*nresln $
;     , bkspace = resln*nresln $
     , maxiter = 10, maxrej = 5, upper = 3, lower = 3 $
     , sensfit = sensfit, sensfunc = sensfunc $
     , wave_min = wave_min, wave_max = wave_max $
     , outmask = outmask) 
   
   if keyword_set(CHK) then x_splot, wave, sensfit, /blo 
   IF KEYWORD_SET(sensfuncfile) THEN mwrfits, mag_set, sensfuncfile, /create
endif else return, 1


RETURN, mag_set
END
