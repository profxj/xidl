;+
; NAME:
;        TSPEC_CLEAN
;
;
; PURPOSE:
;        Removes quadrant-to-quadrant crosstalk and capacitive
;        coupling effects in raw Palomar TSPEC spectra.
;
;
; CALLING SEQUENCE:
;        result = tspec_clean(image[,parameters accepted by readfits])
;
;
; INPUTS:
;        image: A 2048x1024 Palomar TSPEC image OR a Palomar TSPEC
;        fits file.  This way you can use it in the place of readfits
;
;
; OPTIONAL INPUTS:
;        Any named parameters accepted by readfits with be passed.
;
;
; KEYWORD PARAMETERS:
;        sigma: optinally NAN all deviant pixels changed by
;        sigma_filter from the Astro IDL library
;
;
; OUTPUTS:
;        result: A 2048x1024 DOUBLE array containing the cleaned
;        image, with crosstalk and capacitive coupling removed.
;
;
; OPTIONAL OUTPUTS:
;        header: String containing the header of the fits file passed,
;        if a fits file is passed for image.
;
;        variance: A 2048x1024 DOUBLE array containing the calculated
;        photon-noise variances of the image (sigma^2).
;
; RESTRICTIONS
;        Requires the IDL Astronomy User's Library: http://idlastro.gsfc.nasa.gov/contents.html
;
; SIDE EFFECTS:
;        Converts image to double.
;
;
; MODIFICATION HISTORY:
;  Written by Phil Muirhead <muirhead@astro.cornell.edu> with
;       parameters derived by Terry Herter <herter@astro.cornell.edu>.
;
;  Modififed by JXP for RN (April 2011)
;
;       Version 1.0
;-


function tspec_clean, exposure, header, ivar=ivar, sigma=sigma, _extra=extra

; Set the default xtalk factor to 170
 xtalk_factor = 170
  
; If exposure is a filename, read in the file
  if ( (n_params() NE 0) AND (size(exposure, /tname) EQ 'STRING')) then begin
     ifexists = findfile(exposure +'*', count = count)
     if (count EQ 0) then begin
        print, 'ERROR: File not found!'
        return, -1
     endif else begin
        ;;image = readfits(exposure, header, _extra = extra)
        image = xmrdfits(exposure, 0, header,  _extra = extra)
     endelse
  endif else begin
  
; Check for existence of array
     if ( (n_params() NE 0) AND (size(exposure, /tname) NE 'STRING') AND $
          (size(exposure, /tname) EQ 'UNDEFINED')) then begin
        print, 'ERROR: Data array does not exist!'
        return, -1
     endif else begin
        image = exposure
     endelse

  endelse

; Check to make sure it has the right elements
  if (((size(image))[0] NE 2) OR $
      ((size(image))[1] NE 2048) OR $
      ((size(image))[2] NE 1024)) then begin
     print, 'ERROR: Input data must be a TSPEC exposure (2048x1024)!'    
     return, -1
  endif

; Convert to double for algorhythms
  image = double(temporary(image))

; Correct for crosstalk
  for i=0, 1023 do begin
     left_row = mean(image[0:1023,i],/nan,/double) / double(xtalk_factor)
     right_column = mean(image[1024+i,*],/nan,/double) / double(xtalk_factor)
     image[0:1023,i] -= (left_row + right_column)
     image[1024+i,*] -= (left_row + right_column)
  endfor

; Correct for capacitive coupling
; Ok, now run a derivative filter across the image
; Left-to-right derivative for the left quadrant

  filtered = dblarr(2048,1024)
  left_filt = [0d,-1d,1d]
  left_filtered = convol(image, left_filt,/center)
  right_filt = [[-1d],[1d],[0d]]
  right_filtered = convol(image, right_filt,/center)

  deriv_image = [left_filtered[0:1023,*],right_filtered[1024:2047,*]]

; Terry's matricies in % of derivativ to subtract in each channel
  CC_Q1n = $
     [[ 0.00,  1.30,  0.19,  0.13, -0.11, -0.10, -0.15, -0.26], $
      [ 1.51,  0.00,  0.84,  0.29, -0.10, -0.11, -0.15, -0.21], $
      [ 0.34,  0.94,  0.00,  1.34,  0.28,  0.22, -0.13, -0.22], $
      [ 0.52,  0.51,  1.39,  0.00,  0.88,  0.26, -0.16, -0.22], $
      [ 0.00, -0.07,  0.22,  0.82,  0.00,  1.36,  0.28, -0.01], $
      [ 0.05, -0.08,  0.23,  0.16,  1.33,  0.00,  1.15, -0.09], $
      [-0.02, -0.11, -0.08, -0.17,  0.61,  1.25,  0.00,  1.18], $
      [ 0.01,  0.02,  0.01, -0.14,  0.26,  0.07,  1.48,  0.00]]
  
  CC_Q1p = $
     [[ 0.00,  1.08,  0.20,  0.14, -0.08, -0.05, -0.12, -0.25], $
      [ 1.27,  0.00,  0.75,  0.23, -0.05, -0.06, -0.11, -0.21], $
      [ 0.26,  0.76,  0.00,  1.09,  0.22,  0.22, -0.10, -0.20], $
      [ 0.51,  0.42,  1.23,  0.00,  0.77,  0.22, -0.12, -0.21], $
      [ 0.04, -0.05,  0.18,  0.67,  0.00,  1.12,  0.24,  0.00], $
      [ 0.06, -0.03,  0.19,  0.15,  1.14,  0.00,  0.95, -0.08], $
      [-0.01, -0.02, -0.04, -0.15,  0.51,  1.02,  0.00,  1.00], $
      [ 0.05,  0.03, -0.03, -0.08,  0.22,  0.08,  1.21,  0.00]] ;
  
  CC_Q2 = $
     [[  0.0, 1.38, -0.13, 0.17, 0.07, 0.04, 0.17, 0.13], $
      [ 1.12, 0.00,  0.82, 0.36, 0.04, 0.02, 0.10, 0.09], $
      [ 0.01, 1.08,  0.00, 1.05, 0.33, 0.28, 0.16, 0.13], $
      [ 0.11, 0.32,  0.98, 0.00, 0.82, 0.22, 0.14, 0.10], $
      [-0.08,-0.11, -0.05, 0.65, 0.00, 1.30, 0.52, 0.62], $
      [-0.11,-0.09,  0.00, 0.10, 1.25, 0.00, 0.89, 0.39], $
      [-0.05,-0.08, -0.29, 0.00, 0.35, 0.81, 0.00, 1.22], $
      [-0.07,-0.08, -0.20, 0.00, 0.28, 0.28, 1.26, 0.00]] ;
  
;;   channel=4
;;   channel_corrected=5
; Start with right quadrant
  for channel_corrected=0, 7 do begin
     for channel=0, 7 do begin
        image[1024:2043, channel_corrected * 128: channel_corrected * 128 + 127] -= $
           deriv_image[1024:2043, channel * 128: channel * 128 + 127] * $
           CC_Q2[channel,channel_corrected] * 0.01
     endfor
  endfor  


; Now for 1st quadrant
  for channel_corrected=0, 7 do begin
     for channel=0, 7 do begin
        image[channel_corrected * 128: channel_corrected * 128 + 127,0:1023] += $
           deriv_image[channel * 128: channel * 128 + 127,0:1023] * $
           (CC_Q1p[channel,channel_corrected] * 0.01 + $
            CC_Q1n[channel,channel_corrected] * 0.01)/2.
     endfor
  endfor  

; NAN bad pixels
  if keyword_set(sigma) then begin
     print, 'Running sigma_filter on image...'
     filtered = sigma_filter(image,11,/iterate,n_sigma=5)
     filter_pix = where(image ne filtered, filter_count)
     IF filter_count GT 0 THEN image[filter_pix] = 0.0 ;;!values.d_nan
  endif
  
; NAN saturating pixels
  saturating = where(image ge 25000d, saturating_count)
  ;;??????  What should we do with saturated pixels??
  ;;if saturating_count ne 0 then image[saturating] = 0.0 ;;!values.d_nan
  
  ;; Convert to electrons
  gain = double((strsplit(sxpar(header, 'GAIN'), /extract))[0])
;  gain = 3.8d ;ph/ND
  
  ;; This assumes FOWLER sampling
  dnsamp = double((strsplit(sxpar(header, 'FSAMPLE'), /extract))[0])
  rnoise = 14.0D/sqrt(dnsamp)

  image = gain*image
  
  ;; Inverse variance :: The dark current appears low and there are no
  ;; bias counts.  Probably do not need to account for the overscan in
  ;; the statistics.
  ivar = 1.0/(abs(image - sqrt(2.0)*rnoise) +rnoise^2)
  
  ;; JFH 08/04/12 I don't know who wrote down these lines below
  ;; but they are incorrect. Also appears that image is not being
  ;; being multiplied by the gain. 
  ;variance = image/gain + rnoise^2
  ;ivar = 1./(variance)
  IF saturating_count GT 0 THEN ivar[saturating] = 0.0
  IF KEYWORD_SET(SIGMA) THEN BEGIN
     IF filter_count GT 0 THEN ivar[filter_pix] = 0.0
  ENDIF
  
  ;; Transpose
  image = transpose(reverse(image, 1)) 
  ivar = transpose(reverse(ivar, 1)) 

  return, image

end

  
