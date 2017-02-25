;+
; NAME:
;   gmos_transimg
;
; PURPOSE:
;
;   Main program for the Low-redux pipeline.  This set of algorithms
;   runs mainly as a black box.
;
; CALLING SEQUENCE:
;  
;
; INPUTS:
;  utransimg -- File created by long_plan which guides the reduction
;               process
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;               pixels)
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
;
; OUTPUTS:
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;-----------------------------------------------------------------------------
;;flag=0 means from from untransformed coords to transformed coords
;;flag=1 means go from transformed coords to untransformed coords

PRO gmos_trnimg3to1, utranfile, outimg, outivar, hdr = hdr $
                  , gain = gain, rnoise = rnoise $
                  , verbose = verbose, bin = bin

  ;; Determine raw image, binning, and overscan sizes from headers. 
  hdrt = xheadfits(utranfile, exten = 1)
  datasec = strcompress(sxpar(hdrt, 'DATASEC'), /rem)
  data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
  nx_raw = data_arr[1]
  ny_raw = data_arr[3]
  
  xbin = long(2L*1056/nx_raw)
  ybin = long(4L*1152/ny_raw)
  
  IF xbin EQ 2 AND ybin EQ 2 THEN BTAG = '2x2' $
  ELSE IF xbin EQ 1 AND ybin EQ 1 THEN BTAG = '1x1' $
  ELSE IF xbin EQ 2 AND ybin EQ 1 THEN BTAG = '2x1' $
  ELSE message, 'Your binning is not supported'

  nx = 2048/xbin
  ny = ny_raw

  ;; Note that x and y are reversed because we are transposed
  gmos_gaps, ny, nx, ygap1 = ygap1, ygap2 = ygap2
  nmx = 3*nx + ygap1 + ygap2
  
  outimg = fltarr(nmx, ny)
  outivar = fltarr(nmx, ny)
  
  raw1 = fltarr(nx_raw, ny_raw)
  raw2 = fltarr(nx_raw, ny_raw)
  raw3 = fltarr(nx_raw, ny_raw)
  ivar1 = fltarr(nx_raw, ny_raw)
  ivar2 = fltarr(nx_raw, ny_raw)
  ivar3 = fltarr(nx_raw, ny_raw)
  
  ;; For chips 1 and 2 have oscan at the beginning
  datacols = lindgen(nx) + data_arr[0] - 1
  ;; For chips 3 the oscan is at the end
  datacols3 = lindgen(nx)
  long_oscan, utranfile, img_dum, ivar_dum, ccdonly = 1, hdr = hdr $
              , gain = gain, rnoise = rnoise, bin = bin, verbose = verbose
  ;; transform back to native orientation
  img_dum = reverse(transpose(img_dum), 1)
  ivar_dum = reverse(transpose(ivar_dum), 1)
  raw1[datacols, *] = img_dum
  ivar1[datacols, *] = ivar_dum
  ;; Read in image 2
  long_oscan, utranfile, img_dum, ivar_dum, ccdonly = 2, hdr = hdr $
              , gain = gain, rnoise = rnoise, bin = bin, verbose = verbose
  ;; transform back to native orientation
  img_dum = reverse(transpose(img_dum), 1)
  ivar_dum = reverse(transpose(ivar_dum), 1)
  raw2[datacols, *] = img_dum
  ivar2[datacols, *] = ivar_dum
  ;; Read in image 3
  long_oscan, utranfile, img_dum, ivar_dum, ccdonly = 3, hdr = hdr $
              , gain = gain, rnoise = rnoise, bin = bin, verbose = verbose
  ;; transform back to native orientation
  img_dum = reverse(transpose(img_dum), 1)
  ivar_dum = reverse(transpose(ivar_dum), 1)
  raw3[datacols3, *] = img_dum
  ivar3[datacols3, *] = ivar_dum

  ;; transform ccd1
  path = getenv('LONGSLIT_DIR') + '/calib/transform/GMOS/'
  rdfloat, path + 'chip1-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  polywarp, ux, uy, tx, ty, 2, kx, ky
  cutim1 = poly_2d(raw1, kx, ky, 1, nmx, ny, missing = -99999)
  ivarm1 = poly_2d(ivar1, kx, ky, 1, nmx, ny, missing = -99999)
  ;; transform ccd2
  rdfloat, path + 'chip2-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  polywarp, ux, uy, tx, ty, 2, kx, ky
  cutim2 = poly_2d(raw2, kx, ky, 1, nmx, ny, missing = -99999)
  ivarm2 = poly_2d(ivar2, kx, ky, 1, nmx, ny, missing = -99999)
  ;; transform ccd3
  ;; chip3 only have 169 elements?
  rdfloat, path + 'chip3-'+btag+'.coo', tx, ty, ux, uy, /double, /silent
  polywarp, ux, uy, tx, ty, 2, kx, ky
  cutim3 = poly_2d(raw3, kx, ky, 1, nmx, ny, missing = -99999)
  ivarm3 = poly_2d(ivar3, kx, ky, 1, nmx, ny, missing = -99999)
  
  qui = where(cutim1 gt -9999 AND ivarm1 gt -9999)
  outimg[qui] = cutim1[qui]
  outivar[qui] = ivarm1[qui]
  qui = where(cutim2 gt -9999 AND ivarm2 gt -9999)
  outimg[qui] = cutim2[qui]
  outivar[qui] = ivarm2[qui]
  qui = where(cutim3 gt -9999 AND ivarm3 gt -9999)
  outimg[qui] = cutim3[qui]
  outivar[qui] = ivarm3[qui]
  ;; now transform to my reverse transpose convention
  outimg  = transpose(reverse(outimg, 1))
  outivar = transpose(reverse(outivar, 1))
  
  RETURN
END

  
