; NAME:
;   long_oscan
;
; PURPOSE:
;   Overscan subtract a standard CCD.  There are separate subroutines
;   for LRIS, DEIMOS, BCS, GMOS and KAST.
;
; CALLING SEQUENCE:
;   long_oscan, filename, [ rawsub, rawivar, hdr= , $
;    gain= , rnoise= , image=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   verbose    - If set, then verbose output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   rawsub     - Bias subtracted image
;   rawivar    - Inverse variance of image
;   hdr        - Fits header of raw image
;   gain       - vector of gains for each amp
;   rnoise     - vector or readnoise for each amps
;   image      - original data image
;   bin        - Binning of the CCD [x,y]
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   mrdfits
;   djs_avsigclip
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB) and Scott Burles (MIT)
;-
;------------------------------------------------------------------------------

function gmos_badpixelmask, ncol, nrow, ccd_num

  mask = bytarr(ncol, nrow) + 1B

  xbin = long(2048/ncol)
  ybin = long(4608/nrow)

  if ccd_num EQ 1 then begin
     mask[88/xbin:89/xbin, *] = 0
     mask[772/xbin:775/xbin, 708*4/ybin:*] = 0
     mask[986/xbin:987/xbin, 768*4/ybin:*] = 0
     return, mask
  endif

  if ccd_num EQ 2 then begin
     mask[1252/xbin:1255/xbin, 820*4/ybin:*] = 0
     mask[1280/xbin:1281/xbin, 315*4/ybin:*] = 0
     return, mask
  endif

  ;; no bad columns for CCD 3?
  return, mask
end


; add support for the Bok/B&C Spectrograph
pro bok_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)

data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region too close
; to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

; generate the bad-pixel mask - assumes no spectral binning!!!
mask = rawsub*0+1
;mask[405-1:408-1,*] = 0
;mask[507-1:507-1,*] = 0
;mask[538-1:538-1,*] = 0
;mask[650-1:650-1,*] = 0
;mask[681-1:681-1,*] = 0
;mask[706-1:706-1,*] = 0
;mask[767-1:771-1,*] = 0
;mask[897-1:898-1,*] = 0

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar*mask) ; zero out!

sxaddpar, hdr, 'NAXIS1', (size(rawsub,/dim))[0] ; note!
sxaddpar, hdr, 'NAXIS2', (size(rawsub,/dim))[1]
sxdelpar, hdr, 'biassec'
sxdelpar, hdr, 'trimsec'
sxdelpar, hdr, 'datasec'
sxdelpar, hdr, 'ccdsec'
sxdelpar, hdr, 'origsec'
RETURN
END

; add support for the Bok/B&C Spectrograph
pro dupont_oscan, filename, rawsub, rawivar, hdr = hdr $
                  , gain = gain, rnoise = rnoise, image = image $
                  , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]

if sxpar(hdr, 'BINNING') NE '1x1' then stop
bin_col = 1
bin_row = 1
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'EGAIN'))
rnoise = double(sxpar(hdr, 'ENOISE'))

nbias = ncol - 2050 - 1
biascols = lindgen(nbias) + 2050
oscan1 = djs_avsigclip(image[biascols, 1:512], 1)

sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])

oscan = convol(oscan1, kernel, /edge_truncate)
imagecol = 2048L
imagerow = 513L

osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 1:imagerow-1L] - osub)*gain

rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; Rotate
rawsub = rotate(rawsub,3)
rawivar = rotate(rawivar,3) ; zero out!

sxaddpar, hdr, 'NAXIS1', (size(rawsub,/dim))[0] ; note!
sxaddpar, hdr, 'NAXIS2', (size(rawsub,/dim))[1]
;sxdelpar, hdr, 'biassec'
;sxdelpar, hdr, 'trimsec'
;sxdelpar, hdr, 'datasec'
;sxdelpar, hdr, 'ccdsec'
;sxdelpar, hdr, 'origsec'
RETURN
END

pro gmos_oscan, filename, flux, invvar, hdr = hdr, thdr = thdr $
                , gain = gain, rnoise = rnoise_vec, image = image $
                , VERBOSE = VERBOSE, bin = bin, CCDONLY = CCDONLY


ccd = [1, 2, 3]
IF NOT KEYWORD_SET(CCDONLY) THEN i = ccd[0] $
ELSE i = ccd[ccdonly-1]
;  no overscan for rows...
;  look for filename
filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0], exten = 0)
rawdata = xmrdfits(filelist[0], i, hdrt, /fscale)

if size(hdr0, /tname) EQ 'INT' then begin
    print, 'having trouble with ', filename
    return
endif
telescope = strcompress(sxpar(hdr0, 'TELESCOP'), /rem)

detsize = sxpar(hdr0, 'DETSIZE')
hccds   = sxpar(hdr0, 'NCCDS')
hamps   = sxpar(hdr0, 'NAMPS')

if hamps NE 1 then begin
    print, 'Only ready for single amplifier reads'
    return
endif

if detsize NE '[1:6144,1:4644]' OR hccds NE 3 then begin
    print, 'Expecting 3 CCDs of total size: [1:6144,1:4644]'
    return
endif

rnoise = sxpar(hdrt, 'RDNOISE')
gain    = sxpar(hdrt, 'GAIN')

rnoise_vec = fltarr(3)
rnoise_vec[i-1L] = rnoise

IF KEYWORD_SET(VERBOSE) THEN print, "CCD Number: ", i, " Gain:", gain $
  , " Read Noise:", rnoise, format = '(a,i2,a,f6.3,a,f5.2)'

rawcol = (size(rawdata))[1]
rawrow = (size(rawdata))[2]

biassec = strcompress(sxpar(hdrt, 'BIASSEC'), /rem)
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
oscan_pix=bias_arr[1]-bias_arr[0] + 1L
colbin = long(2L*1024/(rawcol-oscan_pix))
rowbin = long(4L*1152/rawrow)
bin = [colbin, rowbin]
bin = reverse(bin)

ncol = 2048L/colbin
nrow = rawrow

colsplit = lindgen(64L/colbin) + long(strmid(sxpar(hdrt, 'BIASSEC'), 1)) -1
datacols = lindgen(ncol) + long(strmid(sxpar(hdrt, 'DATASEC'), 1)) - 1
oscan = djs_avsigclip(rawdata[colsplit[2:*], *], 1)
suboscan = rebin(transpose(oscan), ncol, nrow)

flux1 = (rawdata[datacols, *] - suboscan) * gain

invvar1 = 1.0/(abs(flux1 - sqrt(2.0)*rnoise) + rnoise^2)
IF strcmp(telescope, 'Gemini-North') THEN  invvar1 = invvar1*gmos_badpixelmask(ncol, nrow, i)


n0 = n_elements(hdr0)
nt = n_elements(hdrt)
hdr = strarr(max([n0, nt]), 4)
hdr[0:n0-1L, 0] = hdr0[0:n0-1L]
hdr[0:nt-1L, 1] = hdrt[0:nt-1L]
nccd = n_elements(ccd)
IF KEYWORD_SET(CCDONLY) THEN BEGIN
    flux = transpose(reverse(flux1, 1))
    invvar = transpose(reverse(invvar1, 1))
    RETURN
 ENDIF
flux1 = rebin(flux1, ncol, nrow, nccd)
invvar1 = rebin(invvar1, ncol, nrow, nccd)
for iccd = 1, nccd-1 do begin
    i = ccd[iccd]
    rawdata = xmrdfits(filename, i, hdrt, /silent, /fscale)
    hdrt    = xheadfits(filename, exten = i)
    rnoise = sxpar(hdrt, 'RDNOISE')
    rnoise_vec[i-1] = rnoise
    gain    = sxpar(hdrt, 'GAIN')
    IF KEYWORD_SET(VERBOSE) THEN $
      print, "CCD Number: ", i, " Gain:", gain, " Read Noise:", rnoise $
      , format = '(a,i2,a,f6.3,a,f5.2)'

    colsplit = lindgen(64L/colbin) +  $
      long(strmid(sxpar(hdrt, 'BIASSEC'), 1)) -1
    datacols = lindgen(ncol) + long(strmid(sxpar(hdrt, 'DATASEC'), 1)) -1
    oscan = djs_avsigclip(rawdata[colsplit[2:*], *], 1)
    ;; This kludge deals with the fact that the oscan is impacted by a
    ;; charge-trap for the middle chip on GMOS-South
    IF strcmp(telescope, 'Gemini-South') AND i EQ 2 THEN $
       oscan[0:(1000L/rowbin)] = median(oscan[(1000L/rowbin):*])
    suboscan = rebin(transpose(oscan), ncol, nrow)

    f = (rawdata[datacols, *] - suboscan) *gain

    invv = 1.0/(abs(f-sqrt(2.0)*rnoise) + rnoise^2)

;    badpixels from Adam's pipeline:
    IF strcmp(telescope, 'Gemini-North') THEN invv = invv*gmos_badpixelmask(ncol, nrow, i)

    flux1[*, *, i-1] = f
    invvar1[*, *, i-1] = invv

    hdr[0:nt-1L, iccd+1L] = hdrt
endfor

flux = gmos_mosaic(flux1)
invvar = gmos_mosaic(invvar1)

return
end

pro mmt_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif


;; Parse the header
gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))

      nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

case strtrim(detector,2) of
   '*ccd34*': BEGIN
      rawsub = transpose(rawsub)
      rawivar = transpose(rawivar)
;  different trim for red channel
      rawsub = rawsub[80:610, 160:1029] ; 125:574
      rawivar = rawivar[80:610, 160:1029]
      dims = size(image, /dim)
      bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too (not sure this was done right)
   END
   'mmtredchan': BEGIN
      ;; No trimming
      ;; Binning?
      dims = size(image, /dim)
      bin= [1,round(1032./dims[1])]  ;; This might not be right (JXP: 06/21/2010)
   end
   'mmtbluechan': BEGIN  ;; KHRR added 28Apr2014
      rawsub = reverse(transpose(rawsub),2)  
      rawivar = reverse(transpose(rawivar),2)
      binning = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
      bin1 = strmid(binning,0,1)
      bin2 = strmid(binning,1,2)
      bin = [long(bin2),long(bin1)]
   end 
   else: stop
endcase

;;  different trim for upgraded red channel (using transposed coordinates)
;;  keep commented for now since as it will require changing
;;  wavelength solution etc  24-04-2009
;IF strmatch(detector, 'mmtredchan*') THEN BEGIN
;    rawsub  =  rawsub[*, 7:1006]
;    rawivar = rawivar[*, 7:1006]
;ENDIF


RETURN
END

pro kpno_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar)
RETURN
END

pro caha_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBINX'))
bin_row = long(sxpar(hdr, 'CCDBINY'))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'CCDSENS'))
rnoise = double(sxpar(hdr, 'CCDRON'))
detector = strcompress(sxpar(hdr, 'CCDNAME'), /rem)
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[2]-bias_arr[0]

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = bias_arr[0] - data_arr[0] + 1L
imagerow = data_arr[3]-data_arr[1] + 1L

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; Pitch the first four pixels for both the red and blue side as these
;; are hot pixels.

rawsub  = rawsub[4:*, *]
rawivar = rawivar[4:*, *]
bin = reverse(bin)
IF strcmp(detector, 'SITe#22b_14') THEN BEGIN ;; BLUE SIDE
   rawsub = transpose(rawsub)
   rawivar = transpose(rawivar)
ENDIF ELSE IF strcmp(detector, 'SITe#20b_12') THEN BEGIN ;; RED SIDE
   rawsub  = reverse(transpose(rawsub), 2)
   rawivar = reverse(transpose(rawivar), 2)
ENDIF
RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro caha_22m_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan when windowed!
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBINX'))
bin_row = long(sxpar(hdr, 'CCDBINY'))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'CCDSENS'))
rnoise = double(sxpar(hdr, 'CCDRON'))
detector = strcompress(sxpar(hdr, 'CCDNAME'), /rem)
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
if data_arr[0] LT bias_arr[0] then stop ;; You can do overscan subtraction and need to trim

temp_image = image * gain  ;; Electrons
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2) ;; This is wrong!!

RETURN
END

pro ntt_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin
  
;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
nx = dims[0]
ny = dims[1]
ind_binx = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINX', /bool))
binx = double(strmid(hdr[ind_binx], 30, 14))
ind_biny = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINY', /bool))
biny = double(strmid(hdr[ind_binx], 30, 14))
bin = [binx, biny]

ind_gain = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 GAIN', /bool))
gain = (double(strmid(hdr[ind_gain], 30, 14)))[0]
ind_ron = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 RON', /bool))
rnoise = (double(strmid(hdr[ind_ron], 30, 14)))[0]
ind_det =  WHERE(stregex(hdr, 'HIERARCH ESO DET CHIP1 NAME', /bool))
detector = strcompress(strmid(hdr[ind_det], 30, 14))

ind_prscx = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 PRSCX', /bool))
prscx = (long(strmid(hdr[ind_prscx], 30, 14)))[0]
ind_ovscx = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 OVSCX', /bool))
ovscx = (long(strmid(hdr[ind_ovscx], 30, 14)))[0]

ind_prscy = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 PRSCY', /bool))
prscy = (long(strmid(hdr[ind_prscy], 30, 14)))[0]
ind_ovscy = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 OVSCY', /bool))
ovscy = (long(strmid(hdr[ind_ovscy], 30, 14)))[0]

data_arr = [prscx+1L, prscy+1L, nx-ovscx, ny-ovscy]
bias_arr = [prscx+1L, ny-ovscy+1L, nx-ovscx, ny]
nbias = bias_arr[3]-bias_arr[1] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 1L
imagex = data_arr[2]-data_arr[0] + 1L
imagey = data_arr[3]-data_arr[1] + 1L

rawsub = fltarr(imagex, imagey)
rawivar = fltarr(imagex, imagey)

;biascols = lindgen(nbias-oscan_buffer) + imagey + oscan_buffer + prscy[0]
;oscan1 = djs_avsigclip(image[0:imagex-1L, biascols], 0)
;sig_res = 7
;nhalf =  long(sig_res)*4L
;xkern = dindgen(2*nhalf+1)-nhalf
;kernel = gauss1(xkern, [0.0, sig_res, 1.0])
;oscan = convol(oscan1, kernel, /edge_truncate)
;osub = oscan # replicate (1, imagey)

;; EFOSC oscan is not a true overscan and suffers from bleeding
;; from the image. Must do bias subtraction
temp_image = (image[prscx:prscx + imagex-1L $
                    , prscy: prscy + imagey -1L])*gain 
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

IF bin[0] EQ 2.0 AND bin[1] EQ 2.0 THEN BEGIN
   badpixmaskfile = getenv('XIDL_DIR') + $
                    '/Spec/Longslit/calib/masks/badpix_2x2.fits.gz'
   badmask1 = xmrdfits(badpixmaskfile, 0, /fscale)
   badmask = badmask1[prscx:prscx + imagex-1L, prscy: prscy + imagey -1L]
   ;;rawsub = rawsub*(badmask GT 0.0)
   rawivar = rawivar*(badmask GT 0.0)
ENDIF
;; Now transpose so that spectral direction vertical with blue down
;rawsub = transpose(rawsub)
;rawivar = transpose(rawivar) 
;bin = reverse(bin)

RETURN
END


;;
pro imacs_oscan, filename, rawsub, rawivar, hdr = hdr $
                 , gain = gain, rnoise = rnoise, image = image $
                 , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too

gain = double(sxpar(hdr, 'EGAIN'))
rnoise = double(sxpar(hdr, 'ENOISE'))
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)


; Removed by JXP on 08 Dec 2014
;if keyword_set(LONGSLIT) then begin
;   ;; trim the spectra (Longslit???)
;   rawsub  = rawsub[*, 1000:*]
;   rawivar = rawivar[*, 1000:*]
;   ;; bad pixel mask (red camera)
;   rawsub[1688, 2504:*] = 0
;   rawivar[1688, 2504:*] = 0
;   rawsub[248, 2476:*] = 0
;   rawivar[248, 2476:*] = 0
;   rawsub[977, *] = 0
;   rawivar[977, *] = 0
;endif

;; First row is junk
rawivar[*,0] = 0

;; Transpose?  Will depend on chip#
ipos = strpos(filename, '.fits')
cchip = strmid(filename,ipos-2,2)

case cchip of
   'c1': begin
      rawsub = rotate(rawsub,3)
      rawivar = rotate(rawivar,3)
   end
   'c2': begin
      rawsub = rotate(rawsub,3)
      rawivar = rotate(rawivar,3)
   end
   'c3': begin
      rawsub = rotate(rawsub,3)
      rawivar = rotate(rawivar,3)
   end
   'c4': begin
      rawsub = rotate(rawsub,3)
      rawivar = rotate(rawivar,3)
   end
   else: stop
endcase

RETURN
END



pro dis_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBIN1'))
bin_row = long(sxpar(hdr, 'CCDBIN2'))
bin = [1, 1]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

If strmatch(detector, '*blue*') THEN BEGIN
    rawsub = transpose(rawsub)
    rawivar = transpose(rawivar)
ENDIF ELSE IF strmatch(detector, '*red*') THEN BEGIN
    rawsub  = reverse(transpose(rawsub), 2)
    rawivar = reverse(transpose(rawivar), 2)
    ;; mask bad column on red side
    rawsub[212:*, 773] = 0.0D
    rawsub[212:262, 770:777] = 0.0
    rawivar[212:*, 773] = 0.0D
    rawivar[212:262, 770:777] = 0.0
ENDIF ELSE message, 'Unrecognized detector'

RETURN

END

pro p200_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))

bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RON'))
detector = sxpar(hdr, 'FPA')
ccd = sxpar(hdr, 'DETNAM')
datasec = strcompress(sxpar(hdr, 'TSEC1'), /rem)
biassec = strcompress(sxpar(hdr, 'BSEC1'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]-data_arr[0] + 1L
icol = data_arr[0] -1L
imagerow = data_arr[3]-data_arr[2] + 1L
biascol = bias_arr[0]-1L
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)
biascols = lindgen(nbias-oscan_buffer) + biascol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
oscan = djs_median(oscan1, width = 200, boundary = 'reflect')
osub = replicate(1, imagecol) # oscan
temp_image = (image[icol:icol+imagecol-1L, *] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

If strmatch(detector, 'DBSP_BLUE') THEN BEGIN
   rawsub  = reverse(rawsub, 2)
   rawivar = reverse(rawivar, 2)
   dims = size(rawsub, /dim)
   case ccd of
      'CCD44-82':
      else: begin
         IF dims[1] NE 1880 THEN BEGIN
            IF dims[1] GT 1880 THEN BEGIN
               rawsub1 = fltarr(dims[0], 1880)
               rawivar1 = fltarr(dims[0], 1880)
               rawsub1 = rawsub[*, 550L:(1880+550-1L)]
               rawivar1 = rawivar[*,550L:(1880+550-1L)]
;          rawsub1 = rawsub[*, 788L:(788+1880-1L)]
;          rawivar1 = rawivar[*, 788:(788+1880-1L)]
               rawsub = rawsub1
               rawivar = rawivar1
               splog, 'WARNING: Truncating P200 blue side arrays to 1880'
            ENDIF ELSE BEGIN
               rawsub1 = fltarr(dims[0], 1880)
               rawivar1 = fltarr(dims[0], 1880)
               rawsub1[*, 0:dims[1]-1L] = rawsub
               rawivar1[*, 0:dims[1]-1L] = rawivar
               rawsub = rawsub1
               rawivar = rawivar1
               splog, 'WARNING: Padding P200 blue side arrays to 1880'
            ENDELSE
         endif
      end
   endcase
ENDIF ELSE IF strmatch(detector, 'DBSP_RED') THEN BEGIN
   rawsub  = transpose(rawsub)
   rawivar = transpose(rawivar)
;;MF2013 added this two lines
ENDIF ELSE IF strmatch(detector, 'DBSP_RED2') THEN BEGIN
   rawsub  = transpose(rawsub)
   rawivar = transpose(rawivar)
ENDIF ELSE message, 'Unrecognized detector'

RETURN

END



pro mars_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

bin = [1,1]

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

ncol = (size(image))[1]
nrow = (size(image))[2]

biassec = [long(strmid(sxpar(hdr, 'BIASSEC'), 1))-1L $
           , long(strmid(sxpar(hdr, 'BIASSEC'), 6))-1L]
postcol = biassec[1]-biassec[0] + 1L

rowtrim = [215, 544]
gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
imagecol = ncol-postcol
imagerow = rowtrim[1]-rowtrim[0] + 1L
post_buffer = 4L
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biaspix = 5L
biascols = reverse(ncol-post_buffer-lindgen(biaspix))
goodrows = rowtrim[0]+ lindgen(imagerow)
biasrows = [rowtrim[0]+150L, rowtrim[0]+199L]
oscan = djs_median(image[biascols, biasrows[0]:biasrows[1]])
temp_image = (image[0:imagecol-1L, rowtrim[0]:rowtrim[1]] - oscan)*gain
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar)
bin = reverse(bin)

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deimos_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

bin = [1,1]

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

ncol = (size(image))[1]
nrow = (size(image))[2]

datasec = xregtovec(sxpar(hdr,'DATASEC'))
postcol = ncol-datasec[1]

rowtrim = [datasec[2]-1,datasec[3]-1]
biasrows = rowtrim
;gain = double(sxpar(hdr, 'GAIN'))
;rnoise = double(sxpar(hdr, 'RDNOISE'))
gain = 1.25
rnoise = 2.4

imagecol = datasec[1]-datasec[0]+1L  ; ncol-postcol
imagerow = rowtrim[1]-rowtrim[0] + 1L
post_buffer = 4L

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biaspix = (ncol - datasec[1] - 10)
biascols = reverse(ncol-post_buffer-lindgen(biaspix))
goodrows = rowtrim[0]+ lindgen(imagerow)

;; Get the overscan
oscan = djs_median(image[biascols, biasrows[0]:biasrows[1]],1)
;; Subtract
temp_image = (image[0:imagecol-1L, rowtrim[0]:rowtrim[1]] - $
              (replicate(1.,imagerow)#oscan))*gain
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;rawsub = transpose(rawsub)
;rawivar = transpose(rawivar)
;bin = reverse(bin)

RETURN
END


pro lris_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin
;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
     endif
    ;stop
    ;; Images are packed differently when Red side upgraded
    if strmid(sxpar(hdr,'DATE'),0,10) GT '2009-07-01' then begin
        image = x_readmhdufits(filelist[0], header=hdr, /notrim, $
                               /nobias)
        ;; Flip in y for blue side (back to the old way)
        if strcmp(strcompress(sxpar(hdr, 'INSTRUME'), /rem),'LRISBLUE') then $
          image = rotate(image,7) else begin
            ;; Trim top and bottom
            PRELINE = sxpar(hdr, 'PRELINE')
            POSTLINE = sxpar(hdr, 'POSTLINE')
            szi = size(image, /dimen)
            image = image[*,preline:szi[1]-preline-postline-1]
         endelse
        IF strcmp(strcompress(sxpar(hdr, 'INSTRUME'), /rem), 'LRIS') THEN BEGIN
           dims = size(image, /dim)
           nx = dims[0]
           ny = dims[1]
           IF nx MOD 2048 GT 0 OR ny MOD 2048 GT 0 THEN BEGIN
              ;; This data has been windowed, as a kludge, just use
              ;; the Keck Bias subtraction and an approximate read
              ;; noise. 
              if sxpar(hdr, 'NUMAMPS') NE 4 then stop
              namps = 4
              if strmid(sxpar(hdr, 'DATE'), 0, 10) GT '2010-12-03' then begin
                 ;; Second upgrade
                 gain = [1.25, 1.18, 1.19, 1.17]
                 rnoise_vec = [4.65, 4.7, 4.5, 4.65]
                 print, 'long_oscan: LRISr re-upgrade :: Using gain = ', gain
              endif else begin ;; LRIRs upgrade (the first)
                 gain = [0.955, 1.022, 0.916, 0.877]
                 rnoise_vec = [3.97, 4.57, 3.38, 3.39]
              endelse
              ;; JFH 12-31-2013, as a hack, just use the Keck routine
              ;; for procing for windowed data
              n = n_elements(gain)
              vidinp = ['VidInp1', 'VidInp2', 'VidInp3', 'VidInp4']
              foo = {GAINDATA, vidinp:'', gain:0.}
              gaindata = replicate({GAINDATA}, n)
              for i = 0, n-1 do $
                 gaindata[i] = {GAINDATA, vidinp:vidinp[i], gain:gain[i]}
              rawsub = x_readmhdufits(filelist[0], gaindata = gaindata)
              rnoise_mean = total(rnoise_vec)/4.0d
              rawivar = 1.0/(abs(rawsub - $
                                 sqrt(2.0)*rnoise_mean) +rnoise_mean^2)
              bin = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
              rnoise = rnoise_mean
              image = rawsub
              RETURN
           ENDIF
        ENDIF
     endif else $
           image = xmrdfits(filelist[0], 0, hdr, /fscale)
 ENDIF

ncol = (size(image))[1]
nrow = (size(image))[2]
bin = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
xbin = bin[0]
ybin = bin[1]
postpix = long(sxpar(hdr, 'POSTPIX'))
if (keyword_set(verbose)) then $
  splog, filelist[0], ': Binning is ', xbin, ' x', ybin
; are we red or blue?
instrument = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
IF strcmp(instrument, 'LRISBLUE') THEN BEGIN
;   prepix is in different headers for red versus blue (old red chip)
    prepix = long(sxpar(hdr, 'PRECOL'))
    ;; Sometimes only 2 amps for Longslit data
    namps = sxpar(hdr,'NUMAMPS')
    ;; JXP -- Likely bug in headers for TWOAMPTOP
    if namps EQ 2 then begin
        prepix = 100/xbin
        postpix = 160/xbin
    endif
    gain   = ([1.55,1.54471,1.63,1.65991])[4-namps:*]
    rnoise = ([3.9, 4.2, 3.6, 3.6])[4-namps:*]
ENDIF ELSE IF strcmp(instrument, 'LRIS') THEN BEGIN
    if strmid(sxpar(hdr,'DATE'),0,10) GT '2009-07-01' then begin
       if sxpar(hdr, 'NUMAMPS') NE 4 then stop
       namps = 4
       if strmid(sxpar(hdr,'DATE'),0,10) GT '2010-12-03' then begin
          ;; Second upgrade
          gain = [1.25, 1.18, 1.19, 1.17]
          rnoise = [4.65, 4.7, 4.5, 4.65]
          print, 'long_oscan: LRISr re-upgrade :: Using gain = ', gain
       endif else begin ;; LRIRs upgrade (the first)
          gain = [0.955, 1.022, 0.916, 0.877]
          rnoise = [3.97, 4.57, 3.38, 3.39]
       endelse
       prepix = long(sxpar(hdr, 'PRECOL'))
    endif else begin
        namps = 2
        gain = [1.98, 2.17]
        rnoise = [6.1, 6.3]
        prepix = long(sxpar(hdr, 'PREPIX'))
    endelse
ENDIF ELSE message, 'Cant figure out what instrument this is'

; These are small buffers to avoid using the overscan region to close to data
post_buffer1 = 4L
post_buffer2 = 8L
; The factor of namps is because of the two different amps

;;Begin - MF Jan 2011: Original line
; imagecol = ncol-namps*(prepix+postpix)/xbin
; Updated line:
IF (strcmp(instrument, 'LRIS')) and (strmid(sxpar(hdr, 'DATE'), 0, 10) GT '2010-12-10') THEN BEGIN
   IF xbin EQ 2 THEN imagecol = ncol-namps*(prepix+postpix-1)/xbin $
   ELSE IF xbin EQ 1 THEN imagecol = ncol-namps*(prepix+postpix)/xbin
;;
;; JFH 12-29-2013: the prepix + postpix = 87 for the new chip, so we need this -1
;; for xbin =2 
ENDIF ELSE imagecol = ncol-namps*(prepix+postpix)/xbin
;Note that this might be a problem even before 2010-12-10 for binned data.
;;End - MF Jan 2011
imagerow = nrow
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

datacol =  namps* (prepix/xbin) + lindgen(namps)*1024/xbin
postcol =  datacol[namps-1L] + (1024+post_buffer1)/xbin
FOR iamp = 0, namps-1L DO BEGIN
    biascols = lindgen((postpix-post_buffer2)/xbin) $
      + (iamp*postpix)/xbin + postcol
    oscan = djs_avsigclip(image[biascols, *], 1)
    osub = replicate(1, imagecol/namps) # oscan
    imagecols = lindgen(1024/xbin)+ iamp*1024/xbin
    temp_image = (image[imagecols + namps*(prepix/xbin), *] - osub)*gain[iamp]
    rawsub[imagecols, *] = temp_image
    rawivar[imagecols, *] = 1.0/(abs(temp_image - $
                                     sqrt(2.0)*rnoise[iamp]) +rnoise[iamp]^2)
 ENDFOR

IF strcmp(instrument, 'LRIS') THEN BEGIN 
   IF strmid(sxpar(hdr, 'DATE'), 10) LT '2009-07-01' then begin
      ;; bad pixel mask
      rawsub = transpose(rawsub)
      rawivar = transpose(rawivar)
      bin = reverse(bin)
      IF nrow EQ 2048 THEN BEGIN
         rawivar[0:722, 822] = 0.0
         rawivar[0:1482, 956] = 0.0
         rawivar[0:1485, 960] = 0.0
         rawivar[*, 2047] = 0.0
      ENDIF ELSE IF nrow EQ 1000 THEN BEGIN
         rawivar[0:224, 822] = 0.0
         rawivar[0:444, 956] = 0.0
      ENDIF ELSE message $
         , 'Problem in lris_oscan readout geometry or binning not supported'
   ENDIF ELSE IF strmid(sxpar(hdr, 'DATE'), 10) GT '2010-12-10' THEN BEGIN 
      rawivar[2044:2055, *] = 0.0
   ENDIF
ENDIF

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;
pro kast_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin


;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
endif

;mjd_newcontrol = 2454534L ;; New controller
mjd_newcontrol = 2454892L ;; New controlor (March, 2009)

;; Blue or Red?
if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
  strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then side = 1 else side = 2
;if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' then side = 1 else side = 2

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]


dateobs = sxpar(hdr, 'DATE-OBS')
mjd = x_setjdate(strmid(dateobs,0,10))

if mjd LT mjd_newcontrol then begin
    ;; Old controllers
    bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too

    if side EQ 1 then begin  ;; Blue side
        case strtrim(sxpar(hdr,'CLKSPEED'),2) of
            'Slow': begin
                dobs = sxpar(hdr,'DATE-OBS')
                jdate = x_setjdate(strmid(dobs,0,10))
                if jdate GE 2454526L then begin
                   print, 'long_oscan:  Assuming CCD#9 for the Kast blue CCD'
                   if jdate GE 2454721L then begin
                      gain = 1.2 ;; Dewar#9 on Sep 11, 2008 by E. Gates
                      rnoise = 3.7
                   endif else begin
                      gain = 2.3 ;; Dewar#9 on Feb 29, 2008 as measured by E. Gates on Mar 1, 2008
                      rnoise = 8
                   endelse
                endif else begin
                   gain = 3.8 ;; Original Kast CCD
                   rnoise = 6.
                endelse
             end
            else: stop
        endcase
    endif else begin  ;; Red side
        case strtrim(sxpar(hdr,'CLKSPEED'),2) of
            'Slow': begin
                gain = 3.0
                rnoise = 12.5
            end
            else: stop
        endcase
    endelse
    rawsub = image * gain
    rawivar = 1. / ( (image>3.)*gain + rnoise^2)
    ;; Trim
    if side EQ 1 then begin
        rawsub = rawsub[*,10:*]
        rawivar = rawivar[*,10:*]
    endif else begin
        rawsub = rawsub[*,0:144]
        rawivar = rawivar[*,0:144]
    endelse
endif else begin
  ;; New controller (installed -- )
    bin = [sxpar(hdr,'CBIN'), sxpar(hdr,'RBIN')]
    case round(sxpar(hdr,'READ-SPD')) of
       10: begin ;; Fast
          if side EQ 1 then begin
             gain = 1.3
             rnoise = 6.5
          endif else stop
       end
       80: begin ;; Slow
          if side EQ 1 then begin
             gain = 1.2
             rnoise = 4.
          endif else begin
             gain = 3.
             rnoise = 12.5
          endelse
       end
       ;;following the repair of the red side after Oct 24 2011
       ;;the new value for the read out slow on the red side is
       ;;160. The values for the RN and Gain are from 2008
       160: begin ;; Slow
          splog, "Check the gain/rn value!"
          if side EQ 1 then begin
             gain = 1.2
             rnoise = 4.
          endif else begin
             gain = 3.
             rnoise = 12.5
          endelse
       end
       else: stop
    endcase

    if side EQ 2 then begin ;; Red side (easier)
        ;; Subtract overscan
        ;RdS 2011  
        col_sci = sxpar(hdr, 'DNAXIS1') > (sxpar(hdr, 'NAXIS1')-sxpar(hdr, 'COVER'))
        ;col_sci = sxpar(hdr, 'DNAXIS1')
        col_origin = sxpar(hdr, 'CRVAL1')
        ncol = sxpar(hdr, 'NAXIS1')
        col_oscan = col_sci-col_origin  ;; This might not be right

        ;; Get the overscan
        oscan = djs_median(image[col_oscan:*, *], 1)
        temp_image = (image[0:col_oscan-1L, *] - $
                  (replicate(1.,ncol)#oscan))*gain
        rawsub = temp_image
        rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)
    endif else begin  ;; BLUE SIDE
        ncol_over = sxpar(hdr, 'COVER')
        namps = sxpar(hdr, 'AMPSROW')
        ncol = sxpar(hdr, 'NAXIS1')
        if namps NE 2 then stop

        ;RdS 2011 ;;MF 2012 - attempt to account for binning 2x2
        col_sci = sxpar(hdr, 'DNAXIS1')/bin[0] > (sxpar(hdr, 'NAXIS1')-2*sxpar(hdr, 'COVER'))
        ;col_sci = sxpar(hdr, 'DNAXIS1')
        sci1 = [0L, col_sci/2 - 1]
        sci2 = [sci1[1]+1,col_sci-1]
        bias1 = [col_sci+2,col_sci+ncol_over-2]
        bias2 = [col_sci+ncol_over+2,ncol-1]
        oscan1 = djs_median(image[bias1[0]:bias1[1], *], 1)
        oscan2 = djs_median(image[bias2[0]:bias2[1], *], 1)
        ;; Save
        temp_image = fltarr(col_sci,nrow)
        temp_image[sci1[0]:sci1[1],*] = (image[sci1[0]:sci1[1], *] - $
                  (replicate(1.,sci1[1]-sci1[0]+1)#oscan1))*gain
        temp_image[sci2[0]:sci2[1],*] = (image[sci2[0]:sci2[1], *] - $
                  (replicate(1.,sci2[1]-sci2[0]+1)#oscan2))*gain
        rawsub = temp_image
        rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)
    endelse
endelse

;; Spectra running up and down
rawsub = transpose(rawsub)
rawivar = transpose(rawivar)

RETURN

END

pro fire_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

  tmp = mrdfits(filename, 0, hh)

  hdr = hh
  gain = 1.2
  rnoise = 17.0/sqrt(4)
  bin = 1

  rawsub = (tmp)
  rawivar = 1.0/(abs(tmp - sqrt(2.0)*rnoise) +rnoise^2)
  rawivar = (rawivar)


END

pro fors2_oscan, filename, rawsub, rawivar, hdr = hdr $
                 , gain = gain, rnoise = rnoise, image = image $
                 , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale)
 endif

dims = size(image, /dim)
nx = dims[0]
ny = dims[1]
binx = esopar(hdr, 'HIERARCH ESO DET WIN1 BINX')
biny = esopar(hdr, 'HIERARCH ESO DET WIN1 BINY')
bin = [binx, biny]

;; CONAD = Conversion from ADUs to electrons. ESO uses different
;; convention!! The gain in the header is actually ADU/e^-
gain = esopar(hdr, 'HIERARCH ESO DET OUT1 CONAD')
rnoise = esopar(hdr, 'HIERARCH ESO DET OUT1 RON')
prscx = esopar(hdr, 'HIERARCH ESO DET OUT1 PRSCX')
ovscx = esopar(hdr, 'HIERARCH ESO DET OUT1 OVSCX')
ind_ovscx = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 OVSCX', /bool))
ovscx = (long(strmid(hdr[ind_ovscx], 30, 14)))[0]

prscy = esopar(hdr, 'HIERARCH ESO DET OUT1 PRSCY')
ovscy = esopar(hdr, 'HIERARCH ESO DET OUT1 OVSCY')

data_arr = [prscx+1L, prscy+1L, nx-ovscx, ny-ovscy]
bias_arr = [prscx+1L, ny-ovscy+1L, nx-ovscx, ny]
;datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
;biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
;data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
;bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[3]-bias_arr[1] + 1L
; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 1L
imagex = data_arr[2]-data_arr[0] + 1L
imagey = data_arr[3]-data_arr[1] + 1L

rawsub = fltarr(imagex, imagey)
rawivar = fltarr(imagex, imagey)

biascols = lindgen(nbias-oscan_buffer) + imagey + oscan_buffer + prscy[0]
oscan1 = djs_avsigclip(image[0:imagex-1L, biascols], 0)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = oscan # replicate (1, imagey)
temp_image = (image[prscx:prscx + imagex-1L $
                    , prscy: prscy + imagey -1L] - osub)*gain
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; Now transpose so that spectral direction vertical with blue down
rawsub = transpose(rawsub)
rawivar = transpose(rawivar)
bin = reverse(bin)

RETURN
END

pro mods_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

  ;; Red or Blue?
  chip = strtrim(strsplit(sxpar(hdr,'INSTRUME'),/extract),2)

;  no overscan for rows...
  if NOT keyword_set(image) then begin
;  look for filename
     filelist = lookforgzip(filename)
     if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
     endif
     image = xmrdfits(filelist[0], 0, hdr, /fscale)
  endif

  ;dims = size(image, /dim)
  ;nx = dims[0]
  ;ny = dims[1]

  ;; Translated from modsBias.py
  ind_naxis1 = WHERE(stregex(hdr, 'NAXIS1',/bool))
  naxis1 = long(strmid(hdr[ind_naxis1], 17, 14))
  naxis1 = naxis1[0]
  ind_naxis2 = WHERE(stregex(hdr, 'NAXIS2',/bool))
  naxis2 = long(strmid(hdr[ind_naxis2], 17, 14))
  naxis2 = naxis2[0]
  channel = strtrim(sxpar(hdr,'CHANNEL'),2)

  ind_binx = WHERE(stregex(hdr, 'CCDXBIN',/bool))
  binx = long(strmid(hdr[ind_binx], 17, 14))
  ind_biny = WHERE(stregex(hdr, 'CCDYBIN',/bool))
  biny = long(strmid(hdr[ind_biny], 17, 14))
  bin = [binx, biny]

  overscanx = 48L
  overscany = 0L

  c1 = overscanx 		; first image column counting from *zero*
  c2 = long(0.5*naxis1)-1	; last image column on first half
  c3 = c2+1			; first image column on second half
  c4 = naxis1-overscanx-1 	; last image column
  r1 = overscany 		; first image row
  r2 = long(0.5*naxis2)-1	; last image row on first half
  r3 = r2+1			; first image row on second half
  r4 = naxis2-overscany-1  	; last image row
  outnaxis1 = c4-c1+1		; columns in output, trimmed image
  outnaxis2 = r4-r1+1		; rows in output, trimmed image
  collen = long(0.5*outnaxis1)	; number of rows in an image quadrant
  rowlen = long(0.5*outnaxis2)	; number of rows in an image quadrant

  ; Identify the columns to use to calculate the bias level
  ; Skip the first and last columns of the overscan

  ; Identify the even and odd columns in the overscan

  cols_over_q1e = long(2.0 * findgen((overscanx-8.0)/2.0) + 4.0)
  cols_over_q2e = cols_over_q1e
  cols_over_q1o = long(2.0 * findgen((overscanx-8.0)/2.0) + 5.0)
  cols_over_q2o = cols_over_q1o
  cols_over_q3e = long(2.0 * findgen((overscanx-8.0)/2.0) + (naxis1-overscanx+4.0))
  cols_over_q4e = cols_over_q3e
  cols_over_q3o = long(2.0 * findgen((overscanx-8.0)/2.0) + (naxis1-overscanx+5.0))
  cols_over_q4o = cols_over_q3o

  ; Identify the even and odd columns in each quadrant

  cols_q1e = long(2.0 * findgen((c2-c1+2)/2.0) + c1)
  cols_q2e = cols_q1e
  cols_q1o = long(2.0 * findgen((c2-c1+2)/2.0) + c1 + 1.0)
  cols_q2o = cols_q1o
  cols_q3e = long(2.0 * findgen((c4-c3+2)/2.0) + c3)
  cols_q4e = cols_q3e
  cols_q3o = long(2.0 * findgen((c4-c3+2)/2.0) + c3 + 1.0)
  cols_q4o = cols_q3o

  ; Create arrays with the median overscan vs. row for each amplifier

  bias_q1e = fltarr(rowlen)
  bias_q1o = fltarr(rowlen)
  bias_q2e = fltarr(rowlen)
  bias_q2o = fltarr(rowlen)
  bias_q3e = fltarr(rowlen)
  bias_q3o = fltarr(rowlen)
  bias_q4e = fltarr(rowlen)
  bias_q4o = fltarr(rowlen)

  ; Calculate 1-D bias arrays for each amplifier

  for i=r1,r2 do begin
     bias_q1e[i] = median(image[cols_over_q1e,i])
     bias_q1o[i] = median(image[cols_over_q1o,i])
     bias_q2e[i] = median(image[cols_over_q2e,i+rowlen])
     bias_q2o[i] = median(image[cols_over_q2o,i+rowlen])
     bias_q3e[i] = median(image[cols_over_q3e,i])
     bias_q3o[i] = median(image[cols_over_q3o,i])
     bias_q4e[i] = median(image[cols_over_q4e,i+rowlen])
     bias_q4o[i] = median(image[cols_over_q4o,i+rowlen])
  endfor

  ; Subtract a single bias value for each amplifier

  bq1e = median(bias_q1e)
  bq1o = median(bias_q1o)
  bq2e = median(bias_q2e)
  bq2o = median(bias_q2o)
  bq3e = median(bias_q3e)
  bq3o = median(bias_q3o)
  bq4e = median(bias_q4e)
  bq4o = median(bias_q4o)

  temp_image = 0.0*image

  temp_image[cols_q1e,r1:r2] = image[cols_q1e,r1:r2] - bq1e
  temp_image[cols_q1o,r1:r2] = image[cols_q1o,r1:r2] - bq1o
  temp_image[cols_q2e,r3:r4] = image[cols_q2e,r3:r4] - bq2e
  temp_image[cols_q2o,r3:r4] = image[cols_q2o,r3:r4] - bq2o
  temp_image[cols_q3e,r1:r2] = image[cols_q3e,r1:r2] - bq3e
  temp_image[cols_q3o,r1:r2] = image[cols_q3o,r1:r2] - bq3o
  temp_image[cols_q4e,r3:r4] = image[cols_q4e,r3:r4] - bq4e
  temp_image[cols_q4o,r3:r4] = image[cols_q4o,r3:r4] - bq4o

  gain = 2.0  ;; from MODS ObsTools Preliminary Sensitivity estimate page
  rnoise = 2.5

  temp_rawsub = temp_image[c1:c4,r1:r4]


  gain_im = 0.0*temp_rawsub

  if chip eq 'MODS1R' then begin
     ;; The following values are good for Jun2012, from
     ;; mods1r20120623.0024  -- JXP 17 Aug 2012 :: This won't
     ;;                         work for the flats !!
     gain_im[0:4095,0:1543] = gain * 2.98/2.69
     gain_im[0:4095,1544:3087] = gain * 2.98/2.69 * (15022./13078) ;; JXP fix
     gain_im[4096:8191,0:1543] = gain
     gain_im[4096:8191,1544:3087] = gain ;* 3.21/3.03  ;; JXP fix
  endif else begin
     gain_im[0:4095,0:1543] = gain * 2.98/2.69
     gain_im[4096:8191,0:1543] = gain
     gain_im[4096:8191,1544:3087] = gain * 3.21/3.03
     gain_im[0:4095,1544:3087] = gain * 5.87/5.68
  endelse



  ;rawsub = gain * temp_image[c1:c4,r1:r4]
  rawsub = gain_im * temp_rawsub
  rawivar = 0.0 * rawsub
  rawivar = 1.0/(abs(rawsub - sqrt(2.0)*rnoise) +rnoise^2)

  ;stop

  ;; Now transpose so that spectral direction vertical with blue down
  if chip eq 'MODS1B' then begin
     rawsub = transpose(rawsub)
     rawivar = transpose(rawivar)
     bin = reverse(bin)
  endif else if chip eq 'MODS1R' then begin
     rawsub = reverse(transpose(rawsub),2)
     rawivar = reverse(transpose(rawivar),2)
     bin = reverse(bin)
  endif



RETURN
END


pro osiris_oscan, filename, rawsub, rawivar, hdr = hdr, $
                  gain=gain, rnoise=rnoise, image=image, $
                  verbose=verbose, bin=bin

  gain = sxpar(hdr, 'GAIN')
  ;; Gain/RDnoise combinations taken from
  ;; http://www.gtc.iac.es/en/pages/instrumentation/osiris.php#Detector_Setup
  if size(gain, /tname) eq 'STRING' then begin
     if gain eq 'GAIN_x4_75' then begin ;; guessing that this is what is indicated by 'GAIN_x4_75'
        gain = 1.18
        rnoise = 3.5
     endif else if gain eq 'GAIN_x9_5' then begin ;; this one too
        gain = 1.46
        rnoise = 8.0
     endif else begin
        message, 'Unknown gain/rnoise combination'
        return
     endelse
  endif else begin ;; keyword format changed sometime b/n 2011-05-09 and 2011-10-20
     if gain eq 0.95 then rnoise = 4.5 $
     else if gain eq 1.18 then rnoise = 3.5 $
     else if gain eq 1.46 then rnoise = 8.0 $
     else begin
        message, 'Unknown gain/rnoise combination'
        return
     endelse
  endelse

  ;; Read in left and right chip
  a1=mrdfits(filename, 1, hdr1, /silent)
  a2=mrdfits(filename, 2, hdr2, /silent)

  ;; Rescale
  a1=a1*sxpar(hdr1,'BSCALE')+sxpar(hdr1,'BZERO')
  a2=a2*sxpar(hdr2,'BSCALE')+sxpar(hdr2,'BZERO')

  ncols = (size(a1, /dim))[0]

  biassec1 = sxpar(hdr1, 'BIASSEC')
  biassec1 = long(strsplit(biassec1, '[*:*,*:*]', /extract))-1
  biassec2 = sxpar(hdr2, 'BIASSEC')
  biassec2 = long(strsplit(biassec2, '[*:*,*:*]', /extract))-1

  ;; Looked like a buffer was a good idea in the images...
  biasbuffer = 7

  oscan1 = djs_avsigclip(a1[biassec1[0]+biasbuffer:biassec1[1], *], 1)
  oscan2 = djs_avsigclip(a2[biassec2[0]+biasbuffer:biassec2[1], *], 1)

  osub1 = replicate(1,ncols) # oscan1
  osub2 = replicate(1,ncols) # oscan2

  a1 -= osub1
  a2 -= osub2

  ;; Trim to data section
  datasec1 = sxpar(hdr1, 'DATASEC')
  datasec1 = long(strsplit(datasec1, '[*:*,*:*]', /extract))-1
  datasec2 = sxpar(hdr2, 'DATASEC')
  datasec2 = long(strsplit(datasec2, '[*:*,*:*]', /extract))-1

  a1 = a1[datasec1[0]:datasec1[1], datasec1[2]:datasec1[3]]
  a2 = a2[datasec2[0]:datasec2[1], datasec2[2]:datasec2[3]]

  rawsub = transpose([[transpose(a1)],[transpose(a2)]]*gain)
  rawivar = 1.0/(abs(rawsub - sqrt(2.0)*rnoise)+rnoise^2)

  ;; Hard coded for now, should be determined from header in the
  ;; future
  bin = [2,2]

  return
end


;;Added ISIS/WHT - MF2015 
pro isis_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
       print, 'Could not find file named ', filename
       return
    endif
    hdr=headfits(filelist[0])
    image = xmrdfits(filelist[0], 1, hdr_img, /fscale)
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr_img, 'GAIN'))
rnoise = double(sxpar(hdr_img, 'READNOIS'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr_img, 'TRIMSEC'), /rem)
biassec = strcompress(sxpar(hdr_img, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))

;;copy values of interest in master header for later use
sxaddpar, hdr, 'GAIN', gain
sxaddpar, hdr, 'RON', rnoise

;These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

;Pick a median value for bias - can probably do better 
osub = djs_median(image[bias_arr[0]:bias_arr[1]-1L,bias_arr[2]+oscan_buffer:bias_arr[3]-1L])
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;;Blue/Red images are native with blue down
;;If strmatch(detector, '*blue*') THEN BEGIN
;;    rawsub = transpose(rawsub)
;;    rawivar = transpose(rawivar)
;;ENDIF ELSE IF strmatch(detector, '*red*') THEN BEGIN
;;    rawsub  = reverse(transpose(rawsub), 2)
;;    rawivar = reverse(transpose(rawivar), 2)
;;ENDIF ELSE message, 'Unrecognized detector'

RETURN

END


;;;;;;;;;;;;;;;;;;;;;;;
pro long_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin=bin, CCDONLY = CCDONLY $
                , TRANSFORM = TRANSFORM


hdr1 = xheadfits(filename)

if (size(hdr1, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    flux = 0
    invvar = 0
    return
endif

;  Is this Keck or Gemini?
telescope = strcompress(sxpar(hdr1, 'TELESCOP'), /rem)
instrument = strcompress(sxpar(hdr1, 'INSTRUME'), /rem)
detector = strcompress(sxpar(hdr1, 'DETECTOR'), /rem)
telid =  strcompress(sxpar(hdr1, 'TELID'), /rem)
;stop
IF strcmp(instrument, 'LRIS') OR $
   strcmp(instrument, 'LRISBLUE') THEN $
      lris_oscan, filename, rawsub, rawivar, hdr = hdr $
                  , gain = gain, rnoise = rnoise, image = image $
                  , verbose = verbose, bin = bin $
ELSE IF (strcmp(telescope, 'Gemini-North') OR $
         strcmp(instrument, 'GMOS-N') OR $
         strcmp(telescope, 'Gemini-South')) THEN BEGIN
   IF KEYWORD_SET(TRANSFORM) THEN $
      gmos_trnimg3to1, filename, rawsub, rawivar, hdr = hdr $
                    , gain = gain, rnoise = rnoise $
                    , verbose = verbose, bin = bin $
   ELSE gmos_oscan, filename, rawsub, rawivar, hdr = hdr, thdr = thdr $
                    , gain = gain, rnoise = rnoise, VERBOSE = VERBOSE $
                    , bin = bin, CCDONLY = CCDONLY
ENDIF ELSE IF strcmp(telescope, 'mmt') THEN $
   mmt_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
              , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(detector, 'mars') THEN $
   mars_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(telescope, 'kp4m') THEN $
   kpno_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, 'DEIMOS*') THEN $
  deimos_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF (stregex(instrument, 'KAST*',/bool,/fold_case) EQ 1) OR $
         (stregex(sxpar(hdr1,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then $
  kast_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, 'DIS*') THEN  $
  dis_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, '*ISIS*') THEN  $
  isis_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(telid, '200') THEN $
  p200_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, '*IMACS*') THEN $
  imacs_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF strmatch(instrument, 'FIRE') THEN $
  fire_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF  strcmp(telescope, 'CA-3.5') THEN $
   caha_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(telescope, 'CA-2.2') THEN $
   caha_22m_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                   , rnoise = rnoise, image = image, verbose = verbose $
                   , bin = bin $
ELSE IF  strcmp(telescope, 'DuPont') THEN $
   dupont_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                   , rnoise = rnoise, image = image, verbose = verbose $
                   , bin = bin $
ELSE IF strmatch(instrument, '*bcspeclamps*') THEN $ ; jm11jun08ucsd
   bok_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
              , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF strmatch(instrument, '*FORS2*') THEN $
   fors2_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, '*EFOSC*') THEN $
   ntt_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
              , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strcmp(telescope, 'LBT-SX') THEN $
   mods_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strcmp(instrument, 'OSIRIS') THEN $
   osiris_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                 , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, '*ISAAC*') THEN $
   isaac_proc, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE message, 'Not sure what instrument you want here'

RETURN
END
