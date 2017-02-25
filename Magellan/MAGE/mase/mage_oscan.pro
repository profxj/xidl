pro mage_oscan, filename, rawsub, rawivar, hdr = hdr $
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

ncol = (size(image))[1]
nrow = (size(image))[2]
bin = long(strsplit(sxpar(hdr, 'BINNING'), 'x', /extract))
xbin = bin[0]
ybin = bin[1]
if (keyword_set(verbose)) then $
  splog, filelist[0], ': Binning is ', xbin, ' x', ybin
gain = double(sxpar(hdr, 'EGAIN'))
rnoise = double(sxpar(hdr, 'ENOISE'))
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[3]-bias_arr[2] + 1L

; This is a small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = 2048L
imagerow = 1024L

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagerow + oscan_buffer
oscan = median(image[0:imagecol-1L, biascols])
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - oscan)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(reverse(rawsub, 1))
rawivar = transpose(reverse(rawivar, 1)) 

RETURN
END
