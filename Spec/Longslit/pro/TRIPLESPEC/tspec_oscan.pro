pro tspec_oscan, filename, rawsub, rawivar, hdr = hdr $
                 , gain = gain, rnoise = rnoise, image = image $
                 , verbose = verbose 
  
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

gain = double((strsplit(sxpar(hdr, 'GAIN'), /extract))[0])

;; This assumes FOWLER sampling
dnsamp = double((strsplit(sxpar(hdr, 'FSAMPLE'), /extract))[0])
rnoise = 14.0D/sqrt(dnsamp)
;rnoise = 7.; hack
;biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
;bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
;nbias = bias_arr[3]-bias_arr[2] + 1L
;nbias = 10.; hack, this is not really used.

; This is a small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = ncol
imagerow = nrow

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

;biascols = lindgen(nbias-oscan_buffer) + imagerow + oscan_buffer
;oscan = median(image[0:imagecol-1L, biascols])
oscan = 0. ; The H2RG images are delivered with oscan subtraction already done
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - oscan)*gain 
rawsub[0:imagecol-1L, *] = temp_image

rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)
rawsub = transpose(reverse(rawsub, 1))
rawivar = transpose(reverse(rawivar, 1)) 
RETURN
END
