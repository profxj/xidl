PRO NOT_BIASSUB, filename, outimag, outivar, superbias = superbias $
                 , delta_bias = delta_bias $
                 , instbiasfile = instbiasfile, superdark = superdark $
                 , ihdr = hdr, gain = gain, rnoise = rnoise, hdr_arr = hdr_arr $
                 , mask = mask

  image = not_mosaic(filename, ihdr = hdr, hdr_arr = hdr_arr, mask = mask)
  ndim = size(image, /n_dimen)
  dims = size(image, /dimens)
  outimag = fltarr(dims[0], dims[1])
  outivar = fltarr(dims[0], dims[1])

  IF NOT KEYWORD_SET(SUPERBIAS) OR NOT KEYWORD_SET(DELTA_BIAS) THEN $
     message, 'Currently superbias required'
  if (size(superbias, /n_dimen) NE ndim $
      OR total(size(superbias, /dimens) NE dims) NE 0) then $
         message, 'Dimensions of image and superbias image do not agree'
  IF KEYWORD_SET(instbiasfile) THEN $
     instbias = not_mosaic(instbiasfile, ihdr = ibiashdr $
                           , hdr_arr = ibiashdr_arr)
  IF KEYWORD_SET(instbias) THEN BEGIN
     if (size(instbias, /n_dimen) NE ndim $
         OR total(size(instbias, /dimens) NE dims) NE 0) then $
            message, 'Dimensions of image and bias image do not agree'
  endif else instbias = superbias

  if (keyword_set(superdark)) then begin
     if (size(superdark, /n_dimen) NE ndim $
         OR total(size(superdark, /dimens) NE dims) NE 0) then $
            message, 'Dimensions of image and bias image do not agree'
  endif else superdark = superbias

  darkminbias = superdark - superbias
  FOR chip = 1L, 4L DO BEGIN
     gain = double(sxpar(hdr_arr[*, chip-1], 'GAIN'))
     rnoise = double(sxpar(hdr_arr[*, chip-1], 'RDNOISE'))
     ipix = WHERE(mask EQ chip, npix)
     inst_level = median(instbias[ipix])
     dark_current = median(darkminbias[ipix])
     biasnow = inst_level + dark_current + delta_bias[ipix]
     outimag[ipix] = (image[ipix] - biasnow)*gain
     outivar[ipix] = 1.0/(abs(outimag[ipix] - sqrt(2.0)*rnoise) + rnoise^2)
  ENDFOR
  outimag = outimag*float(mask GT 0)
  outivar = outivar*float(mask GT 0)

  RETURN
END
