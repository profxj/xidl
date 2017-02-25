FUNCTION GMOS_MOSAIC, imag2, VERBOSE = VERBOSE, ygap1 = ygap1, ygap2 = ygap2

; transpose so that spectral direction is row direction
imag1 = transpose(imag2, [1, 0, 2])
; already been transposed by gemini_oscan
rawspat = (size(imag1[*, *, 0]))[1]
rawspec = (size(imag1[*, *, 0]))[2]

spatbin = long(4L*1152/rawspat)
specbin = long(2L*1024/rawspec)

IF specbin EQ 4 THEN BEGIN
    ygap1 = 9L
    ygap2 = 9L
ENDIF ELSE IF specbin EQ 2 THEN BEGIN
    ygap1 = 18L
    ygap2 = 18L
ENDIF ELSE IF specbin EQ 1 THEN BEGIN
    ygap1 = 36L
    ygap2 = 36L
ENDIF ELSE message, 'Gmosaic not supported for your binning'

nspec = 3*rawspec + ygap1 + ygap2
nspat = rawspat

imag = dblarr(nspat, nspec)

;IF KEYWORD_SET(hdr) THEN grating = sxpar(hdr[*, 0], 'GRATING') $
;ELSE grating = ''
;IF grating EQ 'R150+_G5306' THEN BEGIN
;    imag[*, 0:rawspec-1L] = shift(reverse(imag1[*, *, 2], 2), [1, 0])
;    imag[*, rawspec + ygap1:2*rawspec + ygap1-1L] = $
;      shift(reverse(imag1[*, *, 1], 2), [1, 0])
;    imag[*, 2*rawspec+ygap1+ygap2:*] = reverse(imag1[*, *, 0], 2)
;ENDIF ELSE BEGIN
imag[*, 0:rawspec-1L] = reverse(imag1[*, *, 2], 2)
imag[*, rawspec + ygap1:2*rawspec + ygap1-1L] = reverse(imag1[*, *, 1], 2)
imag[*, 2*rawspec+ygap1+ygap2:*] = reverse(imag1[*, *, 0], 2)
;ENDELSE

return, imag
end

       
