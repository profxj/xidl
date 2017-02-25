FUNCTION GMOS_UNMOSAIC, imag2, VERBOSE = VERBOSE

; transpose so that spectral direction is x direction
imag1 = transpose(imag2)
; already been transposed by gemini_oscan
rawspat = (size(imag1[*, *, 0]))[2]
rawspec = (size(imag1[*, *, 0]))[1]

IF       (rawspec-2L*9L) MOD  512 EQ 0 THEN specbin = 4 $
ELSE IF (rawspec-2L*18L) MOD  512 EQ 0 THEN specbin = 2 $
ELSE IF (rawspec-2L*36L) MOD  512 EQ 0 THEN specbin = 1 $
ELSE message, 'gunmosaic not supported for your binning'

spatbin = long(4L*1152/rawspat)

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

nspat = rawspat
nspec = (rawspec-ygap1-ygap2)/3L

imag = dblarr(nspec, nspat, 3)

imag[*, *, 2] = reverse(imag1[0:nspec-1L, *], 1)
imag[*, *, 1] = reverse(imag1[nspec+ygap1:2*nspec+ygap1-1L, *], 1)
imag[*, *, 0] = reverse(imag1[2*nspec+ygap1+ygap2:*, *], 1)

;imag[*, 0:rawspec-1L] = reverse(imag1[*, *, 2], 2)
;imag[*, rawspec + ygap1:2*rawspec + ygap1-1L] = reverse(imag1[*, *, 1], 2)
;imag[*, 2*rawspec+ygap1+ygap2:*] = reverse(imag1[*, *, 0], 2)

return, imag
end

       
