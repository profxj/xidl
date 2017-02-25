FUNCTION GMOS_CASCADE, imag1, VERBOSE = VERBOSE, YGAP1 = YGAP1, YGAP2 = YGAP2

; already been transposed by gemini_oscan
dims = size(imag1, /dim)
rawspat = dims[0]
rawspec = dims[1]

spatbin = long(4L*1152/rawspat)
specbin = long(2L*1024/rawspec)

IF specbin EQ 2 THEN BEGIN
    ygap1 = 18L
    ygap2 = 18L
ENDIF ELSE IF specbin EQ 1 THEN BEGIN
    ygap1 = 36L
    ygap2 = 36L
ENDIF ELSE message, 'Gmosaic not supported for your binning'

nspec = 3*rawspec + ygap1 + ygap2
nspat = rawspat

imag = dblarr(nspat, nspec)

imag[*, 0:rawspec-1L] = imag1[*, *, 2]
imag[*, rawspec + ygap1:2*rawspec + ygap1-1L] = imag1[*, *, 1]
imag[*, 2*rawspec+ygap1+ygap2:*] = imag1[*, *, 0]

return, imag
end

       
