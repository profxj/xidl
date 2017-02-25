PRO GMOS_GAPS, nx, ny, YGAP1 = YGAP1, YGAP2 = YGAP2

spatbin = long(4L*1152/nx)
specbin = long(2L*1024/ny)


IF specbin EQ 4 THEN BEGIN
   ygap1 = 9L
   ygap2 = 9L
ENDIF ELSE IF specbin EQ 2 THEN BEGIN
   ygap1 = 18L
   ygap2 = 18L
ENDIF ELSE IF specbin EQ 1 THEN BEGIN
   ygap1 = 36L
   ygap2 = 36L
ENDIF ELSE message, 'gaps not supported for your binning'

return
end

       
