; wrapper for spherematch using pixels instead of ra and dec
PRO long_pixelmatch, x1t, y1t, x2t, y2t, matchlength, match1, match2 $
                , distance12, maxmatch = maxmatch
  
p1t = 0.01d
p2t = 0.01d
; fake conversion to RA and DEC
scale1 = p1t/3600.D
ra1 = x1t*scale1
dec1 = y1t*scale1
scale2 = p2t/3600.0D
ra2 = x2t*scale2
dec2 = y2t*scale2

spherematch, ra1, dec1, ra2, dec2, (matchlength*scale1) $
             , match1, match2, distance12_temp, maxmatch = maxmatch
distance12 = distance12_temp/scale1

RETURN
END
