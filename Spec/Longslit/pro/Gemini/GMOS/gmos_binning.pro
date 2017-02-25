PRO GMOS_BINNING,  nx, ny, SPECBIN=SPECBIN,SPATBIN=SPATBIN,BTAG=BTAG
  
  IF      (ny-2L*9L)  MOD  512L EQ 0 THEN specbin = 4 $
  ELSE IF (ny-2L*18L) MOD  512L EQ 0 THEN specbin = 2 $
  ELSE IF (ny-2L*36L) MOD  512L EQ 0 THEN specbin = 1 $
  ELSE message, 'Transformation not supported for your binning'
  spatbin = long(4L*1152L/nx)

  ;; Binning is spec x spatial (i.e. X x Y in Gemini convention)
  IF specbin EQ 2 AND spatbin EQ 2 THEN BTAG = '2x2' $
  ELSE IF specbin EQ 1 AND spatbin EQ 1 THEN BTAG = '1x1' $
  ELSE IF specbin EQ 2 AND spatbin EQ 1 THEN BTAG = '2x1' $
  ELSE message, 'Your binning is not supported'

  
  RETURN
  END
