FUNCTION gmos_alignmask,flat,xstart,xend,silent=silent
  
  nslit=n_elements(xstart)
  dim=size(flat,/dim)
  nx = dim[0]
  ny = dim[1]
  slitmask = lonarr(nx,ny)
  xarr = findgen(nx) # replicate(1.0,ny)
  FOR islit=0L,nslit-1L DO BEGIN
     ipix=where(xarr GE xstart[islit] AND xarr LE xend[islit])
     slitmask[ipix]=1L
  ENDFOR
  if (NOT keyword_set(maxlag)) then maxlag = 15L
  lags = lindgen(2L * maxlag + 1L) - 15L
  flat = (flat < 1.0d4) > 0.0
  corr = c2_correlate(slitmask, flat < 1.0d4, lags)
  xshift = long_find_nminima(-corr, lags, nfind = 1, width = 4, minsep = 2)
  IF NOT (keyword_set(silent)) THEN BEGIN
     splog, 'Correlation gives shift of ', xshift
     splog, 'Shifting xstart and xend positions....'
  ENDIF 
     
  RETURN,xshift
END





