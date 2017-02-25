FUNCTION ISAAC_IRWAVE_MASK, wave

  ;; Define clean lineless regions for fitting continuum
  fit_wv = [ [0.8137525, 0.8244728],  $
             [0.8621566, 0.8739599], $
             [0.9137012, 0.9278867], $
             [1.060408, 1.063224], $  
             [1.079467, 1.082174], $
             [1.188598, 1.192172], $
             [1.207332, 1.210255], $
             [1.495439, 1.498038], $
             [1.535289, 1.537021], $
             [1.756367, 1.761998], $
             [1.862055, 1.867686], $
             [2.076419, 2.082700], $
             [2.136626, 2.146156], $ 
             [2.235341, 2.242704], $   
             [2.256565, 2.263712], $ 
             [2.279305, 2.288401], $    
             [2.332366, 2.343627], $    
             [2.369833, 2.376170]]  
  szf = size(fit_wv, /dimen)
  gdi = [0]
  for ii = 0L, szf[1]-1 do begin
     morei = where( wave GT fit_wv[0, ii] AND wave LT fit_wv[1, ii], nmore)
     if nmore GT 0 then gdi = [gdi, morei]
  endfor
  
  RETURN, gdi
END
