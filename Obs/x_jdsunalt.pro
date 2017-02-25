function x_jdsunalt_altit, dec, ha, lat

  DEG_IN_RADIAN = 180.d/!dpi
  HRS_IN_RADIAN = DEG_IN_RADIAN/15.d

  rdec = dec / DEG_IN_RADIAN       ;
  rha = ha / HRS_IN_RADIAN       ;
  rlat = lat / DEG_IN_RADIAN     ;  /* thank heavens for pass-by-value */
  x = DEG_IN_RADIAN * asin(cos(rdec)*cos(rha)*cos(rlat) + sin(rdec)*sin(rlat)) ;

  return, x
end

function x_jdsunalt, alt, jdguess, lat, longit

  del = 0.002                     ;
  i = 0                   ;

;	/* first guess */
  x_sunradec, jdguess, ra, dec       ;
  ha = x_getlst(jdguess,longit) - ra ;
  alt2 = x_jdsunalt_altit(dec,ha,lat)      ;

  jdguess = jdguess + del     
  x_sunradec, jdguess, ra, dec 
  alt3 = x_jdsunalt_altit(dec,(x_getlst(jdguess,longit) - ra),lat) ;
  err = alt3 - alt              ;
  deriv = (alt3 - alt2) / del   ;

  ;; Loop
  while( (abs(err) > 0.1) AND (i < 10)) do begin
      jdguess = jdguess - err/deriv 
      x_sunradec, jdguess, ra, dec   ;
      alt3 = x_jdsunalt_altit(dec,(lst(jdguess,longit) - ra),lat) ;
      err = alt3 - alt          ;
      i = i+1
      if i EQ 9 then return, -1.0e10 ; /* bad status flag */
  endwhile
  jdout = jdguess               ;
  return, jdout                 ;

end
