pro MW_dust, Rv, wav, Awav, wvmin=wvmin, wvmax=wvmax, wvstep=wvstep

; MW_dust -- Calculates A(lambda)/A(V) for MW extinction using
;    the parameterization given in Cardelli, Clayon, & Mathis (1989)

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'MW_dust, Rv, wav, Awav, [WVMIN= , WVMAX= , WVSTEP= ]'
    return
  endif 

  if not keyword_set( WVMIN ) then    wvmin    = 1000.
  if not keyword_set( WVMAX ) then    wvmax    = 10000.
  if not keyword_set( WVSTEP ) then    wvstep    = 50.


  if(wvmin LT 1000) then begin
      print, 'CCM89 does not extend below 1000A'
      return
  endif

; Loop on wavelengths

  nstep = long ( (wvmax-wvmin) / wvstep ) + 1
  wav = fltarr(nstep)
  Awav = fltarr(nstep)
  for i=0,nstep-1 do begin
      wav[i] = wvmin + i * wvstep
      x = 10000/wav[i]

; Far-UV
      if(x GE 8) then begin
          ax = -1.073 - 0.628*(x-8) + 0.137*(x-8)^2 - 0.070*(x-8)^3
          bx = 13.670 + 4.257*(x-8) - 0.420*(x-8)^2 + 0.374*(x-8)^3
      endif
; UV
      if(x GE 3.3 AND x LT 8) then begin
          if(x GT 5.9) then begin
              Fax = -0.04473*(x-5.9)^2 - 0.009779*(x-5.9)^3
              Fbx = 0.2130*(x-5.9)^2 + 0.1207*(x-5.9)^3
          endif else begin
              Fax = 0.0
              Fbx = 0.0
          endelse
          ax = 1.752 - 0.316*x - 0.104/ ( (x-4.67)^2 + 0.341) + Fax
          bx = -3.090 + 1.825*x + 1.206 / ( (x-4.62)^2 + 0.263 ) + Fbx
      endif

; Optical

      if( x LT 3.3 AND x GE 1.1) then begin
          y = x - 1.82
          ax = 1 + 0.17699*y - 0.50447*y^2 - 0.02427*y^3 + 0.72085*y^4 $
            + 0.01979 * y^5 - 0.77530*y^6 + 0.32999*y^7
          bx = 1.41338*y + 2.28305*y^2 + 1.07233*y^3 - 5.38434*y^4 $
            - 0.62251*y^5 + 5.30260*y^6 - 2.09002*y^7
      endif

; IR

      if (x LT 1.1) then begin
          ax = 0.574 * x^(1.61)
          bx = -0.527 * x^(1.61)
      endif

      Awav[i] = ax + bx / Rv
  endfor

return
end
