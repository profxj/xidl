FUNCTION hires_murphy_thar, wave_thar= wave_thar

  ;; Read in Murphy spectrum
  thar_file=getenv('XIDL_DIR') + $
            '/Keck/HIRES/Redux/pro/Arcs/ThAr/thar_spec_MM201006.fits'
  obj=mrdfits(thar_file,0,hdr)
  thar=obj[*,0]
  sz=size(thar,/dim)
  nspec=sz[0]
  wave_thar=x_setwave(hdr,nspec)
  
  RETURN, thar
END
