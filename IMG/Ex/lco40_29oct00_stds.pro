pro lco40_29oct00_stds, struct

  STDSR = where(struct.type EQ 'STD' AND struct.flg_anly NE 0 $
	AND struct.filter EQ 'R')
  xdimg_proc, struct, STDSR, 'Flats/SkyFltN', /DELOV

  ; Didn't have a good B SkyFlat using Twilight

  STDSB = where(struct.type EQ 'STD' AND struct.flg_anly NE 0 $
	AND struct.filter EQ 'B')

  xdimg_proc, struct, STDSB, 'Flats/TwiFltN', /DELOV

  ; Find stars in the image
  xdimg_starid, struct

  ; Do Aperture photometry
  xdimg_stdmag, struct

  ; Do photometric calibration
  xdimg_photcal, struct, XSIZE=1000, YSIZE=800

end
  
