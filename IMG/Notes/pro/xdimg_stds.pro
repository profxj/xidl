pro lris_13jan99_stds, struct

  STDS = where(struct.type EQ 'STD' AND struct.flg_anly NE 0 )
	
  xdimg_proc, struct, STDS, 'Flats/SkyFltN', /DELOV

  ; Find stars in the image
;  xdimg_starid, struct

  ; Do Aperture photometry
;  xdimg_stdmag, struct

  ; Do photometric calibration
;  xdimg_photcal, struct, XSIZE=1000, YSIZE=800

end
  
