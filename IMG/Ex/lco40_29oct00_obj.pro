pro lco40_29oct00_obj, struct

  OBJR = where(struct.type EQ 'OBJ' AND struct.flg_anly NE 0 $
	AND struct.filter EQ 'R')
  xdimg_proc, struct, OBJR, 'Flats/SkyFltN'

  ; Didn't have a good B SkyFlat using Twilight

  OBJB = where(struct.type EQ 'OBJ' AND struct.flg_anly NE 0 $
	AND struct.filter EQ 'B')

  xdimg_proc, struct, OBJB, 'Flats/TwiFltN'

end
  .com 
