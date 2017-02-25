PRO flux,mike
	for i=1L,13 do begin
		mike_flux,mike,1,i,1,fluxfil='Extract/Sens_mb0020.idl',/clobber
		mike_combspec,mike,1,i,1
		mike_1dspec,mike,1,i,1
		mike_flux,mike,1,i,2,fluxfil='Extract/Sens_mr0020.idl',/clobber
		mike_combspec,mike,1,i,2
		mike_1dspec,mike,1,i,2
	endfor
return
END
