function ccd_geom,sys
; see trace.x

; CCDGEOM: specify the geometry of the CCDs in the mosaic
; Probably want to add FCS devices also (note that FDS CCDs are rotated 90deg)
; Order is x(pix), y(pix), theta(deg)


NCCD=!NCCD

a=dblarr(NCCD,3)


xdimeff = !CCDXPIX + 2.*!CCDXEDG/!PIX_SZ + !NOMXGAP/!PIX_SZ	; 2135.200
ydimeff = !CCDYPIX + 2.*!CCDYEDG/!PIX_SZ + !NOMYGAP/!PIX_SZ	; 4112.000

; #	coeff  = pix-off + nom.gap	+ adjustment

a[0,0] = -1.5 * xdimeff		- 20.05
a[0,1] = -0.5 * ydimeff		+ 14.12
a[0,2] = 0.			+ !dtor*(-0.082)

a[1,0] = -0.5 * xdimeff		- 12.64
a[1,1] = -0.5 * ydimeff		+ 7.25
a[1,2] = 0.			+ !dtor*(0.030)

a[2,0] =  0.5 * xdimeff
a[2,1] = -0.5 * ydimeff
a[2,2] = 0.

a[3,0] =  1.5 * xdimeff		- 1.34
a[3,1] = -0.5 * ydimeff		- 19.92
a[3,2] = 0.			+ !dtor*(-0.1206)

a[4,0] = -1.5 * xdimeff		- 19.02
a[4,1] =  0.5 * ydimeff		+ 16.46
a[4,2] = 0.			+ !dtor*(0.136)

a[5,0] = -0.5 * xdimeff		- 9.65
a[5,1] =  0.5 * ydimeff		+ 8.95
a[5,2] = 0.			+ !dtor*(-0.06)

a[6,0] =  0.5 * xdimeff		+ 1.88
a[6,1] =  0.5 * ydimeff		+ 1.02
a[6,2] = 0.			+ !dtor*(-0.019)

a[7,0] =  1.5 * xdimeff		+ 4.81
a[7,1] =  0.5 * ydimeff		- 24.01
a[7,2] = 0.			+ !dtor*(-0.082)

a[*,0] = a[*,0] + sys.CN_XERR
a[*,1] = a[*,1] + sys.CN_YERR
a[*,2] = a[*,2] + sys.CN_RERR

return,a
end
