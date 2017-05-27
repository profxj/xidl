function gsetup,sys

; see qmodel.x
; sets up one grating for qmodel

; ... below assumes phin=0. (ie adopts the roll/yaw approach)



thetan= - sys.MU
xsi= sys.GR_ZERR
rhon = sys.GR_YERR


cost = cos (thetan)
sint = sin (thetan)
cosx = cos (xsi)
sinx = sin (xsi)
cosr = cos (rhon)
sinr = sin (rhon)


a3=dblarr(3,3)

a3[0,0] =  cosx*cosr
a3[0,1] =  sint*sinr + cost*sinx*cosr
a3[0,2] = -cost*sinr + sint*sinx*cosr
a3[1,0] = -sinx
a3[1,1] =  cost*cosx
a3[1,2] =  sint*cosx
a3[2,0] =  cosx*sinr
a3[2,1] = -sint*cosr + cost*sinx*sinr
a3[2,2] =  cost*cosr + sint*sinx*sinr

return,a3
end





