pro mask_to_proj,xc,yc,ac,xp,yp,ap
; see deimos_util.x
; can vectorize


; MASK_TO_PROJ: project slitmask coords (curved surface) onto plane
; Double inputs, outputs
; Note that this is pure geometry -- any empirically determined corrections
; should go elsewhere...

;double	xc, yc			# x,y values on mask surface
;double	ac			# position angle on curved surface
;double	xp, yp			# returned x,y in (focal)-plane system
;double	ap			# returned position angle in planar system

;double	mu, cosm			# mu, cos (mu)
;double	cost, tant			# cos, tan of mask tilt angle
;double	tanpa				# tan PA

;double	rho			# radius from telescope optical axis
;double	hs, hm			# height of image surface, mask above datum
;double	xx, yy		# Work variables corresponding to xp, yp

mu = xc / !M_RCURV
cosm = cos (mu)
cost = cos (!dtor*(!M_ANGLE))
tant = tan (!dtor*(!M_ANGLE))
xx =  !M_RCURV * sin (mu)
yy =  (yc - !M_RCURV * tant * (1. - cosm)) * cost + !ZPT_YM
tanpa = (tan (!dtor*(ac)) - tant * xx / !M_RCURV) * cost / cosm
ap = atan (tanpa)

; What follows is a small correction for the fact that the mask does
; not lie exactly in the spherical image surface (where the distortion
; values and gnomonic projection are calculated) and the rays are arriving
; from the pupil image; thus, the exact locations are moved _slightly_
; wrt the telescope optical axis.  Note also that these corrections are
; only calculated to first order.

; Spherical image surface height:

rho = sqrt (xx * xx + yy * yy)
hs = !R_IMSURF * (1. - sqrt(1. - (rho / !R_IMSURF)^2))

; Mask surface height:
hm = !MASK_HT0 + yc * sin(!dtor*(!M_ANGLE)) + !M_RCURV*(1. - cosm)
yp = yy - (hs - hm) * yy / !PPLDIST 
xp = xx - (hs - hm) * xx / !PPLDIST 


return
end
