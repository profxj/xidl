pro proj_to_mask,xp,yp,ap,xc,yc,ac

;see deimos_util.x

; inputs: xp,yp,ap, outputs: xc,yc,ac




; PROJ_TO_MASK: project planar onto curved slitmask coordinates
; Double inputs, outputs
; Note that this is pure geometry -- any empirically determined corrections
; should go elsewhere...

;double	xp, yp			# x,y in (focal)-plane system
;double	ap			# position angle in planar system (deg)
;double	xc, yc			# returned x,y values on mask surface
;double	ac			# returned angle on curved surface (deg)

;double	mu, cosm			# mu, cos(mu)
;double	cost, tant			# cos, tan of mask tilt angle
;double	tanpa				# tan PA

;double	rho			# radius from telescope optical axis
;double	hs, hm			# height of image surface, mask above datum
;double	xx, yy		# Work variables corresponding to xc, yc

xp=double(xp)
yp=double(yp)
ap=double(ap)

mu = asin (xp / !M_RCURV)
cosm = cos(mu)
cost = cos(!dtor*(!M_ANGLE))
tant = tan(!dtor*(!M_ANGLE))
xx =  !M_RCURV * mu
yy =  (yp - !ZPT_YM) / cost + !M_RCURV * tant * (1. - cosm)

tanpa = tan(!dtor*(ap)) * cosm / cost + tant * xp / !M_RCURV
ac = !radeg*(atan(tanpa))

; What follows is a small correction for the fact that the mask does
; not lie exactly in the spherical image surface (where the distortion
; values and gnomonic projection are calculated) and the rays are arriving
; from the pupil image; thus, the exact locations are moved _slightly_
; wrt the telescope optical axis.  Note also that these corrections are
; only calculated to first order.

; Spherical image surface height:
	rho = sqrt (xp * xp + yp * yp)
	hs = !R_IMSURF * (1. - sqrt (1. - (rho / !R_IMSURF)^2))
; Mask surface height:
	hm = !MASK_HT0 + yy * sin (!dtor*(!M_ANGLE)) + !M_RCURV * (1. - cosm)
; Correction:
	yc = yy + (hs - hm) * yp / !PPLDIST / cost
	xc = xx + (hs - hm) * xp / !PPLDIST / cosm

return
end
