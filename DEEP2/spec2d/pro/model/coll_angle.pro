pro coll_angle,xp, yp, sys, tc, pc
; see refl.x
; can vectorize

; NOTE that returned tc, pc are the downward normal

defsysv,'!R_CURV',2124.71D0
defsysv,'!D_1',20018.4D0
defsysv,'!K_COLL',-0.75D0
defsysv,'!R_COLL',4394.2D0
   
xp=double(xp)
yp=double(yp)

;double	xp, yp			# projected x,y in plane of telescope
;double	sys[NPARAM]		# system parameters
;double	tc, pc			# returned theta, phi of collimator surface

;double	rp, rc			# radius projected, collimator
;double	hm			# height at mask
;double	cott			# cotangent theta
;double	d2			# pupil to collimator dist
;double	d, k			# work factors

pc = atan(yp, xp)		;phi will be same on coll
rp = sqrt (xp*xp + yp*yp)
hm = !R_CURV - sqrt(!R_CURV*!R_CURV - rp*rp)
d2 = !D_1 + sys.COL_DST

cott = (!D_1 - hm) / rp
k = 1. + (1. + !K_COLL) * cott*cott
d = d2 * (1. + !K_COLL)

;The following is general for conic sections.  In practice, a parabola is fine
; Note the switch in sign for quadratic root
	if (!R_COLL - d) gt 0. then begin
		rc = (!R_COLL - d) / k * $
(sqrt(cott*cott + d2*k*(2.*!R_COLL - d) /(!R_COLL - d)^2) - cott)
	endif else begin
		rc = (d - !R_COLL) / k * $
(sqrt (cott*cott + d2*k*(2.*!R_COLL - d)/(!R_COLL - d)^2) + cott)
	endelse

; This is for parabola:
;	rc = !R_COLL * (sqrt (cott*cott + 2. * d2 / !R_COLL ) - cott)
	
; The general conic form (important)
	rc = rc / !R_COLL
	tc = atan(rc / sqrt (1. - (1.+!K_COLL)*rc*rc))

return
end


