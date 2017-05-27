pro cam_distort,thetain,thetaout,apply=apply
; vectorized already

; see deimos_util.x

; CAM_DISTORT: remove camera distortion from distorted angles.
; Currently DEIMOS theoretical distortion curves apply. !!!
; Double inputs, outputs; bool apply

if n_elements(apply) eq 0 then apply=0.


MAX_TAN = 0.07852	; tan (angle) at one side

CAMD_C0		= 1.D0
CAMD_C2		= 0.0457563
CAMD_C4		= -0.3088123
CAMD_C6		= -14.917

if apply then thetaout=thetain* $
	(CAMD_C0 + thetain*thetain * $
	(CAMD_C2 + thetain*thetain * $
	(CAMD_C4 + thetain*thetain * CAMD_C6))) $ 
else thetaout=thetain/ $
			(CAMD_C0 + thetain*thetain * $
			(CAMD_C2 + thetain*thetain * $
			(CAMD_C4 + thetain*thetain * CAMD_C6)))

return
end
