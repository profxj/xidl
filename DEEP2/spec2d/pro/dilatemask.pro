function dilatemask,mask,radius
;+
; NAME:
;   dilatemask
;
; PURPOSE:
;   Given a mask (1 where bad, 0 where good) expand bad regions
;
; CALLING SEQUENCE:
;   newmask=dilatemask( mask, dilation radius)
; 
; INPUTS:
;   mask  - 2d array containing the mask to be expanded
;   radius - radius (in pixels) to expand by (default 1)
;
; OUTPUTS:
;   newmask - expanded version of mask
;
; MODIFICATION HISTORY:
;    15-Aug-2002 JAN - original version
;-


; borrow some code from CR_REJECT, altered to reflect dist not working 
; right on aquila

	if n_elements(radius) eq 0 then radius=1

	dilation=radius
	kdim = 1 + 2*floor(dilation+1.e-4)

	kernel = make_array(kdim, kdim, value=1b)
	half_kern = fix(kdim/2)

	kernind=lindgen(kdim,kdim)
	kernx=kernind MOD kdim
	kerny = kernind / kdim
 
	kerndist=((kernx-half_kern)^2+(kerny-half_kern)^2)
	wkz = where(kerndist GT (dilation+0.0001)^2,ckz)
	IF ckz GT 0 THEN kernel[wkz] = 0b

	; mask out regions around 'dust' for good measure
	newmask=dilate(mask,kernel)

	return,newmask OR mask
end






