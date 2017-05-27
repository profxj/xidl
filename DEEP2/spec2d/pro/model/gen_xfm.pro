pro gen_xfm,r,a,forward=forward

; see refl.x

; GEN_XFM: general transform of r[3] into another CS desribed by a; "forward"
; is set or not to describe if xform is into or out-of CS
; Note that the appropriate operation (eg transmission, reflection) must be
; applied afterward
;

r=double(r)
a=double(a)

if n_elements(forward) eq 0 then forward = 0

if forward then r=matrix_multiply(a,r) $
	 else r=matrix_multiply(a,r,/atranspose)

return
end
