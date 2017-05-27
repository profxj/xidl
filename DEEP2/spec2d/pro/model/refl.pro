pro refl,r,a
; see refl.x

; "Reflection Algorithm"

r=double(r)
a=double(a)

; transform
rp=matrix_multiply(a,r)
; for a reflection, change (zp) --> (-zp)
rp(2)=-rp(2)
; re-transform
r=matrix_multiply(a,rp,/atranspose)

return

end
