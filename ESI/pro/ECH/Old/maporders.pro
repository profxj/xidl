;******************************************

pro maporders,image,outimage,coeff
;+
; routine takes coefficents (4th ordr)
; to remap data in spatial direction, order by order
; specialized for echellete mode of ESI
; 5apr00, MD
; 6 apr 01, AIS
;
; INPUT: image to be remapped
;	 coeff,  [4, 10, 5] matrices for quadratic fits, 4th order, 10 spectral
;	 orders, 5 stars per order....
;
; OUTPUT: outimage: remapped image, with no transformation in Lambda direction
;-


outimage=image*0.

;wset,0

;tv,rebin(image,1024,512)/80

;ypos=findgen(4096)
;ypos=[200,1024,4000]
ypos=(findgen(50)+1)*4095/50
slicelength=170

; begin starting loop
for kk=1,9 do begin 
   xoffset=180*kk  ; displacement of subsequent orders


   w0= (coeff[0,kk,0] + coeff[1,kk,0]*ypos + $
	coeff[2,kk,0]*ypos*ypos + coeff[3,kk,0]*ypos*ypos*ypos)  

;   w0=w0*(w0 gt 0) +0*(w0 lt 0)
 

   w1=( coeff[0,kk,1] + coeff[1,kk,1]*ypos +  $
	 coeff[2,kk,1]*ypos*ypos + coeff[3,kk,1]*ypos*ypos*ypos)
;   w1=w1*(w1 gt 0) +0*(w1 lt 0)

   w2=( coeff[0,kk,2] + coeff[1,kk,2]*ypos +  $
	 coeff[2,kk,2]*ypos*ypos + coeff[3,kk,2]*ypos*ypos*ypos)
;   w2=w2*(w2 gt 0) +0*(w2 lt 0)

   w3=( coeff[0,kk,3] + coeff[1,kk,3]*ypos +  $
	 coeff[2,kk,3]*ypos*ypos + coeff[3,kk,3]*ypos*ypos*ypos)
;   w3=w3*(w3 gt 0) +0*(w3 lt 0)

   w4=( coeff[0,kk,4] + coeff[1,kk,4]*ypos +  $
	 coeff[2,kk,4]*ypos*ypos + coeff[3,kk,4]*ypos*ypos*ypos)
;   w4=w4*(w4 gt 0) +0*(w4 lt 0)
  
;    w4=w1*(w4 gt 0) + 0*(w4 gt 0) ; bound by 0
   
;
; warp each order, all columns at once, 5 points per column
; 

  
    X0=[ypos,ypos,ypos,ypos,ypos]

    X1=X0

    Y0=[w0,w1,w2,w3,w4]

    wf0=replicate(xoffset+5.,n_elements(ypos))
    wf1=replicate(xoffset+45.,n_elements(ypos))
    wf2=replicate(xoffset+85.,n_elements(ypos))
    wf3=replicate(xoffset+125.,n_elements(ypos))
    wf4=replicate(xoffset+165.,n_elements(ypos))

    Y1=[wf0,wf1,wf2,wf3,wf4]

;tv,rebin(image,1024,512)/80
;xyouts,x0/4,y0/4,'-',/device






     polywarp,x0,y0,x1,y1,2,P,Q

    newimage=poly_2d(image,P,Q,cubic=-0.5)	;cubic spline interpolation
;     newimage=poly_2d(image,P,Q,1)		;or blinear interpolation
					;cubic spline is a little better but a lot slower
    outimage(*,xoffset:xoffset+slicelength-1)=$
      newimage(*,xoffset:xoffset+slicelength-1)
;wset,0


;tv,rebin(newimage,1024,512)/80
;xyouts,x1/4,y1/4,'-',/device
;wait,3
;tv,rebin(outimage,1024,512)/80
;xyouts,x1/4,y1/4,'-',/device
;wait,1



endfor

return

end







     
