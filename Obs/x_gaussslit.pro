;+ 
; NAME:
; x_gaussslit   
;    Version 1.1
;
; PURPOSE:
;    Calcualte the slit lost assuming a Gaussian object and a box slit
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function x_gaussslit, w, h, xo, yo

  if  N_params() LT 4 then begin 
      print,'Syntax - ' + $
        'flux = x_gaussslit(width, height, x, y) (v1.0)'
      return, -1
  endif 
  
  psf  = [1.000,   .995,   .985,   .971,   .954,        $
          .933,   .911,   .886,   .860,   .833, $
          .804,   .774,   .743,   .713,   .682, $
          .651,   .620,   .594,   .559,   .529, $
          .500,   .471,   .443,   .417,   .391, $
          .366,   .342,   .319,   .297,   .276, $
          .256,   .237,   .218,   .202,   .187, $
          .172,   .158,   .145,   .132,   .122, $
          .113,   .104,   .097,   .089,   .082, $
          .077,   .072,   .065,   .059,   .057, $
          .052,   .049,   .046,   .042,   .039, $
          .037,   .034,   .032,   .029,   .027, $
          .026,   .024,   .023,   .021,   .019, $
          .018,   .017,   .017,   .016,   .016, $
          .015,   .014,   .013,   .012,   .011, $   
          .010,   .010,   .009,   .009,   .008, $
          .008,   .007,   .007,   .006,   .006, $   
          .005,   .005,   .005,   .004,   .004, $
          .004,   .004,   .003,   .003,   .003, $
          .003,   .003,   .002,   .002,   .002]

  width = 20*w
  height = 20*h 
  xoff = 40*xo
  yoff = 40*yo

  xin=0.
  xout=0.
  for i=1,199 do begin
      y=float(100-i)
      dy=abs(y-yoff)
      for j=1,199 do begin
          x=float(j-100)
          dx=abs(x-xoff)
          radius=sqrt(x*x+y*y)
          if radius GE 99. then flux = 0. else begin
              irad=long(radius)
              drad=radius-long(radius)
              flux=(1.-drad)*psf[irad]+drad*psf[irad+1]
          endelse
          if dy LT height and  dx LT width then xin=xin+flux
          if dy GT height or  dx GT width then xout=xout+flux
          if (dy EQ height and  dx LE width) OR $
            (dx EQ width and  dy LE height) then begin
              xin=xin+.5*flux
              xout=xout+.5*flux
          endif 
      endfor
  endfor
  slit=xin/(xin+xout)
  return, slit

end


