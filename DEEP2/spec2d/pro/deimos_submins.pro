;+
; NAME:
;   deimos_submins
;
; PURPOSE:
;     subtract off a two-plane model of the background light from an image
;
; CALLING SEQUENCE:
;   fixedimage=deimos_submins(image,planeparams)
; 
; INPUTS:
;   image  - flat/data frame to subtract off background from
;
; OUTPUTS:
;   planeparams - parameters defining a two-plane fit to the set of minima in image
;   This is a 3-element x 2 array of the form 
;    [constant, x coefficient, y coefficient (where e.g. x=0...2047)], 
;   where the 0th row is for the bottom half of the chip, and the 1st row is fit for the top half.
;
; MODIFICATION HISTORY:
;   JAN 8/16/02    
;-
function deimos_submins,image,planeparams

	nrows=4096
	ncols=2048

	mask=image eq 0
	
	xarr=lindgen(ncols,nrows)
	yarr = xarr / ncols
	xarr = xarr MOD ncols
	
;	top=where(yarr ge nrows/2)
;	bot=where(yarr lt nrows/2)

	botparams=planeparams[*,0]
	topparams=planeparams[*,1]

	planefitb=botparams[0]+botparams[1]*xarr+botparams[2]*yarr
	planefitt=topparams[0]+topparams[1]*xarr+topparams[2]*yarr

; want to smoothly transition from bottom to top
	fractop= ( (yarr-0.4*nrows)/nrows*5. <1) > 0
	fracbot= ( (-yarr+0.6*nrows)/nrows*5. <1) > 0

; check for finiteness of the fits
        nofix = total(finite(planefitb+planefitt) eq 0) NE 0
	if NOT nofix then begin
          image = image-planefitb*fracbot-planefitt*fractop
          whmask = where(mask, maskct)
          if maskct gt 0 then image[whmask] = 0.
        endif

  return, image
end








