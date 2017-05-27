function deimos_omodel, ccdnumber, bluslits, header
;+
; NAME:
;   deimos_omodel
;
; PURPOSE:
;   generate estimated wavelengths based on optical model of DEIMOS
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   model = deimos_omodel( ccdnumber, bluslits, header )
;
; INPUTS:
;   ccdnumber -- which ccd is this?
;   bluslits -- blueprint slit positions cut for this mask (from FITS header)
;               header   -- FITS header of image detailing grating,
;                           Ggrating angle
; KEYWORDS:
;
; OUTPUTS:
;   model -- lambda as a function of position for each slitlet,
;            returned as a cubic polynomial function for each slitlet 
;               model(4, nslits) 
;
; MODIFICATION HISTORY:
;  md 11feb02
;-

deimos_grating,header,g_rule,grangle,lambda_c

slider=sxpar(header,'GRATEPOS')

if slider eq 3 then tltval=sxpar(header,'G3TLTVAL') $
  else tltval=sxpar(header,'G4TLTVAL') 

grname = floor(g_rule+1E-5)



; determine parameters given grating, slider, and tilt
omodel_params, grname, slider, tltval, roll, o3, mu


; set up system variables
deimos_map=sysinit(mu,gmm=g_rule,o3=o3,roll3=roll)


npoints = 250
nslits = n_elements(bluslits)
lambda = findgen(npoints)*24.  +4000.


;note sign error in mask 3110 -- to be fixed in next mask!!
;qmodel, deimos_map, replicate(1, npoints)#bluslits.xmm, $

; it looks like that error no longer is occurring.

; pick the appropriate amap & bmap for this data
choose_amapbmap,header,amapfile,bmapfile

; now call optical model for an array of points

;xi1/yi1 = mosaic coordinate system
;xp1/yp1 = CCD coordinate system

xmmplus=bluslits.xmm+(bluslits.xmmt-bluslits.xmm)*0.05d0
ymmplus=bluslits.ymm+(bluslits.ymmt-bluslits.ymm)*0.05d0

xmmminus=bluslits.xmm+(bluslits.xmmb-bluslits.xmm)*0.05d0
ymmminus=bluslits.ymm+(bluslits.ymmb-bluslits.ymm)*0.05d0


qmodel, deimos_map, replicate(1, npoints)#bluslits.xmm, $
    replicate(1, npoints)#bluslits.ymm, lambda#replicate(1, nslits),  $
     xi1, yi1,  ccdnum1, xp1, yp1, amapfile, $
      bmapfile, /cubic, goodloc=goodloc

; go plus or minus a bit to establish slope
qmodel, deimos_map, replicate(1, npoints)#xmmplus, $
    replicate(1, npoints)#ymmplus, lambda#replicate(1, nslits),  $
     xi1p, yi1p,  ccdnum1p, xp1p, yp1p, amapfile, $
      bmapfile, /cubic

qmodel, deimos_map, replicate(1, npoints)#xmmminus, $
    replicate(1, npoints)#ymmminus, lambda#replicate(1, nslits),  $
     xi1m, yi1m,  ccdnum1m, xp1m, yp1m, amapfile, $
      bmapfile, /cubic

xp1 = reform(xp1, npoints, nslits, /over)
yp1 = reform(yp1, npoints, nslits, /over)
ccdnum1 = reform(ccdnum1, npoints, nslits, /over)
goodloc = reform(goodloc, npoints, nslits, /over)

xp1p = reform(xp1p, npoints, nslits, /over)
yp1p = reform(yp1p, npoints, nslits, /over)
xp1m = reform(xp1m, npoints, nslits, /over)
yp1m = reform(yp1m, npoints, nslits, /over)


ncoeff=6
ntilt=3

model = { lambda_y:dblarr(ncoeff), $
          lambda_y_top:dblarr(ncoeff),lambda_y_bottom:dblarr(ncoeff) , $
          xb: 0., xt: 0., tiltx: dblarr(ntilt)}
model =  replicate( model,  nslits ) 
npix=4096

; determine wavelength solution for each slit
for i=0, nslits-1 do begin

   onchip = where(ccdnum1[*,i] eq ccdnumber-1 AND (goodloc[*,i] eq 1) ,onct)
   if onct gt 10 then onchip = onchip[1 < (onct-1):(onct-2)>0]
                                ;fetch solution for this chip
   if n_elements(onchip) gt 10 then begin 
      xon = yp1[onchip,i]
      xfit= reform(xon/(npix/2.) -1.)
      yon = lambda[onchip] 

;fit lambda(y) where y is spectral direction

      model[i].lambda_y=svdfit(xfit,yon,ncoeff, /double, /legendre,yfit=yfit)
; set up fit for tilt
      tiltfit=-(yp1p[*,i]-yp1m[*,i])/((xp1p[*,i]-xp1m[*,i]) > 1.E-10)

; do things this way to avoid floating-point error messages
      derivs=tiltfit*0.
      whok=where(yp1[*,i]-shift(yp1[*,i],1) NE 0.,okct)
      derivs[whok]=deriv(yp1[whok,i],lambda[whok])
      tiltfit=tiltfit*derivs

      tiltfit=tiltfit / polyleg(yp1[*,i]/(npix/2.) -1., model[i].lambda_y)


; restrict range to avoid bad regions/derivatives
      xontilt=yp1[onchip[1:n_elements(onchip)-2],i]
      xfit = reform(xontilt/(npix/2.) -1.)

      tiltfit=tiltfit[onchip[1:n_elements(onchip)-2]]
; fir for tilt
      model[i].tiltx=svdfit(xfit,tiltfit,ntilt, /double, /legendre,yfit=yfit)

   endif
endfor



; REPEAT FOR SLIT BOTTOMS...

qmodel, deimos_map, replicate(1, npoints)#bluslits.xmmb, $
    replicate(1, npoints)#bluslits.ymmb, lambda#replicate(1, nslits),  $
     xi1b, yi1b,  ccdnum1b, xp1b, yp1b, amapfile, $
      bmapfile, /cubic, goodloc=goodloc

xp1b = reform(xp1b, npoints, nslits, /over)
xi1b = reform(xi1b, npoints, nslits, /over)
yp1b = reform(yp1b, npoints, nslits, /over)
ccdnum1b = reform(ccdnum1b, npoints, nslits, /over)
goodloc = reform(goodloc, npoints, nslits, /over)

for i=0, nslits-1 do begin
   onchip = where(ccdnum1b[*,i] eq ccdnumber-1 AND (goodloc[*,i] eq 1), onct )
   if onct gt 10 then onchip = onchip[1 < (onct-1):(onct-2)>0]
 
                                ;fetch solution for this chip
   if n_elements(onchip) gt 10 then begin 
      xon = yp1b[onchip,i]
      xfit= reform(xon/(npix/2.) -1.)
      yon = lambda[onchip] 
      model[i].lambda_y_bottom=svdfit(xfit,yon,ncoeff, $
                                      /double, /legendre,yfit=yfit)
      model[i].xb=median(xp1b[onchip,i],/even)
   endif
endfor

; AND FOR SLIT TOPS

qmodel, deimos_map, replicate(1, npoints)#bluslits.xmmt, $
    replicate(1, npoints)#bluslits.ymmt, lambda#replicate(1, nslits),  $
     xi1t, yi1t,  ccdnum1t, xp1t, yp1t, amapfile, $
      bmapfile, /cubic, goodloc=goodloc
xp1t = reform(xp1t, npoints, nslits, /over)
xi1t = reform(xi1t, npoints, nslits, /over)
yp1t = reform(yp1t, npoints, nslits, /over)

ccdnum1t = reform(ccdnum1t, npoints, nslits, /over)
goodloc = reform(goodloc, npoints, nslits, /over)

for i=0, nslits-1 do begin
   onchip = where(ccdnum1t[*,i] eq ccdnumber-1 AND (goodloc[*,i] eq 1) ,onct)
     if onct gt 10 then onchip = onchip[1 < (onct-1):(onct-2)>0]
 
                                ;fetch solution for this chip
   if n_elements(onchip) gt 10 then begin 
      xon = yp1t[onchip,i]
      xfit= reform(xon/(npix/2.) -1.)
      yon = lambda[onchip] 
      model[i].lambda_y_top=svdfit(xfit,yon,ncoeff, $
                                      /double, /legendre,yfit=yfit)
      model[i].xt=median(xp1t[onchip,i],/even)
   endif else begin
       whgood=where(goodloc[*,i] eq 1,goodct)
; check if this chip is red or blue
       if goodct gt 0 AND model[i].xb gt 0 then begin
           if long(ccdnumber)/4 eq 0 then whgood=whgood[0:goodct/2] $
             else whgood=whgood[goodct/2:*]
           model[i].xt=model[i].xb+median( (xi1t-xi1b)[whgood,i],/even)
       endif else model[i].xt = -1
   endelse

   if model[i].xb le 0 then begin
         whgood=where(goodloc[*,i] eq 1,goodct)
          if goodct gt 0 AND model[i].xt gt 0 then begin
           if long(ccdnumber)/4 eq 0 then whgood=whgood[0:goodct/2] $
             else whgood=whgood[goodct/2:*]
           model[i].xb=model[i].xt-median( (xi1t-xi1b)[whgood,i],/even)
       endif else model[i].xb = -1
   endif

endfor

; Note: consider adding 1 extra slit in either direction (off-chip, extrapolated)


return, model
end



