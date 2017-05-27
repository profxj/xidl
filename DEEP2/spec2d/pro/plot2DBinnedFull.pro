;------------------------------------------------------
;Function from book to do histogram clipping
function imclip, image, PERCENT=PERCENT

;Skipping argument checks

if (n_elements(percent) eq 0) then percent = 2.0

;Get image max and min
min_value = min(image, max=max_value)

;Compute histogram 
nbins = 100
binsize = float(max_value-min_value) / float(nbins)
hist = histogram(float(image),binsize=binsize)
bins = lindgen(nbins+1)*binsize + min_value

;Compute normalized cumulative sum
sum = fltarr(n_elements(hist))
sum[0]=hist[0]
for i = 1L,n_elements(hist)-1L do $
   sum[i] = sum[i-1]+hist[i]
sum=100.0*(sum/float(n_elements(image)))

;Find and return range
range = [min_value, max_value]
index = where((sum ge percent) and $
(sum le (100.0-percent)),count)
if (count ge 2) then $
   range = [bins[index[0]],bins[index[count-1]]]
return, range

end

;---------------------------------------------------------------------------------------------------------------------

pro plot2DBinnedFull, bluespecfile,redspecfile, z, xrebin=xrebin, yrebin=yrebin, line_array=line_array, specpos=specpos, TITLE=TITLE, xunbin=xunbin, yunbin=yunbin
  
;  bluespecfile - Blue spectrum structure returned from call to mrdfits, i.e. from bluespec = mrdfits('slit.COS.019B.fits',1)
;  redspecfile - Red spectrum structure returned from call to mrdfits, i.e. from redspec = mrdfits('slit.COS.019R.fits',1)
;  z - Redshift
;  xrebin - Factor by which to bin the x dimension, can be any positive integer. Default is 4.
;  yrebin - Factor by which to bin the y dimension; any positive integer. Default is 1.
;  line_array - Array holding any combination of integer values from 0-11. The presence of a particular integer 
;     indicates that the corresponding spectral feature will be overlaid on the 2-D spectrum. Default is to show all lines.
;     The lines are as follows:
;
;     0 - CIII, 1176
;     1 - Lyalpha, 1215.7
;     2 - NV, 1238
;     3 - SiII, 1260.4
;     4 - OI/SiII, 1303
;     5 - CII, 1334.5
;     6,7 - SiIV, 1393.8, 1402.8
;     8 - SiII, 1526.7
;     9 - CIV, 1549
;     10 - FeII, 1608.5
;     11 - FeII/AlII, 1670.8
;  
;  specpos - Fraction (0-1) giving how far up into the 2D spectrum the source is located, to adjust plot lines. Default of 0.5.
;  TITLE - Optional title from the user.  Default is to display the redshift.
;  xunbin - Factor by which to "unbin" the binned image, for display size.  Default is 0 -> isn't done.
;  yunbin - " "
;
;  Example call: plot2DBinnedFull,bb,rr,4.11,xrebin=4,yrebin=3,line=[2,4,5],Title='Redshift of 4.11'
;  Here an object with redshift 4.11 is plotted, with the lines overlaid from NV, OI/SiII, and CII, and the title 'Redshift of 4.11'
;  It is binned by a factor of 4 in the x dimension and 3 in the y dimension.

bluespec = mrdfits(bluespecfile,1)
redspec = mrdfits(redspecfile,1)
  
if(n_params() ne 3) then $
   message, 'Usage: bluespec, redspec, z'
if(n_elements(bluespec) eq 0) then $
   message, 'Argument bluespec is undefined'
if(n_elements(redspec) eq 0) then $
   message, 'Argument redspec is undefined'
if(n_elements(z) eq 0) then $
   message, 'Argument z is undefined'
if(n_elements(xrebin) eq 0) then $
   xrebin = 4 ;Default is to rebin x by factor of 4
if(n_elements(yrebin) eq 0) then $
   yrebin = 1 ;Default is no binning in y
if(n_elements(line_array) eq 0) then $
   line_array = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] ;Default is to show all lines
if(n_elements(specpos) eq 0) then $
   specpos = 0.5
if(n_elements(TITLE) eq 0) then $
   TITLE = '!5z = ' + string(z,format='(F5.2)')
if(n_elements(xunbin) ne 0) then begin $
   if (n_elements(yunbin) eq 0) then $
      yunbin = 1 ;if xunbin is defined but yunbin isn't, set yunbin to 1
endif else begin $
   if (n_elements(yunbin) ne 0) then $
      xunbin = 1 ;if yunbin is defined by xunbin isn't, set xunbin to 1
endelse
   
if((n_elements(yunbin) eq 0) and (n_elements(xunbin) eq 0)) then begin $
   yunbin = 0
   xunbin = 0
endif

;-------------------------------------------
;Binning bluespec:

s = size(bluespec.flux)
xnumpix_1 = s[1]
ynumpix_1 = s[2]

while (xnumpix_1 mod xrebin) ne 0 do begin
   xnumpix_1 = xnumpix_1+1
endwhile

while (ynumpix_1 mod yrebin) ne 0 do begin
   ynumpix_1 = ynumpix_1+1
endwhile

;Create new variance and flux array for rebinning
flux = fltarr(xnumpix_1,ynumpix_1)
ivar = fltarr(xnumpix_1,ynumpix_1)

flux(0:s[1]-1,0:s[2]-1)=bluespec.flux
ivar(0:s[1]-1,0:s[2]-1)=bluespec.ivar

;Get the binned flux
newxpix_1 = xnumpix_1/xrebin
newypix_1 = ynumpix_1/yrebin
binnedflux_1 = rebin(flux*ivar,newxpix_1,newypix_1)/sqrt(rebin(ivar,newxpix_1,newypix_1))

;-------------------------------------------
;Binning redspec:

s = size(redspec.flux)
xnumpix_2 = s[1]
ynumpix_2 = s[2]

while (xnumpix_2 mod xrebin) ne 0 do begin
   xnumpix_2 = xnumpix_2+1
endwhile

while (ynumpix_2 mod yrebin) ne 0 do begin
   ynumpix_2 = ynumpix_2+1
endwhile

;Create new variance and flux array for rebinning
flux = fltarr(xnumpix_2,ynumpix_2)
ivar = fltarr(xnumpix_2,ynumpix_2)

flux(0:s[1]-1,0:s[2]-1)=redspec.flux
ivar(0:s[1]-1,0:s[2]-1)=redspec.ivar

;Get the binned flux
newxpix_2 = xnumpix_2/xrebin
newypix_2 = ynumpix_2/yrebin
binnedflux_2 = rebin(flux*ivar,newxpix_2,newypix_2)/sqrt(rebin(ivar,newxpix_2,newypix_2))

;-----------------------------------------
;Perform the unbinning, if specified

if (xunbin ne 0) then begin $
   newxpix_1 = newxpix_1 * xunbin
   newypix_1 = newypix_1 * yunbin
   xrebin_1 = xrebin / xunbin
   yrebin_1 = yrebin / yunbin
   binnedflux_1 = rebin(binnedflux_1,newxpix_1,newypix_1)
   newxpix_2 = newxpix_2 * xunbin
   newypix_2 = newypix_2 * yunbin
   xrebin_2 = xrebin / xunbin
   yrebin_2 = yrebin / yunbin
   binnedflux_2 = rebin(binnedflux_2,newxpix_2,newypix_2)
endif

;-----------------------------------------
;Combine the two binned spectra into one big spectrum.  
;The missing pixels in between are set to -1000.

;Determine the number of blank pixels to place between the two sides, 
;since they don't align perfectly
wv_break = redspec.lambda0[0]-bluespec.lambda0[4095]
wv_range = redspec.lambda0[4095]-bluespec.lambda0[0]

;Note - newxpix_1 and newxpix_2 should be the same...just lazy coding
junkpix = round ((newxpix_1+newxpix_2)*wv_break / (wv_range-wv_break)) 

xpix = max([newxpix_1,newxpix_2])
ypix=max([newypix_1,newypix_2])

binnedflux = fltarr(2*xpix+junkpix,ypix)
binnedflux(0:xpix-1,0:newypix_1-1)=binnedflux_1
binnedflux(newxpix_1:(xpix+junkpix)-1,0:ypix-1) = -1000
binnedflux(xpix+junkpix:(2*newxpix_2+junkpix)-1,0:newypix_2-1) = binnedflux_2

newxpix = 2*xpix+junkpix
newypix = ypix

;Create a new lambda vector, for the purpose of mapping wavelengths to pixels:
junkentries =round(8192*(wv_break / wv_range))
lambda0 = fltarr(8192+junkentries)
lambda0(0:4095) = bluespec.lambda0
lambda0(4096:4096+junkentries-1) = 0
lambda0(4096+junkentries:8191+junkentries)=redspec.lambda0


;--------------------------------------------------
;Complete the display of the binned, combined spectrum

;Find the max flux, for setting the brightness
;First check for and remove -NaN values
missing = -1000
index = where(finite(binnedflux) eq 0, count)
if (count gt 0) then binnedflux[index] = missing
maxflux = max(binnedflux)
;index = where(binnedflux eq missing, count)
;if (count gt 0) then binnedflux[index] = maxflux ;Replace NaN vals with the maximum flux
index2 = where(binnedflux gt -1000)
range = imclip(binnedflux(index2))


;Get colors ready, decide on window size
;device, decomposed=0
;tvlct,0B,255B,0B,128 ;green
;tvlct,185B,211B,238B,129 ;slate gray

xwinsize = 1.2*newxpix;
ywinsize = 7*newypix;

window,/free, xsize=xwinsize, ysize=ywinsize
;tv, imscale(binnedflux,range=imclip(binnedflux[index2])), newxpix / 10, 3*newypix, channel=3
tv, bytscl(binnedflux,range[0],range[1], top=(!d.table_size-1)), newxpix / 10, 3*newypix, channel=3

;-------------------------------------------------------------------------
;Moving from plotting the 2D spectrum to plotting lines/text on top of it

;Comute redshifted lines
Lines = fltarr(12)
;CIII
  Lines[0] = 1176*(z+1)
;Lyalpha
  Lines[1] = 1215.7*(z+1)
;NV
  Lines[2] = 1238*(z+1)
;SiII
  Lines[3] = 1260.4*(z+1)
;OI/SiII
  Lines[4] = 1303*(z+1)
;CII
  Lines[5] = 1334.5*(z+1)
;SiIV
  Lines[6] = 1393.8*(z+1)
  Lines[7] = 1402.8*(z+1)
;SiII
  Lines[8] = 1526.7*(z+1)
;CIV
  Lines[9] = 1549*(z+1)
;FeII
  Lines[10] = 1608.5*(z+1)
;FeII/AlII
  Lines[11] = 1670.8*(z+1)

;If xunbin and yunbin were set to zero, set them to 1 so they have no effect on the following calculations.
if (xunbin eq 0 and yunbin eq 0) then begin $
   xunbin = 1
   yunbin = 1
endif

;Plot A Band and B Band
abandwv_diff = abs(lambda0 - 7630)
bbandwv_diff = abs(lambda0 - 6890)
minidx1 = where(abandwv_diff eq (min(abandwv_diff)))
minidx2 = where(bbandwv_diff eq (min(bbandwv_diff)))

if (min(abandwv_diff) lt 3) then begin $
   destpix = fix(float(minidx1) / (xrebin/xunbin))
   plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/6)],thick=.5,color='B8D380'XL,linestyle=2, /device
   plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/6),(3*newypix+2.5*newypix)],thick=.5,color='B8D380'XL,linestyle=2, /device
   xyouts,[newxpix/10+destpix,newxpix/10+destpix],3*newypix+2.6*newypix, 'A Band', alignment=0.5,charsize=1.2,/device
endif

if (min(bbandwv_diff) lt 3) then begin $
   destpix = fix(float(minidx2) / (xrebin/xunbin))
   plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/6)],thick=.5,color='B8D380'XL,linestyle=2, /device
   plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/6),(3*newypix+2.5*newypix)],thick=.5,color='B8D380'XL,linestyle=2, /device
   xyouts, [newxpix/10+destpix,newxpix/10+destpix],3*newypix+2.6*newypix, 'B Band', alignment=0.5,charsize=1.2,/device
endif

;Plot specified lines
for i = 0,11 do begin
   line = where(line_array eq i, count)
   if (count gt 0) then begin $
      ;Print redshift
      ;zstr = string(z,format='(F5.2)')
      xyouts,0.6*newxpix,6.3*newypix, TITLE,font=1,charsize=2.5,charthick=2,alignment=0.5, /device
      ;Find the correct pixel position to place this line
      lambda = Lines[i]
      wv_diff = abs(lambda0 - lambda)
      minidx = where(wv_diff eq min(wv_diff))
      if (min(wv_diff) lt 3) then begin $
         destpix = fix(float(minidx) / (xrebin/xunbin))
         case i of
            0: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6CIII', alignment=0.5,charsize=1.2,/device
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[0],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            1: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6Lyalpha', alignment=0.5,charsize=1.2,/device
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[1],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            2: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6NV', alignment=0.5,charsize=1.2,/device           
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[2],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
               end
            3: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6SiII', alignment=0.5,charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[3],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end

            4: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6OI/SiII', alignment=0.5,charsize=1.2,/device   
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[4],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            5: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6CII', alignment=0.5,charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[5],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            6: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6SiIV', alignment=0,charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[6],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            7: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
              ;xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6SiIV', charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[7],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            8: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6SiII', alignment=0.5,charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[8],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end
            9: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6CIV', alignment=0.5,charsize=1.2,/device 
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[9],format='(F7.1)'),alignment=0.5, charsize=1.2,/device
            end 
            10: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6FeII', alignment=0.5,charsize=1.2,/device  
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[10],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end 
            11: begin 
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix,3*newypix+(specpos*newypix-newypix/8)],thick=.5,color='00ff00'XL, /device
               plots,[newxpix/10+destpix,newxpix/10+destpix],[3*newypix+(specpos*newypix+newypix/8),(3*newypix+2*newypix)],thick=.5,color='00ff00'XL, /device
               xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, '!6FeII/AlII', alignment=0.5,charsize=1.2,/device
              ; xyouts, newxpix/10+destpix,3*newypix+2.1*newypix, string(Lines[11],format='(F7.1)'),alignment=0.5,charsize=1.2, /device
            end 
         endcase
      endif
   endif
endfor
      



end

