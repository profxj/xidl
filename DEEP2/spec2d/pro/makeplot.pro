;------------------------------------------------------
FUNCTION XAXIS, axis, index, value
pix = STRING(value,FORMAT='(I10)')
RETURN, pix
END

;------------------------------------------------------

pro makeplot, file, scale, zp, rebin, sm, z, xwinsize= xwinsize, ywinsize = ywinsize, $
                                  plot2D=plot2D, _EXTRA = EXTRA_KEYWORDS

;
;  Purpose: Make a 1-D binned plot of the spectrum.
;  file - the spec1D file to be plotted  
;  scale - factor by which to multiply the spectrum
;  zp - offset term to the spectrum
;  rebin - factor for rebinning; must divide evenly into 4096
;  z - redshift
;  sm - smoothing window.  Must be > 0 and odd?
;  xwinsize - Number of windox pixels in x dimension.  Default of 1200.
;  ywinsize - Number of window pixels in y dimension.  Default of 700.
; -------------------------------------------------------------------------------------------------
;  plot2D - a flag (default of 0) to determine whether to plot the 2D spectrum, using plot2DBinnedFull.pro
;  If plot2D is set to 1, then bluespec and redspec must be passed. 
;  
;  The EXTRA_KEYWORDS that can be passed to this routine are those that can be used by the routine plot2DBinnedFull to modify the 2D spectrum. 
;  They are as follows:
;
;  xrebin - Factor by which to bin the x dimension; any positive integer. Default is 1.
;  yrebin - Factor by which to bin the y dimension; any positive integer. Default is 1.
;  line_array - Array holding any combination of integer values from 0-11. The presence of a particular integer 
;       indicates that the corresponding spectral feature will be overlaid on the 2-D spectrum. Default is all lines.
;       The lines are as follows:
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
;  specpos - Fraction (0-1) giving where in the 2D spectrum the source is located, to adjust plot lines. Default of 0.5.
;  TITLE - Optional title from the user.  Default is to display the redshift.
;  xunbin - Factor by which to "unbin" the binned image.  Default is 0 -> isn't done.
;  yunbin - " "

if(n_elements(plot2D) eq 0) then $
   plot2D=0 ;flag to specify whether to make 2D spectrum plot
if(n_elements(xwinsize)) eq 0 then $
   xwinsize = 1200
if(n_elements(ywinsize)) eq 0 then $
   ywinsize = 700

blue = mrdfits(file,1)
red = mrdfits(file,2)
print,blue.OBJPOS
print,red.OBJPOS

;Number of pixels in binned image
binsize=4096/rebin

;Perform the desired binning
bspecbin = rebin(blue.spec*blue.ivar,binsize)/rebin(blue.ivar,binsize)*scale
bwavebin = rebin(blue.lambda*blue.ivar,binsize)/rebin(blue.ivar,binsize)
bivarbin = rebin(blue.ivar,binsize)
rspecbin = rebin(red.spec*red.ivar,binsize)/rebin(red.ivar,binsize)*scale
rwavebin = rebin(red.lambda*red.ivar,binsize)/rebin(red.ivar,binsize)
rivarbin = rebin(red.ivar,binsize)

;Put together red.lambda and blue.lambda, to help with plotting the pixel value of the wavelengths
lambdatot = [blue.lambda,red.lambda]
;pix = findgen(8192)
;plot,x,lambdatot

;shortest line viewable, for this redshift
lmin = min(blue.lambda)/(1.0+z)
;longest line viewable, for this redshift
lmax = max(red.lambda)/(1.0+z)

;window, /free, xsize=xwinsize, ysize=ywinsize

;Find the min and max flux values to set plot size
;Note that the presence of -NaN values in the flux seems to require the use of this method of finding max
blueflux = ivarsmooth(bspecbin,bivarbin,sm)+zp
redflux = ivarsmooth(rspecbin,rivarbin,sm)+zp

s=size(blueflux)
maxflux=0
for i=0,s[1]-1 do begin
   if blueflux[i] gt maxflux AND blueflux[i] gt -10000 AND blueflux[i] lt 10000 then begin
      maxflux = blueflux[i]
   end
endfor

r=size(redflux)
for i=0,r[1]-1 do begin
   if redflux[i] gt maxflux AND redflux[i] gt -10000 AND redflux[i] lt 10000 then begin
      maxflux = redflux[i]
   end
endfor

minflux = maxflux
for i=0,s[1]-1 do begin
   if blueflux[i] lt minflux AND blueflux[i] gt -10000 AND blueflux[i] lt 10000 then begin
      minflux = blueflux[i]
   end
endfor

r=size(redflux)
for i=0,r[1]-1 do begin
   if redflux[i] lt minflux AND redflux[i] gt -10000 AND redflux[i] lt 10000 then begin
      minflux = redflux[i]
   end
endfor

;plot blue portion of spectrum
plot,bwavebin/(1.0+z),blueflux,xrange=[lmin,lmax],yrange=[round(minflux)-round((maxflux-minflux)/5),(maxflux+maxflux/10)], ytickinterval=50,ystyle=1, thick=1, xstyle=8, $
XTITLE='Rest-frame Wavelength (Angstroms)', YTITLE='Flux'
;plot red portion of spectrum (oplot adds to plot  
oplot, rwavebin/(1.0+z),redflux,thick=1

AXIS,XAXIS=1, XRANGE = [1,binsize*rebin*2], XTICKFORMAT='XAXIS',XTICKINTERVAL=500, XSTYLE=1, XTITLE='Approximate Pixel Position',XTICKLEN=0.06

;Plot a fit for high-z LBG galaxy
fmt = 'F5,F'
readfmt,'/home/phobos/capak/IDL/lbg_shapley.spec',fmt,lambda,spec
;not sure what this does...have to do with skylines?
oplot,lambda,spec*7e30,color='ff00ff'XL 

;these i don't quite understand, but think they give some idea of the variance on the plot at that wavelength.
oplot,bwavebin/(1.0+z),1.0/sqrt(bivarbin*rebin*sm),color='FFFF00'XL
oplot,rwavebin/(1.0+z),1.0/sqrt(rivarbin*rebin*sm),color='FFFF00'XL

;Display low redshift lines
oplot,[3727,3727],[maxflux/3,-20],color='00ff00'XL
xyouts, 3727, -47, '!6OII!X',orientation=90,charsize=1.3
oplot,[4861,4861],[maxflux/3,-20],color='00ff00'XL
xyouts, 4861, -65, '!6H-beta!X',orientation=90,charsize=1.3
oplot,[4363,4363],[maxflux/3,-20],color='00ff00'XL
xyouts, 4363, -47, '!6OIII!X',orientation=90,charsize=1.3
oplot,[4959,4959],[maxflux/3,-20],color='00ff00'XL
xyouts, 4959, -47, '!6OIII!X',orientation=90,charsize=1.3
oplot,[5007,5007],[maxflux/3,-20],color='00ff00'XL
xyouts, 5007, -47, '!6OIII!X',orientation=90,charsize=1.3
oplot,[6563,6563],[maxflux/3,-20],color='00ff00'XL
xyouts, 6563, -70, '!6H-alpha!X',orientation=90,charsize=1.3

;------------------------------------------------------------------
;mark EM line 

;oplot,[7110/(1.0+z),7110/(1.0+z)],[60,-20],color='ff0000'XL

;------------------------------------------------------------------
;Display high-redshift lines

;FeII/AlII
oplot,[1670.79,1670.79],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1670.79, round(minflux)-round((maxflux-minflux)/5)+3, '!6FeII/AlII!X',orientation=90,charsize=1.3

;FeII
oplot,[1608.45,1608.45],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1 
xyouts, 1608.45,round(minflux)-round((maxflux-minflux)/5)+3, '!6FeII!X',orientation=90,charsize=1.3

;CIV
oplot,[1550.78,1550.78],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
oplot,[1548.20,1548.20],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1549,round(minflux)-round((maxflux-minflux)/5)+3, '!6CIV!X',orientation=90,charsize=1.3

;SiII
oplot,[1526.71,1526.71],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1526.71,round(minflux)-round((maxflux-minflux)/5)+3, '!6SiII!X',orientation=90,charsize=1.3

;SiIV
oplot,[1402.77,1402.77],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
oplot,[1393.76,1393.76],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1401, round(minflux)-round((maxflux-minflux)/5)+3, '!6SiIV!X',orientation=90,charsize=1.3

;CII
oplot,[1334.53,1334.53],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1334.53, round(minflux)-round((maxflux-minflux)/5)+3, '!6CII!X',orientation=90,charsize=1.3

;OI/SiII
oplot,[1304.37,1304.37],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
oplot,[1302.17,1302,17],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1304.4, round(minflux)-round((maxflux-minflux)/5)+3, '!6OI/SiII!X',orientation=90,charsize=1.3

;SiII
oplot,[1260.42,1260.42],[round(minflux)-20,maxflux],color='0000ff'XL, linestyle=1
xyouts, 1260.42, round(minflux)-round((maxflux-minflux)/5)+3, '!6SiII!X',orientation=90,charsize=1.3

;CIII
oplot,[1176,1176],[round(minflux)-20,maxflux],color='00ff00'XL, linestyle=1
xyouts, 1176, round(minflux)-round((maxflux-minflux)/5)+3, '!6CIII!X',orientation=90,charsize=1.3

;?
;oplot,[1343,1343],[-100+zp, maxflux],color='00ff00'XL, linestyle=1
;xyouts, 1343, -25, '!6?!X',orientation=90,charsize=1.3,alignment=0.5
;?
;oplot,[1501,1501],[round(minflux)-20,maxflux],color='00ff00'XL, linestyle=1
;xyouts, 1501, -25, '!6?!X',orientation=90,charsize=1.3,alignment=0.5

;Lyalpha
oplot,[1215.67,1215.67],[round(minflux)-20,maxflux],color='00ff00'XL, linestyle=1
xyouts, 1216, round(minflux)-round((maxflux-minflux)/5)+3, '!6Lyalpha!X',orientation=90,charsize=1.3

;NV
oplot,[1238,1238],[round(minflux)-20,maxflux],color='ffa300'XL, linestyle=1
oplot,[1242,1242],[round(minflux)-20,maxflux],color='ffa300'XL, linestyle=1
xyouts, 1242, round(minflux)-round((maxflux-minflux)/5)+3, '!6NV!X',orientation=90,charsize=1.3

;oplot,[1393,1393],[round(minflux)-20,maxflux+maxflux/10],color='ffa300'XL, linestyle=1
;oplot,[1402,1402],[round(minflux)-20,maxflux+maxflux/10],color='ffa300'XL, linestyle=1

;CIV
oplot,[1548,1548],[round(minflux)-20,maxflux],color='ffa300'XL, linestyle=1
xyouts, 1548, round(minflux)-round((maxflux-minflux)/5)+3, '!6CIV!X',orientation=90,charsize=1.3

;oplot,[1550,1550],[round(minflux)-20,maxflux+maxflux/10],color='ffa300'XL, linestyle=1

;HeII
oplot,[1640,1640],[round(minflux)-20,maxflux],color='ffa300'XL, linestyle=1
xyouts, 1640, round(minflux)-round((maxflux-minflux)/5)+3, '!6HeII!X',orientation=90,charsize=1.3

;oplot,[1215.67,1215.67],[round(minflux)-20,maxflux+maxflux/10],color='ff00ff'XL, linestyle=1
;oplot,[1238.81,1238.81],[round(minflux)-20,maxflux+maxflux/10],color='ff00ff'XL, linestyle=1

aband=7630.0/(z+1.0)
oplot,[aband,aband],[round(minflux)-20,maxflux+maxflux/10],color='ff0000'XL,thick=3, linestyle=2
xyouts, aband, round(minflux)-round((maxflux-minflux)/5)+3, '!6A Band!X',orientation=90,charsize=1.3

bband=6890/(z+1.0)
oplot,[bband,bband],[round(minflux)-20,maxflux+maxflux/10],color='ff0000'XL,thick=3, linestyle=2
xyouts, bband, round(minflux)-round((maxflux-minflux)/5)+3, '!6B Band!X',orientation=90,charsize=1.3

;oplot,[],[60,0],color='ff00ff'XL
;oplot,[],[60,0],color='ff00ff'XL

;------------------------------------------------------------------------
;Make the 2D plot, if specified
;First determine the slit number:
s1 = strmid(file,7)
pos = strpos(s1,'.')
s2 = strmid(s1,pos+1,3) ;This should be the slit number
bluespecfile = 'slit.COS.'+s2+'B.fits'
redspecfile = 'slit.COS.'+s2+'R.fits'

if (plot2D ne 0) then $
   plot2DBinnedFull, bluespecfile, redspecfile, z, _extra = extra_keywords




end ;makeplot

