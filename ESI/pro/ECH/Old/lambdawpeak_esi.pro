;+
; NAME:
;  
;
; PURPOSE:
;   Associate the tabulated arc lines with the central position of
;   spectral lines. Must be run after an initial cross correlation identification
; 
;
; CALLING SEQUENCE:
;   lambdawpeack spec, wset,  lampsold, gap, xstart, lamps, bracket=bracket,  bin=bin,  $
;                   smooth=smooth, parab=parab, show=show  
;
; INPUTS:
;   spec   - 1D spectrum (along row ymid)
;   wset   - structure containing the guessed coefficients of lambda=lambda(x)
;   lampsold  - structure containing the tabulated arc lines given in input
;   gap    - initial guess for the gap across the dispersion
;            direction (integer pixel)
;
; OPTIONAL INPUTS:
;   
;
; REQUIRED KEYWORDS:
;   bin    - pixel range for peak searching right and left of the
;            cross-correlated pixel position tabulated arc lines 
;           
; OPTIONAL KEYWORDS: 
;   bracket - pixel range right and left of the central part of the spectrum over which  
;             tabulated wavelengths  and peaks are associated
;   smooth  - width for boxcar smoothing the spectrum (default = 0
;             i.e. no smoothing) 
;   parab   - number N such that 2N+1 points centered on the  peak maximum are
;             interpolated with a parabola. N=0 indicates  no interpolation  
;   show    - if set to 1, a plot of the final line-peak association
;             is shown 
;
; OUTPUTS:
;   xstart  - spec1D (row ymid) peak centers (in pix units) corresponding to
;             tabulated arclines  
;   lamps   - structure containing the tabulated arc lines (bracketed
;             and not saturated)
;
; OPTIONAL OUTPUTS:
;   
;   
; COMMENTS:
;
;
; EXAMPLES:
;
;
; BUGS:
;
;
; PROCEDURES CALLED:
;   
;
; REVISION HISTORY:
;   19-Mar-2001 Written by Chris Marinoni, Berkeley
;   4-sept-2001 modified by Andrew Sheinis,UCSC, for ESI (4096 x 2048 images)
;   
;-
;--------------------------------------------------------------------------- 

pro  lambdawpeak, spec, wset, lampsold,gap, xstart, lamps, bracket=bracket,  bin=bin,  $
                   smooth=smooth,  parab=parab,  show=show

    
   ny = (n_elements(spec)-gap)
   nny = ny+(ny+gap)*(gap gt 0)

   if NOT keyword_set(bracket) then bracket = ny
   if NOT keyword_set(parab) then parab = 0
   if NOT keyword_set(smooth) then smooth = 0
   if NOT keyword_set(show) then show = 0
   if NOT keyword_set(lev) then lev = 7
   
   xstart = traceset2pix(wset,lampsold.lambda) ;pixel value of lampsold.lambda
   keep = where(xstart le 4097)
   xstart = xstart[keep]

   qtrim = xstart GE ny-bracket AND xstart LT ny+gap+bracket-1 ; AND lamps.good 
   itrim = where(qtrim, ct)

   if (ct EQ 0) then print, 'No arc lines bracketed'
   xstart = xstart(itrim)
   lamps = lampsold(itrim)

   xstart_old = xstart
   lamps_old = lamps 

   if smooth gt 0 then spec = smooth(spec,  smooth)

   n = n_elements(xstart)-1
   

   ; Peak finding
   for i=0,n do begin  
      am = (xstart(i)-bin >  0)
      bm = (xstart(i)+bin < nny-1)
      z = max(spec[am:bm],  subscr)
      gain = 1.2                ;ccd gain from Sandy
      backg = median(spec[(xstart(i)-100 > 0):(xstart(i)+100 < nny-1)])      
      ; select only line with s/N>10
      if (z gt backg+10*sqrt(backg/gain)) then begin 
         xstart[i] = fix(am)+subscr ;pixel value (integer) of the peak center
      endif else begin 
         xstart[i] = 0         
      endelse
   endfor
   j = where(xstart gt 0, count) & if count eq 0 then print,"no lines found!"
   xstart = xstart[j]
   lamps = lamps[j]
  
   ;  Do not fit saturated lines and lines at the
   ;  edges of the chips (they can be deformed)
   ; 
lev=7
   i = where(spec[xstart] le 0.95*max(spec) and ((xstart ge lev and xstart le ny-lev) $ 
       or (xstart ge ny+gap+lev and xstart le  nny-lev)))
   xstart = xstart[i]
   lamps = lamps[i]
   a = where(xstart le 4095, count)
   if show eq 1 then begin 
      print, 'Now fitting ', n_elements(xstart), ' tabulated lines'
      print,  xstart
      print, 'Number of tabulated lines on the blue chip=',  count
      print, 'Number of tabulated lines on the red  chip=', n_elements(xstart)-count
   endif
;   if n_elements(xstart) le 20 or  n_elements(where(xstart le 4096)) le 10 $
;     or n_elements(where(xstart ge 4096)) le 10 then begin 
;      print,  '******WARNING********'
;      print,  'Fitting few lines!!'
;      stop
;   endif
 
print, 'number of lines fitting equals ', n_elements(xstart)

   ;   arclines centers as centers of the best fitting parabola
   ;
   n = n_elements(xstart)-1
   if (parab gt 0 ) then begin  
      binp = parab
      xstartp = fltarr(n+1)
      p = [20000, 0.0001, -10000]
      xp = -binp+findgen(2*binp+1)
      for i=0, n do begin 
         yp = spec[xstart(i)-binp:xstart(i)+binp]
;         R = SVDFIT(xp, yp, 3, A=p, CHISQ=chi, yfit=speci) ;only if binp>=1
         R = POLY_FIT(xp, yp, 2 ,  YFIT=speci)
         xstartp(i)=xstart(i)-R(1)/(2.*R(2))
         xp2 = -binp+(2*binp+1)*findgen(100)/(100-1)
         if show eq 1 then begin 
            plot,  xp,  yp, psym=4
            oplot, xp2, R(0)+R(1)*xp2+R(2)*xp2^2 ;plot each fitted parabola            
         endif
         if abs(xstartp(i)-xstart(i)) gt 1 then begin  
            print,  'PARABOLA CENTER NOT CONSTRAINED! need to increase bin size'
;            stop
         endif 
      endfor
      if show eq 1 then begin
         plot, xstart-xstartp,  ytitle='pixel centers-parabola centers (ARCLINES)',$
           xtitle='arc lines',  ps=4
         stop
      endif
      xstart = xstartp
   endif

   ;  check positions of the centers
   if (show eq 1) then begin
      traceset2xy, wset, xtemp, lambda ; lambda=lambda(x)
      plot, xtemp, spec
      oplot, xstart_old, lamps_old.intensity*2, ps=1 ; (best crosscorelated waves)
      oplot, xstart, lamps.intensity*2, ps=4 ;diamond assigned centers
      stop
   endif
   return
end

