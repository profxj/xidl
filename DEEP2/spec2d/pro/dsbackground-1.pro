;+
; NAME:
;       DSBACKGROUND
;
; PURPOSE:
;       Fit and subtract a row or column background
;
; INPUTS:
;       A 2D image array
;
; KEYWORD PARAMETERS:
;
;       TITLE      - prefix to output fits name; default='s_'
;       AXIS       - axis along which to fit (1=row, 2=column) [default=2]
;       ORDER      - polynomial order to fit
;	ORDERSCALE - factor used to automatically derive ORDER [default=64]
;       HSIG       - high sigma threshold over projected median
;       LSIG       - low sigma threshold over projected median
;       NITERATE   - number of rejection iterations
;       GROW       - distance around rejected pixels to also reject [default=0]
;       DISPLAY    - switch; display resulting image
;       WRITE      - switch; write out resulting fits imager
;       HELP       - switch; returns quick keyword summary
;       VERBOSE    - switch; give extra information
;
; OUTPUTS:
;       Fits a polynomial to each row or column of an image.  This function is 
;       then subtracted to create an output (sky-subtracted) image.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       image_s = dsbackground('image.fits', TITLE = title, $
;             AXIS = axis, ORDER = order, ORDERSCALE = orderscale, $ 
;             HSIG = hsig, LSIG = lsig, NITERATE = niterate, GROW = grow, $
;             /DISPLAY, /WRITE, /HELP, /VERBOSE)
;
; COMMENTS:
;       None.
;
; PROCEDURES USED:
;       LITERSTAT
; 
; MODIFICATION HISTORY:
;       DS '03oct22 - writing
;	   '03dec02 - added ORDERSCALE
;-
function dsbackground, file,                 $
                    TITLE = title,           $
                    AXIS = axis, 	     $
                    ORDER = order,           $
		    ORDERSCALE = orderscale, $
                    HSIG = hsig,             $
                    LSIG = lsig,             $
                    NITERATE = niterate,     $
                    GROW = grow,             $
                    DISPLAY = display,       $
                    WRITE = write,           $
                    HELP = help,             $
		    TEST = test,	     $
                    VERBOSE = verbose

     if n_elements(help) ne 0 then begin
         print,' '
         print,'image_s = dsbackground(''image.fits'', TITLE = title, $'
         print,'      AXIS = axis, ORDER = order, ORDERSCALE = orderscale, $'
         print,'      HSIG = hsig, LSIG = lsig, $'
         print,'      NITERATE = niterate, GROW = grow, $'
         print,'      /DISPLAY, /WRITE, /HELP, /VERBOSE)'
         print,' '
         return,-1
     endif

     if n_elements(file) eq 0 then begin
         print,'DSBACKGROUND: No fits image provided'
         return,-1
     endif

     colors=getcolor(/load)
     if n_elements(TITLE)      eq 0 then title='s_'
     if n_elements(AXIS)       eq 0 then AXIS = 2
     if n_elements(ORDER)      eq 0 then ORDER = 8
     if n_elements(ORDERSCALE) eq 0 then ORDERSCALE = 64
     if n_elements(HSIG)       eq 0 then HSIG = 2.0
     if n_elements(LSIG)       eq 0 then LSIG = 4.0
     if n_elements(NITERATE)   eq 0 then NITERATE = 10
     if n_elements(GROW)       eq 0 then GROW = 0


     ; check if input is FITS file or IDL array:
     if size(file,/type) eq 7 then begin
       print,'% DSBACKGROUND: background subtracting '+file
       image=readfits(file,hd)
     endif else if size(file,/type) eq 4 then begin
       print,'% DSBACKGROUND: background subtracting image'
       image = file
     endif else if size(file,/type) eq 5 then begin
       print,'% DSBACKGROUND: background subtracting image'
       image = file
     endif else begin
       print,'% DSBACKGROUND: Please provide a fits image'
       return,-1
     endelse


     image_s = image
     sz=size(image,/dimen)
     ; implementation of AXIS, part I (setup):
     if (axis eq 2) then begin 
       x = findgen(sz[1])
       lastsub = sz[0] - 1
       lastsub2 = sz[1] - 1
     endif else begin
       x = findgen(sz[0])
       lastsub = sz[1] - 1
       lastsub2 = sz[0] - 1
     endelse

     ; ORDER depends on length of slit:
     ; JAN13 04: cut off at 13th, not 15th - may be too much to calc.
     ; JAN16 04: cut off at 12th, not 15th - may be too much to calc.
     ; JAN22 04: cut off at 11th, not 12th - helped one, but not ms16m2
     order = sz[1] / orderscale
     if order lt 3 then order = 3
     if order gt 11 then order = 11
     print,'% DSBACKGROUND: fitting with polynomial of order',strcompress(order)

     ; DISPLAY original image, if requested:
     if n_elements(display) ne 0 then begin
        literstat,image,imagestat,/silent
        loadct,0
        tvlct,r,g,b,/get
        tvlct,reverse(r),reverse(g),reverse(b)
	window,0
        plotimage,bytscl(image, $
           min=imagestat.median-imagestat.sigma, $
           max=imagestat.median+5.*imagestat.sigma),/preserve
     endif

     if n_elements(test) ne 0 then begin
       lastsub = 100
     endif
     for i = 0, lastsub do begin 

       ; implementation of AXIS, part II (for each row/col):
       if (axis eq 2) then begin
         y = image[i,*]
       endif else begin
         y = image[*,i]
       endelse
       y = y[*]
       y2 = y

       if n_elements(test) ne 0 then begin
	 plot,y,xrange=[0,385],/xst,yrange=[1000,1100],/yst, $
           psym=10,charthick=2 
       endif

       for j = 1, NITERATE do begin
         poof = poly_fit(x,y2,ORDER,yfit=yy)
	 literstat,(y-yy),ystat,/silent
	 hthresh = yy + HSIG*ystat.sigma
	 lthresh = yy - LSIG*ystat.sigma
	 badpix = where(y lt lthresh or y gt hthresh)

	 ; implementaion of GROW parameter:
	 if total(badpix) ne -1 then begin
	   if (grow eq 0) then begin
	     y2(badpix) = yy(badpix)
	   endif else begin
	     sz2 = size(badpix,/dimen)
             badpixbot = (badpix - grow > 0)
	     badpixtop = (badpix + grow < lastsub2)
             for k = 0, sz2[0]-1 do begin
               y2(badpixbot[k]:badpixtop[k]) = $ 
                 yy(badpixbot[k]:badpixtop[k])
             endfor
	   endelse
	 endif

         if n_elements(test) ne 0 then oplot,yy,color=colors.green
       endfor
       if n_elements(test) ne 0 then oplot,yy,color=colors.red

       ; implementation of AXIS, part III (output):
       if (axis eq 2) then begin
         image_s[i,*] = y - yy
       endif else begin
         image_s[*,i] = y - yy
       endelse
     endfor


     ; DISPLAY sky-subtracted image, if requested:
     if n_elements(display) ne 0 then begin
        literstat,image_s,imagestat,/silent
	window,1
        plotimage,bytscl(image_s, $
           min=imagestat.median-imagestat.sigma, $
           max=imagestat.median+5.*imagestat.sigma), /preserve
	beep
     endif

     ; WRITE out result, if requested AND input is a fits image:
     if n_elements(write) ne 0 and size(file,/type) eq 7 then begin
       sxaddhist,'DSBACKGROUND applied - extra sky subtraction',hd
       outfile=strmid(file,0,strpos(file,'/',/reverse_search)+1)+$
         title+strmid(file,strpos(file,'/',/reverse_search)+1,100)
       writefits,outfile,image_s,hd
     endif

     return,image_s
 end
