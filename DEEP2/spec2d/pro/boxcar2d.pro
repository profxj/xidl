function boxcar2d, image, smoothdimensions,  ALL_PIXELS=all_pixels,DOUBLE=DOUBLE,NO_FT=NO_FT
; NAME:
;      BOXCAR2D (derived from FILTER_IMAGE)
;
; PURPOSE:
;       Identical to SMOOTH but handle edges and allow rectangular boxes
; EXPLANATION:
;       Computes the average  of pixels in moving box,
;       replacing center pixel with the computed average,
;       (using the Goddard CONVOLVE() function).
;       The main reason for using this function is the options to
;       also process the pixels at edges and corners of image, and,
;       to use non-square kernels.
;
; CALLING SEQUENCE:
;       Result = boxcar2d( image, smoothdimensions, /ALL_PIXELS)

;
; INPUT:
;       image = 2-D array (matrix)
;       smoothdimensions = scalar or 2 element array of the smoothing
;       size to use in X and Y
;
; OPTIONAL INPUT KEYWORDS:
;   
;       /ALL_PIXELS causes the edges of image to be filtered as well.   This
;               is accomplished by reflecting pixels adjacent to edges outward
;               (similar to the /EDGE_WRAP keyword in CONVOL).
;               Note that this is a different algorithm from the /EDGE_TRUCATE; 
;               keyword to SMOOTH or CONVOL, which duplicates the
;               nearest pixel.   This option is set by default. 
;       /DOUBLE causes the convolution to be performed in double precision.
;

; RESULT:
;       Function returns the smoothed image.
; EXTERNAL CALLS:
;       function convolve
;
; PROCEDURE:
;       If both /ALL_PIXELS (or /ITERATE)  keywords are set then
;       create a larger image by reflecting the edges outward, then call the 
;       IDL MEDIAN() or SMOOTH() function on the larger image, and just return;       the central part (the original size image).
;
;       NAN values are recognized during calls to MEDIAN() or SMOOTH(), but 
;       not for convolution with a Gaussian (FWHM keyword supplied). 
; HISTORY:
;       Written, 1991, Frank Varosi, NASA/GSFC.
;       FV, 1992, added /ITERATE option.
;       FV, 1993, added FWHM_GAUSSIAN= option.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use /EVEN call to median, recognize NAN values in SMOOTH 
;                  W. Landsman   June 2001
;       Added PSF keyword,   Bjorn Heijligers/WL, September 2001
;-

  if N_params() LT 1 then begin
      print,'Syntax - Result = boxcar2d( image, smoothdimensions)'
      return, -1
  endif

        sim = size( image )
        Lx = sim[1]-1
        Ly = sim[2]-1

        if (sim[0] NE 2) OR (sim[4] LE 4) then begin
                message,"input must be an image (a matrix)",/INFO
                return,image
           endif
           
        if n_elements(all_pixels) eq 0 then all_pixels=1

        if n_elements(smoothdimensions) eq 2 then begin
            box_x=smoothdimensions[0]
            box_y=smoothdimensions[1]
        endif else begin
            box_x=smoothdimensions
            box_y=smoothdimensions
        endelse

                npix = box_x;( 3 * box_x[ 0: ( (N_elements( box_X )-1) < 1 ) ] ) > 3
                npix = 2 * fix( npix/2 ) + 1    ;make # pixels odd.
                box_x = box_x > max( [npix] )

                npix = box_y;( 3 * box_y[ 0: ( (N_elements( box_y )-1) < 1 ) ] ) > 3
                npix = 2 * fix( npix/2 ) + 1    ;make # pixels odd.
                box_y = box_y > max( [npix] )



        if (box_x LT 3 AND box_y LT 3) then return, image

        if keyword_set(all_pixels) then begin
                
                box_x = fix( box_x )
                box_y = fix( box_y )

                radiusx = (box_x/2) > 1
                radiusy = (box_y/2) > 1
                radiusx=radiusx
                radiusy=radiusy
                Lxr = Lx+radiusx
                Lyr = Ly+radiusy
                rrx = 2*radiusx
                rry = 2*radiusy
                imf = fltarr( sim[1]+rrx, sim[2]+rry )
                imf[radiusx,radiusy] = image              ; reflect edges outward


                                                        ; to make larger image.
                imf[  0,0] = rotate( imf[radiusx:rrx,*], 5 )      ;Left
                imf[Lxr,0] = rotate( imf[Lx:Lxr,*], 5 )         ;right
                imf[0,  0] = rotate( imf[*,radiusy:rry], 7 )      ;bottom
                imf[0,Lyr] = rotate( imf[*,Ly:Lyr], 7 )         ;top

          endif else begin

                radiusx=0
                radiusy=0
                imf = image
           endelse


           if n_elements(no_ft) eq 0 then no_ft=0
           if keyword_set(double) then begin
               kernel=dblarr(box_x,box_y)+1
               kernel=kernel/total(kernel)
               imf = convolve( double(imf),  NO_FT=no_ft, kernel) 
           endif else begin
               kernel=fltarr(box_x,box_y)+1
               kernel=kernel/total(kernel)
               imf = convolve( (imf),  NO_FT=no_ft, kernel) 
           endelse

    if radiusx GT 0 AND radiusy gt 0 then $
                return, float(imf[ radiusx:(Lx+radiusx), radiusy:(Ly+radiusy) ]) $
           else return, float(imf)
end
