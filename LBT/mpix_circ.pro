;+ 
; NAME:
; mpix_circ
;   Version 1.1
;
; PURPOSE:
;  Finds all pixels within a circle on an image.
;
; CALLING SEQUENCE:
;   pval = mpix_circ(x0, y0, radius, COUNT=, XSIZE=, YSIZE=, ANNULUS=,
;   WIDTH=,/silent)
;
; INPUTS:
;   x0 -- Circle x position
;   y0 -- Circle y position
;   radius -- Radius of the circle
;   XSIZE -- x image size
;   YSIZE -- y image size
;   WIDTH -- the width of the annulus
; RETURNS:
;  pval -- pixels in the circle:  1D array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /ANNULUS   --compute pixels in annulus. Radius becomes inner
;                radius. Oute radius is set by radius+width
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
; Oct,2009 Written by MF
;-    
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mpix_circ, x0, y0, radius, COUNT=count, XSIZE=xsize, YSIZE=ysize, $
                    ANNULUS=ANNULUS, WIDTH=width, silent=silent

; Finds all pixels (integer values >= 0) within a circle

; 
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'arrpix = mpix_circ(x0, y0, radius, COUNT=, XSIZE=, YSIZE=, ANNULUS=, WIDTH=)'
    return, -1
  endif 

  if not keyword_set(XSIZE) or not keyword_set(YSIZE) then begin
     splog, "Image size not specified!"
     return, -1
  endif
  
  if keyword_set(annulus) and not keyword_set(width) then begin
     splog, "Set annulus width!"
     return, -1
  endif 


;test if inside the image
  if not keyword_set(annulus) then begin 
     xmax=x0+radius
     xmin=x0-radius
     ymax=y0+radius
     ymin=y0-radius
  endif else begin
     xmax=x0+radius+width
     xmin=x0-radius-width
     ymax=y0+radius+width
     ymin=y0-radius-width
  endelse

  new_rad=radius

  ;resize if needed 
  if not keyword_set(annulus) then begin 

     if(xmax GT xsize) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        diff=xmax-xsize+1
        new_rad=radius-diff
     endif
     
     if(ymax GT ysize) then begin
       if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
       diff=ymax-ysize+1
        new_rad=radius-diff
     endif
     
     if(xmin LT 0) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        new_rad=radius-abs(xmin)-1
     endif
     
     if(ymin LT 0) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        new_rad=radius-abs(ymin)-1
     endif
   
  endif else begin

     
     if(xmax GT xsize) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        diff=xmax-xsize+1
        new_rad=radius-diff-width
     endif
     
     if(ymax GT ysize) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        diff=ymax-ysize+1
        new_rad=radius-diff-width
     endif
     
     if(xmin LT 0) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        new_rad=radius-abs(xmin)-1-width
     endif
     
     if(ymin LT 0) then begin
        if ~keyword_set(silent) then $
           splog, "Radius exeeds image dimension.. shrink to img size"
        new_rad=radius-abs(ymin)-1-width
     endif
  endelse

  ;finally, if radius becomes negative or 0, set to 1
  if(new_rad LT 1) then new_rad=1


;make index image
  
  ;x
  xarray=make_array(xsize,/index)
  ximage=rebin(xarray,xsize,ysize)
   
  ;y
  yarray=make_array(ysize,/index)
  yimage=rebin(transpose(yarray),xsize,ysize)

  ;square
  xyimage=(ximage-x0)^2+(yimage-y0)^2
  ;atv, xyimage, /block
  
  ;find indexes
  if not keyword_set(annulus) then index=where(xyimage LT new_rad^2, count) $
  else index=where(xyimage LT (new_rad+width)^2 and xyimage GT new_rad^2, count)
  
  
  ;xyimage[index]=xyimage[index]+1D20
  ;atv, xyimage, /block

  return, index
  
  
  
end
          
