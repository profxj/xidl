;+ 
; NAME:
; x_invertarc   
;    Version 1.0
;
; PURPOSE:
;    Converts an Arc image in the flat frame to the original
;     Uses a C program, is a memory hog
;
; CALLING SEQUENCE:
;   
;   arcimg = x_invertarc( arcimg, map )
;
; INPUTS:
;   arcimg  - Arc image in the flat frame (data or fits file)
;   map     - Map image (data or fits file)
;
; RETURNS:
;   newarc  - Arc image in the Original frame
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  YSTRT  - Row defining the arcfit (default = middle)
;  LINELIST  - Arc line list
;  NSIG  - Sig of RMS from the fit that line must match
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   arcimg = x_invertarc( trcstrct, arcfit, imgsz)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_invertarc, arcimg, map, DBL=dbl

;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'newarc = x_invertarc( arcimg, map, /DBL ) [v1.0]'
    return, -1
  endif 

;  Optional Keywords

; Read in the Arc image + convert to double
  daimg = x_readimg(arcimg, /dscale)
  sz_aimg = size(daimg, /dimensions)

; Read in the map
  dmap = x_readimg(map, /dscale)

; Create the output image
  newarc = dblarr(sz_aimg[0], sz_aimg[1])

; Call the C program

  ndim = 2L
  soname = filepath('libxmath.so', $
                    root_dir=getenv('XIDL_DIR'), subdirectory='/lib')
  retval = call_external(soname, 'invert_arc', $
                         ndim, long(sz_aimg), newarc, daimg, dmap)
  
; Free memory
  delvarx, daimg, dmap

; Return
  if keyword_set( DBL ) then return, newarc else return, float(newarc)
end

