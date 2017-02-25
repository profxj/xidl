;+ 
; NAME:
; x_invertarc   
;    Version 1.1
;
; PURPOSE:
;    Converts an Arc image in the flat frame to the original
;     Uses a C program, is a memory hog
;
; CALLING SEQUENCE:
;   
;   newarc = x_invertarc( arcimg, map )
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
;  /DBL  - Use double precision
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   newarc = x_invertarc( arcimg, map, /DBL)
;
; PROCEDURES/FUNCTIONS CALLED:
;  C code: invert_arc
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
             'newarc = x_invertarc( arcimg, map, /DBL ) [v1.1]'
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
  soname = filepath('libxmath.' +idlutils_so_ext(), $
                    root_dir=getenv('XIDL_DIR'), subdirectory='/lib')
  retval = call_external(soname, 'invert_arc', $
                         ndim, long(sz_aimg), newarc, daimg, dmap)
  
; Free memory
  delvarx, daimg, dmap

; Return
  if keyword_set( DBL ) then return, newarc else return, float(newarc)
end

