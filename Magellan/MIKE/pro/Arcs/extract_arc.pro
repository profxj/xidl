;+ 
; NAME:
; extract_arc
;     Version 1.1
;
; PURPOSE:
;  Straighten each order indivdiually (mike_rectify) and then extract
;  a boxcar down the center of each order by taking the average flux
;  in two regions (each side of the center) of width 1/3 the order
;  width.
;
; CALLING SEQUENCE:
;   
;  flux = extract_arc( arc_img, ordr_str )
;
; INPUTS:
;   arc_img  -  2D Arc image
;   ordr_str -  Order structure
;
; RETURNS:
;   flux     -  1D spectrum down the center of each order
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  flux = extract_arc( arc_img, ordr_str )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_rectify
;
; REVISION HISTORY:
;   ??--2004 Written by SB
;-
;------------------------------------------------------------------------------

function extract_arc, arc_img, ordr_str

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'flux = extract_arc(arc_img, ordr_str ) [v1.1] ' 
      return, -1
  endif

  nord = n_elements(ordr_str)
  sz_arc = size(arc_img,/dimen)
  arc_flux = fltarr(sz_arc[1], nord)
  
  for i=0, nord-1 do begin
      rect_arc = mike_rectify(arc_img, ordr_str[i].lhedg, $
                              ordr_str[i].rhedg, /nocorrect)
      ncols = (size(rect_arc))[1]
      n3 = fix(ncols/3)
      n33 = fix(n3/3)
      
;      lower = djs_median(rect_arc[n33:n33+n3-1,*],1)
;      upper = djs_median(rect_arc[ncols-n33-n3:ncols-n33-1,*],1)
;      arc_flux[*,i] = 0.5*(upper+lower) 

      arc_flux[*,i] = rect_arc[ncols/2,*]

  endfor
  return, arc_flux
end
       
      

      
