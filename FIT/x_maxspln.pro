;+ 
; NAME:
; x_maxspln   
;   Version 1.1
;
; PURPOSE:
;    Finds the maximum in a set of data using a spline.  The program
;    simply splines the data and then finds the max.
;
; CALLING SEQUENCE:
;   
;   maxsp = x_maxspln( xval, yval, EDGES=, EDGVAL=, SPMX=)
;
; INPUTS:
;   xval       - x pos
;   yval       - y pos (values to be maximized)
;
; RETURNS:
;   maxsp      - Maximum of the values
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  EDGES   -- x-value Edges of peak corresponding to 0.333
;  EDGVAL  -- y-value Edges of peak corresponding to 0.333
;  SPMX    -- Maximum in the spline
;
; COMMENTS:
;
; EXAMPLES:
;   maxsp = x_maxspln(x, y)
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_golden
;
; REVISION HISTORY:
;   29-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_maxspln_spline, x

common x_maxspln, splin, spl_x, spl_y

  val = spl_interp(spl_x, spl_y, splin, x, /double)
  return, val[0]

end

;;;;


function x_maxspln, xval, yval, EDGES=edges, EDGVAL=edgval, SPMX=spmx

common x_maxspln

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'maxsp = x_maxspln(xval, yval, EDGES=, EDGVAL=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords

;  Deal with endpoints

  nval = n_elements(yval)

  mx = max(yval, imx)

  ; Return endpoint as necessary
  if imx EQ 0 OR imx EQ nval-1 then return, xval[imx]

;  Spline

  spl_x = xval
  spl_y = -yval  ; Take negative to find minimum

  splin = spl_init(spl_x, spl_y, /double)

  if not finite(splin[0]) then begin
      ;; Spline crashed!  Not sure why..
      print, 'x_maxspln: Bad spline!  Returning min'
      spmx = xval[imx]
      return, mx
  endif


;  Use Golden for the Max

  spmx = x_golden('x_maxspln_spline', xval[imx-1],xval[imx],xval[imx+1])


; EDGES

  if arg_present( EDGES ) then begin
      edges = dblarr(2)

      ; mxval
      if not keyword_set( mxval ) then mxval = -x_maxspln_spline(spmx)
      ; edgval
      if not keyword_set( EDGVAL ) then edgval = 0.333

      ; value
      val = -mxval*edgval

      ; Left edge
      edges[0] = x_fndspln(spl_x, spl_y, val, splin, IPX=imx, /neg)
      ; Right
      edges[1] = x_fndspln(spl_x, spl_y, val, splin, IPX=imx)

 endif

  ; Release memory
;  delvarx, spl_x, spl_y, splin
      
;  if arg_present( MXVAL ) then mxval = -x_maxspln_spline(spmx)
  return, -1.*x_maxspln_spline(spmx)
  
end
