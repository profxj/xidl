;+ 
; NAME:
; x_mkspln   
;   Version 1.0
;
; PURPOSE:
;    Finds the center of a splined 'line'
;
; CALLING SEQUENCE:
;   
;   xcen = x_mkspln(xval, yval, [fracpk])
;
; INPUTS:
;   xval       - x pos
;   yval       - y pos 
;
; RETURNS:
;   xcen      - Center
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
;   xcen = x_mkspln(data)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_mkspln_spline, x

common x_mkspln, splin, spl_x, spl_y

  val = spl_interp(spl_x, spl_y, splin, x, /double)
  return, val[0]

end

;;;;


function x_mkspln, xval, yval, fracpk, SILENT=silent, FORCE=force, EDGES=edges

common x_mkspln

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'xcen = x_mkspln(xval, yval, [fracpk]) [v1.0]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( FRACPK ) then fracpk = 0.3333

;  Deal with endpoints

  nval = n_elements(yval)
  mx = max(yval, imx)
  if mx LT 0. then stop ; shouldn be here

  mn = min(yval) > 0.

  ; Return endpoint as necessary
  if imx EQ 0 OR imx EQ nval-1 then return, xval[imx]

;  Spline

  spl_x = xval
  spl_y = -yval  ; Take negative to find minimum

  splin = spl_init(spl_x, spl_y, /double)

;  Use Golden for the Max

  spmx = x_golden('x_mkspln_spline', xval(imx-1),xval(imx),xval(imx+1))

  mxval = x_mkspln_spline(spmx)

; EDGES

  edges = dblarr(2)

  ; value
  val = -mn + (mxval+mn)*fracpk  ; + because of the sign flip

  ; Left edge
  edges[0] = x_fndspln(spl_x, spl_y, val, splin, IPX=imx, /neg, SILENT=silent)
  ; Right
  edges[1] = x_fndspln(spl_x, spl_y, val, splin, IPX=imx, SILENT=silent)

  ; Catch on edges ;; 
  if min(edges) EQ -1 then begin
      if keyword_set( FORCE ) then return, spmx  $ ; Take max
      else return, -1
  endif

  ; Release memory
  delvarx, spl_x, spl_y, splin
      
  return, (edges[0]+edges[1])/2.
  
end
