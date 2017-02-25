;+ 
; NAME:
; x_centspln   
;   Version 1.1
;
; PURPOSE:
;  Given a 'peak', this routine will find the center of that peak in a
;  non-parameteric manner.  The routine first spliens the peak.  It
;  then steps in from both sides
;  until it hits FRACPK of the peak value of the line.  Finally, the
;  centroid is the midpoint of these two spots.
;
; CALLING SEQUENCE:
;   
;   xcen = x_centspln(xval, yval, [fracpk])
;
; INPUTS:
;   xval       - x values of the peak
;   yval       - y values of the peak
;   [FRACPK]   - Fraction of the peak for centroiding [default=0.3333]
;
; RETURNS:
;   xcen      - Center
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SILENT -- Turn off warnings
;  /FORCE  -- Calculate a centroid even if the edges are non-sensical
;
; OPTIONAL OUTPUTS:
;  EDGES=  -- Values of the spots on the peak corresponding to FRACPK
;
; COMMENTS:
;
; EXAMPLES:
;   xcen = x_centspln(data)
;
;
; PROCEDURES/FUNCTIONS CALLED:
; x_fndspln
; x_golden
;
; REVISION HISTORY:
;   07-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_centspln_spline, x

common x_centspln, splin, spl_x, spl_y

  val = spl_interp(spl_x, spl_y, splin, x, /double)
  return, val[0]

end

;;;;


function x_centspln, xval, yval, fracpk, SILENT=silent, $
                     FORCE=force, EDGES=edges

common x_centspln

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'xcen = x_centspln(xval, yval, [fracpk], /SILENT, ' + $
      '/FORCE, EDGES=) [v1.1]'
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

  spmx = x_golden('x_centspln_spline', xval(imx-1),xval(imx),xval(imx+1))

  mxval = x_centspln_spline(spmx)

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
