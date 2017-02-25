;+ 
; NAME:
; x_calc2dfit   
;   Version 1.1
;
; PURPOSE:
;    Calculates a 2d fit-function at a set of xy pairs.  
;    Restricted to a 2D polynomial for now.
;
; CALLING SEQUENCE:
;   
;   fit = x_calc2dfit(xyval, func, ffit, nx, ny, NRM=, FITSTR=)
;
; INPUTS:
;   xyval      - 2D array [npix, 2] of the x,y values
;   [func]     - String for Fitting function (POLY)
;   [ffit]     - Output from the fitting stuff
;   [nx]       - order of the 2d surface in x
;   [ny]       - order of the 2d surface in y
;
; RETURNS:
;   fit        - Values at each xyval [1D array]
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FITSTR= -- IDL structure for 2D fitting (recommended)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_calc2dfit(xy, 'POLY', ffit, 3, 3)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_calc2dfit, xyval, func, ffit, nx, ny, NRM=nrm, FITSTR=fitstr, $
                      MMEM=mmem

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'fit = x_calc2dfit(xyval, [func, ffit, nx, ny], NRM=, FITSTR=) (V1.1)'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( MMEM ) then mmem = 200.

  ; FIT structure
  if keyword_set( FITSTR ) then begin
      nrm = fitstr.nrm
      ffit = *fitstr.ffit
      nx = fitstr.nx
      ny = fitstr.ny
      func = fitstr.func
  endif

; Normalization
  if keyword_set( NRM ) then begin
      xnrm = 2. * (xyval[*,0] - nrm[0,0])/nrm[0,1]
      ynrm = 2. * (xyval[*,1] - nrm[1,0])/nrm[1,1]
  endif else begin
      xnrm=xyval[*,0]
      ynrm=xyval[*,1]
  endelse
  ; Set to float if xyval is float
  if size(xyval,/type) LE 4 then begin
      xnrm = float(xnrm)
      ynrm = float(ynrm)
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Here we go

  ; Loop for Memory issues  (assume doubles!)

  nval = n_elements(xyval[*,0])
  mval = long((mmem * 1024.^2) / (8. * 2. * nx * ny)) 

  ; Return val as double if xval is double or float otherwise
  val = make_array(nval, type=size(xyval,/type))

  i = 0L
  while(i LE (nval-1)) do begin
      tval = (mval < (nval - i)) < nval
      case func of 
          'POLY': begin
              if nx EQ ny then flg = 0 else flg = nx
              tmp = f2dpoly(0,-1,$
                            XVAL=xnrm[i:i+tval-1],$
                            YVAL=ynrm[i:i+tval-1],$
                            FLG=flg)
              
              ;; Basis vectors
              if keyword_set( BASISARR ) then delvarx, basisarr
              basisarr = f2dpoly(lindgen(tval),nx*ny)
              
              ;; Nice 1 shot
              if tval GT 1 then val[i:i+tval-1] = basisarr # ffit $
                else val[i:i+tval-1] = total(basisarr * ffit)
          end
          else: begin
              print, 'x_calc2dfit: Invalid function declaration! '
              return, -1
          end
      endcase
      i = i+mval
  endwhile

;;;;;;;;;;
; Clean up
  case func of
      'POLY' : begin
          tmp = f2dpoly(0,-2) 
          delvarx, basisarr
      end
      else: 
  endcase
  

  delvarx, xnrm, ynrm
  return, val

end
