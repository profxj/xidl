;+ 
; NAME:
; x_calcfit   
;   Version 1.2
;
; PURPOSE:
;    Calculates a fit-function at a set of x points
;
; CALLING SEQUENCE:
;   
;   fit = x_calcfit(xval,func,ffit, NORD=)
;
; INPUTS:
;   xval       - Value(s) along one dimension
;   [func]       - String for Fitting function (POLY, LEGEND, BSPLIN,
;                GAUSS)
;   [ffit]       - Output from the fitting stuff
;
; RETURNS:
;   fit        - Values at each xval
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   nord=      - Required for LEGEND
;   FITSTR=    - Fit structure (OVERRIDES input values of func, ffit)
;
; OPTIONAL OUTPUTS:
;   NRM=       - Normalization parameters (-1 to 1)
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_calcfit(x, 'POLY', 5)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  POLY_FIT
;  POLY
;  SVDFIT
;  SVLEG
;  BSPLIN
;  GAUSS
;
; REVISION HISTORY:
;   20-Nov-2001 Written by JXP
;   31-Jan-2002 Added NRM, CHEBY func, FITSTR
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_calcfit, xval, func, ffit, NORD=nord, NRM=nrm, FITSTR=fitstr

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'fit = x_calcfit(xval, [func, ffit], NORD=, NRM=, FITSTR=) (V1.2)'
    return, -1
  endif 

;  Optional Keywords

  ; FIT structure
  if keyword_set( FITSTR ) then begin
      nrm = fitstr.nrm
      ffit = *fitstr.ffit
      nord = fitstr.nord
      func = fitstr.func
  endif
  if keyword_set( NRM ) then xnrm = 2. * (xval - nrm[0])/nrm[1] $
  else xnrm=xval

  ; Error checking

  if (func EQ 'LEGEND' OR func EQ 'CHEBY' AND not keyword_set( NORD )) $
    then begin
      print, 'Nord must be set for LEGEND or CHEBY'
      return, -1
  endif


; Function
  case func of 
      'POLY': begin 
          if n_elements(ffit) NE 1 then val = poly(xnrm,ffit) $
          else val = replicate(ffit[0], n_elements(xnrm))
      end
      'LEGEND': val = flegendre(xnrm,nord) # ffit
      'CHEBY': val = fchebyshev(xnrm,nord) # ffit
      'BSPLIN': val = bspline_valu(xnrm, ffit)
      'GAUSS': begin
          case n_elements(ffit) of 
              3: val = ffit[0]*exp(- ((ffit[1]-xnrm)/ffit[2])^2 / 2)
              4: val = ffit[0]*exp(- ((ffit[1]-xnrm)/ffit[2])^2 / 2) $
                + ffit[3]
              5: val = ffit[0]*exp(- ((ffit[1]-xnrm)/ffit[2])^2 / 2) $
                + ffit[3] + ffit[4]*xnrm
              6: val = ffit[0]*exp(- ((ffit[1]-xnrm)/ffit[2])^2 / 2) $
                + ffit[3] + ffit[4]*xnrm + ffit[5]*xnrm^2
              else: begin
                  print, 'In correct number of ffit elements'
                  return, -1
              end
          endcase
      end
;  No function set
      else: begin
          print, 'Invalid function declaration! '
          return, -1
      end
  endcase

  delvarx, xnrm
  ; Return val as double if xval is double or float otherwise
  case size(xval,/type) of
      5 : return, double(val)
      else : return, float(val)
  endcase

end
