;+ 
; NAME:
; x_fit   
;   Version 1.2
;
; PURPOSE:
;    Fits a function to a set of x,y data
;
; CALLING SEQUENCE:
;   
;   fit = x_fit(xdat,ydat,func,nord,sig=,reg=) 
;
; INPUTS:
;   xdat       - Values along one dimension
;   ydat       - Values along the other
;
; RETURNS:
;   fit        - Values of the fit at each xdat value
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   sig=       - Errors in the yval points
;   IVAR=      - Inverse variance of the ydat values
;   func=      - String for Fitting function (POLY, LEGEND, BSPLIN,
;                GAUSS)
;   nord=      - Order of the fit
;   reg=       - Regions of data to fit
;   msk=       - Mask  (0 = Do NOT include)
;   guess=     - First guess for GAUSS routine
;   FITSTR=    - 1D Fit structure used to set many of the above
;                parameters (recommended)
;   FLG_BSP=   - Flag controlling BSPLINE options (1: nord = everyn)
;   /NONRM     - Do not normalize the xdat from -1 to 1
;
; OPTIONAL OUTPUTS:
;   ffit     - Functional form
;   NRM      - Normalization numbers for xdat :: dblarr(2)
;   RMS      - RMS of the fit
;
; COMMENTS:
;
; EXAMPLES:
;   fit = x_fit(x, y, 'POLY', nord=5)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  POLY_FIT
;  POLY
;  SVDFIT
;  SVLEG
;  FUNC_FIT
;
	; REVISION HISTORY:
;   20-Nov-2001 Written by JXP
;   31-Jan-2002 Revised significantly.  Added fit structure, CHEBY
;   27-Feb-2002 Added ivar option, am using func_fit for legendre
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fit, xdat, ydat, FUNC=func, NORD=nord, SIG=sig, REG=reg, FFIT=ffit, $
                MSK=msk, GUESS=guess, NRM=nrm, FITSTR=fitstr, NONRM=nonrm, $
                RMS=rms, IVAR=ivar, FLG_BSP=flg_bsp, SILENT=silent

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'fit = x_fit(xdat, ydat, FUNC=func, NORD=, '
    print, '        SIG=, REG=, FFIT=, MSK=msk, GUESS=, NRM=, FITSTR=)'
    print, '        /NONRM, IVAR= [v1.2]'
    return, -1
  endif 

  igood = lindgen(n_elements(xdat))

;  Optional Keywords

  ; Fit structure
  if keyword_set(FITSTR) then begin
      func = fitstr.func
      if keyword_set( FLG_BSP ) then begin
          case flg_bsp of
              0: nord = fitstr.nord
              else: begin
                  nord = flg_bsp
                  everyn = fitstr.nord
              end
          endcase
      endif else nord = fitstr.nord
  endif

  ; func, nord
  if not keyword_set( FUNC ) then func = 'POLY'
  if not keyword_set( nord ) then nord = 0

  ; mask
  if keyword_set( MSK ) then igood = where(msk NE 0)

; Normalize xdat

  if not keyword_set( NONRM ) then begin
      nrm = dblarr(2)
      mnx = min(xdat, MAX=mxx)
      nrm[0] = 0.5 * (mnx + mxx)
      nrm[1] = mxx - mnx
      xnrm = 2. * (xdat - nrm[0])/nrm[1]
  endif else xnrm = xdat

;   Pre-set Regions
  if keyword_set( REG ) then begin
      regpt = x_fndreg(xdat[igood], reg, NPNT=npnt)  ; Find the regions 
      if npnt EQ 0 then return, -1
      igood = temporary(igood[regpt])
  endif 

; IVAR
  if keyword_set( IVAR ) then begin
      a = where(ivar[igood] GT 0.)
      igood = igood[a]
  endif


  xfit = xnrm[igood] ; Normalized
  yfit = ydat[igood]
                                
; Function

  case func of 
      'POLY': begin
          ; sigma array
          if not keyword_set( SIG ) then begin
              if keyword_set( IVAR ) then begin
                  a = where( IVAR GT 0.)
                  sig = fltarr(n_elements(ydat))
                  sig[a] = 1./ sqrt(ivar[a])
               endif            ; else sig = fltarr(n_elements(ydat))+1.  ;; JXP
          endif
          if keyword_set(SIG) then sfit = sig[igood]
          ; Fit
         ffit = poly_fit(xfit,yfit,nord, measure_errors=sfit,$
                          /double, status=status)
;         if status EQ 0 OR status EQ 2 then begin
          if status EQ 0 then begin
              if nord EQ 0 then fit = replicate(ffit, n_elements(xnrm)) $
              else fit = poly(xnrm,ffit) 
          endif else begin
              print, 'POLYFIT Failed!', status
              return, -1
          endelse
      end
      'LEGEND': begin
          if keyword_set(IVAR) then ivfit = ivar[igood]
          if not keyword_set( IVAR ) AND keyword_set(SIG) then begin
              gdsig = where( sig[igood] GT 0. )
              ivfit = dblarr( n_elements(igood) )
              ivfit[gdsig] = 1. / (sig[igood[gdsig]])^2
          endif
          ffit = func_fit( xfit, yfit, nord, func='flegendre', invvar=ivfit)
          fit = flegendre(xnrm, nord) # ffit
      end
      'CHEBY': begin
          if not keyword_set( REG ) AND NOT keyword_set( MSK ) then begin
              ffit = svdfit(xfit,yfit,nord, measure_errors=sfit,$
                            function_name='fchebyshev', $
                            /double, yfit=fit)
          endif else begin
              ffit = svdfit(xfit, yfit, nord, measure_errors=sfit, $
                            function_name='fchebyshev', $
                            /double, singular=singular)
              if(singular EQ 0) then fit = fchebyshev(xnrm,nord) # ffit $
              else begin
                  print, 'Problem with svdfit: Order too high?', singular
                  return, -1
              endelse
          endelse
      end
      'BSPLIN': begin
          ffit= bspline_iterfit(xfit, yfit, $
                                nbkpts=nord, yfit = bsplinefit, $
                                everyn=everyn, /SILENT)
          if (size(ffit,/tname) EQ "INT") then begin
              print, 'Problem with bspline_fit: Order is probably too high'
              return, -1 
          endif else fit = bspline_valu(xnrm, ffit)
      end
      'GAUSS': begin
          if keyword_set( NORD ) then nord = 3 > nord < 6
          fit = gaussfit(xfit, yfit, ffit, ESTIMATE=guess, $
                                  NTERMS=nord)
      end
;  No function set
      else: begin
          print, 'x_fit: Invalid function declaration! '
          return, -1
      end
  endcase

; RMS

  rms = sqrt(total((fit-yfit)^2)/(n_elements(xfit)-1.))
  if arg_present(fitstr) then fitstr.rms = rms

  ; Release normalized data
  delvarx, xnrm, xfit, yfit, sfit

; Write to fit structure
  if arg_present(fitstr) then begin
      if not keyword_set( NONRM ) then fitstr.nrm = nrm
      if ptr_valid(fitstr.ffit) EQ 1 then ptr_free, fitstr.ffit
      fitstr.ffit = ptr_new(ffit)
  endif


  ; Return val as double if xval is double or float otherwise
  if size(ydat,/type) EQ 5 OR  size(xdat,/type) EQ 5 then $
    return, double(fit) $
  else return, float(fit)

end
