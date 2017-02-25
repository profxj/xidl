;+ 
; NAME:
; x_centerpk   
;    Version 1.1
;
; PURPOSE:
;    Finds the center of a peak
;
; CALLING SEQUENCE:
;   center = x_centerpk(xdat, ydat, [fracpk], VAR=, FUNC=, NORD=, HSIG=, LSIG=,
;   IPX=, MINVAL=, FFIT=, /INTER)
;
; INPUTS:
;   xdat       - x Values 
;   ydat       - y Values 
;   [fracpk]     - Fraction of peak relative to the minimum 
;                      (default=1/3 max + min)
;
; RETURNS:
;   center     - center of the peak
;
; OUTPUTS:
;   ffit - Fit parameters
;
; OPTIONAL KEYWORDS:
;   func - Function 
;   VAR  - variance in ydat
;   nord - order number
;   HSIG - upper sigma
;   LSIG - lower sigma
;   IPX  - Allows user to give best guess
;   INTER - Interative fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   center = x_centerpk(xdat, ydat, /inter)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_centerpk, xdat, ydat, fracpk, VAR=var, FUNC=func, $
                     NORD=nord, HSIG=hsig, $
                     LSIG=lsig, IPX=ipx, MINVAL=minval, $
                     FFIT=ffit, INTER=inter, SILENT=silent


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'center = x_centerpk(xdat, ydat, [fracpk], FUNC=, NORD=, HSIG=, LSIG='
    print, '          FFIT=, IPX=, MINVAL=, /INTER) [V1.1]'
    return, -1
  endif 

  if n_elements(xdat) NE n_elements(ydat) then message, 'Wrong size for xdat, ydat'

;  Optional Keywords

  if not keyword_set(fracpk) then fracpk = 0.33333
  if not keyword_set( FUNC ) then func = 'SPLIN'
  if not keyword_set( NORD ) then nord = long( n_elements(xdat) - 1)

; xmin/max

  xmin = min(xdat, MAX=xmax)


  if func EQ 'SPLIN' then begin
     ; Low order fit to take off 'bias' level
      fit = x1dfit(xdat, ydat, FUNC='POLY', NORD=1, HSIG=2., LSIG=3., $
                   INTER=INTER, NRM=nrm)
      yfit = ydat - fit

      ; Center with a spline
      xcen = x_centspln(xdat, yfit, fracpk)
      return, xcen

  endif else begin
; FIT
      fitstr = { fitstrct }
      fitstr.func = func
      fitstr.nord = nord
      if keyword_set( HSIG ) then fitstr.hsig = hsig
      if keyword_set( LSIG ) then fitstr.lsig = lsig

      fit = x1dfit(xdat, ydat, FITSTR=fitstr, INTER=INTER)

;  Find min
      if not keyword_set(MINVAL) then minval = min(fit) > min(ydat)

;  Calculate new fit
      if n_elements(ydat) LT 1000 then begin
          xfit = findgen(1000)*(xmax-xmin)/1000. + xmin
          fit = x_calcfit(xfit, FITSTR=fitstr)
      endif

; Max
      if keyword_set( IPX ) then begin
          step = round(1000./n_elements(ydat))
          sep = abs(xfit-ipx)
          minsep = min(temporary(sep), itmp)
          maxval = max(fit[(itmp-step)>0:(itmp+step)<(n_elements(ydat)>1000)], imax)
          imax = round(imax+itmp-step)
      endif else maxval = max(fit, imax)
      
; Step to the left of peak
      
      val = (maxval - minval)*fracpk + minval
      lft = x_fndfitval(val, fitstr, xfit, fit, IPX=imax, /NEG, SILENT=silent)
      rght = x_fndfitval(val, fitstr, xfit, fit, IPX=imax, SILENT=silent)
      
  ; Catch the error
      if lft EQ -1 OR rght EQ -1 then return, -1 else $
        return, (lft+rght)/2.
  endelse
      
end
