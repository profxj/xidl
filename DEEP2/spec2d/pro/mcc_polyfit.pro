;+
;
; NAME
;      mcc_polyfit.pro
;
; PURPOSE
;      The mcc_polyfit procedure performs an ordinary least-squares
;      fit to a polynomial of the form y = a0 + a1*x + a2*x^2 + a3*x^3
;      + ... + an*x^n where any of the orders can be omitted (taken to
;      be zero). That is, the procedure will fit a polynomial of the
;      form y = a0 + a3*x^3. Note that all of the independent measures
;      (the x's) are taken to be identical. The procedure will not fit
;      a function y = a0*x^2 + a1*v^3 + a2*t^(1/2). However,
;      mcc_polyfit.pro will perform multivariate chi-squared fitting
;      by use of the keyword "disp" which allows each measurement of
;      ydata to be weighted differently. The xdata values are assumed
;      to have no uncertainties in their values.
;
; SYNTAX
;      mcc_polyfit, xdata, ydata, pow, [a=a, cov=cov, disp=disp, 
;                    yfit=yfit, /plot, /verb]
;
; INPUTS
;	xdata = a vector containing M independent values of the quantity x. 
;	ydata = a vector containing M measurements of the observed 
;		quantity y corresponding (one-to-one) to the xdata values.
;	pow = a vector containing the powers of the polynomial. More 
;	      specifically, if the fit polynomial is y = a0 + a1*x^1 + a2*x^3
;	      then the vector pow is given by [0, 1, 3].
;
; KEYWORDS
;	a = a variable that will be set to the vector containing the
;	    fit values of the parameters a0, a1, a2, ... Note that a 
;           will have the same length as pow. 
;	cov = a variable that will be set to an array of dimension np
;	      x np where np is the number of parameters to be fit. The 
;             cov array will be filled with the covariance matrix as 
;             defined by the fit.
;	disp = a vector containing the uncertainties in the measurement
;	       of each y point. Note that disp can be given as a scalar
;	       in which case the uncertainties in the y points are taken 
;	       to be the same for all y points. If this keyword is not set,
;	       then each point is simply weighted evenly (disp =
;	       1.0). Recall that in multivariate chi-squared fitting
;	       the weight w is given by w = 1/(disp^2). 
;       red_chisq = a variable that will be set to the reduced
;                   chi-squared value for the fit.
;       yfit = a variable that will be set to a vector containing the
;              fit. That is, the variable will filled with the y values
;              corresponding to the polynomial fit.
;       plot = if this keyword is set, then a plot will be generated
;              showing the data (ydata vs. xdata) and the polynomial
;              fit will be overplotted. Also the residuals will be
;              plotted.
;       verb = if this keyword is set, then output is printed to the
;              terminal in which IDL is running. The value for the
;              parameters along with the errors in each. 
;
; OUTPUTS
;	The inputs (a), (cov), (red_chisq), and (yfit) are filled
;	according to the least-squares fit routine. 
;
; PROCEDURES CALLED
;	None
;
; EXAMPLES
;	(1) Given a data set of x values and a y data set of noisey values.
;		IDL> x = findgen(1000)
;		IDL> y = 10. + 2.*x 
;		IDL> y = y + (100. * randomn(seed, 1000))
;	Then fit the data with linear polynomial
;		IDL> a = findgen(2)
;		IDL> cov = fltarr(2,2)
;		IDL> mcc_polyfit, x, y, [0,1], a, cov
;		IDL> yfit = a[0] + a[1]*x
;		IDL> PLOT, x, y
;		IDL> OPLOT, x, yfit
;       (2) Create a data set with noise according to a double-sided
;       exponential.
;               IDL> x = findgen(4999) + 1.0
;               IDL> y = 5.0 + 2.5*x^(0.25) + 3.25*x^(-0.5)
;               IDL> y = y + mcc_noisegen(4999, /dbl)
;               IDL> mcc_polyfit, x, y, [0,0.25,-0.5], pars, var, $
;               IDL> yfit=outvec, /plot, /verb 
;
; HISTORY
;	Created March 15, 2002 by mcc.
;       Revised July 27, 2002 by mcc - revised documentation.
;
;-

PRO mcc_polyfit, xdata, ydata, pow, a=a, cov=cov, disp=disp, yfit=yfit, $
                 plot=plot, verb=verb, red_chisq=red_chisq

;;; CHECK THAT ENOUGH PARAMETERS WERE SUPPLIED 
  IF N_PARAMS(0) LT 3 THEN BEGIN
    PRINT, 'CALLING SEQUENCE: mcc_polyfit, xdata, ydata, pow'
    PRINT, 'KEYWORD PARAMETERS: a= , cov= , disp= , yfit= , /plot, /verb'
    RETALL
  ENDIF

;;; CHECK TO SEE THAT LENGTHS OF ARGUMENTS ARE OKAY
  IF N_ELEMENTS(xdata) NE N_ELEMENTS(ydata) THEN BEGIN
    PRINT, 'The arguments (xdata) and (ydata) must have the same length, stupid!'
    RETALL
  ENDIF

;;; CHECK TO SEE IF KEYWORD disp IS SET
  IF N_ELEMENTS(disp) LT 1 THEN disp = 1.0

  IF N_ELEMENTS(disp) GT 1 AND N_ELEMENTS(disp) NE N_ELEMENTS(xdata) THEN BEGIN
    PRINT, 'IF YOU SET THE disp INPUT, THEN IT MUST HAVE THE ' + $
      'SAME LENGTH AS xdata AND ydata!'
    RETALL
  ENDIF

  IF n_elements(disp) EQ 1 THEN disp=disp[0]

;;; DETERMINE THE NUMBER OF MEASUREMENTS (M) 
  M = n_elements(xdata)

;;; DETERMINE THE NUMBER OF PARAMETERS TO FIT (np)
  np = n_elements(pow)

;;; ALTERNATE WAY TO POPULATE MATRICES - IT IS FASTER. 
  Xtx = fltarr(np, M)
  FOR i=0, np-1 DO BEGIN
    Xtx[i,*] = xdata^pow[i] / disp
  ENDFOR
  Ytx = fltarr(1,M) + ydata / disp

;;; CALCULATE alpha AND beta MATRICES
  alpha = transpose(Xtx) ## Xtx
  beta = transpose(Xtx) ## Ytx

;;; DETERMINE COVARIANCE MATRIX 
  cov = invert(alpha)

;;; DERIVE PARAMETER MATRIX 
  a = cov ## beta

;;; DERIVE THE NORMALIZED COVARIANCE MATRIX
  diag = cov[(np+1) * indgen(np)]
  ncov = cov / sqrt(diag##diag)

;;; GENERATE A VECTOR yfit CONTAINING THE FIT 
;;; VALUES. THAT IS, A VECTOR CONTAINING THE 
;;; VALUES OF THE POLYNOMIAL FIT. 
  yfit = findgen(M) * 0.0
  FOR i=0, np-1 DO BEGIN
    yfit = yfit + a[i]*xdata^(pow[i]) 
  ENDFOR

;;; CALCULATE THE REDUCED CHI-SQUARED VALUE
;delYtx = Ytx-(Xtx##a)
;red_chisq = transpose(delYtx) ## delYtx/(M-np)
;;; ALTERNATE WAY TO CALCULATE REDUCED CHI-SQUARED
  red_chisq = total( (ydata - yfit)^2 / disp^2 ) / ( M - np)

;;; CHECK IF KEYWORD /plot IS SET. IF SET, THEN PLOT 
;;; THE DATA (ydata vs. xdata) AND OVERPLOT THE FIT. 
  IF KEYWORD_SET(plot) THEN BEGIN
    !P.MULTI=[0,1,2]
    WINDOW, 7
    PLOT, xdata, ydata
    OPLOT, xdata, yfit, LINESTYLE=2, THICK=2, COLOR=getcolor('red')
    PLOT, ydata-yfit, YTITLE="Residuals", PSYM=4
    OPLOT, [0,M], [0,0], LINESTYLE=2, THICK=2, COLOR=getcolor('red')
  ENDIF


;;; CHECK IF KEYWORD /verb IS SET. IF SET, THEN PRINT 
;;; OUTPUT TO SCREEN. PRINT PARAMETER VALUES, UNCERTAINTIES 
;;; IN PARAMETERS, NORMALIZED COVARIANCE MATRIX, REDUCED 
;;; CHI-SQUARED VALUE.
  IF KEYWORD_SET(verb) THEN BEGIN
    print
    print, 'FIT RESULTS FROM mcc_polyfit.pro:'
    print, '-----------------------------------'
    FOR i=0,np-1 DO BEGIN    
      print
      print, 'Parameter a' + STRING(i, FORMAT='(I1)') + ' = ' + STRING(a[i])
      print, 'Uncertainty in a' + STRING(i, FORMAT='(I1)') + ' = ' + $
        STRING( SQRT(cov[i,i]) )
      print, '-----------------------------------'
    ENDFOR
    print, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print, 'NOTE: the uncertainties in the fit parameters quoted above are '
    print, 'accurate for an ordinary least-squares fit (no weighting used - '
    print, '"disp" argument not specified by user). If performing a multivariate '
    print, 'chi-squared fit, then these uncertainties should be accurate - but '
    print, 'only if the reduced chi-squared value (see below) is close to unity!!!'
    print, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print, '-----------------------------------'
    print
    print, 'The Normalized Covariance Matrix: '
    print, ncov
    print, '-----------------------------------'
    print 
    IF N_ELEMENTS(disp) GT 1 OR (disp[0] NE 1.0) THEN BEGIN 
      print, 'The Reduced Chi-Squared Value: '
      print, red_chisq
      print, '-----------------------------------'
      print, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print, 'NOTE: if reduced chi-square value is NOT ~ 1, then multiply ' 
      print, 'estimated measurment errors (disp parameter) by the reduced '
      print, 'chi-squared value and redo fit!!!' 
      print, 'Note that when the reduced chi-squared is roughly unity, then '
      print, 'variances quoted above should agree with those below.'
      print, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print, '-----------------------------------'
      print
      print, 'Parameter Uncertainties (multivariate chi-squared fit): '
      print, SQRT(red_chisq*diag)
      print, '-----------------------------------'
      print
    ENDIF

  ENDIF



END


