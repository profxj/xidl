;+
;
; NAME
;      mcc_gauss1dfit.pro
;
; PURPOSE
;      The mcc_gauss1dfit procedure fits a single 1-dimensional
;      Gaussian of the general form ydata = C + D*xdata +
;      E*e^{-(xdata-u)^2/2*s^2} to the data set (xdata, ydata). 
;
;
; SYNTAX
;	mcc_gauss1dfit, xdata, ydata, [a=a, cov=cov, yfit=yfit, 
;                       weight=weight, swidth=swidth, inum=inum, 
;                       pdif=pdif, /slope, /cont, guesscont=guesscont, 
;                       /hanning, nhan=nhan, /verbose, /movie, /plot,
;                       /integra] 
;
; INPUTS
;	xdata = a vector containing the data values for the
;	        independent variable x in the Gaussian function 
;               y = C + D*x + E*e^{-(x-u)^2/2*s^2}. 
;       ydata = a vector containing the measurements of the dependent 
;               variable y. 
;
; KEYWORDS
;       a = if this keyword is set, then the variable passed as a will
;           be set as a vector containing the values of the fit
;           parameters in the order a = [C, D, E, u, s] or 
;           a = [C, E, u, s], or a = [E, u, s] depending on the form
;           of the fitting function.
;       cov = if this keyword is set, then the variable passed as cov will
;             be set as an array containing covariance of the fit
;             parameters. That is, the variable will be set as the
;             covariance matrix. 
;	inum = the number of iterations to perform. Default value is 100.
;	swidth = the width of the smoothing window to be sent to the 
;		 boxcar smoothing precedure SMOOTH. The boxcar
;		 smoothing function is the default smoothing kernel
;		 used to smooth the ydata and find the Gaussian peak
;		 in the ydata values. Note there is an option for a
;		 hanning kernel.
;	yfit = a vector to be filled with the fit y values.
;	slope = if this keyword is set, then the Gaussian fit includes a 
;		fit to the continuum (w/ slope) C + D*xdata as in the full 
;		form ydata = C + D*xdata + E*e^{-(xdata-u)^2/2*s^2}. This 
;		keyword and the keyword \cont should NOT be both set in a 
;		call to mcc_gauss1dfit. If this keyword is set, then 
;		the keyword guesscont MUST be supplied.
;	cont = if this keyword is set, then the Gaussian fit includes a 
;	       fit to the continuum (w/o slope) C as in the full 
;	       form ydata = C +  E*e^{-(xdata-u)^2/2*s^2}. If this 
;	       keyword is NOT set and the /slope keyword is also not set, 
;	       then the data is simply fit to a function of form 
;	       ydata = E*e^{-(xdata-u)^2/2*s^2}. If this keyword is set, then 
;	       the keyword guesscont MUST be supplied.
;       movie = if this keyword is set, then a plotting window (number
;               4 and 5) is opened up on your display. In that window,
;               the iterative fits (as the fit converges) will be
;               plotted as a movie. 
;       plot = if this keyword is set, then a plotting window (number
;              6) will be opened on your display. After fitting the data
;              with the Gaussian of desired form, the data (ydata
;              vs. xdata) will be plotted with the Gaussian fit
;              overlayed. In a second plot (in the same window), the
;              residuals will also be plotted.
;       verbose = if this keyword is set, then the value of each fit
;                 parameter is printed to the screen (the terminal
;                 running IDL) along with the uncertainty (sigma) in
;                 the fit of the parameter.
;       weight = the weight keyword allows the user to apply weights
;                to the ydata values. If the weight keyword is not 
;                passed to the procedure, then a uniform weighting is
;                assumed. If specified the weight should be a vector
;                of the same length as xdata and ydata.
;       hanning = if this keyword is set, then a hanning kernel is
;                 applied to the ydata to remove noise and find the
;                 peak in the Gaussian. If this keyword is not set,
;                 then the default smoothing function is a boxcar.   
;       nhan = the number of times the user wishes to smooth with the
;              hanning kernel. Note that in order to use the hanning
;              kernel the /hanning keyword must be set. The default
;              value for nhan is 1. 
;       pdif = this keyword allows the user to specify a threshold at
;              which to halt the fitting routine. As the fitting
;              routine iterates (getting closer and closer to the fit
;              value), the program checks if the changes in the
;              parameters from iteration to iteration are small enough
;              that the fit can be considered to have converged. The
;              pdif value is the percentage of the current (at that
;              time in the program) parameter values that the change
;              in the parameter values must be less than in order to
;              consider the fit to have converged. For example if the
;              parameter is D, then we will consider D to have been
;              fit if (D - D_guess) LT pdif*D_guess. The default value
;              is 1% (pdif = 0.01).
;       integra = if set, the integrated area under the Gaussian
;                 profile will be printed to the display. The
;                 integration is performed via the Simpson & Romberg
;                 methods and using IDL's int_tabulated function. Note
;                 that in using this option, a system variable called
;                 !mcc_gauss is defined as a 5 element vector. 
;
; OUTPUTS
;       See Keywords (a, cov, yfit).
;
; PROCEDURES CALLED
;       mcc_pclip
;       mcc_hanning
;       getcolor
;
; EXAMPLES
;	(1) Fit a Gaussian and flat continuum to this fake data. 
;             IDL> x = findgen(1000)/10 - 100
;             IDL> y = 1.0 - 0.01*x + 4.0 * EXP((-(x + 40)^2)/(2*2.8^2))
;             IDL> y = y + randomn(seed, 1000)
;             IDL> mcc_gauss1dfit, x, y, /slope, guesscont=2.0, $
;                  /movie, /plot, /verbose, pdif=0.0001
;
; HISTORY
;       Created March 19, 2002 by mcc.
;
;-

;------------
FUNCTION mcc_gauss_slope, xdata
   RETURN, !mcc_gauss[0] + !mcc_gauss[1] * xdata + $
           !mcc_gauss[2] * EXP( -(xdata - !mcc_gauss[3])^2 / (2. * !mcc_gauss[4]^2) )
END

FUNCTION mcc_gauss_cont, xdata
   RETURN, !mcc_gauss[0] + !mcc_gauss[2] * $
           EXP( -(xdata - !mcc_gauss[3])^2 / (2. * !mcc_gauss[4]^2) )
END

FUNCTION mcc_gauss, xdata
   RETURN, !mcc_gauss[2] * EXP( -(xdata - !mcc_gauss[3])^2 / (2. * !mcc_gauss[4]^2) )
END
;------------

PRO mcc_gauss1dfit, xdata, ydata, a=a, cov=cov, yfit=yfit, weight=weight, $
                    swidth=swidth, inum=inum, pdif=pdif, slope=slope, $
                    cont=cont, guesscont=guesscont, hanning=hanning, $
                    nhan=nhan, movie=movie, plot=plot, verbose=verbose, $
                    integra=integra

;;; SPECIFY mcc_pclip AND mcc_hanning AS 
;;; FUNCTIONS RATHER THAN VARIABLES
  FORWARD_FUNCTION mcc_pclip
  FORWARD_FUNCTION mcc_hanning

;;; CHECK THAT ENOUGH PARAMETERS WERE SUPPLIED 
  IF N_PARAMS(0) LT 2 THEN BEGIN
    PRINT, 'CALLING SEQUENCE:  mcc_guassfit, xdata, ydata'
    PRINT, 'KEYWORD PARAMETERS: a= , cov=, swidth= , pdif= , yfit= , ' + $
      'guesscont= , nhan= , /cont, /slope, /movie, /plot, ' + $
      '/verbose, /weight, /hanning'
    RETALL
  ENDIF

;;; CHECK TO SEE IF BOTH /slope AND /cont HAVE BEEN SET.
  IF KEYWORD_SET(slope) AND KEYWORD_SET(cont) THEN BEGIN
    PRINT, 'Keywords /slope and /cont canNOT be passed to' + $
      'mcc_guass1dfit in the same call to procedure!'
    RETALL
  ENDIF

;;; CHECK TO SEE THAT IF /cont OR /slope KEYWORDS ARE GIVEN THAT 
;;; guesscont KEYWORD IS ALSO GIVEN.
  IF (KEYWORD_SET(slope) OR KEYWORD_SET(cont)) AND $ 
    (N_ELEMENTS(guesscont) LT 1) THEN BEGIN
    PRINT, 'If keyword /slope or /cont is passed to procedure, then ' + $
      'guesscont keyword must also be passed!'
    RETALL
  ENDIF

;;; CHECK TO SEE THAT LENGTHS OF ARGUMENTS ARE OKAY
  IF n_elements(xdata) NE n_elements(ydata) THEN BEGIN
    PRINT, 'The arguments xdata and ydata must have the same length, stupid!'
    RETALL
  ENDIF

  IF KEYWORD_SET(weight) AND $ 
    ( n_elements(weight) NE n_elements(ydata) ) THEN BEGIN 
    PRINT, 'The argument weight must have the same length as xdata and ydata!' 
    RETALL
  ENDIF

  IF KEYWORD_SET(hanning) AND NOT(KEYWORD_SET(nhan)) THEN nhan = 1
  IF n_elements(swidth) LT 1 THEN swidth = 5
  IF n_elements(pdif) LT 1 THEN pdif = 0.01
  IF n_elements(inum) LT 1 THEN inum = 100
  IF NOT(KEYWORD_SET(weight)) THEN weight = 1.0

;;; DETERMINE THE NUMBER OF MEASUREMENTS (M)
  M = N_ELEMENTS(xdata)

;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; IF KEYWORD slope IS SET 
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; CHECK IF slope KEYWORD IS SET. AND IF slope KEYWORD 
;;; IS SET, PERFORM FIT USING FITTING FUNCTION
;;; ydata = C + D*xdata + E*e^{-(xdata-u)^2/2*s^2}.
  IF KEYWORD_SET(slope) THEN BEGIN
    N=5.0

;;; GUESS VALUES FOR PARAMETERS
;;; SEE WHAT SMOOTHING KERNEL TO USE (hanning VS boxcar)
    IF KEYWORD_SET(hanning) THEN BEGIN
      ysmooth = mcc_hanning(ydata, nhan)
    ENDIF ELSE BEGIN
;;; ELSE CLIP PEAKS THEN SMOOTH ydata WITH BOXCAR
      ysmooth = SMOOTH(mcc_pclip(ydata), swidth)
    ENDELSE
;;; THEN SUBTRACT guesscont FROM ysmooth
    ycont = ysmooth - guesscont
;;; FIND EXTREMA 
    peak = max( ABS(ycont), sub )
;;; DETERMINE GUESS AMPLITUDE OF GAUSSIAN
    E_guess = ycont[sub]
    u_guess = xdata[sub]
    C_guess = guesscont
    D_guess = (ysmooth[M-1] - ysmooth[0]) / (xdata[M-1] - xdata[0])
    index = where(ABS(ycont) GT ABS(E_guess)/2.0)
    s_guess = ABS( u_guess - MIN(xdata[index]) )

;;; CHECK IF /movie KEYWORD IS SET. IF SET, THEN MAKE MOVIE.
    IF KEYWORD_SET(movie) THEN BEGIN
      window, 4, /pixmap 
      window, 5
      wset, 4
      y_guess = C_guess + D_guess*xdata + E_guess * $
        EXP( -(xdata - u_guess)^2 / (2 * s_guess^2) )
      PLOT, xdata, ydata
      OPLOT, xdata, y_guess, linestyle=2, thick=2, color=getcolor('red')
      wset, 5
      device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
      wait, 0.25
    ENDIF

;;; KEEP TRACK OF NUMBER OF TIMES THROUGH LOOP
    cnt = 0
    JUMP1: cnt = cnt + 1

;;; CALCULATE y VALUES BASED ON GAUSSIAN PROFILE AND 
;;; GUESSED PARAMETERS
    y_guess = C_guess + D_guess*xdata + E_guess * $
      EXP( -(xdata - u_guess)^2 / (2 * s_guess^2) )

;;; DETERMINE DIFFERENCE BETWEEN GUESS y VALUES AND
;;; ACTUAL y VALUES
    del_y = ydata - y_guess

;;; EXPAND del_y ABOUT THE GUESS VALUES AND THEN SOLVE 
;;; THE LINEAR LEAST-SQUARE PROBLEM FOR THE PARAMETERS 
;;; a0 = (C - C_guess), a1 = (D - D_guess), 
;;; a2 = (E - E_guess), a3 = (u - u_guess), AND 
;;; a4 = (s - s_guess).
;;; SO NOW OUR DEPENDENT MEASUREMENTS ARE THE del_y 
;;; VALUES AND THE INDEPENDENT MEASUREMENTS (x's) ARE 
;;; THE PARTIAL DERIVATIVES IN THE TAYLOR EXPANSION AND 
;;; AS STATED BEFORE THE PARAMETERS TO FIT FOR ARE 
;;; (C - C_guess), (D - D_guess), (E - E_guess), 
;;; (u - u_guess), AND (s - s_guess).

    Xtx = fltarr(N, M)
    Xtx[0,*] = 1.0 / weight
    Xtx[1,*] = xdata / weight
    Xtx[2,*] = EXP( -(xdata - u_guess)^2 / (2 * s_guess^2) ) / weight
    Xtx[3,*] = (xdata - u_guess) / s_guess^2 * $
      E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight
    Xtx[4,*] = (xdata - u_guess)^2 / s_guess^3 * $
      E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight
    Ytx = fltarr(1, M)
    Ytx[0, *] = del_y / weight

;;; COMPUTE alpha AND beta MATRICES
    alpha = TRANSPOSE(Xtx) ## Xtx
    beta = TRANSPOSE(Xtx) ## Ytx

;;; NOW COMPUTE THE COEFFICIENTS a
    a = INVERT(alpha) ## beta

;;; CHECK LOOP NUMBER - IF OVER inum, STOP 
;;; AND RETURN ERROR MESSAGE
    IF cnt GT inum-1 THEN BEGIN
      PRINT, 'ERROR: FIT DID NOT CONVERGE!'
      STOP
    ENDIF

    IF KEYWORD_SET(movie) THEN BEGIN
      wset, 4
      PLOT, xdata, ydata
      OPLOT, xdata, y_guess, linestyle=2, thick=2, color=getcolor('red')
      wset, 5
      device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
      wait, 0.25
    ENDIF


;;; CHECK IF FIT WAS CLOSE OR NOT
    a_guess = [C_guess, D_guess, E_guess, u_guess, s_guess]
    junk = where(ABS(a) GT ABS(pdif*a_guess), num)
    IF num GT 0 THEN BEGIN

;;; COMPUTE NEW/CORRECTED GUESS VALUES 
      C_guess = C_guess + a[0]
      D_guess = D_guess + a[1]
      E_guess = E_guess + a[2]
      u_guess = u_guess + a[3]
      s_guess = s_guess + a[4]

;;; LOOP THROUGH AGAIN
      GOTO, JUMP1
    ENDIF ELSE BEGIN

;;; COMPUTE FINAL VALUES FOR PARAMETERS AND FILL a VECTOR
      C = C_guess + a[0]
      D = D_guess + a[1]
      E = E_guess + a[2]
      u = u_guess + a[3]
      s = s_guess + a[4]

      a = [C, D, E, u, s]

;;; NOW COMPUTE THE UNCERTAINTIES IN THE COEFFICIENTS
      y_bar = C + D*xdata + E * EXP(-(xdata - u)^2/(2*s^2))
      resids = ydata - y_bar
      s_avg = total(resids^2) / (M - N)
      cov = s_avg * INVERT(alpha)
      ;diag = cov[(num_a+1) * indgen(num_a)]
      ;ncov = cov / sqrt(diag##diag)

;;; DEFINE yfit TO BE THE y VALUES OF THE GAUSSIAN FIT.
;;; RECALL THAT IF yfit IS DEFINED, THEN FIT IS RETURNED.
      yfit = y_bar

;;; CALCULATE THE CHI-SQUARED VALUE FOR THE FIT.
      chisq = total( (ydata - y_bar)^2 / weight )

;;; CHECK IF KEYWORD /plot IS SET. IF SET, PLOT RESULTS.
      IF KEYWORD_SET(plot) THEN BEGIN
        !P.MULTI = [0,1,2]
        WINDOW, 6
        PLOT, xdata, ydata
        OPLOT, xdata, y_bar, linestyle=2, thick=2, color=getcolor('red')
        PLOT, xdata, resids, psym=4
        OPLOT, [xdata[0], xdata[M-1]], [0, 0], linestyle=2
        !P.MULTI = [0,1,1]
      ENDIF

      IF KEYWORD_SET(verbose) THEN BEGIN
        print
        print, 'Number of iterations: ' + STRING(cnt, FORMAT='(f3.0)')
        print
        print, '    PARAMETER      ' + '   SIGMA   '
        print, '-------------------------------'
        print, 'C = ', a[0], SQRT(cov[0,0])
        print, 'D = ', a[1], SQRT(cov[1,1])
        print, 'E = ', a[2], SQRT(cov[2,2])
        print, 'u = ', a[3], SQRT(cov[3,3])
        print, 's = ', a[4], SQRT(cov[4,4])
        print
      ENDIF

;;; Integrate to find the area under the Gaussian profile. 
;;; This is intended for use as a photometry tool.
      IF KEYWORD_SET(integra) THEN BEGIN
        xmin = min(xdata)
        xmax = max(xdata)
        DEFSYSV, '!mcc_gauss', [C, D, E, u, s]
        int_simp = QSIMP('mcc_gauss_slope', xmin, xmax)
        int_romb = QROMB('mcc_gauss_slope', xmin, xmax)
        int_tab = INT_TABULATED(xdata, y_bar - C - D*xdata, /SORT)
        print
        print, 'Integration Results: '
        print, 'Simpson Method = ', int_simp
        print, 'Romberg Method = ', int_romb
        print, 'Int_tabulate Result = ', int_tab
        print
      ENDIF

    ENDELSE
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; IF KEYWORD cont IS SET 
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; NOW IF THE slope KEYWORD IS NOT SET THEN CHECK IF 
;;; THE KEYWORD cont IS SET. IF cont IS SET, THEN PERFORM 
;;; FIT USING FITTING FUNCTION ydata = C + E*e^{-(xdata-u)^2/2*s^2}.
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(cont) THEN BEGIN 
      N = 4.0
      
;;; GUESS VALUES FOR PARAMETERS
;;; SEE WHAT SMOOTHING KERNEL TO USE (hanning VS boxcar)
      IF KEYWORD_SET(hanning) THEN BEGIN
        ysmooth = mcc_hanning(ydata, nhan)
      ENDIF ELSE BEGIN
;;; ELSE CLIP PEAKS THEN SMOOTH ydata WITH BOXCAR
        ysmooth = SMOOTH(mcc_pclip(ydata), swidth)
      ENDELSE
  
;;; THEN SUBTRACT guesscont FROM ysmooth
      ycont = ysmooth - guesscont
;;; FIND EXTREMA 
      peak = max( ABS(ycont), sub )
;;; DETERMINE GUESS AMPLITUDE OF GAUSSIAN
      E_guess = ycont[sub]
      u_guess = xdata[sub]
      C_guess = guesscont
      index = where(ABS(ycont) GT ABS(E_guess)/2.71828)
      s_guess = ABS( u_guess - MIN(xdata[index]) )

;;; CHECK IF movie KEYWORD IS SET. IF SET MAKE MOVIE.
      IF KEYWORD_SET(movie) THEN BEGIN
        window, 4, /pixmap
        window, 5
        wset, 4
        y_guess = C_guess + E_guess * EXP( -(xdata - u_guess)^2 / $
                                           (2 * s_guess^2) )
        PLOT, xdata, ydata
        OPLOT, xdata, y_guess, linestyle=2, thick=2, color=getcolor('red')
        wset, 5
        device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
        wait, 0.25
      ENDIF

;;; KEEP TRACK OF NUMBER OF TIMES THROUGH LOOP
      cnt = 0
      JUMP2: cnt = cnt + 1

;;; CALCULATE y VALUES BASED ON GAUSSIAN PROFILE AND 
;;; GUESSED PARAMETERS
      y_guess = C_guess + E_guess * EXP( -(xdata - u_guess)^2 / $
                                         (2 * s_guess^2) )

;;; DETERMINE DIFFERENCE BETWEEN GUESS y VALUES AND
;;; ACTUAL y VALUES
      del_y = ydata - y_guess

;;; EXPAND del_y ABOUT THE GUESS VALUES AND THEN SOLVE 
;;; THE LINEAR LEAST-SQUARE PROBLEM FOR THE PARAMETERS 
;;; a0 = (C - C_guess), a1 = (E - E_guess), 
;;; a2 = (u - u_guess), AND a3 = (s - s_guess).
;;; SO NOW OUR DEPENDENT MEASUREMENTS ARE THE del_y 
;;; VALUES AND THE INDEPENDENT MEASUREMENTS (x's) ARE 
;;; THE PARTIAL DERIVATIVES IN THE TAYLOR EXPANSION AND 
;;; AS STATED BEFORE THE PARAMETERS TO FIT FOR ARE 
;;; (C - C_guess), (E - E_guess), (u - u_guess), AND (s - s_guess).

      Xtx = fltarr(N, M)
      Xtx[0,*] = 1.0 / weight
      Xtx[1,*] = EXP( -(xdata -u_guess)^2 / (2 * s_guess^2) ) / weight
      Xtx[2,*] = (xdata - u_guess) / s_guess^2 * $
        E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight
      Xtx[3,*] = (xdata - u_guess)^2 / s_guess^3 * $
        E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight

      Ytx = fltarr(1, M)
      Ytx[0, *] = del_y / weight
      
;;; COMPUTE alpha AND beta MATRICES
      alpha = TRANSPOSE(Xtx) ## Xtx
      beta = TRANSPOSE(Xtx) ## Ytx

;;; NOW COMPUTE THE COEFFICIENTS a
      a = INVERT(alpha) ## beta

;;; CHECK LOOP NUMBER - IF OVER inum, STOP 
;;; AND RETURN ERROR MESSAGE
      IF cnt GT inum-1 THEN BEGIN
        PRINT, 'ERROR: FIT DID NOT CONVERGE!'
        STOP
      ENDIF

;;; CHECK IF movie KEYWORD IS SET. IF SET MAKE MOVIE.
      IF KEYWORD_SET(movie) THEN BEGIN
        wset, 4
        PLOT, xdata, ydata
        OPLOT, xdata, y_guess, linestyle=2, thick=2, color=getcolor('red')
        wset, 5
        device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
        wait, 0.25
      ENDIF

;;; CHECK IF FIT WAS CLOSE OR NOT
      a_guess = [C_guess, E_guess, u_guess, s_guess]
      junk = where(ABS(a) GT ABS(pdif*a_guess), num)
      IF num GT 0 THEN BEGIN

;;; COMPUTE NEW/CORRECTED GUESS VALUES 
        C_guess = C_guess + a[0]
        E_guess = E_guess + a[1]
        u_guess = u_guess + a[2]
        s_guess = s_guess + a[3]

;;; LOOP THROUGH AGAIN
        GOTO, JUMP2
      ENDIF ELSE BEGIN

;;; COMPUTE FINAL VALUES FOR PARAMETERS
        C = C_guess + a[0]
        E = E_guess + a[1]
        u = u_guess + a[2]
        s = s_guess + a[3]

        a = [C, E, u, s]

;;; NOW COMPUTE THE UNCERTAINTIES IN THE COEFFICIENTS
        y_bar = C + E * EXP(-(xdata - u)^2/(2*s^2))
        resids = ydata - y_bar
        s_avg = total(resids^2) / (M - N)
        cov = s_avg * INVERT(alpha)
         ;diag = cov[(num_a+1) * indgen(num_a)]
         ;ncov = cov / sqrt(diag##diag)
         
;;; DEFINE yfit TO BE THE y VALUES OF THE GAUSSIAN FIT.
;;; RECALL THAT IF yfit IS DEFINED, THEN FIT IS RETURNED.
         yfit = y_bar

;;; CHECK IF KEYWORD /plot IS SET. IF SET, PLOT RESULTS.
         IF KEYWORD_SET(plot) THEN BEGIN
           !P.MULTI = [0,1,2]
           WINDOW, 6
           PLOT, xdata, ydata
           OPLOT, xdata, y_bar, linestyle=2
           PLOT, xdata, resids, psym=4
           OPLOT, [xdata[0], xdata[M-1]], [0, 0], linestyle=2
           !P.MULTI=[0,1,1]
         ENDIF

         IF KEYWORD_SET(verbose) THEN BEGIN
           print 
           print, 'Number of iterations: ' + STRING(cnt, FORMAT='(f3.0)')
           print
           print, '     PARAMETER      ' + '   SIGMA   '
           print, '-------------------------------'
           print, 'C = ', a[0], SQRT(cov[0,0])
           print, 'E = ', a[1], SQRT(cov[1,1])
           print, 'u = ', a[2], SQRT(cov[2,2])
           print, 's = ', a[3], SQRT(cov[3,3])
           print
         ENDIF

;;; Integrate to find the area under the Gaussian profile. 
;;; This is intended for use as a photometry tool.
         IF KEYWORD_SET(integra) THEN BEGIN
           xmin = min(xdata)
           xmax = max(xdata)
           DEFSYSV, '!mcc_gauss', [C, 0., E, u, s]
           int_simp = QSIMP('mcc_gauss_cont', xmin, xmax)
           int_romb = QROMB('mcc_gauss_cont', xmin, xmax)
           int_tab = INT_TABULATED(xdata, y_bar - C, /SORT)
           print
           print, 'Integration Results: '
           print, 'Simpson Method = ', int_simp
           print, 'Romberg Method = ', int_romb
           print, 'Int_tabulate Result = ', int_tab
           print
         ENDIF

      ENDELSE
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; IF KEYWORD cont AND slope ARE NOT SET 
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; NOW IF THE slope AND cont KEYWORD ARE NOT SET, 
;;; THEN PERFORM FIT USING FITTING FUNCTION 
;;; ydata = E*e^{-(xdata-u)^2/2*s^2}.
    ENDIF ELSE BEGIN
      N = 3.0

;;; GUESS VALUES FOR PARAMETERS
;;; SEE WHAT SMOOTHING KERNEL TO USE (hanning VS boxcar)
      IF KEYWORD_SET(hanning) THEN BEGIN
        ysmooth = mcc_hanning(ydata, nhan)
      ENDIF ELSE BEGIN
;;; ELSE CLIP PEAKS THEN SMOOTH ydata WITH BOXCAR
        ysmooth = SMOOTH(mcc_pclip(ydata), swidth)
      ENDELSE

;;; FIND EXTREMA 
      peak = max( ysmooth, sub )
;;; DETERMINE GUESS AMPLITUDE OF GAUSSIAN
      E_guess = ysmooth[sub]
      u_guess = xdata[sub]
      index = where(ABS(ysmooth) GT ABS(E_guess)/2.71828)
      s_guess = ABS( u_guess - MIN(xdata[index]) )

;;; CHECK IF movie KEYWORD IS SET. IF SET MAKE MOVIE.
      IF KEYWORD_SET(movie) THEN BEGIN
        window, 4, /pixmap
        window, 5
        wset, 4
        y_guess = E_guess * EXP( -(xdata - u_guess)^2 / $
                                 (2 * s_guess^2) )
        PLOT, xdata, ydata
        OPLOT, xdata, y_guess, linestyle=2, thick=2
        wset, 5
        device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
      ENDIF

;;; KEEP TRACK OF NUMBER OF TIMES THROUGH LOOP
      cnt = 0
      JUMP3: cnt = cnt + 1

;;; CALCULATE y VALUES BASED ON GAUSSIAN PROFILE AND 
;;; GUESSED PARAMETERS
      y_guess = E_guess * EXP( -(xdata - u_guess)^2 / (2 * s_guess^2) )

;;; DETERMINE DIFFERENCE BETWEEN GUESS y VALUES AND
;;; ACTUAL y VALUES
      del_y = ydata - y_guess

;;; EXPAND del_y ABOUT THE GUESS VALUES AND THEN SOLVE 
;;; THE LINEAR LEAST-SQUARE PROBLEM FOR THE PARAMETERS 
;;; a0 = (E - E_guess), a1 = (u - u_guess), AND a2 = (s - s_guess).
;;; SO NOW OUR DEPENDENT MEASUREMENTS ARE THE del_y 
;;; VALUES AND THE INDEPENDENT MEASUREMENTS (x's) ARE 
;;; THE PARTIAL DERIVATIVES IN THE TAYLOR EXPANSION AND 
;;; AS STATED BEFORE THE PARAMETERS TO FIT FOR ARE 
;;; (E - E_guess), (u - u_guess), AND (s - s_guess).

      Xtx = fltarr(N, M)
      Xtx[0,*] = EXP( -(xdata -u_guess)^2 / (2 * s_guess^2) ) / weight
      Xtx[1,*] = (xdata - u_guess) / s_guess^2 * $
        E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight
      Xtx[2,*] = (xdata - u_guess)^2 / s_guess^3 * $
        E_guess * EXP(-(xdata - u_guess)^2/(2*s_guess^2)) / weight
      
      Ytx = fltarr(1, M)
      Ytx[0, *] = del_y / weight
      
;;; COMPUTE alpha AND beta MATRICES
      alpha = TRANSPOSE(Xtx) ## Xtx
      beta = TRANSPOSE(Xtx) ## Ytx

;;; NOW COMPUTE THE COEFFICIENTS a
      a = INVERT(alpha) ## beta

;;; CHECK LOOP NUMBER - IF OVER inum, STOP 
;;; AND RETURN ERROR MESSAGE
      IF cnt GT inum-1 THEN BEGIN
        PRINT, 'ERROR: FIT DID NOT CONVERGE!'
        STOP
      ENDIF

      IF KEYWORD_SET(movie) THEN BEGIN
        wset, 4
        PLOT, xdata, ydata
        OPLOT, xdata, y_guess, linestyle=2, thick=2
        wset, 5
        device, copy = [0,0,!d.x_size,!d.y_size,0,0,4]
        wait, 0.25
      ENDIF

;;; CHECK IF FIT WAS CLOSE OR NOT
      a_guess = [E_guess, u_guess, s_guess]
      junk = where(ABS(a) GT ABS(pdif*a_guess), num)
      IF num GT 0 THEN BEGIN
;;; COMPUTE NEW/CORRECTED GUESS VALUES 
        E_guess = E_guess + a[0]
        u_guess = u_guess + a[1]
        s_guess = s_guess + a[2]

;;; LOOP THROUGH AGAIN
        GOTO, JUMP3
      ENDIF ELSE BEGIN

;;; COMPUTE FINAL VALUES FOR PARAMETERS
        E = E_guess + a[0]
        u = u_guess + a[1]
        s = s_guess + a[2]
        a = [E, u, s]

;;; NOW COMPUTE THE UNCERTAINTIES IN THE COEFFICIENTS
        y_bar = E * EXP(-(xdata - u)^2/(2*s^2))
        resids = ydata - y_bar
        s_avg = total(resids^2) / (M - N)
        cov = s_avg * INVERT(alpha)
         ;diag = cov[(num_a+1) * indgen(num_a)]
         ;ncov = cov / sqrt(diag##diag)
         
;;; DEFINE yfit TO BE THE y VALUES OF THE GAUSSIAN FIT.
;;; RECALL THAT IF yfit IS DEFINED, THEN FIT IS RETURNED.
        yfit = y_bar

;;; CHECK IF KEYWORD /plot IS SET. IF SET, PLOT RESULTS.
        IF KEYWORD_SET(plot) THEN BEGIN
          !P.MULTI = [0,1,2]
          WINDOW, 6
          PLOT, xdata, ydata
          OPLOT, xdata, y_bar, linestyle=2
          PLOT, xdata, resids, psym=4
          OPLOT, [xdata[0], xdata[M-1]], [0, 0], linestyle=2
          !P.MULTI = [0,1,1]
        ENDIF

        IF KEYWORD_SET(verbose) THEN BEGIN
          print 
          print, 'Number of iterations: ' + STRING(cnt, FORMAT='(f3.0)')
          print
          print, '    PARAMETER      ' + '   SIGMA   '
          print, '-------------------------------'
          print, 'E = ', a[0], SQRT(cov[0,0])
          print, 'u = ', a[1], SQRT(cov[1,1])
          print, 's = ', a[2], SQRT(cov[2,2])
          print
        ENDIF

;;; Integrate to find the area under the Gaussian profile. 
;;; This is intended for use as a photometry tool.
        IF KEYWORD_SET(integra) THEN BEGIN
          xmin = min(xdata)
          xmax = max(xdata)
          DEFSYSV, '!mcc_gauss', [0., 0., E, u, s]
          int_simp = QSIMP('mcc_gauss', xmin, xmax)
          int_romb = QROMB('mcc_gauss', xmin, xmax)
          int_tab = INT_TABULATED(xdata, y_bar, /SORT)
          print
          print, 'Integration Results: '
          print, 'Simpson Method = ', int_simp
          print, 'Romberg Method = ', int_romb
          print, 'Int_tabulate Result = ', int_tab
          print
        ENDIF

      ENDELSE

    ENDELSE

  ENDELSE


END






