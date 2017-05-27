;+
;
; NAME
;      mcc_pclip.pro
;
; PURPOSE
;      The mcc_pclip function takes a vector of data values (ydata) as
;      its only required input and removes all outlying data points
;      and replaces those values with the average of the surrounding
;      data points. The peak-removed vector is then returned as the
;      output from the function. The number of surrounding data values
;      included in the average is specified via the swidth
;      keyword. The criteria for whether or not a point is an outlying
;      peak is specified using the /percen and thresh keywords.
;
; SYNTAX
;	yout = mcc_pclip(ydata, [swidth=swdith, thresh=thresh, /percen])
;
; INPUTS
;	ydata = the set of data values to be examined for peaks and have 
;		those peaks removed.	
;
; KEYWORDS
;	swidth = this keyword specifies the number of surrounding data values
;		 that are used in determining the average of the data around 
;		 each data point. The size of the window over which the 
; 		 average is taken is given by 2swidth+1 spanning the indices
;		 in the ydata vector [(i-swidth):(i+swidth)] such that the 
;		 the point in question (the one we are deciding to clip 
;		 or not) is included in the average.
;	thresh = the number of standard deviations or the percentage used 
;		 in determining whether a given data value is an outlier or 
;		 not. If the /percen keyword is not set, then thresh is 
;		 taken to be the number of standard deviations away from the 
; 		 mean such that a given data point will be clipped if it 
;		 lies more than thresh standard deviations away from the 
;		 surrounding average. The default value is 3.0 standard 
;		 deviations. If \percen is set, then thresh must be set too! 
;		 In this case, thresh specifies the percentage of the mean 
;		 that a point can lie away from the mean and not be clipped.
;	\percen = if this keyword is set, then the peak clipping criteria 
;		  will be based on a percentage of the surrounding mean 
;		  value (see thresh). If \percen is set, then thresh must 
;		  be set too! If \percen is not set, then the clipping 
;		  criteria is by default based on standard deviations away 
;		  from the mean.
;
; OUTPUTS
;	yout = a vector of the same original length as ydata containing 
;	       the original ydata values except where peaks were clipped. In
;	       those cases, the peak values have been replaced by the average
;	       of the surrounding data points.
;
; PROCEDURES CALLED 
;	None
;
; EXAMPLES
;	(1) Create the following data set with peaks.
;		IDL> y = randomn(seed,500)
;		IDL> i = findgen(25)*13 + findgen(25)*3
;		IDL> y[i] = y[i]*randomn(seed,25)*50.0
;		IDL> yout = mcc_pclip(y)
;		IDL> PLOT, y
;		IDL> OPLOT, yout, color=getcolor('red')
;
; HISTORY
;	Created March 19, 2002 by mcc.
;
;-

FUNCTION mcc_pclip, ydata, swidth=swidth, thresh=thresh, percen=percen

;;; CHECK THAT ENOUGH PARAMETERS WERE SUPPLIED.
  IF N_PARAMS(0) LT 1 THEN BEGIN
    PRINT, 'CALLING SEQUENCE:  yout = mcc_pclip(ydata)'
    PRINT, 'KEYWORD PARAMETERS: swidth= , thresh= , /percen'
    RETALL
  ENDIF

;;; CHECK IF swidth KEYWORD IS SET.
  IF N_ELEMENTS(swidth) LT 1 THEN swidth = 10.0

;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; CHECK IF /percen KEYWORD HAS BEEN SET. IF SO, THEN PROCEED TO CLIP
;;; PEAKS USING PERCENTAGE GIVEN IN thresh
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IF KEYWORD_SET(percen) THEN BEGIN 
    IF n_elements(thresh) LT 1 THEN BEGIN 
      PRINT, 'IF KEYWORD /percen IS SET, THEN thresh MUST ALSO BE SET!'
    ENDIF 

;;; DETERMINE NUMBER OF ydata VALUES.
    len = n_elements(ydata)
;;; STEP THRU ydata VECTOR CHECKING TO SEE IF EACH POINT IS WITHIN THE
;;; THRESHOLD LEVEL WITH RESPECT TO THE AVERAGE OF THE YDATA VALUES IN
;;; THE WINDOW OF WIDTH swidth.
    FOR i=0,len-1 DO BEGIN
      IF i LT swidth THEN BEGIN 
        avg = MEAN(ydata[0:(i+swidth)])
      ENDIF ELSE BEGIN
        IF len-i LE swidth THEN BEGIN 
          avg = MEAN(ydata[(i-swidth):(len-1)])
        ENDIF ELSE BEGIN
          avg = MEAN(ydata[(i-swidth):(i+swidth)])
        ENDELSE
      ENDELSE
;;; IF ydata[i] VALUE IS TO FAR AWAY FROM MEAN, THEN RESET ITS VALUE
;;; IN ydata VECTOR TO THE AVERAGE VALUE.
      IF ABS(ydata[i] - avg) GT (thresh*avg) THEN ydata[i] = avg
;;; CLOSE FOR-DO LOOP.
    ENDFOR

;;; RETURN THE VECTOR OF ydata VALUES WITH PEAKS REMOVED.
    RETURN, ydata

;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;; NOW IF \percen KEYWORD IS NOT SET, THEN PROCEED TO CLIP PEAKS
;;; USING NUMBER OF STANDARD DEVIATIONS GIVEN BY thresh.
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF ELSE BEGIN
    IF n_elements(thresh) LT 1 THEN thresh = 3.0

;;; DETERMINE NUMBER OF ydata VALUES. 
    len = n_elements(ydata)

;;; STEP THRU ydata VECTOR CHECKING TO SEE IF EACH POINT IS WITHIN THE
;;; THRESHOLD LEVEL WITH RESPECT TO THE AVERAGE OF THE YDATA VALUES IN
;;; THE WINDOW OF WIDTH swidth.
    FOR i=0,len-1 DO BEGIN
      IF i LT swidth THEN BEGIN 
        dset = ydata[0:(i+swidth)]
        avg = mean(dset)
      ENDIF ELSE BEGIN 
        IF len-i LE swidth THEN BEGIN 
          dset = ydata[(i-swidth):(len-1)]
          avg = mean(dset)
        ENDIF ELSE BEGIN
          dset = ydata[(i-swidth):(i+swidth)]
          avg = mean(dset)
        ENDELSE
      ENDELSE

;;; DETERMINE THE DISPERSION OF THE POINTS INSIDE THE WINDOW OF WIDTH
;;; swidth.
      sigma = SQRT(total(dset^2) / (2.0 * swidth + 1.0) - avg^2)

;;; IF ydata[i] VALUE IS TO FAR AWAY FROM MEAN, THEN RESET ITS VALUE
;;; IN ydata VECTOR TO THE AVERAGE VALUE.
      IF ABS(ydata[i] - avg) GT thresh*sigma THEN ydata[i] = avg

;;; CLOSE FOR-DO LOOP.
   ENDFOR

;;; RETURN THE VECTOR OF ydata VALUES WITH PEAKS REMOVED.
   RETURN, ydata

ENDELSE

;;; END OF FUNCTION.
END


