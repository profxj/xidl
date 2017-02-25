 pro ksone_wght, data, func_name, d, prob, WGHT = WGHT, PLOT = plot, TITLE = title, _EXTRA = extra
;+
; NAME:
;       KSONE_WGHT
; PURPOSE:
;       Compute the one-sided Kolmogorov-Smirnov statistic
; EXPLANATION:
;       Returns the Kolmogorov-Smirnov statistic and associated probability for 
;       for an array of data values and a user-supplied cumulative distribution
;       function (CDF) of a single variable.   Algorithm from the procedure of
;       the same name in "Numerical Recipes" by Press et al. 2nd edition (1992)
;
; CALLING SEQUENCE:
;       ksone_wght, data, func_name, D, prob, [ /PLOT ]
;
; INPUT PARAMETERS:
;       data -  vector of data values, must contain at least 4 elements for the
;               K-S statistic to be meaningful 
;       func_name - scalar string giving the name of the cumulative distribution
;               function.    The function must be defined to accept the data
;               vector as its only input (see example), though keywords may be
;               passed via the _EXTRA facility.
;
; OUTPUT PARAMETERS:
;       D - floating scalar giving the Kolmogorov-Smirnov statistic.   It 
;               specified the maximum deviation between the cumulative 
;               distribution of the data and the supplied function 
;       prob - floating scalar between 0 and 1 giving the significance level of
;               the K-S statistic.   Small values of PROB show that the 
;               cumulative distribution function of DATA is significantly 
;               different from FUNC_NAME.
;
; OPTIONAL INPUT KEYWORD:
;       PLOT - If this keyword is set and non-zero, then KSONE_WGHT will display a
;               plot of the CDF of the data with the supplied function 
;               superposed.   The data value where the K-S statistic is 
;               computed (i.e. at the maximum difference between the data CDF 
;               and the function) is indicated by a vertical line.
;               KSONE_WGHT accepts the _EXTRA keyword, so that most plot keywords
;               (e.g. TITLE, XTITLE, XSTYLE) can also be passed to KSONE_WGHT.
;
; EXAMPLE:
;       Determine if a vector created by the RANDOMN function is really 
;       consistent with a Gaussian distribution with unit variance.
;       The CDF of a Gaussian is the error function except that a factor
;       of 2 is included in the error function.   So we must create a special
;       function:
;
;       function gauss_cdf, x
;       return, errorf( x/sqrt(2) )
;       end
;
;       IDL> data = randomn(seed, 50)          ;create data array to be tested
;       IDL> ksone_wght, abs(data), 'gauss_cdf', D, prob, /PLOT     ;Use K-S test
;      
;       PROB gives the probability that the null hypothesis (DATA came from a 
;       Gaussian distribution with unit variance) is correct.
;
; NOTES:
;       The code for PROB_KS is from the 2nd (1992) edition of Numerical 
;       Recipes which includes a more accurate computation of the K-S 
;       significance for small values of N than the first edition.
;
; PROCEDURE CALLS
;       procedure PROB_KS - computes significance of K-S distribution
;
; REVISION HISTORY:
;       Written     W. Landsman                   August, 1992
;       Accept _EXTRA keywords   W. Landsman      September, 1995          
;       Fixed possible bug in plot display showing position maximum difference
;       in histogram   M. Fardal/ W. Landsman      March, 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Documentation updates   W. Landsman   June 2003
;       Pass _EXTRA to func_name  M. Fitzgerald    April, 2005
;       Adopted from IDLUTILS/goddard/pro/math/ksone   K. Cooksey   Dec, 2008
;-
 On_error, 2

 if ( N_params() LT 3 ) then begin
    print,'Syntax - ksone_wght, data, func_name, D, [prob ,/PLOT]'
    return
 endif

 N = N_elements( data )
 if N LT 3 then message, $
   'ERROR - Input data values (first param) must contain at least 3 values'
 if not keyword_set(wght) then wght = replicate(1.,N)

 sortwght = wght[ sort(data) ]
 sortdata = data[ sort( data ) ]                                   

 fn = total(SORTWGHT,/cumulative)/total(SORTWGHT)
 f0 = (total(SORTWGHT,/cumulative) - SORTWGHT[0])/total(SORTWGHT)
 ff = call_function( func_name, sortdata) ; , _EXTRA = extra)

 D = max( [ max( abs(f0-ff), sub0 ), max( abs(fn-ff), subn ) ], msub )

 PROB_KS, D, N, prob           ;Compute significance of K-S statistic

 if keyword_set(plot) then begin
     clr = getcolor(/load)
     if msub EQ 0 then begin 
        yone = f0
        xthr = [sortdata[sub0], sortdata[sub0]]
     endif else begin
        yone = fn
        xthr = [sortdata[subn], sortdata[subn]]
     endelse 

     if not keyword_set(title) then title = ''
     if prob lt 1e-5 then $
        fmt = '("KS1: D = ",f7.5,"; N = ",i5,"; P = ",e9.2)' $
     else fmt = '("KS1: D = ",f7.5,"; N = ",i5,"; P = ",f7.5)' 
     plot, sortdata, yone, psym=10, color=clr.black, background=clr.white,$
           title=title,xtitle='Data',ytitle='Cumulative Distribution',$
           thick=2,charsize=2,yrange=[0.,1.],xrange=[sortdata[0],sortdata[n-1]],$
           /xstyle,/ystyle, _extra=extra
     oplot,sortdata,ff,psym=10, color=clr.red, thick=2
     oplot,xthr,[0,1],psym=10, color=clr.blue, thick=2
     xyouts,0.5,0.15,string(D,N,prob, format=fmt),/normal,$
            alignment=0.5,charsize=2

  endif

 return
 end
