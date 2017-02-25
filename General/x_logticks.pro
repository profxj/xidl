;+
;
; (C) Cooperative Institute for Meteorological Satellite Studies, 2000
;
; NAME:
;       x_logticks
;
; PURPOSE:
;       Function to print logarithmic axis tickmarks with exponential
;       output.
;
; CATEGORY:
;       Graphics
;
; LANGUAGE:
;       IDL v5
;
; CALLING SEQUENCE:
;       PLOT, x, y, /YLOG, YICKFORMAT = 'logticks_exp'
;
; INPUTS:
;       axis:  the axis number. 0 for X axis, 1 for Y axis, 2 for Z axis.
;       index: the tick mark index which starts at 0.
;       value: the default tick mark value (a floating-point number);
;
; OUTPUTS:
;       Function returns a string containing the tick mark labels.
;
; CREATION HISTORY:
;       Written by:     Paul van Delst, CIMSS/SSEC, 08-Nov-2000
;                       paul.vande...@ssec.wisc.edu
;
;-
FUNCTION x_logticks, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
   tickmark = '10!u' + STRTRIM( STRING( exponent ), 2 ) + '!N'
   ; Return the formatted tickmark string
   RETURN, tickmark
END 
