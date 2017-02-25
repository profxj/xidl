;+ 
; NAME:
; x_cumulative
;    Version 1.1
;
; PURPOSE:
;;    Inputs an array and genreates the cumulative arrays for it.
;  
;
; CALLING SEQUENCE:
;  x_cumulative, array, xsrt, ycum
;
; INPUTS:
;   Array
;
; RETURNS:
;
; OUTPUTS:
;  xsrt == Sorted values of the array
;  ycum == Cumulative array (0,1)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Oct-2012 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_cumulative, array, xsrt, ycum

  ; 
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_cumulative, array, xsrt, ycum [v1.0]'
      return
  endif 

  ;; Sort
  srt = sort(array)
  xsrt = array[srt]
  
  ;; Cumulative
  narray = n_elements(xsrt)
  ycum = findgen(narray) / (narray-1)
  
  return

end
