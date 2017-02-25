
;+
; NAME:
;   automated_oscan
;
; PURPOSE:
; Procedure that computes the bias using the overscan
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; chips - the full cdd (image+oscan)
; 
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
; dataimage - return the data only (oscan trimmed)
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Written and revised by MF
;
;-
;------------------------------------------------------------------------------
;; 

pro automated_oscan, chips, dataimage

STOP 

  ;;empty so far
  dataimage=make_array(2048,2048,/float)
  
  ;;run over the chip
  
  ;;this is the image minus the oscan
  data=reform(chips[50:2097,1:2048])
  
  ;;get bias info, this is the overscan array, reform apparently does nothing
  overscan=reform(chips[2105:2195,1:2048])
  
  ;;djs_median collapses an array along the specified axis by median combining.
  ;;so a [91,2048] array collapsed along 1 gives a [2048] array combined across
  ;; the rows
  bia=djs_median(overscan,1)
  ;;undefined until now, transpose makes bia 2-d, [1,2048]. Rebin makes it [2048,2048]
  bias=rebin(transpose(bia),2048,2048)
  
  dataimage=data-bias
  
end
