;+ 
; NAME:
; xdimg_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first file in the directory with 'xdimg*fits'
;
; CALLING SEQUENCE:
;  xdimg = xdimg_ar()
;
; INPUTS:
;
; RETURNS:
;    xdimg -  WFCCD structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg = xdimg_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xdimg_ar

;  Optional Keywords
  
; Find the file
  a = findfile('./xdimg*fits', count=count)
  if count EQ 0 then return, -1 $
  else begin
      print, 'Reading from '+a[0]
      return, mrdfits(a[0],1,STRUCTYP='xdimgstrct')
  endelse

end
      
