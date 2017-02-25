;+
; NAME:
;   automated_makebias
;
; PURPOSE:
;   Make a bias (oscan and image)
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; str     - the structure with relevant information
; index   - the index where the bias is
; 
; OPTIONAL INPUTS:
;  
; plan         - string to save the bias
;
; OUTPUTS: 
;
; Creates a fits of the median bias and passes back in biasmedian
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
  


pro automated_makebias, str, biasmedian, index, plan=plan
 
  if(n_elements(index) gt 10) then begin
      index=index[0:10]
      splog, "Considering only 10 bias frames!"
  endif
  
  ;;make storage for the data chip
  bias1=make_array(N_ELEMENTS(index),2048,2048,/float)
  
  for pos=0, n_elements(index)-1 do begin
      
      ;;open files, tells what object looking at
      splog, "Working on  bias ", str.name[index[pos]]
      ;;chip is 2200x2200, hea is where the header lives
      chip=mrdfits(strtrim(str.name[index[pos]],2),0,hea,/silent,/fscale)
      
      ;;go for oscan subtraction
      automated_oscan, chip, biasimage
      
      ;;store each chip
      bias1[pos,*,*]=biasimage[*,*]
      
  endfor

  ;;combine the final bias
  splog, "Creating median bias... "  
  
  biasmedian=djs_median(bias1,1) 
  
  ;;save median bias
  mwrfits, biasmedian, plan+"medbias.fits", /create

  splog, "All done with median bias!"
 
end
