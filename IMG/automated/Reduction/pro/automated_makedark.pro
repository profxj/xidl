;+
;
; NAME:
;   automated_makedark
;
; PURPOSE:
;   Estimates the dark current across the image
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
;
;
; OUTPUTS: 
;
; Find the dark current value as function of time
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

PRO automated_makedark, str, index,  BIAS=bias, STATUS=status

  
  ;;start here
  ndark=n_elements(index)
  ;;make storage
  dark=make_array(ndark,2048,2048,/float)
    
  valuedark=fltarr(ndark)
  valuetime=fltarr(ndark)
     
  
  ;;now loop over all the images
  for img=0, ndark-1 do begin
      
      ;;open
      splog, "Working on  ", str.name[index[img]]  
      chips=mrdfits(strtrim(str.name[index[img]],2),0,hea,/silent,/fscale)
      valuetime[img]=str.extime[index[img]]
      
      ;;subtract the bias out of the darks
      ;;go for oscan subtraction
      supl_oscan, chips, imageout
  
      ;;remove zero
      if keyword_set(bias) then imageout=imageout-bias
      
      ;;get median dark
      valuedark[img]=djs_median(imageout)
      ;;splog, "Dark current for time ",valuetime[img], " is ", valuedark[img]
     

  endfor
  
  
  ;;plot and fit
  ;;Comment out so extra windows are not interfering with the cron job
  ;;window, 1, xsize=400, ysize=400  
  ;;plot, valuetime, valuedark, psym=1, yrange=[-1,10], xtitle="time
  ;;(s)", ytitle="cnt", title='Dark current'
  res = linfit(valuetime,valuedark,yfit=fit) 
  ;;oplot, valuetime, fit
  
  
  ;;store fit
  ;;m
  status.dark[0]=res[1]
  ;;q
  status.dark[1]=res[0]

  
END
