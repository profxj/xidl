;+
; NAME:
;   automated_cron_wrapper
;
; PURPOSE:
;
;    Program acts as a wrapper for the pipeline. 
;    Called by cron every morning to automatically download, reduce, and analyze 
;    all blazar images taken with Super-LOTIS so I can be lazy
; 
;    1)logs all events
;    2)gets date
;    3)do science
;    4)update website
;
;
;
; CALLING SEQUENCE:
;   automated_cron_wrapper, date=
;
; INPUTS:
;    
; OPTIONAL INPUTS:
;
;   date -- String DDMMMYYYY
;              
; OUTPUTS: 
;  Returns a structure describing the instrument configuration which
;  is used to guide the reduction steps.
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;   Run the wrapper for 04 July 2012 
;       automated_cron_wrapper, date='120704'
;
;   Run the wrapper for last night's observations
;       automated_cron_wrapper
; BUGS:
;
; PROCEDURES CALLED:
;   automated_work, automated_replot
;
; REVISION HISTORY:
;                Written by KK
;   05-Jun-2012  Revised by MF
;-
;------------------------------------------------------------------------------
;; 

pro automated_cron_wrapper,date=date,FLWO=FLWO


  if (~keyword_set(date)) then begin 
      ;;Grab date from last night using the linux date command
      spawn, "date --date='yesterday' +%y%m%d", date
      date = date[0] ;;de-arrayify date so that it's just a string
  endif      
  
  ;;set the journal for the log

  journal, './log/'+'FLWO_'+date+'_reduction.log'
  PRINT,date

  ;;open z window to save output
  set_plot, 'Z'
  connection = automated_canconnect() ;Test if window connection works
  device = !D.NAME
  splog, 'Can we connect to a window buffer? (1=y, 0=n) ', connection
  splog, 'The device is ', device
  
  ;;now run the pipeline to reduce data and do photometry
  ;automated_work, date, REQUEST=getenv("AUTOM_DIR")+'targets/master_target_list.txt',/FLWO, /clean 
   automated_work, date, REQUEST=getenv("AUTOM_DIR")+'queue/FLWO_target_list.txt',/FLWO, /clean,/clobber 
  ;;finally replot the webpage
  ;automated_replot, /UPDATE_WEBPAGE ;Output light-curves and results
  
  ;;close journal
  journal
  
end
