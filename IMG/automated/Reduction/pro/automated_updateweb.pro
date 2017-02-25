
;+
; NAME:
;   automated_phot
;
; PURPOSE:
; 
; Acts as a wrapper for a cron job that is called every morning to 
; automatically update the website, even if no images were
; actually taken by Super-LOTIS
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
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
;           2011  Written by KK
;    05-Jun-2012  Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 

PRO automated_updateweb
  
  JOURNAL, getenv('AUTOM_DIR')+'/log/cron.log'

  ;;Try this in cron job
  
  Set_Plot, 'Z'
  connection = CANCONNECT()     ;Test if window connection works
  device = !D.NAME
  PRINT, 'Can we connect to a window buffer? (1=y, 0=n) ', connection
  PRINT, 'The device is ', device
  
  AUTOMATED_REPLOT, /UPDATE_WEBPAGE, /EMAIL_ALERT ;Output light-curves and results
  ;;Need to automate website updates and emailing out triggers
  
  EXIT
END
