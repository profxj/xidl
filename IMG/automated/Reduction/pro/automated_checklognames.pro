;+
; NAME:
;   automated_checklognames
;
; PURPOSE:
; 
; Allow multiple names for one object. 
; Format new_old_names.txt
; #NewName  OldName
; M87       NGC_4486
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; logfile - the observing log
;             
; OPTIONAL INPUTS:
;
; OUTPUTS: 
;
; A new log file with updated name
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; Current version supports only one new name.
;
; EXAMPLES:
;
; Open file file.log
; automated_checklognames, 'path/to/file.log' 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


pro automated_checklognames, logfile
  
  ;;Read in log file
  readcol, logfile, name, date, void, object, filter, texp, airm, FORMAT='A,A,F,A,A,F,F', COMMENT='#' 

  ;;Read in new/old names catalog
  spawn, 'cp '+logfile+' orig_'+logfile
  readcol, getenv("AUTOM_DIR")+'targets/new_old_names.txt', new_names, old_names, COMMENT='#', FORMAT='A,A' 

  ;;Replace old names with new names
  n_names = n_elements(object)
  for i=0, n_names-1 do begin
     ;;find index of old name for name currently being checked in the log
     j = where(old_names eq object[i]) 
     if j ne [-1] then object[i] = new_names[j]
  endfor
  
  ;;Output log file with new names
  writecol, logfile, name, date, void, object, filter,$
            texp, airm, FMT='(A,X,A,X,F,X,A,X,A,X,F,X,F)' 
 
end
