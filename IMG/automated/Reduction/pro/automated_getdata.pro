;+
; NAME:
;   automated_getdata
;
; PURPOSE:
; Procedure that grab the data and log from a remote site
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; date - day to process 
;        FORMAT YYMMDD
; 
; OPTIONAL INPUTS:
;  
;
; OUTPUTS: 
;
; Untar the observed images
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

pro automated_getdata, date
  
 ;;get the tar files (the elegance of programming...)
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'A.tar.gz'
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'B.tar.gz'
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'C.tar.gz'
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'D.tar.gz'
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'.tar.gz'
  spawn, 'wget http://slotis.kpno.noao.edu/~slotis/data/slotis_'+date+'.tar' ;;Incase super-lotis forgets to gzip it
  
  ;;get the log file
  spawn, 'wget http://slotis.kpno.noao.edu/logs/observations/obs_'+date+'.log'


  ;;untar archive file
  splog, 'Extract images...'
  spawn, 'tar -xf *.tar'
  spawn, 'gunzip *.gz'
  spawn, 'ls slotis*.tar', tarfile
  for i=0, n_elements(tarfile)-1 do begin
     splog, 'Untar '+tarfile[i]
     spawn, 'tar -xf '+tarfile[i]
  endfor

end
