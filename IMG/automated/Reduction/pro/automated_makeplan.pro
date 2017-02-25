;+
; NAME:
;   automated_makeplan
;
; PURPOSE:
; 
; Store relevant info in data structure 
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; logfile - the observing log
;             
; 
; OPTIONAL INPUTS:
;
; path -- path to work in (not active)
;
; OUTPUTS: 
;
; A data structure with image types
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; automated_makeplan, 'path/to/file.log'
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Written and Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 

pro automated_makeplan, logfile,  path=path
  
  if not keyword_set(path) then path='./'
  strname=logfile+'.str'
  
;;get a list of images
  readcol, logfile, name, date, void, object, filter, texp,$
           airm, format='A,A,F,A,A,F,F', comment='#'
  num=n_elements(name)
  
;;get gain and readnoise 
  header=headfits(name[0],exten=0)
  gain=fxpar(header,"GAIN")
  rn=fxpar(header,"RDNOISE")
  
  gain=replicate(gain,num)
  rn=replicate(rn,num)
  
  splog, "Gain is ", gain[0]
  splog, "RN is ", rn[0]
  
;;create type from names
  type=strtrim(STRMID(name,0,3),2)
  date=fxpar(header,"DATE")
  yy=STRMID(date,2,2) 
  mm=STRMID(date,5,2) 
  dd=STRMID(date,8,2)
  date=yy+mm+dd
  
  if(date GT 091201) then begin
     Sci_frame=STRMATCH(type,'p*')
     indx=where(Sci_frame EQ 1,nn)
     if(nn GT 0) then type[indx]='ima'
  endif
    
;;prepare structure 
  str={NAME:name,OBJECT:object,FILTER:filter,EXTIME:texp,AIRMASS:airm,GAIN:gain,$
       RN:rn,TYPE:type}
  mwrfits, str, strname, /create
  splog, "Info structure written to  ", strname
  
end
