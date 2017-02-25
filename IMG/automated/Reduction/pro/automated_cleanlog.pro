;+
; NAME:
;   automated_cleanlog
;
; PURPOSE:
; Clean a log file from all the unwanted images
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; date - day to process 
; request - the queue file                
; 
; OPTIONAL INPUTS:
;  
; OUTPUTS: 
;
; Updated the log 
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
; 
; Run for 04 July 2012 with request file file.txt
; automated_cleanlog, '120704', 'path/to/file.txt'
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


PRO automated_cleanlog, date, request


  ;;open the log file
  readcol, 'obs_'+date+'.log', name, info, flag, type,$
           filt, texp, airmass, format='A,A,F,A,A,F,F' 

  ;;open the request
  readcol, request, junk, junk, junk, junk, junk, junk, $
           junk, filt_req, junk, name_req, format='F,F,F,F,F,F,F,A,A,A'

  ;;make type image. Looks through name and changes sci images to
  ;;'ima' bias to bia etc. Stored in shorttype.
  shorttype=strtrim(STRMID(name,0,3),2)
  if(date GT 091201) then begin
     sci_frame=STRMATCH(shorttype,'p*')
     indx=where(sci_frame EQ 1,nn)
     if(nn GT 0) then shorttype[indx]='ima'
  endif
 
  ;;look for bias and dark
  bias=where(shorttype EQ 'bia', nbias)
  dark=where(shorttype EQ 'dar', ndark)


  ;;correct filter after 100211
  if(date gt 100210) then begin
      which_filter=where(shorttype eq 'fla',qua_f)
      if(qua_f gt 0) then $
        tmpfil=filt[which_filter]
        filt[which_filter]=strtrim(STRMID(type[which_filter],5),2)
  endif


  ;;look for different filters
  for i=0, n_elements(filt_req)-1 do begin
     if(i EQ 0) then filt_used=strcompress(filt_req[i],/remove_all)
     if(i GT 0) then begin
        ;;check if new filter
        res=where(filt_used EQ strcompress(filt_req[i],/remove_all),nummatch)
        ;;if new, append
        if(nummatch EQ 0) then filt_used=[filt_used,strcompress(filt_req[i],/remove_all)]
     endif
  endfor
 

  ;;now get flat for that filter
  for i=0, n_elements(filt_used)-1 do begin
     if(i EQ 0) then flat=where(shorttype EQ 'fla' and filt EQ filt_used[i], nflats)
     if(i GT 0) then begin
        ;;check if new filter
        res=where(shorttype EQ 'fla' and filt EQ filt_used[i],nummatch)
        ;;if new, append
        if(nummatch GT 0) then begin
           flat=[flat,res]
           nflats=nflats+nummatch
        endif
     endif
  endfor

  ;;now get science (check names and the filters to avoid duplicates)
  for i=0, n_elements(name_req)-1 do begin
      if(i EQ 0) then begin
         science=where(shorttype EQ 'ima' and type EQ name_req[i] and filt EQ filt_req[i], nsci)
      endif ELSE BEGIN
         ;;if(i GT 0) then begin
         ;;check if new image
         res=where(shorttype EQ 'ima' and type EQ name_req[i] and filt EQ filt_req[i], nummatch)
         ;;if new, append
         if(nummatch GT 0) then begin              
            science=[science,res]
            nsci=nsci+nummatch
            ;print, 'nummatch greater than 0'
         endif
         ;;endif
      endelse
   endfor
  ;;the first (if i eq 0) is quite tricky since you can get a -1 in the
  ;;first place. Remove it
  if(science[0] eq -1) then science=science[1:n_elements(science)-1]

  
  splog, 'Found BIAS', nbias
  splog, 'Found DARK', ndark
  splog, 'Found FLAT', nflats 
  splog, 'Found SCIENCE', nsci
  
  ;;append all indexes
  all=[bias,dark,flat,science]

  
  ;;save all
  forprint, name[all], info[all], flag[all], type[all], filt[all], texp[all], airmass[all],$ 
            textout='obs_'+date+'.log', format='A25,A25,A10,A20,A5,A20,A10'

  
end
