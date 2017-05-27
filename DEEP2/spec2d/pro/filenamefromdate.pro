;+
; NAME:
;   FILENAMEFROMDATE
;
;
; PURPOSE:
;   given a date, e.g. '2002-01-03', give a filename/directory
;   path/other string containing a date
;   encoded as YYYYmmmDD, where YYYY is 2001-2009, mmm is a string
;   name for the month, and DD is the day within the month, return the
;   corresponding MJD
;
;
; CATEGORY:
;    UTILITIES
;
; CALLING SEQUENCE:
;     filedate = filenamefromdate(date)
;
; INPUTS:
;     date - a string containing the date to be translated
;
; MODIFICATION HISTORY:
;   created 2005mar02
;-


function filenamefromdate,date

if n_elements(date) gt 1 then begin
    array=strarr(n_elements(date))
    for i=0L,n_elements(date)-1 do $
        array[i]=filenamefromdate(date[i])

    return,array
endif

        searchstring = strsplit(date, '-',/extract)
        hasdate = strpos(searchstring, '200') eq 0
        whdate = where(hasdate, datect)
        if datect gt 0 then begin
;            datename=(searchstring[max(whdate)])
;            dateyear=strmid(datename,0,4)
;            datemon=strmid(datename,4,3)
;            dateday=strmid(datename,7,2)
;;            datestring=dateday+'-'+datemon+'-'+dateyear

            dateyear=searchstring[0]
            datemon=searchstring[1]
            dateday=searchstring[2]

            MONTHS = strlowcase(['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG', $
          'SEP','OCT','NOV','DEC'])

            stringout=dateyear+months[fix(datemon)-1]+dateday
            

;            L_MONTH = long(where(strupcase(datemon) eq MONTHS))
;            L_MONTH = L_MONTH[0] + 1 ; Scalarize it...
;            juldate,[long(dateyear),L_month,long(dateday)],juliandate
            return,stringout

        endif else return,'1999jan01'
        
        

end
