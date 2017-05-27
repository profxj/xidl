;+
; NAME:
;   JULDAYFROMFILENAME
;
;
; PURPOSE:
;   given a filename/directory path/other string containing a date
;   encoded as YYYYmmmDD, where YYYY is 2001-2009, mmm is a string
;   name for the month, and DD is the day within the month, return the
;   corresponding MJD
;
;
; CATEGORY:
;    UTILITIES
;
; CALLING SEQUENCE:
;     MJD = juldayfromfilename(filename)
;
; INPUTS:
;     filename - a filename containing the date to be translated
;               (either in its directory path or within the filename).
;               the date may be isolated within the string by either
;               slashes ('/') or periods ('.').
;
;
; MODIFICATION HISTORY:
;   created 2003feb10
;-


function juldayfromfilename,filename

        searchstring = strsplit(filename[0], '/.',/extract)
        hasdate = strpos(searchstring, '200') eq 0
        whdate = where(hasdate, datect)
        if datect gt 0 then begin
            datename=(searchstring[max(whdate)])
            dateyear=strmid(datename,0,4)
            datemon=strmid(datename,4,3)
            dateday=strmid(datename,7,2)
;            datestring=dateday+'-'+datemon+'-'+dateyear

            MONTHS = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG', $
          'SEP','OCT','NOV','DEC']


            L_MONTH = long(where(strupcase(datemon) eq MONTHS))
            L_MONTH = L_MONTH[0] + 1 ; Scalarize it...
            juldate,[long(dateyear),L_month,long(dateday)],juliandate
            return,juliandate  ;+2400000d0

        endif else return,1E10
        
        

end
