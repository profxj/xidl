FUNCTION caha_hdr, header

;+
; NAME:
;     HGREP
;
; PURPOSE:
;       Find a substring in a FITS header (or any other string array)
;
; CALLING SEQUENCE:
;       HGREP, header, substring, [/KEEPCASE, /LINENUM ]
;
; INPUTS: 
;       header -  FITS header or other string array
;       substring - scalar string to find in header; if a numeric value is 
;                 supplied, it will be converted to type string
;
; OPTIONAL INPUT KEYWORDS:
;       /KEEPCASE: if set, then look for an exact match of the input substring 
;                 Default is to ignore case .
;       /LINENUM: if set, prints line number of header in which
;                substring appears 
;
; OUTPUTS:
;       None, results are printed to screen
;
; EXAMPLE: 
;       Find every place in a FITS header that the word 'aperture'
;       appears in lower case letters and print the element number 
;       of the header array:
;       
;       IDL> hgrep, header, 'aperture', /keepcase, /linenum
;
; HISTORY: 
;       Written, Wayne Landsman (Raytheon ITSS)      August 1998
;       Adapted from STIS version by Phil Plait/ ACC November 14, 1997
;       Remove trailing spaces if a non-string is supplied W. Landsman Jun 2002
;-

   if (N_params() LT 1) then begin
      print,'Syntax - HGREP, header'
      return, 0
   endif

   if N_elements(header) eq 0 then begin
      print,'first parameter not defined. Returning...'
      return, 0
   endif
   ;;hh = header
   hh = strtrim(header,2)
   substring = 'HIERARCH'

   if keyword_set(keepcase) then $
         flag = strpos(hh,substring) $
   else  flag = strpos(strlowcase(hh),strlowcase(substring))
     

   g = where(flag NE -1, Ng)
   ;; remove comment
   hh_new = hh
   FOR i = 0L, Ng-1L DO BEGIN
      temp = strsplit(hh_new[g[i]], '/', /extract)
      temp1 = strjoin(strsplit(temp[0], "'", /extract))
      comment = temp[1]
      ind =  strsplit(temp1, '=')
      value = strmid(temp1, ind[1])
      temp2 =  strmid(temp1,ind[0],ind[1]-1L)
      ;;temp3 = strsplit(temp2, ' ', /extract)
      ;;temp4 = temp3[n_elements(temp3)-1L]
      temp3 = strmid(temp2, 18)
      temp4 = STRJOIN(STRSPLIT(temp3, /EXTRACT), '_') 
      nam = strcompress(temp4, /rem)
      ;;while strlen(nam) LT 8 do nam = nam + ' ' ;Make 8 chars long
      sxaddpar, hh_new, string(nam), value, comment
      ;;hh_new = nam + value + '/ ' + comment
   ENDFOR
   
   return, hh_new
   end
