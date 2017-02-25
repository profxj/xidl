
;+
; NAME:
;   automated_makelog
;
; PURPOSE:
; Procedure that reads the FLWO fits headers and constructs a log file
;
;    
; CALLING SEQUENCE:
;
; automated_makelog, date,
;
;
; INPUTS:
;   
; date        -- string day to process  
;                
; 
; OPTIONAL INPUTS:
;  
;
; OUTPUTS: 
;
; Creates a synthetic log file similar to the one use on superlotis  
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; 
; EXAMPLES:
;
;
; BUGS:
;
; PROCEDURES CALLED:
;
;
; REVISION HISTORY:
;    09-Sept-2012 written by Z. Hughes
;
;-
;------------------------------------------------------------------------------
;; 



PRO automated_makelog,date

 
   SPAWN,'ls | less',name

   indx = WHERE(STRMID(name,0,1) EQ '0')
   name = name[indx]
   n = N_ELEMENTS(indx)

   info = STRARR(n)
   flag = STRARR(n)
   type = STRARR(n)
   filt = STRARR(n)
   texp = STRARR(n)
   space = STRARR(n)
   airmass  = STRARR(n)

   FOR i=0,n-1 DO BEGIN
      header = HEADFITS(name[i])
      info[i] = FXPAR(header,"ST")
      flag[i] = '0.00'
      type[i] = FXPAR(header,"OBJECT")
      filt[i] = FXPAR(header,"FILTER")
      texp[i] = FXPAR(header,"EXPTIME")
      space[i] = '   '
      airmass[i] = FXPAR(header,"AIR")      
      
   ENDFOR

   sub = LINDGEN(n)

  forprint, name[sub], info[sub], flag[sub], type[sub], filt[sub], texp[sub],space[sub],airmass[sub],textout='FLWO_obs_'+date+'.log', format='A23,A20,A5,A20,A5,A10,A10,A10'

   


END
