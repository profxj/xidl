;+ 
; NAME:
; esi_echthruput
;     Version 1.0
;
; PURPOSE:
;    Produce the web page for ESI througput measurements.
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with rdxstdered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echrdxstd, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echthruput, std_fil, CLOBBER=clobber

;
;  if  N_params() LT 3  then begin 
;      print,'Syntax - ' + $
;        'esi_echtrcstd, esi, indx, guide [v1.0]'
;      return
;  endif 

  ;; Reduce + Extract standard

  if keyword_set(STD_FIL) then begin
      if x_chkfil(STD_FIL+'*') NE 1 then begin
          print, 'No ESI file '+std_fil
          return
      endif
      ;; Parse header to create output file name
      head = xheadfits(std_fil)
      date = sxpar(head, 'DATE-OBS')
      date = strmid(date,0,4)+x_getmonth(long(strmid(date,5,2)))+$
             strmid(date,8,2)
      frame = sxpar(head,'FRAMENO')
      outfil = getenv('ESI_CALIBS') + '/std_'+date+'_'+strtrim(frame,2)+'.fits'

      ;; Check
      if x_chkfil(outfil+'*') and not keyword_set(CLOBBER) then begin
          print, 'esi_echthruput: File '+outfil+' exists.  Use /CLOBB to overwrite.'
      endif else begin
          ;; Process
          esi_echqckrdx, std_fil, /STD
          ;; Clean up
      endelse
  endif

end
              
      
      
