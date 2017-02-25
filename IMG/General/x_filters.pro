;+ 
; NAME:
; x_filters   
;    Version 1.1
;
; PURPOSE:
;    Grabs out a list of unique filters from a full string
;  array list of filters.
;
; CALLING SEQUENCE:
;   
; x_filters, allfilt, filter, [nfilt], NOSORT=nosort
;
; INPUTS:
;   allfilt - List of all images
;
; RETURNS:
;   filter - Filter list
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOSORT -- Dont bother sorting the unique filter list
;
; OPTIONAL OUTPUTS:
;  nfilt - Number of unique filters
;
; COMMENTS:
;
; EXAMPLES:
;   x_filters, img, filter, nfilt
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-July-2001 Written by JXP
;   08-Aug-2001 Revised to sort
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_filters, allfilt, filter, nfilt, NOSORT=nosort

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'x_filters, allfilt, filter, [nfilt], NOSORT= (v1.1)'
      return
  endif 


  ;; I should replace the following with one uniq command!
  nimg = n_elements(allfilt)
  tmp = strarr(100)
  tmp[1] = allfilt[0]
  nfilt = 1
  for i=0,nimg-1 do begin
      flgfilt = 0
      for j=1,nfilt do begin
          if allfilt[i] EQ tmp[j] then begin
              flgfilt = 1
              break
          endif
      endfor
      if flgfilt NE 1 then begin
          nfilt = nfilt + 1
          tmp[nfilt] = allfilt[i]
      endif
  endfor
  filter = strarr(nfilt)
  for i=0,nfilt-1 do filter[i] = tmp[i+1]

  ; Sort

  if not keyword_set( NOSORT ) then begin
      isort = intarr(nfilt)
      for i=0,nfilt-1 do begin
          case strtrim(filter[i],2) of 
              'U' : isort[i] = 10
              'B' : isort[i] = 20
              'V' : isort[i] = 30
              'Rs' : isort[i] = 35
              'R' : isort[i] = 40
              'I' : isort[i] = 50
              'J' : isort[i] = 60
              'H' : isort[i] = 70 ;Added for Nirc2 images
              'Kp' : isort[i] = 80 ;Added for Nirc2 images
              else : begin
                  print, 'x_filters: Cant handle filter ', filter[i], $
                    ' for sorting!'
                  return
              end
          endcase
      endfor
      finsort = sort(isort)
      filter = strtrim(filter[finsort],2)
  endif
      
  return
end
  
      
