;+ 
; NAME:
; x_filters   
;    Version 1.0
;
; PURPOSE:
;    Grabs out a list of filters from a list of images
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
        'x_filters, allfilt, filter, [nfilt], NOSORT= (v1.0)'
      return
  endif 


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
  
      
