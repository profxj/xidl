;+
; NAME:
;   long_imacs_reduce
;
; PURPOSE:
;   Reduce the IMACS mosaic by looping over chip by chip
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
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
;   2014 Nov 19 Written by JXP 
;-
;------------------------------------------------------------------------------

pro long_imacs_reduce, planfile, CHIPS=chips, _EXTRA=extra, ICHIP=ichip

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'long_imacs_reduce, plan_fil  [v1.0]'
      return
  endif 

  if not keyword_set(CHIPS) then chips = lindgen(8) + 1
  nchip = n_elements(chips)

  ;; Read single plan file
  planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
  nfile = n_elements(planstr)

  ;; Loop on the chips
  if not keyword_set(ICHIP) then ichip = 0
  for qq=ichip,nchip-1 do begin
     cchip = 'c'+strtrim(chips[qq],2)
     ;; Loop on files
     newplan = planstr
     for ii=0L,nfile-1 do begin
        ;; Update filename
        ipos = strpos(newplan[ii].filename, '.fits')
        newplan[ii].filename = strmid(newplan[ii].filename, 0, ipos-2)+cchip+'.fits'
     endfor
     ;; Giddy up
     if cchip EQ 'c3' or cchip EQ 'c8'  or cchip eq 'c4' or cchip eq 'c7' then skytrace = 0L

     if cchip EQ 'c4' or cchip EQ 'c7' then begin
        NO_HEFIXES=1
        CORR_PEAK = 0.35
        MED_ERR = 0.17
     endif else begin
        NO_HEFIXES=0
        CORR_PEAK = 0.8
        MED_ERR = 0.1
     endelse
     long_reduce, planfile, planstr=newplan, _EXTRA=extra, skytrace=skytrace, NO_HEFIXES=no_hefixes, $
                  CORR_PEAK=corr_peak, MED_ERR=med_err
     stop
  endfor
  
  return
END
