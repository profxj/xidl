;+ 
; NAME:
; mike_shifti
;     Version 1.1
;
; PURPOSE:
;  Shift the order position given the values of mike.arc_xyoff
;  This is used to deal with the thermal variations within MIKE.
;
; CALLING SEQUENCE:
;  shft = mike_shifti(arc_xyoff, OSTR=ostr, ORDRS=ordrs)
;
; INPUTS:
;   arc_xyoff - 2 parameter fit to the shift (linear fit)
;
; RETURNS:
;   shft - The shift for the set of orders
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ORDRS --  Physical order numbers
;  OSTR= -- If provided, the routine will update lhedge, rhedge
;              tags
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  shft = mike_shifti(75L, mike[23].arc_xyoff)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-Jun-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function mike_shifti, fitprm, OSTR=ostr, ORDRS=ordrs

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'shft = mike_shifti(fitprm, OSTR=, ORDRS=) [v1.1] ' 
      return, -1
  endif

  if not keyword_set( OSTR ) and not keyword_set( ORDRS ) then stop
  if not keyword_set( ORDRS ) then ordrs = ostr.order 

  ;; Assume linear fit
  shft = poly(ordrs, fitprm)

  ;; Offset OSTR if necessary
  if keyword_set( OSTR ) then begin
      nordr = n_elements(ordrs)
      for ii=0L,nordr-1 do begin
          
          ;; Find the order
          q = where(ostr.order EQ ordrs[ii],nq)
          if nq EQ 0 then stop
          
          ;; Offset edge
          ostr[q].lhedg = ostr[q].lhedg + shft[ii]
          ostr[q].rhedg = ostr[q].rhedg + shft[ii]

          ;; Offset coefficient
          ostr[q].lhc[0] = ostr[q].lhc[0] + shft[ii]
          ostr[q].rhc[0] = ostr[q].rhc[0] + shft[ii]
      endfor
  endif

  return, shft
end

