;+ 
; NAME:
; hires_setgainrn
;     Version 1.1
;
; PURPOSE:
;    Sets the gain and read noise for the HIRES Mosaic using some
;    standard, archived values.
;
; CALLING SEQUENCE:
;
; INPUTS:
;  HIRES -- HIRES structure
;  IDX   -- Index value in the structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP  (These values are only good for after
;   10/27/2004)
;   22-Aug-2005 Updated gain for pre-Oct 2004 based on my own
;   measurements
;-
;------------------------------------------------------------------------------

pro hires_setgainrn, hires, idx, ccdgain

  ;;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_setgainrn, hires, idx, ccdgain [v1.1]'
      return
  endif 

  ;; Date
  if hires[idx[0]].date GT 2453303.6d then flg_date = 1 else flg_date = 0

  ;; 
  case ccdgain of
      'low': begin
          ;; Blue
          blu = where(hires[idx].chip EQ 1, nblu)
          if nblu NE 0 then begin
              case flg_date of
                  0: hires[idx[blu]].gain = 1.1
                  1: hires[idx[blu]].gain = 2.0
              endcase
              hires[idx[blu]].readno = 2.9
          endif
          ;; Green
          gre = where(hires[idx].chip EQ 2, ngre)
          if ngre NE 0 then begin
              case flg_date of
                  0: hires[idx[gre]].gain = 1.24
                  1: hires[idx[gre]].gain = 2.2 
              endcase
              hires[idx[gre]].readno = 2.9
          endif
          ;; Red
          red = where(hires[idx].chip EQ 3, nred)
          if nred NE 0 then begin
              case flg_date of
                  0: hires[idx[red]].gain = 1.15
                  1: hires[idx[red]].gain = 2.9 
              endcase
              hires[idx[red]].readno = 4.2
          endif
      end
      'high': begin
          ;; Blue
          blu = where(hires[idx].chip EQ 1, nblu)
          if nblu NE 0 then begin
              hires[idx[blu]].gain = 0.78
              hires[idx[blu]].readno = 2.3
          endif
          ;; Green
          gre = where(hires[idx].chip EQ 2, ngre)
          if ngre NE 0 then begin
              hires[idx[gre]].gain = 0.84
              hires[idx[gre]].readno = 2.4
          endif
          ;; Red
          red = where(hires[idx].chip EQ 3, nred)
          if nred NE 0 then begin
              hires[idx[red]].gain = 0.89
              hires[idx[red]].readno = 2.8
          endif
      end
      'SINGLE': begin
         hires[idx].gain = 2.45
         hires[idx].readno = 6.0
      end
      else: stop
  endcase
  return

end
    
