;+ 
; NAME:
; dla_calcclm
;  V1.2
;
; PURPOSE:
;    Sets the elemental abundance in the DLA structure
;
; CALLING SEQUENCE:
;   dla_calcclm, dla, nn, Z, Ni, sig
;
; INPUTS:
;  dla -- DLA structure
;  nn  -- Index of the structure
;  Z  -- Atomic number of the ion
;  Ni  -- Column density now measured [assumed linear value]
;  sig  -- Error in the column [assumed linear value]
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_calcclm, dla, nn, Z, Ni, sig

  if (N_params() LT 5) then begin 
    print,'Syntax - ' + $
             'dla_calclcm, dla, nn, Z, Ni, sig [v1.1]'
    return
  endif 

  case dla[nn].elm[Z].flgclm of 
      0: begin  ;; No previous value
          dla[nn].elm[Z].clm = Ni
          dla[nn].elm[Z].sigclm = sig
      end
      1: begin ;; Create weighted mean
         wold = 1. / dla[nn].elm[Z].sigclm^2
         muold = dla[nn].elm[Z].clm*wold
	 
         wnew = wold + 1./sig^2
         
          dla[nn].elm[Z].clm = (muold + Ni/sig^2) / wnew
          dla[nn].elm[Z].sigclm = sqrt(1./wnew)
      end
      2: begin ;; Lower limit
          if dla[nn].elm[Z].clm LT Ni then begin
              dla[nn].elm[Z].clm = Ni
              dla[nn].elm[Z].sigclm = sig
          endif
      end
      3: begin;; Upper limit 
          if dla[nn].elm[Z].clm GT Ni then begin
              dla[nn].elm[Z].clm = Ni
              dla[nn].elm[Z].sigclm = sig
          endif
      end
      else: stop
  endcase

  return
end
