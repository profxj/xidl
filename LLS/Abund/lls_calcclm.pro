;+ 
; NAME:
; lls_calcclm
;  V1.2
;
; PURPOSE:
;    Sets the elemental abundance in the LLS structure
;
; CALLING SEQUENCE:
;   lls_calcclm, lls, nn, sys, Z, Ni, sig
;
; INPUTS:
;  lls -- LLS structure
;  nn  -- Index of the structure
;  sys -- Identifies the sub-system of the LLS
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
;   lls_allabd, lls
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Jul-2008 Written by JXP
;- 
;------------------------------------------------------------------------------
pro lls_calcclm, lls, nn, sys, Z, Ni, sig

  if (N_params() LT 6) then begin 
    print,'Syntax - ' + $
             'lls_calclcm, lls, nn, Z, Ni, sig [v1.1]'
    return
  endif 

  case lls[nn].systems[sys].elm[Z].flgclm of 
      0: begin  ;; No previous value
          lls[nn].systems[sys].elm[Z].clm = Ni
          lls[nn].systems[sys].elm[Z].sigclm = sig
      end
      1: begin ;; Create weighted mean
         wold = 1.d / lls[nn].systems[sys].elm[Z].sigclm^2
         muold = lls[nn].systems[sys].elm[Z].clm*wold
	 
         wnew = wold + 1.d/sig^2
         
          lls[nn].systems[sys].elm[Z].clm = (muold + Ni/sig^2) / wnew
          lls[nn].systems[sys].elm[Z].sigclm = sqrt(1./wnew)
      end
      2: begin ;; Lower limit
          if lls[nn].systems[sys].elm[Z].clm LT Ni then begin
              lls[nn].systems[sys].elm[Z].clm = Ni
              lls[nn].systems[sys].elm[Z].sigclm = sig
          endif
      end
      3: begin;; Upper limit 
          if lls[nn].systems[sys].elm[Z].clm GT Ni then begin
              lls[nn].systems[sys].elm[Z].clm = Ni
              lls[nn].systems[sys].elm[Z].sigclm = sig
          endif
      end
      else: stop
  endcase

  return
end
