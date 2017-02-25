;+ 
; NAME:
; lls_allclm
;  V1.2
;
; PURPOSE:
;  Combines multiple column density measurements of a single ion, e.g.
;  NiII 1741, NiII 1751
;
; CALLING SEQUENCE:
;   lls_allclm, dla, nn, sys, Z, ion, Ni, [sig]
;
; INPUTS:
;  lls -- LLS structure array
;  nn --  Identifies which DLA in the structure array
;  sys -- Identifies the sub-system to analyze
;  Z  -- Atomic number of the ion
;  ion  -- Ionic number (e.g. 1 for I)
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
;   04-Apr-2006 Written by JXP
;- 
;------------------------------------------------------------------------------
pro lls_allclm, lls, nn, sys, Z, ion, Ni, sig, LOG=log

  if (N_params() LT 6) then begin 
    print,'Syntax - ' + $
             'lls_alllcm, lls, nn, sys, Z, Ni, sig, /LOG [v1.1]'
    return
  endif 

  if keyword_set(LOG) then begin
      x_logclm, tmpN, tmpS, Ni, sig, /reverse
      Ni = tmpN
      sig = tmpS
  endif

  case lls[nn].systems[sys].ion[Z].state[ion].flgclm of 
      0: begin  ;; No previous value
          lls[nn].systems[sys].ion[Z].state[ion].clm = Ni
          lls[nn].systems[sys].ion[Z].state[ion].sigclm = sig
      end
      1: begin ;; Create weighted mean
         wold = 1. / lls[nn].systems[sys].ion[Z].state[ion].sigclm^2
         muold = lls[nn].systems[sys].ion[Z].state[ion].clm*wold
	 
         wnew = wold + 1./sig^2
         
          lls[nn].systems[sys].ion[Z].state[ion].clm = (muold + Ni/sig^2) / wnew
          lls[nn].systems[sys].ion[Z].state[ion].sigclm = sqrt(1./wnew)
      end
      2: begin ;; Lower limit
          if lls[nn].systems[sys].ion[Z].state[ion].clm LT Ni then begin
              lls[nn].systems[sys].ion[Z].state[ion].clm = Ni
              lls[nn].systems[sys].ion[Z].state[ion].sigclm = sig
          endif
      end
      else: stop
  endcase

  return
end
