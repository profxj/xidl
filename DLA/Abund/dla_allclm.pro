;+ 
; NAME:
; dla_allclm
;  V1.1
;
; PURPOSE:
;  Combines multiple column density measurements of a single ion, e.g.
;  NiII 1741, NiII 1751
;
; CALLING SEQUENCE:
;   dla_allclm, dla, nn, Z, ion, Ni, [sig]
;
; INPUTS:
;  dla -- DLA structure array
;  nn --  Identifies which DLA in the structure array
;  Z  -- Atomic number of the ion
;  ion  -- Ionic number (e.g. 1 for I)
;  Ni  -- Column density now measured [assumed linear value]
;  sig  -- Error in the column [assumed linear value]
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_allclm, dla, nn, Z, ion, Ni, sig, LOG=log

  if (N_params() LT 6) then begin 
    print,'Syntax - ' + $
             'dla_alllcm, dla, nn, Z, ion, Ni, sig, /LOG [v1.1]'
    return
  endif 

  if keyword_set(LOG) then begin
      x_logclm, tmpN, tmpS, Ni, sig, /reverse
      Ni = tmpN
      sig = tmpS
  endif

  case dla[nn].ion[Z].state[ion].flgclm of 
      0: begin  ;; No previous value
          dla[nn].ion[Z].state[ion].clm = Ni
          dla[nn].ion[Z].state[ion].sigclm = sig
      end
      1: begin ;; Create weighted mean
         wold = 1. / dla[nn].ion[Z].state[ion].sigclm^2
         muold = dla[nn].ion[Z].state[ion].clm*wold
	 
         wnew = wold + 1./sig^2
         
          dla[nn].ion[Z].state[ion].clm = (muold + Ni/sig^2) / wnew
          dla[nn].ion[Z].state[ion].sigclm = sqrt(1./wnew)
      end
      2: begin ;; Lower limit
          if dla[nn].ion[Z].state[ion].clm LT Ni then begin
              dla[nn].ion[Z].state[ion].clm = Ni
              dla[nn].ion[Z].state[ion].sigclm = sig
          endif
      end
      else: stop
  endcase

  return
end
