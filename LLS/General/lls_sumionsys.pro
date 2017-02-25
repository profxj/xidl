;+ 
; NAME:
; lls_sumionsys
;  V1.1
;
; PURPOSE:
;    Sums up the ionic columns for the subsystems
;
; CALLING SEQUENCE:
;   
;   lls_sumionsys, lls, nn
;
; INPUTS:
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-July-2013 Written by JXP
;-
;------------------------------------------------------------------------------
pro lls_sumionsys, lls, nn 

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'lls_sumionsys, lls, nn [v1.0]' 
    return
  endif 

  sz = size(lls[nn].systems.ion.state, /dim)

  ;; Find those with a measurement
  val = where( lls[nn].systems.ion.state.flgclm GT 0 and $
               lls[nn].systems.ion.state.lambda LT 1.,nval) 

  ;; Parse ion, state
  allZ = val / (sz[0]*sz[1]) MOD sz[2]
  allion = val MOD (sz[0]*sz[1])
  
  ;; Loop on ions
  for qq=0L,nval-1 do begin
     clm = 0.
     flg = 0
     Z = allZ[qq]
     ion = allion[qq]

     ;; At least one measurement or lower limit
     detect = where(lls[nn].systems.ion[Z].state[ion].flgclm EQ 1 OR $ 
                    lls[nn].systems.ion[Z].state[ion].flgclm EQ 2, ndet)

     if ndet GT 0 then begin
        for jj=0L,ndet-1 do begin
           ii = detect[jj]
           ;; Flag
           lls[nn].ion[Z].state[ion].flgclm = lls[nn].ion[Z].state[ion].flgclm > $
                                              lls[nn].systems[ii].ion[Z].state[ion].flgclm 
           ;; Sum
           x_logclm, clm, sig, $
                     lls[nn].systems[ii].ion[Z].state[ion].clm, lls[nn].systems[ii].ion[Z].state[ion].sigclm, /rever
           if jj EQ 0 then begin
              lls[nn].ion[Z].state[ion].clm = clm
              lls[nn].ion[Z].state[ion].sigclm = sig
           endif else begin
              wold = 1.d / lls[nn].ion[Z].state[ion].sigclm^2
              muold = lls[nn].ion[Z].state[ion].clm*wold
              wnew = wold + 1.d/sig^2
              
              lls[nn].ion[Z].state[ion].clm = (muold + clm/sig^2) / wnew
              lls[nn].ion[Z].state[ion].sigclm = sqrt(1./wnew)
           endelse
        endfor
        ;; Log and flag
        x_logclm, lls[nn].ion[Z].state[ion].clm, lls[nn].ion[Z].state[ion].sigclm, logN, logsig
        lls[nn].ion[Z].state[ion].clm = logN
        lls[nn].ion[Z].state[ion].sigclm = logsig
     endif else begin ;; Only upper limits
        upper = where(lls[nn].systems.ion[Z].state[ion].flgclm EQ 3, nupper)
        if nupper EQ 0 then stop
        lls[nn].ion[Z].state[ion].flgclm = 3
        ;; Summing over the values / sqrt(N)  [A bit kludgy]
        lls[nn].ion[Z].state[ion].clm = alog10( total(10.^lls[nn].systems[upper].ion[Z].state[ion].clm ) $
                                                /sqrt(nupper))
        lls[nn].ion[Z].state[ion].sigclm = 0.
     endelse
  endfor
  
  return
end
