;+ 
; NAME:
; dla_elmlmts
;  V1.0
;
; PURPOSE:
;    Given a DLA structre and ion, figure out the limits
;
; CALLING SEQUENCE:
;   dla_elmlmts, dla, nn
;
; INPUTS:
;   dla -- DLA structure
;   nn  -- Index of the structure
;   i   -- Atomic number
;   j   -- Ion State
;
; RETURNS:
;
; OUTPUTS:
;  llmt -- Lower limit
;  ulmt -- Upper limit
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  flg -- Describes the limit
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
;   06-Oct-2004 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_elmlmts, dla, nn, i, j, llmt, ulmt, FLG=flg

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'dla_elmlmts, dla, nn, i, j, llmt, ulmt, FLG= [v1.0]'
    return
  endif 

  ;; Init
  flg = 0
  llmt = -999.
  ulmt = 999.

  ;;	Abund
  getabnd, nm, i, abnd, flag=1

  for k=1,dla[nn].ion[i].indx[j] do begin
      ;; Lower
      if dla[nn].ion[i].state[j,k].flgclm EQ 3 then begin
          if (flg MOD 2) NE 1 then flg=flg + 1
          x_logclm, dla[nn].ion[i].state[j,k].clm, $
            dla[nn].ion[i].state[j,k].sigclm, ans, sig 
          ans = ans + 12. - abnd - dla[nn].NHI
          llmt = llmt > ans
      endif
      ;; Upper
      if dla[nn].ion[i].state[j,k].flgclm EQ 5 then begin
          if (flg MOD 4) LT 2 then flg=flg + 2
          x_logclm, dla[nn].ion[i].state[j,k].clm, $
            dla[nn].ion[i].state[j,k].sigclm, ans, sig 
          ans = ans + 12. - abnd - dla[nn].NHI
          ulmt = ulmt < ans
      endif
  endfor
  return
end
