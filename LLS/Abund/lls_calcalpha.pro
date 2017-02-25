;+ 
; NAME:
; lls_calcalpha
;  V1.1
;
; PURPOSE:
;    Determine the alpha abundance [A/H] by parsing through the
;     elemental abundances in a DLA strucure.  The user can control
;     which element is used by setting the tag flgalpha in advance
;     to calling this routine.
;
; CALLING SEQUENCE:
;   lls_calcalpha, lls, nn
;
; INPUTS:
;   lls -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOHSIG -- Do not include uncertainties in NHI in the calculation
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Values of flgalpha --   
;                       1 -- Si abundance
;                       2 -- Si lower limit
;                       3 -- Si upper limit
;                       4 -- S abundance
;                       5 -- O abundance
;                       6 -- Ar abundance
;                       -4 -- Use S abundance
;                       -5 -- Use O abundance
;                       -6 -- Use Ar abundance
;
; EXAMPLES:
;   lls_allabd, lls
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   21-Jul-2008 Written by JXP
;- 
;------------------------------------------------------------------------------
pro lls_calcalpha, lls, nn, sys, NOHSIG=nohsig

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lls_calcalpha, lls, nn, sys [v1.1]'
    return
  endif 

  llmt = -999.
  ulmt = 999.
  lmts = fltarr(3,3)
  lmts[*,1] = -999.
  lmts[*,2] = -999.
  flg = intarr(3)

  ;; AUTOMATIC
  nelc = lls[nn].systems[sys].XH.flgclm
  if lls[nn].systems[sys].flg_Alpha GE 0 then begin
      ;; O value
      if nelc[8] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[8].clm
          sig = lls[nn].systems[sys].XH[8].sigclm
          lls[nn].systems[sys].flg_alpha = 5
          goto, jump1
      endif 
      ;; Si value
      if nelc[14] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[14].clm
          sig = lls[nn].systems[sys].XH[14].sigclm
          lls[nn].systems[sys].flg_alpha = 1
          goto, jump1
      endif 
      ;; S value
      if nelc[16] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[16].clm
          sig = lls[nn].systems[sys].XH[16].sigclm
          lls[nn].systems[sys].flg_alpha = 4
          goto, jump1
      endif 
      ;; Si+S limits
;      lls_elmlmts,lls,nn,14,2, d1, d2, FLG=f1
;      lmts[1,1] = d1
;      lmts[1,2] = d2
;      flg[1] = f1
;      lls_elmlmts,lls,nn,16,2, d1, d2, FLG=f1
;      lmts[2,1] = d1
;      lmts[2,2] = d2
;      flg[2] = f1
      ;; Combined Limits
;      totflg = total(flg)  ; Si + S
;      if (totflg GE 4 AND flg[1] NE 2) OR totflg EQ 3 then begin
;          llmt = lmts[1,1] > lmts[2,1]
;          ulmt = lmts[1,2] < lmts[2,2]
;          lls[nn].Alpha = (llmt+ulmt)/2.
;          lls[nn].sigAlpha = 0.15 >  (ulmt-lls[nn].Alpha)
;          lls[nn].flgAlpha = 13
;          if lls[nn].sigAlpha LT 0.3 then return
;      endif
      ;; Ar 
;      if nelc[18] EQ 1 then begin
;          lls_calcabnd, lls, nn, 18, 1, ans, sig, NOSIGY=NOHSIG
;          lls[nn].flgalpha = 6
;          goto, jump1
;      endif 
      ;; Si limits [JXP; August 2011]
      ans = lls[nn].systems[sys].XH[14].clm
      sig = lls[nn].systems[sys].XH[14].sigclm
      lls[nn].systems[sys].flg_alpha = nelc[14]
;      lls[nn].systems[sys].flg_Alpha = nelc[14]
;      if lls[nn].systems[sys].flg_Alpha NE 0 then begin
;          lls_calcabnd, lls, nn, sys, 14, 1, ans, sig, NOSIGY=NOHSIG
;          goto, jump1
;      endif
      return
  endif else begin ;;  By choice
      stop
      ;; Si value
      if lls[nn].systems[sys].flg_Alpha GT (-4) then begin
          lls_calcabnd, lls, nn, l4, 1, ans, sig, NOSIGY=NOHSIG
      endif else begin
          case lls[nn].systems[sys].flg_Alpha of
              -4: lls_calcabnd, lls, nn, 16, 1, ans, sig, NOSIGY=NOHSIG ; S
              -5: lls_calcabnd, lls, nn, 8, 1, ans, sig, NOSIGY=NOHSIG ; O
              -6: lls_calcabnd, lls, nn, 18, 1, ans, sig, NOSIGY=NOHSIG ; Ar
              else: stop
          endcase
      endelse
;      lls[nn].Alpha = ans
;      lls[nn].sigAlpha = sig
  endelse

  jump1:
  lls[nn].systems[sys].AlphaH = ans 
  lls[nn].systems[sys].sig_AlphaH = sig

  return
end

