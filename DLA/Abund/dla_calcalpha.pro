;+ 
; NAME:
; dla_calcalpha
;  V1.1
;
; PURPOSE:
;    Determine the alpha abundance [A/H] by parsing through the
;     elemental abundances in a DLA strucure.  The user can control
;     which element is used by setting the tag flgalpha in advance
;     to calling this routine.
;
; CALLING SEQUENCE:
;   dla_calcalpha, dla, nn
;
; INPUTS:
;   dla -- DLA structure
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
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_calcalpha, dla, nn, NOHSIG=nohsig

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dla_calcalpha, dla, nn [v1.1]'
    return
  endif 

  llmt = -999.
  ulmt = 999.
  lmts = fltarr(3,3)
  lmts[*,1] = -999.
  lmts[*,2] = -999.
  flg = intarr(3)

  ;; AUTOMATIC
  nelc = dla[nn].elm.flgclm
  if dla[nn].flgAlpha GE 0 then begin
      ;; O value
      if nelc[8] EQ 1 then begin
          dla_calcabnd, dla, nn, 8, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgalpha = 5
          goto, jump1
      endif 
      ;; S value
      if nelc[16] EQ 1 then begin
          dla_calcabnd, dla, nn, 16, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgalpha = 4
          goto, jump1
      endif 
      ;; Si value
      if nelc[14] EQ 1 then begin
          dla_calcabnd, dla, nn, 14, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgalpha = 1
          goto, jump1
      endif 
      ;; Si limits
      dla_elmlmts,dla,nn,14,2, d1, d2, FLG=f1
      lmts[1,1] = d1
      lmts[1,2] = d2
      flg[1] = f1
      dla_elmlmts,dla,nn,16,2, d1, d2, FLG=f1
      lmts[2,1] = d1
      lmts[2,2] = d2
      flg[2] = f1
      ;; Combined Limits
      totflg = total(flg)  ; Si + S
      if (totflg GE 4 AND flg[1] NE 2) OR totflg EQ 3 then begin
          llmt = lmts[1,1] > lmts[2,1]
          ulmt = lmts[1,2] < lmts[2,2]
          dla[nn].Alpha = (llmt+ulmt)/2.
          dla[nn].sigAlpha = 0.15 >  (ulmt-dla[nn].Alpha)
          dla[nn].flgAlpha = 13
          if dla[nn].sigAlpha LT 0.3 then return
      endif
      ;; Ar 
      if nelc[18] EQ 1 then begin
          dla_calcabnd, dla, nn, 18, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgalpha = 6
          goto, jump1
      endif 
      ;; Si limits
      dla[nn].flgAlpha = nelc[14]
      if dla[nn].flgAlpha NE 0 then begin
          dla_calcabnd, dla, nn, 14, 1, ans, sig, NOSIGY=NOHSIG
          goto, jump1
      endif
      return
  endif else begin ;;  By choice
      ;; Si value
      if dla[nn].flgAlpha GT (-4) then begin
          dla_calcabnd, dla, nn, l4, 1, ans, sig, NOSIGY=NOHSIG
      endif else begin
          case dla[nn].flgAlpha of
              -4: dla_calcabnd, dla, nn, 16, 1, ans, sig, NOSIGY=NOHSIG ; S
              -5: dla_calcabnd, dla, nn, 8, 1, ans, sig, NOSIGY=NOHSIG ; O
              -6: dla_calcabnd, dla, nn, 18, 1, ans, sig, NOSIGY=NOHSIG ; Ar
              else: stop
          endcase
      endelse
;      dla[nn].Alpha = ans
;      dla[nn].sigAlpha = sig
  endelse

  jump1:
  dla[nn].Alpha = ans 
  dla[nn].sigAlpha = sig

  return
end

