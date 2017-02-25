;+ 
; NAME:
; lls_calcfeh
;  V1.2
;
; PURPOSE:
;    Determine the Fe-peak abundance [Fe/H] by parsing through the
;     elemental abundances in a DLA structure.  The user can control
;     which element is used by setting the tag flgFe in advance
;     to calling this routine.
;
; CALLING SEQUENCE:
;   lls_calcfeh, lls, nn
;
; INPUTS:
;   lls -- DLA structure
;   nn  -- Index of the structure
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
;  Values of flgFe --   
;                       1 -- Fe abundance
;                       2 -- Fe lower limit
;                       3 -- Fe upper limit
;                       4 -- Ni abundance offset by -0.1
;                       5 -- Cr abundance offset by -0.2
;                       6 -- Al abundance 
;                       11 -- Fe limits from a pair of transitions
;                       13 -- Limit from Fe+Ni
;                       -4 -- Use the Ni abundance
;                       -5 -- Use the Cr abundance
;                       -6 -- Use the Al abundance
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
pro lls_calcfeh, lls, nn, sys, NOHSIG=nohsig, FINE=fine

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lls_calcfeh, lls, nn, sys, /NOHSIG, /FINE [v1.1]'
    return
  endif 

  llmt = -999.
  ulmt = 999.
  lmts = fltarr(5,3)
  lmts[*,1] = -999.
  lmts[*,2] = -999.
  flg = intarr(5)

  nelc = lls[nn].systems[sys].XH.flgclm
  ;; AUTOMATIC
  if lls[nn].systems[sys].flgFe GE 0 then begin
      ;; Fe value
      if nelc[26] EQ 1 then begin
          if keyword_set(FINE) then begin
              feidx = [41 + lindgen(4)]
              gdfe = where(lls[nn].ion[26].state[feidx].flgclm EQ 1, ngdf)
              if ngdf GT 1 then begin
                  lls[nn].elm[26].clm = lls[nn].elm[26].clm + $
                    total(lls[nn].ion[26].state[feidx[gdfe]].clm)
              endif 
          endif
          ans = lls[nn].systems[sys].XH[26].clm
          sig = lls[nn].systems[sys].XH[26].sigclm
          lls[nn].systems[sys].flgFe = 1
          goto, jump1
      endif 
      ;; Ni value
      if nelc[28] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[28].clm
          sig = lls[nn].systems[sys].XH[28].sigclm
          lls[nn].systems[sys].flgFe = 4
          goto, jump1
      endif 
      ;; Fe limits
;      lls_elmlmts, lls, nn, 26, 2, d1, d2, FLG=f1
;      lmts[1,1] = d1
;      lmts[1,2] = d2
;      flg[1] = f1
;      if flg[1] EQ 3 then begin
;          lls[nn].FeH = (lmts[1,1] + lmts[1,2])/2.
;          lls[nn].sigFeH = 0.1 >  (lmts[1,2]-lls[nn].FeH)
;          lls[nn].flgFe = 11
;          if lls[nn].sigFeH LT 0.3 then return
;      endif
      ;; Fe+Ni
;      lls_elmlmts, lls, nn, 28, 2, d1, d2, FLG=f1 ;; Ni
;      lmts[2,1] = d1
;      lmts[2,2] = d2
;      flg[2] = f1
      ;; check limits
;      if(flg[1]+flg[2] GT 2) AND ( (flg[1]+flg[2] MOD 2) EQ 1) then begin 
;          llmt = lmts[1,1] > (lmts[2,1]-0.1)
;          ulmt = lmts[1,2] < (lmts[2,2]-0.1)
;          lls[nn].FeH = (llmt+ulmt)/2.
;          lls[nn].sigFeH = 0.15 >  (ulmt-lls[nn].FeH)
;          lls[nn].flgFe = 13
;          if lls[nn].sigFeH LT 0.3 then return
;      endif
      ;; Cr 
      if nelc[24] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[24].clm
          sig = lls[nn].systems[sys].XH[24].sigclm
          lls[nn].systems[sys].flgFe = 5
          goto, jump1
      endif
      ;; Al
      if nelc[13] EQ 1 then begin
          ans = lls[nn].systems[sys].XH[13].clm
          sig = lls[nn].systems[sys].XH[13].sigclm
          lls[nn].systems[sys].flgFe = 6
          goto, jump1
      endif
      ;; All Limits
;      lls_elmlmts,lls, nn,24,2, d1, d2, FLG=f1 ; Cr
;      lmts[3,1] = d1
;      lmts[3,2] = d2
;      flg[3] = f1
;      lls_elmlmts,lls, nn,13,2, d1, d2, FLG=f1 ; Al
;      lmts[4,1] = d1
;      lmts[4,2] = d2
;      flg[4] = f1
;      if total(flg MOD 2) GE 1 AND total(flg GT 1) GT 1 then begin
;          llmt = lmts[1,1] > (lmts[2,1]-0.1)
;          llmt = llmt > (lmts[3,1]-0.2)
;          llmt = llmt > lmts[4,1]
;          ulmt = lmts[1,2] < (lmts[2,2]-0.1)
;          ulmt = ulmt < (lmts[3,2]-0.2)
;          ulmt = ulmt < lmts[4,2]
;          ;; values
;          lls[nn].FeH = (llmt+ulmt)/2.
;          lls[nn].sigFeH = 0.15 >  (ulmt-lls[nn].FeH)
;          lls[nn].flgFe = 25
;          if lls[nn].sigFeH LT 0.3 then return
;      endif
;      ;; Fe limits
      lls[nn].systems[sys].flgFe = nelc[26]
      if lls[nn].systems[sys].flgFe NE 0 then begin
          ans = lls[nn].systems[sys].XH[26].clm
          sig = lls[nn].systems[sys].XH[26].sigclm
          goto, jump1
      endif
      return
  endif else begin
      stop
      ;; By choice
      ;; Fe value
      if lls[nn].flgFe GT (-4) then begin
          lls_calcabnd, lls, nn, 26, 1, ans, sig, NOSIGY=NOHSIG
          goto, jump1
      endif else begin
      ;; Ni value
          case lls[nn].flgFe of
              -4: begin 
                  lls_calcabnd, lls, nn, 28, 1, ans, sig, NOSIGY=NOHSIG
                  ans = ans - 0.1
              end
              -5: begin
                  lls_calcabnd, lls, nn, 24, 1, ans, sig, NOSIGY=NOHSIG
                  ;; offset
                  ans = ans - 0.2 
              end
              -6: lls_calcabnd, lls, nn, 13, 1, ans, sig, NOSIGY=NOHSIG ;; Al value
              else: stop
          endcase
      endelse
  endelse	

  jump1:
  lls[nn].systems[sys].FeH = ans 
  lls[nn].systems[sys].sig_FeH = sig

  return
end
