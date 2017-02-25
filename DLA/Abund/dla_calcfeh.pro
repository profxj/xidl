;+ 
; NAME:
; dla_calcfeh
;  V1.2
;
; PURPOSE:
;    Determine the Fe-peak abundance [Fe/H] by parsing through the
;     elemental abundances in a DLA structure.  The user can control
;     which element is used by setting the tag flgFe in advance
;     to calling this routine.
;
; CALLING SEQUENCE:
;   dla_calcfeh, dla, nn
;
; INPUTS:
;   dla -- DLA structure
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
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_calcfeh, dla, nn, NOHSIG=nohsig, FINE=fine

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dla_calcfeh, dla, nn, /NOHSIG, /FINE [v1.1]'
    return
  endif 

  llmt = -999.
  ulmt = 999.
  lmts = fltarr(5,3)
  lmts[*,1] = -999.
  lmts[*,2] = -999.
  flg = intarr(5)

  nelc = dla[nn].elm.flgclm
  ;; AUTOMATIC
  if dla[nn].flgFe GE 0 then begin
      ;; Fe value
      if nelc[26] EQ 1 then begin
          if keyword_set(FINE) then begin
              feidx = [41 + lindgen(4)]
              gdfe = where(dla[nn].ion[26].state[feidx].flgclm EQ 1, ngdf)
              if ngdf GT 1 then begin
                  dla[nn].elm[26].clm = dla[nn].elm[26].clm + $
                    total(dla[nn].ion[26].state[feidx[gdfe]].clm)
              endif 
          endif
          dla_calcabnd, dla, nn, 26, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgFe = 1
          goto, jump1
      endif 
      ;; Ni value
      if nelc[28] EQ 1 then begin
          dla_calcabnd, dla, nn, 28, 1, ans, sig, NOSIGY=NOHSIG
          ans = ans - 0.1  ;; offset, NOSIGY=NOHSIG
          dla[nn].flgFe = 4
          goto, jump1
      endif 
      ;; Fe limits
      dla_elmlmts, dla, nn, 26, 2, d1, d2, FLG=f1
      lmts[1,1] = d1
      lmts[1,2] = d2
      flg[1] = f1
      if flg[1] EQ 3 then begin
          dla[nn].FeH = (lmts[1,1] + lmts[1,2])/2.
          dla[nn].sigFeH = 0.1 >  (lmts[1,2]-dla[nn].FeH)
          dla[nn].flgFe = 11
          if dla[nn].sigFeH LT 0.3 then return
      endif
      ;; Fe+Ni
      dla_elmlmts, dla, nn, 28, 2, d1, d2, FLG=f1 ;; Ni
      lmts[2,1] = d1
      lmts[2,2] = d2
      flg[2] = f1
      ;; check limits
      if(flg[1]+flg[2] GT 2) AND ( (flg[1]+flg[2] MOD 2) EQ 1) then begin 
          llmt = lmts[1,1] > (lmts[2,1]-0.1)
          ulmt = lmts[1,2] < (lmts[2,2]-0.1)
          dla[nn].FeH = (llmt+ulmt)/2.
          dla[nn].sigFeH = 0.15 >  (ulmt-dla[nn].FeH)
          dla[nn].flgFe = 13
          if dla[nn].sigFeH LT 0.3 then return
      endif
      ;; Cr 
      if nelc[24] EQ 1 then begin
          dla_calcabnd, dla, nn, 24, 1, ans, sig, NOSIGY=NOHSIG
          ans = ans - 0.2  ; offset
          dla[nn].flgFe = 5
          goto, jump1
      endif
      ;; Al
      if nelc[13] EQ 1 then begin
          dla_calcabnd, dla, nn, 13, 1, ans, sig, NOSIGY=NOHSIG
          dla[nn].flgFe = 6
          goto, jump1
      endif
      ;; All Limits
      dla_elmlmts,dla, nn,24,2, d1, d2, FLG=f1 ; Cr
      lmts[3,1] = d1
      lmts[3,2] = d2
      flg[3] = f1
      dla_elmlmts,dla, nn,13,2, d1, d2, FLG=f1 ; Al
      lmts[4,1] = d1
      lmts[4,2] = d2
      flg[4] = f1
      if total(flg MOD 2) GE 1 AND total(flg GT 1) GT 1 then begin
;      if((mod(flg(1),2).eq.1.or.mod(flg(2),2).eq.1.or.
;     x        mod(flg(3),2).eq.1.or.mod(flg(4),2).eq.1).and.
;     x      (flg(1).gt.1.or.flg(2).gt.1.or.
;     x        flg(3).gt.1.or.flg(4).gt.1)) then
          llmt = lmts[1,1] > (lmts[2,1]-0.1)
          llmt = llmt > (lmts[3,1]-0.2)
          llmt = llmt > lmts[4,1]
          ulmt = lmts[1,2] < (lmts[2,2]-0.1)
          ulmt = ulmt < (lmts[3,2]-0.2)
          ulmt = ulmt < lmts[4,2]
          ;; values
          dla[nn].FeH = (llmt+ulmt)/2.
          dla[nn].sigFeH = 0.15 >  (ulmt-dla[nn].FeH)
          dla[nn].flgFe = 25
          if dla[nn].sigFeH LT 0.3 then return
      endif
      ;; Fe limits
      dla[nn].flgFe = nelc[26]
      if dla[nn].flgFe NE 0 then begin
          dla_calcabnd, dla, nn, 26, 1, ans, sig, NOSIGY=NOHSIG
          goto, jump1
      endif
      return
  endif else begin
      ;; By choice
      ;; Fe value
      if dla[nn].flgFe GT (-4) then begin
          dla_calcabnd, dla, nn, 26, 1, ans, sig, NOSIGY=NOHSIG
          goto, jump1
      endif else begin
      ;; Ni value
          case dla[nn].flgFe of
              -4: begin 
                  dla_calcabnd, dla, nn, 28, 1, ans, sig, NOSIGY=NOHSIG
                  ans = ans - 0.1
              end
              -5: begin
                  dla_calcabnd, dla, nn, 24, 1, ans, sig, NOSIGY=NOHSIG
                  ;; offset
                  ans = ans - 0.2 
              end
              -6: dla_calcabnd, dla, nn, 13, 1, ans, sig, NOSIGY=NOHSIG ;; Al value
              else: stop
          endcase
      endelse
  endelse	

  jump1:
  dla[nn].FeH = ans 
  dla[nn].sigFeH = sig

  return
end
