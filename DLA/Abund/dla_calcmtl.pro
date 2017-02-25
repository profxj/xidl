;+ 
; NAME:
; dla_calcmtl
;  V1.2
;
; PURPOSE:
;    Determine the metallicity [M/H] by parsing through the
;     elemental abundances in a DLA strucure.  The user can control
;     which element is used by setting the tag flgmtl in advance
;     to calling this routine.
;
; CALLING SEQUENCE:
;   dla_calcmtl, dla, nn
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
;  Values of flgalpha --   
;                       1,4,5,6 -- Si, S, O, Ar
;                       2 -- Zn
;                       13 -- Bracketed by two values.  Takes mean.
;                       14,15 -- Fe value (adds +0.3dex) or limit
;                       -1 -- Use Alpha
;                       -2 -- Use Zn
;                       -3 -- Combine Alpha and Zn limits
;                       -4 -- Use Fe
;                       Subtract 10 for lower limit
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   06-Oct-2004 Written by JXP
;   07-Jul-2011 Changed mtl = Feh + 0.3 (instead of 0.4) MAR
;   03-Aug-2011 Changed sigmtl to add 0.16 in quadrature when using Feh MAR
;- 
;------------------------------------------------------------------------------
pro dla_calcmtl, dla, nn

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dla_calcmtl, dla, nn [v1.1]'
    return
  endif 

  llmt = -999.
  ulmt = 999.
  lmts = fltarr(5,3)
  lmts[*,1] = -999.
  lmts[*,2] = -999.
  flg = intarr(5)
  totflg = intarr(3)

  ;; AUTOMATIC
  if dla[nn].flgmtl GE 0 then begin
      ;; Alpha
     if dla[nn].flgAlpha EQ 1 OR dla[nn].flgAlpha EQ 4 OR dla[nn].flgAlpha EQ 5 $
        OR abs(dla[nn].flgAlpha-5) LE 1 then begin
          dla[nn].mtl = dla[nn].alpha
          dla[nn].sigmtl = sqrt(dla[nn].sigalpha^2 + dla[nn].signhi[0]^2)
          dla[nn].flgmtl = dla[nn].flgalpha
          return
      endif
      ;; Zn
      if dla[nn].flgZn EQ 1 then begin
          dla[nn].mtl = dla[nn].ZnH
          dla[nn].sigmtl = sqrt(dla[nn].sigZnH^2 + dla[nn].signhi[0]^2)
          dla[nn].flgmtl = 2
          return
      endif
      ;; Si, S, O, Zn limits
      dla_elmlmts, dla, nn,14,2, d1, d2, FLG=f1
      lmts[1,1] = d1
      lmts[1,2] = d2
      flg[1] = f1
      dla_elmlmts, dla, nn,16,2, d1, d2, FLG=f1
      lmts[2,1] = d1
      lmts[2,2] = d2
      flg[2] = f1
      dla_elmlmts, dla, nn,8,2, d1, d2, FLG=f1
      lmts[3,1] = d1
      lmts[3,2] = d2
      flg[3] = f1
      dla_elmlmts, dla, nn,30,2, d1, d2, FLG=f1
      lmts[4,1] = d1
      lmts[4,2] = d2
      flg[4] = f1
      for i=1,4 do begin
         if (flg[i] MOD 2) EQ 1 then totflg[1] = 1
         if (flg[i] MOD 4) GT 1 then totflg[2] = 1
         llmt = llmt > lmts[i,1]
         ulmt = ulmt < lmts[i,2]
      endfor
      ;;
      if total(totflg) EQ 2 then begin
          dla[nn].mtl = (llmt+ulmt)/2.
          dla[nn].sigmtl = 0.15 > (ulmt-dla[nn].mtl)
          dla[nn].flgmtl = 13
          if dla[nn].sigmtl LT 0.3 then return
      endif
      ;; Fe
      if dla[nn].flgFe EQ 1 OR abs(dla[nn].flgFe-5) LE 1 then begin
          dla[nn].mtl = dla[nn].FeH + 0.3
;          dla[nn].sigmtl = dla[nn].sigFeH + 0.1
          dla[nn].sigmtl = sqrt(dla[nn].sigFeH^2 + 0.16^2 + dla[nn].signhi[0]^2)
          dla[nn].flgmtl = 14
          return
      endif
      if (dla[nn].flgFe EQ 11 OR dla[nn].flgFe EQ 13) AND $
        dla[nn].sigFeH LT 0.25 then begin
          dla[nn].mtl = dla[nn].FeH + 0.3
;          dla[nn].sigmtl = dla[nn].sigFeH + 0.1
          dla[nn].sigmtl = sqrt(dla[nn].sigFeH^2 + 0.16^2 + dla[nn].signhi[0]^2)
          dla[nn].flgmtl = 15
          return
      endif
      dla[nn].mtl = 0.
      dla[nn].sigmtl = 0.
      dla[nn].flgmtl = 0
  endif else begin
      case (dla[nn].flgmtl mod 10) of 
          -1: begin ; Alpha
              dla[nn].mtl = dla[nn].alpha
              dla[nn].sigmtl = dla[nn].sigalpha
          end
          -2: begin ; Zn
              dla[nn].mtl = dla[nn].ZnH
              dla[nn].sigmtl = dla[nn].sigZnH
          end
          -3: begin ; Combination of alpha/Zn limits
              dla[nn].mtl = (dla[nn].ZnH + dla[nn].alpha) / 2.
              dla[nn].sigmtl = dla[nn].sigFeH + 0.1
          end
          -4: begin ; Fe
              dla[nn].mtl = dla[nn].FeH + 0.3
              ;dla[nn].sigmtl = dla[nn].sigFeH + 0.1
              dla[nn].sigmtl = sqrt(dla[nn].sigFeH^2 + 0.16^2)
          end
          -5: begin ; Fe
              stop  ;; Juse use -4 if you are trying to do Fe
          end
          else: stop
      endcase
  endelse

  return
end
