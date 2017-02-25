;+ 
; NAME:
; mk_colines
;  (V1.0)
;
; PURPOSE:
;    Generate an ASCII file of CO lines for the A-X bandheads
;
; CALLING SEQUENCE:
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Sep-2008 Written by JXP with guidance from Y Sheffer
;-
;------------------------------------------------------------------------------
pro mk_colines

  if not keyword_set(OUTFIL) then outfil = 'co_levels.dat'
 ;;
 c = x_constants()

; Fudget factors:: Tuned to match Morton 1994 values
 fudge = fltarr(20)
 fudge[0] = 0.088
 fudge[1] = -0.0728
 fudge[2] = 0.009
 fudge[3] = -0.0085
 fudge[4] = -0.0208
 fudge[5] = 0.002
 fudge[6] = 0.0203

 ;; Loop on v' (v'' = 0)
 nlin = 0L
 svlbl = strarr(10000L)
 svlbl2 = strarr(10000L)
 svjp = lonarr(10000L)
 svjpp = lonarr(10000L)
 svvp = lonarr(10000L)
 svvpp = lonarr(10000L)
 svflam = strarr(10000L)
 svwv = dblarr(10000L)
 vpp = 0L

 for vp = 0L,8L do begin
     ;; Loop on J''
     for Jpp = 0L, 22L do begin
         
         ;; Loop on dJ
         for dJ = -1, 1L do begin
             Jp = Jpp + dJ
             if Jp LT 0 then continue
             if Jpp EQ 0 and Jp EQ 0 then continue  ;  no 0-0 transition
             if Jpp EQ 1 and Jp EQ 0 then continue  ;  no P(1)

             ;; Save
             svvp[nlin] = vp
             svjp[nlin] = jp 
             svjpp[nlin] = jpp

             ;; Label
             case dJ of
                 -1: lbl='P('
                 0: lbl='Q('
                 1: lbl='R('
                 else: stop
             endcase
             lbl = lbl + strtrim(round(Jpp),2)+')'

             ;; Energy of X level
             EX = co_xlevel(vpp,jpp)  ;; Relative to potential minimum;  cm^-1
             EA = co_alevel(vp,jp) ; Same
             dE = abs(EA-EX)
             wave = 1. / dE * 1e8 + fudge[vp] ; Ang
             svwv[nlin] = wave

             ;; flambda
             svflam[nlin] = alog10(co_flambda(vp, jp, jpp))

             ;; Save + Print
             svlbl[nlin] = ' A-X'+strtrim(vp,2)+'-0'+lbl
             svlbl2[nlin] = ' A'+strtrim(vp,2)+'-0'+lbl
             if strmatch(lbl,'R(0)') then print, wave, svlbl[nlin]

             ;; Increment
             nlin = nlin + 1
         endfor
     endfor
 endfor
 

 srt = sort(svwv[0:nlin-1])
 fval = svflam[srt] - alog10(svwv[srt])

 ;; Morton style
 writecol, outfil, svwv[srt], svlbl[srt], svflam[srt]
 ;; COsort.dat
 close, /all
 openw, 1, 'COsort.dat'
 writecol, 'dum', lindgen(nlin)+1, replicate('12  1',nlin), $
           svvp[srt], replicate('0',nlin), svjp[srt], svjpp[srt], $
           svwv[srt], fval, svflam, replicate(0,nlin), $
           replicate('8.000e+08',nlin), svlbl2[srt], $
           FMT='(i3,2x,a5,2x,i2,3x,a1,2x,i2,2x,i2,3x,f8.3,4x,'+$
               'f6.3,6x,f6.3,4x,i1,2x,a9,3x,a10)', FILNUM=1
 close, /all

 return
end

 ;; Rotation constants from Tilford & Simmons 1972
; Bcalc = [1.6002d, $
;          1.5787, $
;          1.5571, $
;          1.5348, $
;          1.5115, $
;          1.4876, $
;          1.4633, $
;          1.4389, $
;          1.4146, $
;          1.3903, $
;          1.3660, $
;          1.3413, $
;          1.3162, $
;          1.2905, $
;          1.2639, $
;          1.2364, $
;          1.2078, $
;          1.1779, $
;          1.1460, $
;          1.1108, $
;          1.0700, $
;          1.0195, $
;          0.9529, $
;          0.8608 $
;         ]  ; cm^-1

; v0 = [ 64748.53d, $
;         66229.86, $
;         67676.19, $
;         69088.21, $
;         70465.92, $
;         71809.13, $
;         73117.76, $
;         74391.85, $
;         75631.59, $
;         76837.20, $
;         78008.89, $
;         79146.70, $
;         80250.53, $
;         81320.00, $
;         82354.52, $
;         83353.21, $
;         84314.82, $
;         85237.63, $
;         86119.15, $
;         86955.60, $
;         87741.17, $
;         88466.84, $
;         89118.77, $
;         89676.07 $
;       ]  ; cm^-1

 ;; Energies:  dE = dE_el + dE_vib + B' J'(J'+1) - B'' J''(J''+1)
