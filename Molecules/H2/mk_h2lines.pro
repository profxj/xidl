;+ 
; NAME:
; mk_h2lines
;  (V1.0)
;
; PURPOSE:
;    Generate an ASCII file of H2 lines for the Lyman and Werner
;    bandheads
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
;   Oct-2008 Written by JXP with guidance from Y Sheffer
;-
;------------------------------------------------------------------------------
pro mk_h2lines 

  if not keyword_set(OUTFIL) then outfil = 'h2_levels.dat'
  if not keyword_set(FITSFIL) then fitsfil = 'h2_lines.fits'
 ;;
 c = x_constants()

 ;; Original line list (for gamma values)
 h2orig = read_h2lin(/ORIG)
 norig = n_elements(h2orig)

 ;; Lyman first
 readcol, 'Lyman.AAS93.dat', vu, ju, vl, jl, A, E, format='L,L,L,L,D,D'
 gdX = where(vl EQ 0, nlyman)
; vu = vu[gd]
; ju = ju[gd]
; vl = vl[gd]
; jl = jl[gd]
; A = A[gd]
; E=E[gd]  ; cm^-1

 dJ = ju-jl
 rlvl = where(dJ EQ 1, complement=plvl)
 ;; Flambda
 flam = 1.49919d4 * (2*ju+1) * (0.01/E)^2 * A / (2*jl+1)

 ;; Label
 lbl = strarr(n_elements(A))
 branch = strarr(n_elements(A))
 branch[rlvl] = 'R'
 branch[plvl] = 'P'

 for ii=0L,nlyman-1 do begin
     lbl[gdx[ii]] = 'B'+strtrim(vu[gdX[ii]],2)+'-0'+branch[gdx[ii]]+'('+$
               strtrim(jl[gdX[ii]],2)+')'
 endfor

 wave = 1./E * 1d8 ;; Ang

 ;; Gamma
 gamma = replicate(1e9,nlyman)
 for ii=0L,max(vu) do begin
     gdj = where(vu EQ ii)
     for jj=0L,max(ju[gdj]) do begin
         gdvj = where(ju[gdj] EQ jj, ngdvj)
         if ngdvj EQ 0 then continue
         gamma[gdj[gdvj]] = total(A[gdj[gdvj]]) > 1e9
     endfor
 endfor
 

 tmp = {h2linstrct}
 lymanlin = replicate(tmp, nlyman)

 lymanlin.wrest = wave[gdX]
 lymanlin.f = flam[gdX]
 lymanlin.gamma = gamma[gdX]
 lymanlin.el = 2
 lymanlin.np = vu[gdX]
 lymanlin.npp = vl[gdX]
 lymanlin.Jp = Ju[gdX]
 lymanlin.Jpp = Jl[gdX]
 lymanlin.label = lbl[gdX]

 h2lin = lymanlin

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;; Werner 
 readcol, 'Werner_m.AAS93.dat',  mvu, mju, mvl, mjl, mA, mE, format='L,L,L,L,D,D'
 readcol, 'Werner_p.AAS93.dat',  pvu, pju, pvl, pjl, pA, pE, format='L,L,L,L,D,D'
 vu = [mvu,pvu]
 vl = [mvl,pvl]
 ju = [mju,pju]
 jl = [mjl,pjl]
 A  = [mA,pA]
 E  = [mE,pE]
 gdX = where(vl EQ 0, nwerner)

 dJ = ju-jl

 ;; Flambda
 flam = 1.49919d4 * (2*ju+1) * (0.01/E)^2 * A / (2*jl+1)
 
 ;; Label
 lbl = strarr(n_elements(A))
 branch = strarr(n_elements(A))
 qlvl = where(dJ EQ 0, complement=olvl)
 branch[qlvl] = 'Q'
 rlvl = where(dJ[olvl] EQ 1, complement=plvl)
 branch[olvl[rlvl]] = 'R'
 branch[olvl[plvl]] = 'P'
 
 for ii=0L,nwerner-1 do begin
     lbl[gdx[ii]] = 'C'+strtrim(vu[gdX[ii]],2)+'-0'+branch[gdx[ii]]+'('+$
               strtrim(jl[gdX[ii]],2)+')'
 endfor
 
 wave = 1./E * 1d8 ;; Ang
 
 ;; Gamma
 gamma = replicate(1e9,nwerner)
 for ii=0L,max(vu) do begin
     gdj = where(vu EQ ii)
     for jj=0L,max(ju[gdj]) do begin
         gdvj = where(ju[gdj] EQ jj, ngdvj)
         if ngdvj EQ 0 then continue
         gamma[gdj[gdvj]] = total(A[gdj[gdvj]]) > 1e9
     endfor
 endfor
 
 tmp = {h2linstrct}
 wernerlin = replicate(tmp, nwerner)
 
 wernerlin.wrest = wave[gdX]
 wernerlin.f = flam[gdX]
 wernerlin.gamma = gamma[gdX]
 wernerlin.el = 3
 wernerlin.np = vu[gdX]
 wernerlin.npp = vl[gdX]
 wernerlin.Jp = Ju[gdX]
 wernerlin.Jpp = Jl[gdX]
 wernerlin.label = lbl[gdX]

 ;; Save
 h2lin = [h2lin, wernerlin]
 srt = sort(h2lin.wrest)
 h2lin = h2lin[srt]

 ;; Toss out lines below the Lyman limit
 gd = where(h2lin.wrest GT 911.)

 ;; Use the gamma values from the original line list
 for jj=0L,norig-1 do begin
     ;; Match
     mt = where(strmatch(strtrim(h2lin.label,2),strtrim(h2orig[jj].label,2)),nmt)
     if nmt NE 1 then stop
     ;; Set gamma
     h2lin[mt].gamma = h2orig[jj].gamma
 endfor

 ;; Out
 mwrfits, h2lin[gd], fitsfil, /create
 spawn, 'gzip -f '+fitsfil

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
