;+ 
; NAME:
; dla_fndtran   
;   Version 1.1
;
; PURPOSE:
;    Given a DLA structure, calculates the rest EW of weak, rare
;    transitions like OI 1355, BII, etc.
;
; CALLING SEQUENCE:
;   dla_fndtran, dla, [fil], OUTFIL=, LMT=
;
; INPUTS:
;   dla --  IDL DLA structure
;   [fil] -- List of weak line transitions [default:
;            '/u/xavier/DLA/Abund/weak_lin.dat']
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL= -- Writes the values to OUTFIL [default: 'fort.23']
;
; OPTIONAL KEYWORDS:
;   LMT= -- Minimum EW to print the transition [default: 0.5 mA]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_fndtran, dla, outfil='weak_lin.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_fndtran, dla, fil, OUTFIL=outfil, LMT=lmt

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dla_fndtran, dla, [fil], OUTFIL=, LMT= [v1.1]'
    return
  endif 

; Optional Keywords
  if not keyword_set( fil ) then fil = '~/DLA/Abund/weak_lin.dat'
  if not keyword_set( outfil ) then outfil = 'fort.23'
  if not keyword_set( lmt ) then lmt = 0.5

; FILE
  close, 23
  openw, 23, outfil
  
; DLA
  ndla = n_elements(dla) 

; Grab the transitions
  readcol, fil, wv
  nwv = n_elements(wv)

  ;; Grab abnd
  all_abnd = fltarr(nwv)
  all_f = dblarr(nwv)
  all_nm = strarr(nwv)
  for i=0L, nwv-1 do begin
      ;; Grab name
      getfnam, wv[i], f, nm
      all_f[i] = f
      all_nm[i] = strtrim(nm,2)
      
      ;; Get abnd
      getabnd, nm, Znum, abnd
      all_abnd[i] = abnd
  endfor

; Create structure
  tmp = { $
          NHI: 0., $
          Z: 0., $
          wv: dblarr(nwv), $
          nm: strarr(nwv), $
          ew: dblarr(nwv) }

; Replicate
  weakstr = replicate(tmp, ndla)
  
; LOOP on DLA
  for kk=0L,ndla-1 do begin

      ;; Get metallicity
      if dla[kk].flgAlpha EQ 1 or dla[kk].flgAlpha GT 3 $
        then Z = dla[kk].alpha $
      else begin
          if dla[kk].flgZn EQ 1 then Z = dla[kk].ZnH $
          else begin
              if dla[kk].flgFe EQ 1 or dla[kk].flgFe GT 3 $
                then Z = dla[kk].FeH + 0.3$
              else continue
          endelse
      endelse

      weakstr[kk].NHI = dla[kk].NHI
      weakstr[kk].Z = Z
      weakstr[kk].nm[0:nwv-1] = all_nm
              
      ;; Calculate M column
      M = dla[kk].NHI + Z

      ;; Calculate expected N
      N = M - 12 + all_abnd

      ;; Save
      weakstr[kk].wv[0:nwv-1] = wv

      ;; Calculate rest EW
      EW = 8.85249d-13 * 10^N * all_f * wv * 1.d-8 * wv * 1.e3 ;; mA

      ;; Save
      weakstr[kk].ew[0:nwv-1] = EW

      ;; Print valid
      gd = where(weakstr[kk].ew[*] GT lmt, ngd)
      if ngd NE 0 then begin
          print, '-----------------------------------------------'
          printf, 23, '-----------------------------------------------'
          print, 'DLA :: ', dla[kk].qso, dla[kk].zabs, dla[kk].qso_mag, $
            weakstr[kk].Z, dla[kk].NHI
          printf, 23, 'DLA :: ', dla[kk].qso, dla[kk].zabs, $
            dla[kk].qso_mag, $
            weakstr[kk].Z, dla[kk].NHI
          print, 'Forest ends at ', (dla[kk].qso_zem+1.)*1215.6701
          printf, 23, 'Forest ends at ', (dla[kk].qso_zem+1.)*1215.6701
          printcol, weakstr[kk].nm[gd], weakstr[kk].wv[gd], $
            weakstr[kk].wv[gd]*(1.+dla[kk].zabs), $
            weakstr[kk].ew[gd]
          writecol, 'tmp', weakstr[kk].nm[gd], weakstr[kk].wv[gd], $
            weakstr[kk].wv[gd]*(1.+dla[kk].zabs), $
            weakstr[kk].ew[gd], FILNUM=23
      endif
  endfor

  close, /all
end
