;+ 
; NAME:
; x_photcal   
;   Version 1.1
;
; PURPOSE:
;    Derives a photometric solution for standard stars
;
; CALLING SEQUENCE:
;   
; x_photcal, obs, landolt, ans, sig, MINSIG=, SETAM=, NCORR=, CHISQ=
;
; INPUTS:
;   obs - Standard star observations  (stdstruct)
;   landolt - Landolt info (lndltstr)
;
; RETURNS:
;   ans - Fit
;   sig - Error on the fit parameters
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  MINSIG - Minimum error for obs+landolt        (default = 0.01)
;  MAXSIG - Maximum error to include in analysis (default = 0.1)
;  SETAM   - Value to assume for the airmass term
;  /NOCLR  - No color term (e.g. only one filter solution)
;  MIN_NOBS  - Minimum number of n epochs to include (default = 4)
;  MIN_MOBS  - Minimum number of m epochs to include (default = 2)
;
; OPTIONAL OUTPUTS:
;  NCORR - Normalized correlation matrix
;  CHISQ - Chisq of the calibrated fit
;
; COMMENTS:
;
; EXAMPLES:
;   x_photcal, obs, landolt, ans, sig
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_photcal, obs, landolt, ans, sig, MINSIG=minsig, MAXSIG=maxsig, $
               SETAM=setam, NCORR=NCORR, NOCLR=noclr, CHISQ=chisq, $
               MIN_NOBS=min_nobs, MIN_MOBS=min_mobs

;
  if  N_params() LT 4  then begin 
      print, 'Syntax - ' +$
        'x_photcal, obs, landolt, ans, sig, MINSIG=, MAXSIG=, SETAM=, NCORR='
      print, '        CHISQ=, /NOCLR (v1.1)'
      return
  endif 

; Optional Keywords

  if not keyword_set( MINSIG ) then minsig = 0.03
  if not keyword_set( MAXSIG ) then maxsig = 0.1
  if n_elements( MIN_NOBS ) EQ 0 then min_nobs = 4
  if n_elements( MIN_MOBS ) EQ 0 then min_mobs = 2

; Begin the bookkeeping

  nobs = n_elements(obs)

  nldlt = n_elements(landolt)


; Big loop

  nstrs = 0
  mobs = fltarr(nobs)
  mtru = fltarr(nobs)
  error = fltarr(nobs)
  airmass = fltarr(nobs)
  color = fltarr(nobs)

  for i=0,nobs-1 do begin
      ; Analysis check
      if obs[i].flg_anly EQ 0 then continue

      ; Restrict on sigma
      if obs[i].sig_Mag GT maxsig then begin
          obs[i].flg_anly = 0
          continue
      endif

      ; Search Landolt catalog
      istr = where(obs[i].Name EQ landolt.Name, cnt)
      if cnt EQ 0 then begin
          print, 'Star not found in Landolt list!: ', obs[i].Name
          obs[i].flg_anly = 0
          continue
      endif
      if cnt GT 1 then begin
          print, 'Found multiple entries in Landolt list!: ', obs[i].Name
          obs[i].flg_anly = 0
          continue
      endif

      ; Restrict on Number of Landolt obs
      if landolt[istr].nobs LT min_nobs OR landolt[istr].mobs LT min_mobs then begin
          obs[i].flg_anly = 0
          continue
      endif

  
      ; Magnitudes
      mobs[nstrs] = obs[i].Mag
      mtru[nstrs] = x_lndltMag(obs[i].filter, landolt[istr], SIGMag=sig_lndlt)

      ; Error
      sig = sqrt(obs[i].sig_Mag^2 + sig_lndlt^2)
      if sig GT maxsig then begin
          obs[i].flg_anly = 0
         continue
      endif
      error[nstrs] = sig > minsig

      ; Airmass
      airmass[nstrs] = obs[i].AM

      ; Color
      if not keyword_set( NOCLR ) then $
        color[nstrs] = x_lndltclr(obs[i].sCLR, landolt[istr])

      ; Index
      nstrs = nstrs + 1
  endfor

  ; Run x_photsol

  if not keyword_set( SETAM ) then begin
      if not keyword_set( NOCLR ) then begin
          ; All 3 parameters are included
          if keyword_set( NCORR ) then $
            x_photsol3, mobs[0:nstrs-1], mtru[0:nstrs-1], error[0:nstrs-1], $
            airmass[0:nstrs-1], color[0:nstrs-1], ans, sig, $
            NCORR=ncorr, CHISQ=tchisq  $
          else $  ; No correlation matrix
            x_photsol3, mobs[0:nstrs-1], mtru[0:nstrs-1], error[0:nstrs-1], $
            airmass[0:nstrs-1], color[0:nstrs-1], ans, sig, CHISQ=tchisq
      endif else $       ; AM, Zero, but NO COLOR
        x_photsol2, mobs[0:nstrs-1], mtru[0:nstrs-1], error[0:nstrs-1], $
        airmass[0:nstrs-1], ans, sig, CHISQ=tchisq
  endif else begin ; AM term is set
      if not keyword_set( NOCLR ) then $
        x_photsol2, mobs[0:nstrs-1], mtru[0:nstrs-1], error[0:nstrs-1], $
        airmass[0:nstrs-1], ans, sig, SETAM=setam, COLOR=color[0:nstrs-1], $
        CHISQ=tchisq $
      else $
        x_photsol1, mobs[0:nstrs-1], mtru[0:nstrs-1], $
        error[0:nstrs-1], airmass[0:nstrs-1], setam, ans, sig, CHISQ=tchisq
  endelse

  if arg_present(CHISQ) then chisq = tchisq
        
      
; Clean up
  delvarx, mobs, mtru, airmass, color, error

end

