;+ 
; NAME:
; x_fndpeaks   
;    Version 1.0
;
; PURPOSE:
;    Finds a set of peaks in a spectrum
;
; CALLING SEQUENCE:
;   
;   x_fndpeaks, spec, center, NSIG=, PEAK=
;
; INPUTS:
;   spec       - Input spectrum (1D)
;
; RETURNS:
;
; OUTPUTS:
;   center     - Array of peak centers
;
; OPTIONAL KEYWORDS:
;  PKWDTH      - Width of peak to center on (pixels)
;
; OPTIONAL OUTPUTS:
;   peak       - Integer values of centeroids
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndpeaks, spec, center
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fndpeaks, spec, center, NSIG=nsig, PEAK=peak, ALL=all, SILENT=silent, $
                PKWDTH=pkwdth, THIN=thin, NORECENT=norecent, FORCE=force, $
                FRACPK=fracpk, EDGES=edges, NORDB=nordb, NEDG=nedg, AFUNC=afunc,$
                MSK=msk, AUTOFIT=autofit


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_fndpeaks, spec, center, NSIG=, /THIN, /NORECENT [v1.0]'
    return
  endif 

  npix = n_elements(spec)

; Optional Keywords

  if not keyword_set( NORDB ) then nordb = 31L
  if not keyword_set( NSIG ) then nsig=5.
  if not keyword_set( PKWDTH ) then pkwdth=3L
  if not keyword_set( NEDG ) then nedg=5L
  if not keyword_set( AFUNC ) then afunc = 'BSPLIN'
  if not keyword_set( MSK ) then msk = bytarr(npix)+1B

  xdat = findgen(npix)
  ; Fit the spectrum with a crude BSPLIN
  if not keyword_set( AUTOFIT ) then begin
      autofit = x_fitrej(xdat, spec, afunc, nordb, MSK=msk, $
                         hsigma=2., lsigma=4., rms=rms)
  endif else begin
      ;; Calculate rms
      djs_iterstat, spec-autofit, sigma=rms, sigrej=2.0
  endelse

  ; Find all nsig sigma peaks and avoid edges
  gdpix = where( spec GT autofit+nsig*rms AND $
                 xdat GE nedg AND MSK EQ 1B AND $
                 xdat LE npix-nedg-1, ngpix)

  ; Error catch
  if ngpix LE 0 then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndpeaks: No', nsig, ' sig features!'
      center = -1
      return
  endif

  ; Include only the inflection points (min of 5 pts)
  npk = 0L
  peak = lonarr(npix)
  if not keyword_set( THIN ) then begin
      for i=0L,ngpix-1 do begin
          if ( spec[gdpix[i]-1] GE spec[gdpix[i]-2] AND $
               spec[gdpix[i]] GE spec[gdpix[i]-1] AND $
               spec[gdpix[i]] GE spec[gdpix[i]+1] AND $
               spec[gdpix[i]+1] GE spec[gdpix[i]+2] ) then begin
              peak[npk] = gdpix[i]
              npk = npk + 1
          endif
      endfor
  endif else begin
      for i=0L,ngpix-1 do begin
          if ( spec[gdpix[i]] GE spec[gdpix[i]-1] AND $
               spec[gdpix[i]] GE spec[gdpix[i]+1] ) then begin
              peak[npk] = gdpix[i]
              npk = npk + 1
          endif
      endfor
      ;; Contract all peaks within 2 of one another
      pkmsk = bytarr(npk) 
      spk = peak[sort(peak[0:npk-1])]
      shft2 = shift(spk,-1)
      max = abs(shft2 - spk)
      for i=0L,npk-2 do begin
          if max[i] LE 4 then begin
              spk[i] = (spk[i]+spk[i+1])/2.
              pkmsk[i+1] = 1B
          endif
      endfor
      
      ;; Reset
      gd = where(pkmsk EQ 0B, npk)
      peak = spk[gd]
  endelse
  ; Error catch
  if npk LE 0 then begin
      if not keyword_set(silent) then print, 'x_fndpeaks: No peaks found!'
      center = -1
      return
  endif

  peak = temporary(peak[0:npk-1])

  dume = fltarr(2)
  ; EDGES
  if arg_present(EDGES) then edges = fltarr(npk,2)

  if keyword_set( NORECENT ) then begin
      center = double(peak)
      return
  endif else begin
      ; Center up the Peaks 
      center = dblarr(npk)
      for i=0L,npk-1 do begin
          dume[*] = 0.
          ; Worry about limits
          pmn = (peak[i]-pkwdth) > 0
          pmx = (peak[i]+pkwdth) < (npix-1)
          ; Center up on the peak
          center[i] = x_centspln(xdat[pmn:pmx], $
                                 spec[pmn:pmx]-autofit[pmn:pmx], $
                                 fracpk, SILENT=silent, FORCE=force, $
                                 EDGES=dume)
          ; Edges
          if arg_present(edges) then edges[i,*] = dume
          ; Catch bad centroids 
          if abs(center[i]-xdat[peak[i]]) GT 1. then begin
              center[i] = xdat[peak[i]]
              if arg_present(edges) then edges[i,*] = -1.
          endif
      endfor
  endelse

  ; Just pass back the good ones
  if not keyword_set( ALL ) then begin
      gdpk = where(center NE -1)
      center = center[gdpk]
      if keyword_set(edges) then edges = temporary(edges[gdpk,*])
      peak = peak[gdpk]
  endif

  return

end
