 ;+ 
; NAME:
; x_fndobj   
;    Version 1.1
;
; PURPOSE:
;    Finds a set of obj in a spectrum.  I think this is a precursor
;  to x_fndpeaks.  I doubt this routine is worth using.
;
; CALLING SEQUENCE:
;  x_fndobj, spec, center, NSIG=, /SILENT, /FORCE, $
;             PKWDTH=, EDGES=, NORDB=, NEDG=, AFUNC=, PEAK=, /NOSMOOTH
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
;   x_fndobj, spec, center
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fndobj, spec, center, NSIG=nsig, SILENT=silent, FORCE=force $
              , PKWDTH = pkwdth, EDGES = edges, NORDB = nordb, NEDG = nedg $
              , AFUNC = afunc, PEAK = peak, NOSMOOTH = nosmooth $
              , FRACPK = FRACPK


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_fndobj, spec, center, NSIG=, /NORECENT, /FORCE, PKWDTH=, EDGES='
     print, '   NORDB=, NEDG=, AFUNC=, PEAK=, /NOSMOOTH [v1.1]'
    return
  endif 

; Optional Keywords

  if not keyword_set( NORDB ) then nordb = 3L
  if not keyword_set( NSIG ) then nsig=3.
  if not keyword_set( PKWDTH ) then pkwdth=3L
  if not keyword_set( NEDG ) then nedg=5L
  if not keyword_set( AFUNC ) then afunc = 'POLY'
  if not keyword_set( SMTHWID ) then smthwid = 5L

; SMOOTH
;  if not keyword_set( NOSMOOTH ) then spec = median(inspec, smthwid) $
;  else spec = inspec
  
; Fit the spectrum with a low-order POLY
  npix = n_elements(spec)
  xdat = findgen(npix)
  autofit = x_fitrej(xdat, spec, afunc, nordb, $
                     hsigma=2., lsigma=4., rms=rms)

  ; Find all nsig sigma peaks and avoid edges
  gdpix = where( spec GT autofit+nsig*rms AND $
                 xdat GE nedg AND $
                 xdat LE npix-nedg-1, ngpix)


  ; Error catch
  if ngpix LE 0 then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndobj: No', nsig, ' sig features!'
      center = -1
      return
  endif

  msk = bytarr(npix)
  msk[gdpix] = 1B

; Find all peaks with groupings of gdpix 

  lft = where(msk EQ 1 AND shift(msk,1) EQ 0, nlft)
  rgt = where(msk EQ 1 AND shift(msk,-1) EQ 0, nrgt)
  if nlft NE nrgt then stop

  szpeak = rgt-lft

  ;; Good peaks
  gdpeak = where(szpeak GE PKWDTH, npk)
  if npk EQ 0 then begin
      center = -1
      return
  endif

; EDGES
  if arg_present(EDGES) then edges = fltarr(npk,2)
; Center?
  if keyword_set( NORECENT ) then begin
      center = double(npk)
      for i=0L,npk-1 do begin
          mx = max( spec[lft[i]:rgt[i]],imn )
          center[i] = double(lft[i]+imn)
      endfor
      if arg_present(peak) then peak = round(center)
      return
  endif else begin
      dume = dblarr(2)
      ;; Center up the Peaks 
      center = dblarr(npk)
      for i=0L,npk-1 do begin
          dume[*] = 0.
          ;; Center up on the peak
          pmn = lft[gdpeak[i]]
          pmx = rgt[gdpeak[i]]
          center[i] = x_centspln(xdat[pmn:pmx], $
                                 spec[pmn:pmx]-autofit[pmn:pmx], $
                                 fracpk, SILENT=silent, FORCE=force, $
                                 EDGES=dume)
          ;; Edges
          if arg_present(edges) then edges[i,*] = dume
          ; Catch bad centroids 
;          if abs(center[i]-xdat[peak[i]]) GT 1. then begin
;              center[i] = xdat[peak[i]]
;              if arg_present(edges) then edges[i,*] = -1.
;          endif
      endfor
      if arg_present(peak) then peak = round(center)
  endelse

  ; Just pass back the good ones
  if not keyword_set( ALL ) then begin
      gdpk = where(center NE -1, ngd)
      if ngd NE 0 then begin
          center = center[gdpk] 
          if keyword_set(edges) then edges = temporary(edges[gdpk,*])
          peak = peak[gdpk]
      endif else begin
          center = -1
          edges = [-1,-1]
          peak = -1
      endelse
  endif

  return

end
