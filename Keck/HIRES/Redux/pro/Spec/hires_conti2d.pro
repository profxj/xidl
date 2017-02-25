;+ 
; NAME:
; hires_conti2d
;     Version 1.1
;
; PURPOSE:
;   Enables fitting of the continuum, one echelle order at a time
;   using x_continuum
;
; CALLING SEQUENCE:
;   hires_conti2d, file, outfil, NEWFIL=, INFIL=, ORDRS=, _EXTRA=extra
;
; INPUTS:
;    file - Filename to continuum fit
;    outfil - Filename of continuum
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ORDRS=  --  Specific orders to fit
;  INFIL=  -- Filename of previously existing continuum fit 
;  ZRO_FIL= -- Filename of zero-level FITS file.  If input, this is
;              subtracted off prior to the continuum fitting
;
; OPTIONAL OUTPUTS:
;  NEWFIL= -- Name of file to hold normalized spectrum
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_continuum
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_conti2d, file, outfil, NEWFIL=newfil, INFIL=infil, ZRO_FIL=zro_fil, $
                   ORDRS=ordrs, _EXTRA=extra, SHOW_PREV=show_prev

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_conti2d, file, outfil, NEWFIL=, INFIL=, ORDRS=, ZRO_FIL=, /SHOW_PREV [v1.1]'
      return
  endif 

  ;; Optional
  if not keyword_set(DUMMF) then dummf = 'dumf.fits'

  ;; Open the file
  spec2d = xmrdfits(file,1,/silent)
  sumv = total(spec2d.var, 1)
  if not keyword_set(ORDRS) then begin
      gdo = where(spec2d.phys_ordr GT 0 AND $
                  sumv GT 0., nord)
  endif else begin
      nord = n_elements(ORDRS)
      gdo = lonarr(nord)
      for jj=0L,nord-1 do begin
          gdo[jj] = where(spec2d.phys_ordr EQ ordrs[jj])
      endfor
  endelse
          
  sz = size(spec2d.fx,/dimensions)

  ;; Zero level
  if keyword_set(ZRO_FIL) then begin
     print, 'hires_conti2d: Correcting zero level with ', zro_fil
     zro_lvl = xmrdfits(ZRO_FIL)
  endif

  ;; Create conti array
  npx = sz[0]
  if not keyword_set(INFIL) then begin
      conti = fltarr(npx, sz[1])
      conti[*] = 1.
  endif else conti = xmrdfits(infil,/silent)

  ;; Loop
  for jj=0L,nord-1 do begin
      qq = gdo[jj]
      ;; Grab good data
      gd = where(spec2d.wave[*, qq] GT 0. AND $
                 spec2d.var[*, qq] GT 0., ngd)
      ;; Start/end pix
      i0 = gd[0]
      i1 = gd[ngd-1]
      np = i1-i0+1
      gdp = i0 + lindgen(np)

      ;; zero out large pixels
      badpix = where(abs(spec2d.fx[*, qq]) GT 1.0d6 OR $
                     finite(spec2d.fx[*, qq]) NE 1, nbad)
      IF nbad GT 0 THEN spec2d.fx[badpix, qq] = 0.0
      
      ;; Zero lvl
      if keyword_set(ZRO_LVL) then zlvl = zro_lvl[gdp,qq] $
      else zlvl = replicate(0., np)
         

      ;; Create dummy file
      mwrfits, spec2d.fx[gdp,qq]-zlvl, dummf, /create
      mwrfits, replicate(1.,np), dummf
      mwrfits, spec2d.wave[gdp,qq], dummf

      ;; INFIL
      if keyword_set(INFIL) then begin
          cfil = 'dumc.fits'
          mwrfits, conti[gdp,qq], cfil, /create
      endif

      ;; Run conti
      if keyword_set(SHOW_PREV) and jj GT 0 then show_conti=conti[gdp,qq-1]
      x_continuum, dummf, inflg = 2, outfil = 'dumc.fits', CONTI = cfil, $
                   _EXTRA=extra, SHOW_CONTI=show_conti, FORDR=7L
      if keyword_set(CFIL) then delvarx, cfil

      ;; Read and save
      sub_cont = xmrdfits('dumc.fits',/silent)
      conti[gdp,qq] = sub_cont

      ;; Save
      mwrfits, conti, outfil, /create
  endfor

  mwrfits, conti, outfil, /create

  ;; Normalize as desired
  if keyword_set(NEWFIL) then begin
      spec2d.fx = spec2d.fx / conti
      spec2d.var = spec2d.var / conti / conti
      IF TAG_EXIST(spec2d, 'NOVAR') THEN $
        spec2d.novar = spec2d.novar/conti/conti
      mwrfits, spec2d, newfil, /create
  endif

  return 

end
      
