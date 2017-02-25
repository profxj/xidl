;+ 
; NAME:
; hires_zerolvl2d
;     Version 1.1
;
; PURPOSE:
;   Enables fitting of the continuum, one echelle order at a time
;   using x_continuum
;
; CALLING SEQUENCE:
;   hires_zerolvl2d, file, outfil, NEWFIL=, INFIL=, ORDRS=, _EXTRA=extra
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
;  /USE_PREV  -- Use previous order as a guess
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

pro hires_zerolvl2d, file, outfil, NEWFIL=newfil, USE_PREV=use_prev, INFIL=infil, $
                   ORDRS=ordrs, _EXTRA=extra

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_zerolvl2d, file, outfil, NEWFIL=, INFIL=, /USE_PREV, ORDRS= [v1.1]'
      return
  endif 

  ;; Optional
  if not keyword_set(DUMMF) then dummf = 'dumfz.fits'

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

  ;; Create conti array
  npx = sz[0]
  if not keyword_set(INFIL) then begin
      zerolvl = fltarr(npx, sz[1])
      zerolvl[*] = 1.
  endif else zerolvl = xmrdfits(infil,/silent)

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

      ;; Create dummy file
      mwrfits, spec2d.fx[gdp,qq], dummf, /create
      mwrfits, replicate(1.,np), dummf
      mwrfits, spec2d.wave[gdp,qq], dummf

      ;; INFIL
      if keyword_set(USE_PREV) AND jj GT 0 then $
         cfil = 'dumz.fits'

      ;; Run conti
      x_spec_zerolvl, dummf, inflg = 2, outfil = 'dumz.fits', ZRO_LVL = cfil, $
                     _EXTRA=extra
      if keyword_set(CFIL) then delvarx, cfil

      ;; Read and save
      sub_cont = xmrdfits('dumz.fits',/silent)
      zerolvl[gdp,qq] = sub_cont

      ;; Save
      mwrfits, zerolvl, outfil, /create
  endfor

  mwrfits, zerolvl, outfil, /create
  spawn, 'gzip -f '+outfil

  ;; Normalize as desired
  if keyword_set(NEWFIL) then begin
      spec2d.fx = spec2d.fx - zerolvl
      spec2d.var = spec2d.var
      IF TAG_EXIST(spec2d, 'NOVAR') THEN $
        spec2d.novar = spec2d.novar
      mwrfits, spec2d, newfil, /create
  endif

  return 

end
      
