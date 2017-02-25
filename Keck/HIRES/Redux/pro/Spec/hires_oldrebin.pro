;+ 
; NAME:
; hires_oldrebin
;    Version 1.1
;
; PURPOSE:
;   Rebin the MAKEE extracted onto a common wavlength scale
;
; CALLING SEQUENCE:
;   hires_combspec, hires, setup, obj_id, chip, exp_id
;
; INPUTS:
;   fil      - MAKEE 2D file
;
; RETURNS:
;
; OUTPUTS:
;   outfil     -  HIRES structure of the orders, rebinned
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hires_oldrebin, '3C196-2D_f.fits', '3C196_rebf.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Sep-2006 Written by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hires_oldrebin, fil, sigfil, contifil, outfil, FSPEC=fspec, $
                    PHYSORDR=physordr, TRIM=trim

;
  if  N_params() LT 4 and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'hires_oldrebin, fil, sigfil, contifil, outfil,   LIST= [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(CRVAL1) then crval1 = 3000.d
  if not keyword_set(CDELT) then cdelt = 2.897300918105E-06
  if not keyword_set(NPIX) then npix = 500000L
  if not keyword_set(PHYSORDR) then physordr = 90L

  ;; Read in the data
  fx = xmrdfits(fil,0,head,/sil)
  sig = xmrdfits(sigfil,/sil)
  conti = xmrdfits(contifil,/sil)


  sz = size(fx,/dime)
  nopix = sz[0]
  nordr = sz[1]

  ;; Set the wavelengths
  wv = hires_oldwave(head)
  tot_wave = 10^(alog10(CRVAL1) + dindgen(npix)*cdelt)

  ;; Create the final sturcture
  if not keyword_set(FSPEC) then begin
      objstr = replicate({hiresobjstrct}, nordr)
      objstr.order = physordr - lindgen(nordr)
      objstr.flg_anly = 1
      objstr.UT = ' '
      objstr.field = ' '
      objstr.img_fil = ' '
      objstr.arc_fil = ' '
      objstr.obj_id = 'a'
  endif else begin
      objstr = { hiresfspecstrct }
  endelse

  ;; Rebin onto a common spectrum (starts at 3000 with 2.1 km/s
       ;; pixels)

  for qq=0L,nordr-1 do begin
      ;; Find starting/ending wave
      mn = min(abs(tot_wave-wv[0,qq]),imi)
      mn = min(abs(tot_wave-wv[nopix-1,qq]),imf)
      npo = imf-imi + 3

      ;; TRIM
      if keyword_set(TRIM) then begin
          sig[0:45,qq] = 0.
      endif

      ;; Rebin flux + var
      x_specrebin, wv[*,qq], fx[*,qq], tot_wave[imi-1:imf+1], nwfx, $
                   VAR=(sig[*,qq]>0.)^2, NWVAR=nwvar

      ;; Rebin continuum
      x_specrebin, wv[*,qq], conti[*,qq], tot_wave[imi-1:imf+1], nwcnt

      ;; Save
      if keyword_set(FSPEC) then begin
          kk = physordr - qq
          objstr.wave[0:npo-1,kk] = tot_wave[imi-1:imf+1]
          objstr.fx[0:npo-1,kk] = nwfx/nwcnt
          objstr.var[0:npo-1,kk] = nwvar/(nwcnt^2)
          objstr.npix[kk] = npo
          objstr.phys_ordr[kk] = kk
      endif else begin
          objstr[qq].wave[0:npo-1] = tot_wave[imi-1:imf+1]
          objstr[qq].fx[0:npo-1] = nwfx/nwcnt
          objstr[qq].var[0:npo-1] = nwvar/(nwcnt^2)
          objstr[qq].npix = npo
      endelse
  endfor

  if keyword_set(FSPEC) then begin
      hires_wrfspec, objstr, outfil 
  endif else begin
      mwrfits, objstr, outfil, /create
  endelse
      
  spawn, 'gzip -f '+outfil

  print, 'hires_oldrebin:  Output is in ', outfil
  print, 'hires_oldrebin:  All done!'


  return
end
  

