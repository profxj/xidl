;+ 
; NAME:
; imacsls_combspec
;    Version 1.1
;
; PURPOSE:
;   Coadd IMACS long slit spectra.  The routine finds a scaling
;  and the relative S/N between multiple exposures.  It then calls
;  x_coadd to combine the spectra accordingly.
;
; CALLING SEQUENCE:
;   imacsls_combspec, imacsls, setup, side, obj_id, [exp_id]
;
; INPUTS:
;   imacsls  -  IMACS structure
;   setup    -  Setup ID value
;   side     -  Specific CCD (1='blue', 2='red')
;   obj_id   -  Object ID value
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
;   imacsls_combspec, imacsls, 0, 1, 2
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by GEP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro imacsls_combspec, imacsls, setup, obj_id, side, exp_id

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'imacsls_combspec, imacsls, setup, obj_id, side, [exp_id] [v1.1]'
    return
  endif 

; Optional Keywords

; Set exp
  allexp = where(imacsls.type EQ 'OBJ' AND imacsls.flg_anly NE 0 AND $
                 imacsls.side EQ side AND imacsls.obj_id EQ obj_id $
                 AND imacsls.mode EQ 1 AND imacsls.setup EQ setup)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
  nexp = n_elements(exp)

  ;; Set x1, x2
  x1 = round(2000./imacsls[exp[0]].rbin)
  x2 = round(4000./imacsls[exp[0]].rbin)

; Obj Structure
  for i=0L,nexp-1 do begin
    tmp = xmrdfits(imacsls[exp[i]].obj_fil, 1, /silent)
    if i EQ 0 then imacslsobj = tmp else imacslsobj = [imacslsobj, tmp]
  endfor
  allobj = where(imacslsobj.exp NE 0)

; OUTPUT
  outfil = 'Extract/Fspec_'+imacsls[exp[0]].img_root

; CREATE FINAL 2D 

  imacslsspec = { lwdfspecstrct }

  ; Copy
  imacslsspec.nexp = nexp

  ; Set texp
  for i=0L,nexp-1 do imacslsspec.texp[i] = imacslsobj[i].exp

  ; other tags
  copy_struct, imacslsobj[allobj[0]], imacslsspec, EXCEPT=["wave","fx","var","npix"]
  for i=0L,nexp-1 do begin
    ipos = strpos(imacslsobj[allobj[i]].spec2d_fil, 'imacsls')
;    obj_fil = 'FSpec/Obj_'+strmid(imacslsobj[allobj[i]].spec2d_fil, ipos)
    obj_fil = 'Extract/Obj_'+strmid(imacslsobj[allobj[i]].spec2d_fil, ipos)
    imacslsspec.obj_fil[i] = obj_fil
  endfor

  svmedian = fltarr(nexp)
  svsnr = fltarr(nexp)

  if nexp EQ 1 then begin
    npix = imacslsobj[allobj].npix
    imacslsspec.wave[0:npix-1] = imacslsobj[allobj].wave[0:npix-1]
    imacslsspec.fx[0:npix-1] = imacslsobj[allobj].flux[0:npix-1]
    imacslsspec.var[0:npix-1] = imacslsobj[allobj].sig[0:npix-1]^2
    imacslsspec.npix = npix
  endif else begin
      npix = imacslsobj[allobj[0]].npix
                                ; Find median ratio of flux and SNR
      med_rto = fltarr(nexp)
      snr_rto = fltarr(nexp)
      med_rto[0] = 1.
      dumi = lindgen(npix)
      gdpix = where(imacslsobj[allobj].sig GT 0. AND $
                   dumi GE x1 AND dumi LE x2)
      fx0 = imacslsobj[allobj[0]].flux[gdpix]
      var0 = imacslsobj[allobj[0]].sig[gdpix]^2
      snr_rto[0] = median(fx0 / sqrt(var0))
      for kk=1L,nexp-1 do begin
          tmpfx = imacslsobj[allobj[kk]].flux[gdpix]
          tmpvar = imacslsobj[allobj[kk]].sig[gdpix]^2
          med_rto[kk] = median(fx0/tmpfx)
          snr_rto[kk] = median(tmpfx/sqrt(tmpvar))
      endfor
      ;; Save
      svmedian[*] = med_rto
      ;; Check this  12/9/03
      svsnr[*] = snr_rto
      print, 'imacsls_combspec: SNR', svsnr
      
      ;; set bad pix to -1
      var = double(imacslsobj[allobj].sig[0:npix-1]^2)
;      for kk=0L,nexp-1 do begin
;          a = where(imacslsobj[allobj[kk]].flux[0:npix-1] LE 0., na)
;          if na GT 0 then var[a,kk] = -1.
;      endfor
      x_combspec, imacslsobj[allobj].flux[0:npix-1], var, fflux, fvar, $
        SCALE=med_rto, SNR=snr_rto
      imacslsspec.wave[0:npix-1] = imacslsobj[allobj[0]].wave[0:npix-1]
      imacslsspec.fx[0:npix-1] = temporary(fflux[0:npix-1])
      imacslsspec.var[0:npix-1] = temporary(fvar[0:npix-1])
      imacslsspec.npix = npix
  endelse 

  x_splot, imacslsspec.wave[0:npix-1], $
    imacslsspec.fx[0:npix-1], YTWO=sqrt(imacslsspec.var[0:npix-1]), /block
  print, 'imacsls_comspec: Updating: ', outfil
  imacsls_wrspec, imacslsspec, outfil
  print, 'imacsls_combspec: All done!'

  return

end


