;+ 
; NAME:
; hires_mk_makee_template
;     Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  hires_mktempl, fitfil, ordrs, outfil
;
; INPUTS:
;  fitfil - IDL save file output by hires_fitarc
;  ordrs  - Orders to be archived
;  outfil - Archive file
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Aug-2005 Created by JXP 
;-
;; hires_mk_makee_template, 'new-Arc-022.fits', 'hires_tmpl2x1S0.idl'
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_mk_makee_template, specfil, outfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_mk_makee_template, idsfil, specfil, iordr, outfil [v1.0]'
      return
  endif 
  
  echnumlist = getenv('XIDL_DIR')+'/Keck/HIRES/Redux/pro/Arcs/Templates/hires_ech_num.dat'
  readcol,echnumlist,ewcen,eordr,format='f,l'

  ;; Read file
  spec = xmrdfits(specfil, 0, head)
  sz = size(spec, /dimen)
  nordr = sz[1]
  npix = long(sz[0])  ;; As generated
  rbin = round(2048./npix)
  if rbin NE 1 then stop

  sv_aspec = fltarr(2048L, 50) ;; This may need to stay 50 (JXP)

  ;; Genreate wavelength array
  hires_rdmakeewave, head, wave, offset=0.

  ;; Orders
  sz=size(wave,/dimensions)
  wmin=wave[floor(sz[0]/2),0]
  wmax=wave[floor(sz[0]/2),sz[1]-1]
  wcen=wave[floor(sz[0]/2),floor(sz[1]/2)]
  wmn=min(abs(wmin-ewcen),indx)
  iordr=eordr[indx]

  guess_ordr = iordr - lindgen(nordr) ;; This is the bluest order

  ;; Dummy structure
  rejtmp = { $
           ngdf: 0L, $
           gdfpt: lonarr(100), $
           gdfpx: dblarr(100), $
           nrej: 0L, $
           rejpt: lonarr(100), $
           rejwv: dblarr(100) $
           }
  rejstr = replicate(rejtmp, 50)
  lintmp = { $
           pix: dblarr(90), $
           wv: dblarr(90), $
           nlin: 0 $
           }
  sv_lines = replicate(lintmp, 50)

  ;; Loop
  all_arcfit = replicate({fitstrct}, nordr)
  for ii=0L,nordr-1 do begin
     bck_i = nordr-ii-1
     ;; Create fit
     fin_fit = {fitstrct}
     fin_fit.func = 'LEGEND'
     fin_fit.nord = 5L
     fin_fit.flg_rej = 1 
     fin_fit.niter = 3 
     fin_fit.maxrej = 10 
     fin_fit.minpt = 5
     fin_fit.hsig = 3.
     fin_fit.lsig = 3.

     ;; Final auto fitting in log10 space!
     fit = x_fitrej(findgen(npix), alog10(wave[*,ii]), FITSTR=fin_fit)
     all_arcfit[bck_i] = fin_fit
     ;newwv = x_calcfit(findgen(npix), fitstr=fin_fit)
     ;x_splot, 10.^newwv, spec[*,ii], xtwo=wave[*,ii], ytwo=spec[*,ii], /bloc, psym1=10, psym2=10
     ;stop
     
     ;; sv_aspec
     sv_aspec[0:npix-1,bck_i] = spec[*,ii]

  endfor

  ;; Save
  guess_ordr = reverse(guess_ordr)
  save, guess_ordr, rejstr, all_arcfit, sv_lines, sv_aspec, $
    filename=outfil, /compress
  print, 'hires_mk_makee_template: Generated ', outfil

  return

end

