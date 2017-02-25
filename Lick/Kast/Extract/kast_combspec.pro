;+ 
; NAME:
; kast_combspec
;    Version 1.1
;
; PURPOSE:
;   Coadd kast spectra using x_combspec
;
; CALLING SEQUENCE:
;   kast_combspec, kast, setup, side, obj_id, [exp_id]
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;   side  --  Specific camera [blue (1) vs. red (2)]
; obj_id  --  Object value
;  [exp_id]  --  Exposure indices
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Show the final combined spectrum
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_combspec, kast, 0, 1, 2
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_combspec
;  kast_wrspec
;
; REVISION HISTORY:
;   28-August-2003 Written by GEP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_combspec, kast, setup, side, obj_id, exp_id, CHK=chk, STD=std, $
                   alrdyfx=alrdyfx

;
if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'kast_combspec, kast, setup, side, obj_id, [exp_id], /CHK [v1.1]'
    return
endif 

; Optional Keywords

; Set exp
if not keyword_Set(std) then $
  allexp = where(kast.type EQ 'OBJ' AND kast.flg_anly NE 0 AND $
                 kast.side EQ side AND kast.obj_id EQ obj_id $
                 AND kast.mode EQ 1 AND kast.setup EQ setup) $
else allexp = where(kast.type EQ 'STD' AND kast.flg_anly NE 0 AND $
                    kast.side EQ side AND kast.obj_id EQ obj_id $
                    AND kast.mode EQ 1 AND kast.setup EQ setup) 

                                ;if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
if arg_present(exp_id) then exp = allexp[exp_id] else exp=allexp
nexp = n_elements(exp)

; Obj Structure
for i=0L,nexp-1 do begin
    tmp = xmrdfits(kast[exp[i]].obj_fil, 1, STRUCTYP='specobjstrct', /silent)
    if i EQ 0 then kastobj = tmp else kastobj = [kastobj, tmp]
endfor
allobj = where(kastobj.exp NE 0)

; OUTPUT
outfil = 'Extract/Fspec_'+kast[exp[0]].img_root

; CREATE FINAL 2D 

kastspec = { kastspecstrct }

                                ; Copy
kastspec.nexp = nexp

                                ; Set texp
for i=0L,nexp-1 do kastspec.texp[i] = kastobj[i].exp

                                ; other tags
copy_struct, kastobj[allobj[0]], kastspec, $
  EXCEPT=["wave","fx","var","npix", "flux","sig"]
for i=0L,nexp-1 do begin
    ipos = strpos(kastobj[allobj[i]].spec2d_fil, 'kast')
;    obj_fil = 'FSpec/Obj_'+strmid(kastobj[allobj[i]].spec2d_fil, ipos)
    obj_fil = 'Extract/Obj_'+strmid(kastobj[allobj[i]].spec2d_fil, ipos)
    kastspec.obj_fil[i] = obj_fil
endfor

if not keyword_set(alrdyfx) then begin
    svmedian = fltarr(nexp)
    svsnr = fltarr(nexp)

    if nexp EQ 1 then begin
        npix = kastobj[allobj].npix
        kastspec.wave[0:npix-1] = kastobj[allobj].wave[0:npix-1]
        kastspec.fx[0:npix-1] = kastobj[allobj].fx[0:npix-1]
        kastspec.var[0:npix-1] = kastobj[allobj].var[0:npix-1]
        kastspec.npix = npix
    endif else begin
        npix = kastobj[allobj[0]].npix
                                ; Find median ratio of flux and SNR
        med_rto = fltarr(nexp)
        snr_rto = fltarr(nexp)
        med_rto[0] = 1.
        snr_rto[0] = 1.
        gdpix = where(kastobj[allobj].var GT 0.)
        fx0 = kastobj[allobj[0]].fx[gdpix]
        var0 = kastobj[allobj[0]].var[gdpix]
        fx0 = fx0[300:1000]
        var0 = var0[300:1000]
        for kk=1L,nexp-1 do begin
            tmpfx = kastobj[allobj[kk]].fx[gdpix]
            tmpvar = kastobj[allobj[kk]].var[gdpix]
            tmpfx = tmpfx[300:1000]
            tmpvar = tmpvar[300:1000]
            med_rto[kk] = median(fx0/tmpfx)
            snr_rto[kk] = median((fx0*sqrt(tmpvar))/(tmpfx*sqrt(var0)))
        endfor
        ;; Save
        svmedian[*] = med_rto
        svsnr[*] = snr_rto

        ;; set bad pix to -1
        var = double(kastobj[allobj].var[0:npix-1])
        for kk=0L,nexp-1 do begin
            a = where(kastobj[allobj[kk]].fx[0:npix-1] LE 0., na)
            if na GT 0 then var[a,kk] = -1.
        endfor
        x_combspec, kastobj[allobj].fx[0:npix-1], var, fflux, fvar, $
          SCALE=med_rto, SNR=snr_rto
        kastspec.wave[0:npix-1] = kastobj[allobj[0]].wave[0:npix-1]
        kastspec.fx[0:npix-1] = temporary(fflux[0:npix-1])
        kastspec.var[0:npix-1] = temporary(fvar[0:npix-1])
        kastspec.npix = npix
    endelse 

endif else begin
;; Combine the flux
;; Added by KLC, 17 Feb 2007 (out of necessity)
    sffil = 'Extract/Fspec_S'+kast[exp[0]].img_root
    spawn,'\cp '+outfil+' '+sffil
    print,'kast_combspec: updating combined flux'
    kast_wrspec, kastspec, outfil, /READ

    svmedian = fltarr(nexp)
    svsnr = fltarr(nexp)

    if nexp EQ 1 then begin
        npix = kastobj[allobj].npix
        kastspec.flux[0:npix-1] = kastobj[allobj].flux[0:npix-1]
        kastspec.var[0:npix-1] = (kastobj[allobj].sig[0:npix-1])^2
    endif else begin
        npix = kastobj[allobj[0]].npix
        ;; Find median ratio of flux and SNR
        med_rto = fltarr(nexp)
        snr_rto = fltarr(nexp)
        med_rto[0] = 1.
        snr_rto[0] = 1.
        gdpix = where(kastobj[allobj].var GT 0.)
        flux0 = kastobj[allobj[0]].flux[gdpix]
        var0 = (kastobj[allobj[0]].sig[gdpix])^2
        flux0 = flux0[300:1000]
        var0 = var0[300:1000]
        for kk=1L,nexp-1 do begin
            tmpflux = kastobj[allobj[kk]].flux[gdpix]
            tmpvar = (kastobj[allobj[kk]].sig[gdpix])^2
            tmpflux = tmpflux[300:1000]
            tmpvar = tmpvar[300:1000]
            med_rto[kk] = median(flux0/tmpflux)
            snr_rto[kk] = median((flux0*sqrt(tmpvar))/(tmpflux*sqrt(var0)))
        endfor
        ;; Save
        svmedian[*] = med_rto
        svsnr[*] = snr_rto

        ;; set bad pix to -1
        var = double((kastobj[allobj].sig[0:npix-1])^2)
        for kk=0L,nexp-1 do begin
            a = where(kastobj[allobj[kk]].flux[0:npix-1] LE 0., na)
            if na GT 0 then var[a,kk] = -1.
        endfor
        x_combspec, kastobj[allobj].flux[0:npix-1], var, fflux, fvar, $
          SCALE=med_rto, SNR=snr_rto
        kastspec.flux[0:npix-1] = temporary(fflux[0:npix-1])
        kastspec.sig[0:npix-1] = temporary(sqrt(fvar[0:npix-1]))
    endelse 
endelse


if keyword_set( CHK ) then $
  x_splot, kastspec.fx[0:npix-1], $
  YTWO=sqrt(kastspec.var[0:npix-1]), /block
print, 'kast_comspec: Updating: ', outfil
kast_wrspec, kastspec, outfil
print, 'kast_combspec: All done!'

return

end


