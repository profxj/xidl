;+ 
; NAME:
; civ_readspec
;    Version 1.0
;
; PURPOSE:
;   Concatenate/splice the spectra from different
;              instruments (STIS and FUSE) for one object into
;              one set of flux, error, wavelength; if \nosplice
;              NOT set, then based on S/N splice spectra
;
; CALLING SEQUENCE:
;   
;   civ_readspec,'lists/3C273instr.lst',strct=strct
;
; INPUTS:
;    instrfil -- ASCII list of spectra in typical format
;
; RETURNS:
;    wave, flux, error -- arrays
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;    nosplice -- return just concatenated flux, error, wave,...
;
; OPTIONAL OUTPUTS:
;    instr -- array of instrument flags corresponding to each
;             piece in concatenated/spliced spectrum
;    outfil -- name of file to write structure (FUSE format)
;    strct -- structure of format {wave,flux,sigma,instr}
;
; COMMENTS:
;
; EXAMPLES:
;   civ_readspec, x, maskid, expsr
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Apr-2008 Written by JXP
;   24-Apr-2008 Modified by KLC to splice data based on S/N
;-
;------------------------------------------------------------------------------
@coadSTIS ; compile measureSN()

pro civ_readspec, file_list, wave, flux, sigma, nosplice=nosplice, $
                  outfil=outfil, instr=instr, strct=strct

readcol, file_list, spec, dw, w0, format='a,f,f'
spec = getenv('MLSS_DIR')+'/'+strtrim(spec,2)
nspec = n_elements(spec)

snr = fltarr(nspec)
wvlim = fltarr(nspec,2)

tmpwv = -1
for ii=0,nspec-1 do begin
    test = file_search(spec[ii],count=ntest)
    if ntest eq 0 then continue ; no spectrum

    if ii le 6 then begin
        ;; FUSE
        fx = x_readspec(spec[ii],wav=wv,sig=err,inflg=3,npix=npix)
    endif else begin            
        ;; HST
;        if strmatch(spec[ii],'F.fits') then $
;          efil = strmid(spec[ii],0,strpos(spec[ii],'F.fits',/reverse_search))+$
;          'E.fits' $
;        else efil = strmid(spec[ii],0,strpos(spec[ii],'f.fits',$
;                                             /reverse_search))+'e.fits' 
;        fx = x_readspec(spec[ii],wav=wv,sig=err,fil_sig=efil,npix=npix)
        fx = x_readspec(spec[ii],wav=wv,sig=err,npix=npix, /AUTO)
    endelse 
    
    ;; Shift
    wv = wv + dw[ii]*w0[ii]/wv

    ;; Measure S/N
    rslt = measureSN(fx,err=err,inflg=4,wave=wv,/silent)

    ;; Store
    snr[ii] = rslt[2]
    wvlim[ii,*] = [min(wv,max=mx),mx]
    if tmpwv[0] eq -1 then begin
        tmpwv = wv
        tmpfx = fx
        tmperr = err
        tmpinst = replicate(2^ii,npix)
        
    endif else begin
        tmpwv = [tmpwv,wv]
        tmpfx = [tmpfx,fx]
        tmperr = [tmperr,err]
        tmpinst = [tmpinst,replicate(2^ii,npix)]
        
    endelse 
endfor                          ;loop spectra

srt = sort(tmpwv)
wave = tmpwv[srt]
flux = tmpfx[srt]
sigma = tmperr[srt]
instr = tmpinst[srt]
if keyword_set(strct) then $
  strct = {wave:wave, flux:flux, sigma:sigma, instr:instr}

;; Return full concatenated arrays without splicing
if keyword_set(nosplice) then begin
    if keyword_set(outfil) then begin
        if not keyword_set(strct) then $
          strct = {wave:wave, flux:flux, sigma:sigma, instr:instr}
        mwrfits,strct,outfil,/create,/silent
        print,'civ_readspec: created ',outfil
    endif 
    
    return
endif 


;; Splice the spectra based on S/N
;; (Can't do much about edges being bad...)
srt = reverse(sort(snr))
tmpwv = -1
mask = replicate(0,n_elements(wave)) ;mask out already considered pixels
for ii=0,nspec-1 do begin
    jj = srt[ii]
    gd = where(instr eq 2^jj and mask eq 0,ngd)
    if ngd eq 0 then continue

    if tmpwv[0] eq -1 then begin
        tmpwv = wave[gd]
        tmpfx = flux[gd]
        tmperr = sigma[gd]
        tmpinst = replicate(2^jj,ngd)
    endif else begin
        tmpwv = [tmpwv,wave[gd]]
        tmpfx = [tmpfx,flux[gd]]
        tmperr = [tmperr,sigma[gd]]
        tmpinst = [tmpinst,replicate(2^jj,ngd)]
    endelse 

    ;; No longer valid is current ranked S/N spec or anything within
    ;; its bounds
    used = where(instr eq 2^jj or $
                 ( wave ge wvlim[jj,0] and wave le wvlim[jj,1] ),nused)
    mask[used] = 1
endfor                          ;loop spectra to splice

srt = sort(tmpwv)
wave = tmpwv[srt]
flux = tmpfx[srt]
sigma = tmperr[srt]
instr = tmpinst[srt]
if keyword_set(strct) then $
  strct = {wave:wave, flux:flux, sigma:sigma, instr:instr}

if keyword_set(outfil) then begin
    if not keyword_set(strct) then $
      strct = {wave:wave, flux:flux, sigma:sigma, instr:instr}
    mwrfits,strct,outfil,/create,/silent
    print,'civ_readspec: created ',outfil
endif 

return
end


