;+ 
; NAME:
; x_readspec
;   Version 1.1
;
; PURPOSE:
;    Input the data spectrum from a FITS file.
;
; CALLING SEQUENCE:
;   dat = x_readspec(spec, dataunit, INFLG=, /FSCALE, HEAD=, /DSCALE, SIG=
;     WAV=, NPIX=, FIL_SIG=, STRUCT=)
;
; INPUTS:
;   spec       - Fits file or data
;
; RETURNS:
;   dat       - Data in fits file or the input data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STRUCT  -- Return a data sturcture containing the key info
;  FSCALE      - Data is float
;  DSCALE      - Data is float
;  INFLG: 0  - Data, sig separate fits files  dat =
;              x_readspec(data_fil, FIL_SIG=sig_fil, SIG=sig, WAV=wav)
;         1  - Flux and Sig in one fits file
;         2  - Flux, sig, wave in one fits file
;         3  - FUSE, WFC3 format (Binary table)
;         4  - Flux in first argument, sigma in SIG (data arrays)
;         5  - SDSS format
;         6  - ??
;         7  - DEIMOS spec1d format
;         8  - ASCII  wv, fx, ivar
;         9  - (old) IRAF 4 'channel' format
;        10  - COS output
;        11  - Flux, ivar, wave in one fits file
;        12  - Howk UV
;        13  - Murphy UVES
;        14  - Meiring COS (wave, flux, sig_up, sig_down)
;        15  - IRAM 30 mt class sepctra
;        16  - IRAF 4 'channel' format
;        17  - COS low-dispersion (Howk LLS)
;        18  - APO TripleSpect merged order file
;        19  - SDSS/DR9/DR10 structure  (BOSS)
;  BIN=    -- Factor to bin data by (odd integers only!; COS data only)
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;  ASIG1, ASIG2=  -- 1-sigma CL values for asymmetric error arrays
;                  (e.g. Poisson for COS).  ASIG1 is the UPPER offset.
;
; COMMENTS:
;
; EXAMPLES:
;   dat = x_readspec('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_readspec, spec, dataunit, INFLG=inflg, FSCALE=fscale, $
                     HEAD=head, DSCALE=dscale, SIG=sig, WAV=wav, NPIX=npix, $
                     FIL_SIG=fil_sig, STRUCT=struct, AUTOFSIG=autofsig, $
                     FLG_FIL=flg_fil, VAC=vac, ASIG1=asig1, ASIG2=asig2, BIN=bin $
                     , NONORM = NONORM, CONT = CONT, NOSTOP=nostop

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dat = x_readspec(spec, [dataunit], INFLG= /fscale, /dscale, HEAD=, SIG='
     print, '   WAV=, NPIX=, FIL_SIG= , /STRUCT, /AUTOFSIG, ASIG1=, ASIG2= [v1.1.]'
     print, 'INFLG: 0  - Data, sig separate fits files  dat ='
     print, '         1  - Flux and Sig in one fits file'
     print, '         2  - Flux, sig, wave in one fits file'
     print, '         3  - FUSE, WFC3 format (Binary table)'
     print, '         4  - Flux in first argument, sigma in SIG (data arrays)'
     print, '         5  - SDSS format'
     print, '         6  - ??'
     print, '         7  - DEIMOS spec1d format'
     print, '         8  - ASCII  wv, fx, ivar'
     print, '         9  - (old) IRAF 4 channel format'
     print, '        10  - COS output'
     print, '        11  - Flux, ivar, wave in one fits file'
     print, '        12  - Howk UV'
     print, '        13  - Murphy UVES'
     print, '        14  - Meiring COS (wave, flux, sig_up, sig_down)'
     print, '        15  - IRAM 30 mt class sepctra'
     print, '        16  - IRAF 4 channel format'
     print, '        17  - COS low-dispersion (Howk LLS)'
     print, '        18  - APO TripleSpect merged order file'
     print, '        19  - SDSS/DR9/DR10 structure  (BOSS)'
     return, -1
  endif 

;  Optional Keywords
  if not keyword_set( dataunit ) then dataunit = 0L
  if not keyword_set( inflg ) then inflg = 0 
  if keyword_set(DSCALE) and keyword_set(FSCALE) then fscale=0
  if keyword_set(AUTOFSIG) and not keyword_set(FIL_SIG) then begin
      pos = strpos(spec, '_f') 
      if pos GT 0 then fil_sig = strmid(spec,0,pos+1)+'e.fits' else begin
          pos = strpos(spec, 'F.fits') 
          fil_sig = strmid(spec,0,pos)+'E.fits' 
      endelse
  endif

  flg_fil = 1

; String or Image?
  if size(spec, /type) EQ 7 then begin
      spec = strtrim(spec,2)
      a = findfile(spec+'*', count=count)
      if count EQ 0 then begin
          print, 'x_readspec: '+spec+' does not exist'
          flg_fil = 0
          if not keyword_set(NOSTOP) then stop
          return, -1
      endif
      flg_fil = 1
  endif

  case inflg of
      0: begin  ; NORMAL -- data, sig separate
          dat = xmrdfits(spec, dataunit, head, fscale=fscale, $
                         dscale=dscale, /silent, FLG_FIL=flg_fil)
          npix = n_elements(dat)
          ;; SIG
          if keyword_set( FIL_SIG ) then sig = xmrdfits(fil_sig, /silent, $
                                                       FLG_FIL=flg_fil) $
          else sig = fltarr(npix)
          ;; Wavelengths
          if arg_present(WAV) or keyword_set(STRUCT) then $
             wav = x_setwave(head, npix) 
      end
      1: begin ; One fits file (flux, sig)
          dat = x_readspec(spec, 0L, /fscale, head=head, FLG_FIL=flg_fil)
          npix = n_elements(dat)
          if arg_present(SIG) or keyword_set(STRUCT) then $
            sig = x_readspec(spec, 1L, fscale=fscale)
          ;; Wavelengths
          if arg_present(WAV) or keyword_set(STRUCT) then $
            wav = x_setwave(head, npix) 
      end
      2: begin ; FLUX, SIG, WAVE (one fits file)
          dat = xmrdfits(spec, 0L, head, FSCALE=fscale, DSCALE=dscale, $
                        FLG_FIL=flg_fil, /silent)
          npix = n_elements(dat)
          if arg_present(SIG) or keyword_set(STRUCT) then $
            sig = xmrdfits(spec, 1L, /silent)
          if arg_present(WAV) or keyword_set(STRUCT) then $
            wav = xmrdfits(spec, 2L, /silent)
      end
      3: begin ; FUSE format (and others)
          if not keyword_set( dataunit) then dataunit = 1
          fuse = xmrdfits(spec, dataunit, head, /silent)
          if arg_present(WAV) or keyword_set(STRUCT) then begin
              if tag_exist(fuse, 'WAVE') then wav = fuse.wave $
              else wav = fuse.wavelength
          endif
          ;; Flux
          if tag_exist(fuse, 'FLUX') then dat = fuse.flux
          if tag_exist(fuse, 'FLAM') then dat = fuse.flam
          ;; Error
          if arg_present(SIG) or keyword_set(STRUCT) then begin
             if tag_exist(fuse, 'ERROR') then sig = fuse.error
             if tag_exist(fuse, 'ERR') then sig = fuse.err
             if tag_exist(fuse, 'SIGMA_FLUX') then sig = fuse.sigma_flux
             if tag_exist(fuse, 'FLAM_SIG') then sig = fuse.flam_sig
          endif
          npix = n_elements(dat)
      end
      4: begin ; 2 arrays fx, sig!
          if keyword_set( FIL_SIG ) then sig = fil_sig
          dat = spec
          npix = n_elements(dat)
          if arg_present(WAV) then wav = findgen(npix)
      end
      5: begin ; SDSS
          ;; Parse
         parse_sdss, spec, dat, wav, SIG = sig, HEAD = head, NPIX = npix
      end
      6: begin ; Low d structure
          ;; Parse
          strct = xmrdfits(spec, 1L,/silent)
          npix = strct.npix
          dat = strct.fx[0:npix-1]
          if arg_present(SIG) then sig = sqrt(strct.var[0:npix-1] > 0.)
          if arg_present(WAV) then wav = strct.wave[0:npix-1]
      end
      7: begin ; Low d structure
          ;; Parse
          if not keyword_set( dataunit) then dataunit = 1
          strct = xmrdfits(spec, dataunit, /silent)
          npix = n_elements(strct.spec)
          dat = strct.spec
          if arg_present(SIG) then begin
             gdi = where(strct.ivar GT 0.)
             sig = fltarr(npix)
             sig[gdi] = sqrt(1./strct.ivar[gdi])
          endif
          if arg_present(WAV) then wav = strct.lambda
      end
      8: begin ; ASCII file
          ;; Parse
          readcol, spec, wav, dat, ivar, /silent
          if keyword_set(VAC) then airtovac, wav
          npix = n_elements(wav)
          if arg_present(SIG) then sig = sqrt(1./ivar)
      end
      9: begin ; IRAF format
          ;; Parse
          idata = xmrdfits(spec, dataunit, head)
          sz = size(idata, /dimension)
          npix = sz[0]
          wav = x_setwave(head, npix) 
;          if keyword_set(VAC) then airtovac, wav
          dat = idata[*,0,0]
          sig = idata[*,0,1]
;          if arg_present(SIG) then sig = sqrt(1./ivar)
      end
      10: begin ; COS output (count rate, counts, etc.)
          ;; Parse
          if not keyword_set( dataunit) then dataunit = 1
          strct = xmrdfits(spec, dataunit)
          if n_elements(strct) GT 1 then stop  ;; Not prepared for this

          ;; Prepare for binning
          counts = strct.counts
          eff_t = strct.eff_t
          dat = strct.cnt_rate
          npix_native = n_elements(dat)
          wave = strct.wavelength

          ;; Background
          bkgd = fltarr(npix_native)
          gd = where(strct.eff_t GT 0.)
          bkgd[gd] = strct.counts[gd]/strct.eff_t[gd] - strct.cnt_rate[gd]  ;; Bkgd count rate

          ;; BIN?
          if keyword_set(BIN) then begin
             ;; Odd integer?
             if (bin MOD 2) NE 1 then stop

             idx = bin*lindgen(npix_native/bin - 1)

             ;; Wavelength and background
             wav = (smooth(wave,bin))[idx]
             npix = n_elements(wav)
             bkgd = (smooth(bkgd,bin))[idx]
             
             ;; Add up counts, teff
             counts = (bin*smooth(counts,bin))[idx]
             eff_t = (bin*smooth(eff_t,bin))[idx]
             
             ;; Count rate
             dat = fltarr(npix)
             gd = where(eff_t GT 0)
             dat[gd] = counts[gd]/eff_t[gd] - bkgd

;          if keyword_set(VAC) then airtovac, wav
          endif else begin
             wav = wave
          endelse

          ;; Proceed
          npix = n_elements(wav)

          ;; Fill up a Poisson CL array (should do a look-up table)
          mx_count = max(counts)
          xpcl = fltarr(mx_count+1, 2)
          for ii=0L,mx_count do xpcl[ii,*] = x_poisscl(ii,sigma=1)
          
          ;; Set the error arrays
          gd = where(eff_t GT 0)
          asig1=fltarr(npix)
          asig2=fltarr(npix)
          asig1 = xpcl[counts[gd],0]/eff_t[gd] - bkgd[gd] - dat[gd]
          asig2 = dat[gd] - (xpcl[counts[gd],1]/eff_t[gd] - bkgd[gd])
          
          sig = (asig1+asig2)/2.
       end
      11: begin                 ; FLUX, IVAR, WAVE (one fits file)
         dat = xmrdfits(spec, 0L, head, FSCALE = fscale, DSCALE = dscale, $
                        FLG_FIL = flg_fil, /silent)
         npix = n_elements(dat)
         if arg_present(SIG) or keyword_set(STRUCT) then BEGIN
            ivar = xmrdfits(spec, 1L, /silent)
            sig = (ivar GT 0.0)/(sqrt(ivar > 0) + (ivar LE 0.0))
         ENDIF
         if arg_present(WAV) or keyword_set(STRUCT) then BEGIN
            loglam = xmrdfits(spec, 2L, /silent)
            wav = 10.0d^loglam
         ENDIF
      end
      12: begin  ;; Howk UV spectra
         ;; Parse
         if not keyword_set(dataunit) then dataunit = 1
         strct = xmrdfits(spec, dataunit, /silent)
         npix = n_elements(strct.flux)
         dat = strct.flux
         IF arg_present(SIG) THEN BEGIN
            sig = strct.err
            imask = WHERE(finite(dat) EQ 0 OR finite(sig) EQ 0, nmask)
            IF nmask GT 0 THEN BEGIN
               dat[imask] = 0.0d
               sig[imask] = 0.0d
            ENDIF
         ENDIF ELSE BEGIN
             imask = WHERE(finite(dat) EQ 0, nmask)
             IF nmask GT 0 THEN dat[imask] = 0.0d
          ENDELSE
         if arg_present(WAV) then wav = strct.wave
      end
      13: begin            ;; Michael Murphy UVES reductions
         array = xmrdfits(spec, 0L, /fscale, head, /silen)
         dims = size(array, /dim)
         npix = dims[0]
         dat = array[*, 0]
         if arg_present(SIG) or keyword_set(STRUCT) then sig = array[*, 1]
         ;; Wavelengths
         if arg_present(WAV) or keyword_set(STRUCT) then wav = x_setwave(head, npix, /silent) 
         CONT = array[*, 3]
         IF KEYWORD_SET(NONORM) THEN BEGIN
            dat = CONT*dat
            IF KEYWORD_SET(SIG) THEN sig = CONT*sig
         ENDIF
      end
      14: begin ; WAVE, FLUX, SIG_UP, SIG_DOWN
          wav = xmrdfits(spec, 0L, head, FSCALE=fscale, DSCALE=dscale, $
                        FLG_FIL=flg_fil, /silent)
          npix = n_elements(wav)
          dat = xmrdfits(spec, 1L, /silent)
          asig1 = xmrdfits(spec, 2L, /silent)
          asig2 = xmrdfits(spec, 3L, /silent)
          sig = (asig1+asig2)/2.
       end
      15: begin  ;IRAM DATA
         dat=mrdfits(spec,0L,head,/fsc,/sil)
         npix = n_elements(dat)
         if arg_present(SIG) or keyword_set(STRUCT) then $
            sig = dat-dat
         if arg_present(WAV) or keyword_set(STRUCT) then begin
            RESTFREQ=fxpar(head,"RESTFREQ")
            CRVAL1=fxpar(head,"CRVAL1")
            CRPIX1=fxpar(head,"CRPIX1")
            CDELT1=fxpar(head,"CDELT1")
            INDX=findgen(n_elements(dat))
            freq=(RESTFREQ+CRVAL1+(INDX-CRPIX1)*CDELT1)/1d9 ;GHz
            wav=2.99792458000D11/(freq*1d9)                 ;mm
         endif
      end
      16: begin ; IRAF format
          ;; Parse
          idata = xmrdfits(spec, dataunit, head)
          sz = size(idata, /dimension)
          npix = sz[0]
          wav = x_setwave(head, npix) 
;          if keyword_set(VAC) then airtovac, wav
          dat = idata[*,0,0]
          sig = idata[*,0,3]
;          if arg_present(SIG) then sig = sqrt(1./ivar)
      end
      17: begin ; COS 
          cos = xmrdfits(spec, 1, head, /silent)
          cosA = cos[0] ;; segments
          cosB = cos[1]
          ;; Wavelength
          wav = [cosA.wavelength, cosB.wavelength]
          ;; Flux
          dat = [cosA.flux, cosB.flux]
          ;; Error
          sig = [cosA.error, cosB.error]
          ;; Sort
          srt = sort(wav)
          wav = wav[srt]
          dat = dat[srt]
          sig = sig[srt]
          ;; CUT
          good  = where(wav GT 1000, npix)
          wav=wav[good]
          dat=dat[good]
          sig=sig[good]
       end
      18: begin                  ; FLUX, SIG, WAVE (one fits file)
         tmp=mrdfits(spec,0L,head,/sil)
         dat = reform(tmp[*,1])
         npix = n_elements(dat)
         sig = reform(tmp[*,2])
         wav = reform(tmp[*,0])*1d4
      end
      19: begin                 ; sdss dr9/dr10 structure
         tmp=xmrdfits(spec,1, /sile)
         dat = tmp.flux
         npix = n_elements(dat)
         sig = sqrt(1/tmp.IVAR)  ;;MF not sure about this line
         wav = 10.^tmp.LOGLAM
      end

      
      ;; ELSE
      else: begin
         print, 'x_readspec: inflg = ', inflg, ' not allowed!'
         stop
         return, -1
      end
  endcase
  if not keyword_set( STRUCT ) then return, dat $
  else begin
      if not keyword_set(SIG) then sig = replicate(1., npix)
      if not keyword_set(WAV) then wav = dindgen(npix)
      str = { $
              npix: npix, $
              fx: dat, $
              sig: sig, $
              velo: dblarr(npix), $
              wv: wav $
            }
      return, str
  endelse
              
end
