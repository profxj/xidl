;+ 
; NAME:
; x_readspec
;   Version 1.0
;
; PURPOSE:
;    Convert input to data whether it is a fits file or not
;
; CALLING SEQUENCE:
;   
;   dat = x_readspec(spec)
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
;  FSCALE      - Data is float
;  DSCALE      - Data is float
;  INFLG: 0  - Data, sig separate fits files  dat =
;              x_readspec(data_fil, FIL_SIG=sig_fil, SIG=sig, WAV=wav)
;         1  - Flux and Sig in one fits file
;         2  - Flux, sig, wave in one fits file
;         3  - FUSE format
;         4  - Flux in first argument, sigma in SIG (data arrays)
;         5  - SDSS format
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
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
                     FIL_SIG=fil_sig

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dat = x_readspec(spec, [dataunit], /fscale, /dscale, HEAD=) [V1.0]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( dataunit ) then dataunit = 0L
  if not keyword_set( inflg ) then inflg = 0 
  if keyword_set(DSCALE) and keyword_set(FSCALE) then fscale=0

  flg_fil = 0

; String or Image?
  if size(spec, /type) EQ 7 then begin
      a = findfile(spec, count=count)
      if count EQ 0 then begin
          print, spec+' does not exist'
          return, -1
      endif
      flg_fil = 1
  endif

  case inflg of
      0: begin  ; NORMAL -- data, sig separate
          dat = xmrdfits(spec, dataunit, head, fscale=fscale, $
                         dscale=dscale, /silent)
          npix = n_elements(dat)
          ;; SIG
          if keyword_set( FIL_SIG ) then sig = xmrdfits(fil_sig, /silent) $
            else sig = fltarr(npix)
          ;; Wavelengths
          if arg_present(WAV) then wav = x_setwave(head, npix) 
      end
      1: begin ; One fits file (flux, sig)
          stop
          ydat = x_readspec(yin, 0L, /fscale, head=head)
          npix = n_elements(ydat)
          if arg_present(SIG) then sig = x_readspec(yin, 1L, fscale=fscale)
          ;; Wavelengths
          if arg_present(WAV) then wav = x_setwave(head, npix) 
      end
      2: begin ; FLUX, SIG, WAVE (one fits file)
          dat = xmrdfits(spec, 0L, head, FSCALE=fscale, DSCALE=dscale)
          npix = n_elements(dat)
          if arg_present(SIG) then sig = xmrdfits(spec, 1L, /silent)
          if arg_present(WAV) then wav = xmrdfits(spec, 2L, /silent)
      end
      3: begin ; FUSE format
          fuse = mrdfits(spec, 1, /silent)
          if arg_present(WAV) then wav = fuse.wave
          dat = fuse.flux
          if arg_present(SIG) then sig = fuse.error
          npix = n_elements(dat)
          delvarx, fuse
      end
      4: begin ; 2 arrays fx, sig!
          if keyword_set( FIL_SIG ) then sig = fil_sig
          dat = spec
          npix = n_elements(dat)
          if arg_present(WAV) then wav = findgen(npix)
      end
      5: begin ; SDSS
          ;; Parse
          parse_sdss, spec, dat, wav, SIG=sig, HEAD=head, NPIX=npix
      end
          
      else: begin
          print, 'x_readspec: inflg = ', inflg, ' not allowed!'
          stop
          return, -1
      end
  endcase
  return, dat
end
