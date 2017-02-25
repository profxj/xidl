;+ 
; NAME:
; x_nrmspec
;   Version 1.1
;
; PURPOSE:
;  Applies the continuum to unnormalized data
;
; CALLING SEQUENCE:
;   x_nrmspec, fluxfil, contfil, outfil, EFIL=, OUTE=
;
; INPUTS:
;  fluxfil  -- Name of fluxed spectrum (or non-normalized)
;  contfil  -- Name of file containing the continuum
;  outfil   -- Name of output file to contain normalized spectrum
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  INFLG= -- Keyword for x_readspec for reading the spectrum [default:
;            0]
;  EFIL=  -- Filename of sigma array
;  OUTE=  -- Filename for normalized sigma array
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
;   Fall-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_nrmspec, ffil, cfil, outfil, EFIL=efil, OUTE=oute, INFLG=inflg

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_nrmspec, ffil, cfil, ofil, EFIL=, OUTE=, INFLG= [V1.1]'
    return
  endif 
  if not keyword_set(INFLG) then inflg = 0

  conti = xmrdfits(cfil, /silent)
  if INFLG LE 1 then begin
      ;; Read 1
      fx = xmrdfits(ffil, 0, head, /silent)
  endif else begin
      fx = x_readspec(ffil,sig=sig,wav=wav,inflg=inflg)
  endelse

  npx = n_elements(fx)

  ;; Norm
  gd = where(conti NE 0.)
  newfx = fltarr(npx)
  newfx[gd] = fx[gd] / conti[gd]

  ;; Error
  if keyword_set(EFIL) then sig = xmrdfits(efil, /silent)
  if keyword_set(SIG) then begin
      newsig = fltarr(npx)
      newsig[gd] = sig[gd] / conti[gd]
  endif

  ;; Output
  if INFLG LE 1 then begin
      mwrfits, newfx, outfil, head, /create
      if keyword_set(OUTE) then mwrfits, newsig, oute, head, /create
  endif else begin
      case inflg of
          2: begin
              mwrfits, newfx, outfil, head, /create
              mwrfits, newsig, outfil
              mwrfits, wav, outfil
           end
          8: begin
             forprint, wav, newfx, 1/newsig^2, textout=outfil
          end
          else: stop
      endcase
  endelse

  return

end
