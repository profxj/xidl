;+ 
; NAME:
; x_specbin
;   Version 1.0
;
; PURPOSE:
;  Simply bins data up in integer pixels    
;
; CALLING SEQUENCE:
;   
;   bin = x_specbin(fx, nbin)
;
; INPUTS:
;   fx       - Flux
;   nbin     - Number of pixels to bin on
;
; RETURNS:
;   bin       - Structure of data
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
;   bin = x_specbin(fx, 3)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_specbin, fx, nbin, VAR=var, WAV=wav

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'bindat = x_specbin(fx, nbin, VAR=, WAV=) [V1.0]'
    return, -1
  endif 


;  Optional Keywords

;  npix
  norg = n_elements(fx)
  npix = norg/nbin + (norg MOD nbin NE 0)
  lst = norg MOD nbin
  if lst EQ 0 then lst = nbin

; Make the arrays
  nwfx = fltarr(npix)
  nwvar = fltarr(npix)
  nwwav = dblarr(npix)

; Bin

  for qq=0L,npix-1 do begin
      ii = qq*nbin

      if qq NE npix-1 then begin
          ; Add em up
          nwfx[qq] = total(fx[ii:ii+nbin-1])/float(nbin)
          if keyword_set( VAR ) then $
            nwvar[qq] = total(var[ii:ii+nbin-1])/float(nbin)
          if keyword_set( WAV ) then $
            nwwav[qq] = total(wav[ii:ii+nbin-1])/float(nbin)
      endif else begin ; Last pixel
          nwfx[qq] = total(fx[ii:ii+lst-1])/float(lst)
          if keyword_set( VAR ) then $
            nwvar[qq] = total(var[ii:ii+lst-1])/float(lst)
          if keyword_set( WAV ) then $
            nwwav[qq] = total(wav[ii:ii+lst-1])/float(lst)
      endelse
  endfor
          
  return, { fx: nwfx, var: nwvar, wave: nwwav}

end
