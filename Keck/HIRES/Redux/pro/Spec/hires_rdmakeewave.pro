;+ 
; NAME:
; hires_rdmakeewave
;    Version 1.1
;
; PURPOSE:
;   Parses the rather painful MAKEE header format (IRAF?)
;
; CALLING SEQUENCE:
;
; INPUTS:
;   header
;
; RETURNS:
;
; OUTPUTS:
;  wave_array  -- Array of wavelengths, order by order   
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hires_combmakee, ['Flux-079.fits', 'Flux-080.fits']
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Sep-2011 Written by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hires_rdmakeewave, header, wave_array, OFFSET=offset

  npix = sxpar(header, 'NAXIS1')
  nordr = sxpar(header, 'NAXIS2')
  crval1 = sxpar(header, 'CRVAL1')

  wave_array = dblarr(npix, nordr)

  if n_elements(offset) EQ 0 then offset = 1.

  ;; Loop
  ncoeff = 8L
  for qq=0l,nordr-1 do begin
     ;; Parse
     c1 = 'WV_0_'+x_padstr(strtrim(qq+1,2), 2, '0', /rev)
     w1 = sxpar(header, c1)
     wvp = double( strsplit(w1, ' ', /extract) )

     c2 = 'WV_4_'+x_padstr(strtrim(qq+1,2), 2, '0', /rev)
     w2 = sxpar(header, c2)
     wvp = [wvp, double( strsplit(w2, ' ', /extract) )]

     ;; Calculate wavelengths
     xv = dindgen(npix) + CRVAL1 + OFFSET
     wvo = dblarr(npix)
     for ii=0L,ncoeff-1 do $
        wvo = wvo + wvp[ii] * (xv^ii)

     ;; Save
     wave_array[*,qq] = wvo
  endfor
        
  return
end
  

