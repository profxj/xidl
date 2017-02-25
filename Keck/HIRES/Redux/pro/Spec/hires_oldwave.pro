;+ 
; NAME:
; hires_oldwave
;    Version 1.1
;
; PURPOSE:
;   Create a wavelength 'image' from the MAKEE header
;
; CALLING SEQUENCE:
;   wv = hires_oldwave(head)
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

function hires_oldwave, head

;
  if  N_params() LT 1 and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'wv = hires_oldwave(head) [v1.1]'
    return, -1
  endif 

  ;; Norder
  nordr = sxpar(head,'NAXIS2')
  if nordr EQ -1 then return, -1

  npix = sxpar(head,'NAXIS1')

  ;; Wavelength array
  wv = dblarr(npix,nordr)
  
  offset = sxpar(head,'CRVAL1')
  pix = dindgen(npix)+offset
  pixarr = pix # replicate(1.,7)

  for jj=1L,7 do pixarr[*,jj-1] = pixarr[*,jj-1]^jj

  ;; Loop
  for qq=1L,nordr do begin

      if qq LT 10 then pad = '0' else pad = ''

      ;; Cards
      crd_a = 'WV_0_'+pad+strtrim(qq,2)
      crd_b = 'WV_4_'+pad+strtrim(qq,2)

      ;; Read
      wv0s = sxpar(head,crd_a)
      wv4s = sxpar(head,crd_b)

      ;; Poly
      polyc = fltarr(8)
      polyc[0:3] = double(strsplit(wv0s,' ',/extract))
      polyc[4:7] = double(strsplit(wv4s,' ',/extract))

      polyarr = replicate(1.,npix) # polyc[1:7]

      ;; Wavelengths
      dum = total(pixarr * polyarr,2)
      wv[*,qq-1] = total(pixarr * polyarr,2) + polyc[0]
  endfor
  

  return,wv
end
  

