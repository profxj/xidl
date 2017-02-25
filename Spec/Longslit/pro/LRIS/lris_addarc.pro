;+
; NAME:
;   lris_addarc
;
; PURPOSE:
;   Combines multiple arcs used when one has to turn 
;   certain lamps on/off  
;
; CALLING SEQUENCE:
; lris_addarc, fil1, fil2, outfil
;
; INPUTS:
;  fil1 -- Filename of Raw image1
;  fil2 -- Filename of Raw image2
;
; OUTPUTS:
;  outfil -- Output filename for the combined arc
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;   18-May-2009  Modified by JXP
;-
;------------------------------------------------------------------------------

pro lris_addarc, fils, outfil, PROC_FIl=proc_fil

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'lris_addarc, fils, outfil [v1.0]' 
      return
  endif 


  raw1 = xmrdfits(fils[0],0,hdr_raw)
  nfil = n_elements(fils)

  sz_raw = size(raw1, /dimens)
  for kk=0L,nfil-1 do begin
      long_proc, fils[kk], arc, arcivar, hdr = head, gain = gain
      if kk EQ 0 then imag = arc else imag = imag + arc
      if kk EQ 0 then ivar = arcivar else ivar = ivar + arcivar
  endfor
  
  if not keyword_set(PROC_FIL) then begin
      sz_proc = size(arc, /dimen)
      arc_fake = fltarr(sz_raw)
      
      ;; Red side
      arc_fake[0:sz_proc[1]/2 -1, *] = $
        transpose(imag[*,0:sz_proc[1]/2.-1])/double(gain[0])
      arc_fake[sz_proc[1]/2:sz_proc[1]-1, *] = $
        transpose(imag[*,sz_proc[1]/2:sz_proc[1]-1])/double(gain[1])
      
      mwrfits, arc_fake, outfil, hdr_raw, /create
      spawn, 'gzip -f '+outfil
  endif else begin
      ;; Header
      szimg = size(imag, /dimen)
      sxaddpar, head, 'NAXIS1', szimg[0]
      sxaddpar, head, 'NAXIS2', szimg[1]
      arcimg = imag
      arcivar = ivar
      hdr = head
;      save, arcimg, arcivar, hdr, filena=proc_fil
;      stop ;; FIX THIS SOMEDAY
      mwrfits, arcimg, PROC_FIL, /create
      mwrfits, arcivar, PROC_FIL, hdr
  endelse

END
