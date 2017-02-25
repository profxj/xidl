;+ 
; NAME:
; hires_1dcoadd
;    Version 1.1
;
; PURPOSE:
;   Combines multiple 1D spectra from HIREDUX
;      Assumes the data were observed with the same slit
;
; CALLING SEQUENCE:
;   hires_1dcoadd, files, outfil
;
; INPUTS:
;   files - Array of normalized flux files
;
; RETURNS:
;
; OUTPUTS:
;   outfil   -  Name for the output file
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Oct-2014 Written by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hires_1dcoadd, files, outfil, CHECK=check

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'hires_1dcoadd, files, outfil'
    print, '       [v1.0]'
    return
  endif 

  ;; 
  nfil = n_elements(files)
  if nfil LE 1 then return

  svhdr = xheadfits(files[0])

  ;; Generate
  maxpix = 200000L
  influx = fltarr(maxpix,nfil)
  inivar = fltarr(maxpix,nfil)
  loglam = dblarr(maxpix,nfil)

  ;; Load them up
  for kk=0L,nfil-1 do begin
     ;; Read
     dat = x_readspec(files[kk], /auto, /struct,head=head)
     cdelt1=sxpar(head,'CDELT1')
     if kk GT 0 then begin
        if abs(cdelt1-sv_cdelt1) GT 1e-8 then stop ;; Pixel size needs to be the same!
     endif else sv_cdelt1 = cdelt1
     ;; 
     influx[0:dat.npix-1, kk] = dat.fx
     gdp = where(dat.sig GT 0.)
     inivar[gdp, kk] = 1./dat.sig[gdp]^2
     loglam[0:dat.npix-1, kk] = alog10(dat.wv)
     ;; 
     if kk EQ 0 then newloglam = alog10(dat.wv) else begin
        tmp = alog10(dat.wv)
        mnwv = min(newloglam, max=mxwv)
        newwv = where(tmp LT (mnwv-1d-7) OR $ 
                      tmp GT (mxwv+1d-7), ntmp)
        if ntmp GT 0 then begin
           newloglam = [newloglam,tmp[newwv]]
        endif
     endelse
  endfor

  ;; Sort
  srt  = sort(newloglam)
  newloglam = newloglam[srt]

  ;; Generate final file
  long_combspec, influx, inivar, loglam, NEWLOGLAM=newloglam, $
                 NEWFLUX=newflux, NEWIVAR=newivar, check=check, $
                 IN_NPOLY=0, FMIN=-5., fmax=5.
  npix = n_elements(newloglam)
  newsig = fltarr(npix)
  gdp = where( newivar GT 0.)
  newsig[gdp] = 1./ sqrt(newivar[gdp])

  ;; Header
  mkhdr, head, newflux
  sxaddpar, 'RA', sxpar(svhdr,'RA')
  sxaddpar, 'DEC', sxpar(svhdr,'DEC')
  sxaddpar, 'TARGNAME', sxpar(svhdr,'TARGNAME')

  ;; Write
  mwrfits, newflux, outfil, /create
  mwrfits, newsig, outfil
  mwrfits, 10.d^newloglam, outfil

  print, 'hires_combspec:  Output is in ', outfil
  print, 'hires_combspec:  All done!'


  return
end
  

