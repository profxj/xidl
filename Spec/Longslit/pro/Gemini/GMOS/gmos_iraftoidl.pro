;+
; NAME:
;   gmos_iraftoidl
;
; PURPOSE:
;  Convert GMOS gs files to IDL type
;  Also produces a variance image
;
; CALLING SEQUENCE:
;   gmos_iraftoidl, 'IRaw/gsN*fits.gz'
;
; INPUTS:
;   infil  -- Filname or list of files
;
; OPTIONAL INPUTS:
;  DAT_EXTEN=  -- Extension of the raw data frame [Default = 2]
;
; OPTIONAL KEYWORDS:
;  /GSKLUDGE -- Some GMOS-S data has a bad illuimation pattern. This
;               sets the variance to zero in that region.  [Not
;               implemented]
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   15-Jun-2009  Written by JXP
;------------------------------------------------------------------------------
; gmos_iraftoidl, 'IRaw/gsN*.fits.gz'
pro gmos_iraftoidl, infil, rawdir, CLOBBER=clobber, DAT_EXTEN=dat_exten;, GSKLUDGE=gskludge

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'gmos_iraftoidl, infil, rawdir (v1.0)'
    return
  endif 

  if not keyword_set(DAT_EXTEN) then dat_exten = 2

  files = findfile(infil, count=nfil)
  if nfil EQ 0 then begin
      print, 'gmos_iraftoidl: No files!', infil
      return
  endif

  cards = ['TELESCOP', 'RDNOISE', 'CCDSUM', 'INSTRUME', 'OBJECT', $
           'EPOCH', 'RA', 'GRTILT', 'EQUINOX', 'ELAPSED', 'OBSERVAT', $
           'DEC', 'AIRMASS', 'UT', 'DATE', 'GAIN', 'OBSMODE', $
           'OBSTYPE', 'GRATING', 'GRWLEN', 'EXPTIME', $
           'UTSTART', 'UTEND', 'EXPOSURE']
  ncrds = n_elements(cards)

  for qq=0L,nfil-1 do begin
      pos = strpos(files[qq],'/',/reverse_sea)
      posf = strpos(files[qq],'.fits',/reverse_sea)

      ;; Outfil
      if pos GT 0 then root = strmid(files[qq],0,pos) else root = ''
      outfil = root+'/idl_'+strmid(files[qq],pos+1,(posf-pos-1))+'.fits'
      b = findfile(outfil+'*',count=nb)
      if nb NE 0 and not keyword_set(CLOBBER) then begin
          print, 'gmos_iraftoidl: Skipping ', outfil
          continue
      endif

          
      ;; Read
      fullhead = xheadfits(files[qq])

      ;; Rawfil
      rawfil = rawdir+strmid(files[qq],pos+3,(posf-pos-3))+'.fits'
      rawhead = xheadfits(rawfil, exten=2)
      dat = xmrdfits(files[qq], DAT_EXTEN, head)
      dat = transpose(dat)
      dat = rotate(dat,7)
;  xatv, dat, /bloc
      sz = size(dat,/dimen)
      
      ;; Create variance image
      rnoise = sxpar(head, 'RDNOISE')
      
      invv = 1.0/(abs(dat-sqrt(2.0)*rnoise) + rnoise^2)
      
      ;; Mask gaps
      cbin = strtrim(strsplit(sxpar(head,'CCDSUM'),/extract),2)
      bin = reverse(long(cbin))
      specbin = bin[1]
      IF specbin EQ 4 THEN BEGIN
          ygap1 = 9L
          ygap2 = 9L
      ENDIF ELSE IF specbin EQ 2 THEN BEGIN
          ygap1 = 18L
          ygap2 = 18L
      ENDIF ELSE IF specbin EQ 1 THEN BEGIN
          ygap1 = 36L
          ygap2 = 36L
      ENDIF 
      rawspec = (sz[1]-(ygap1+ygap2))/3
      invv[*,rawspec-1:rawspec+ygap1+1] = -1.
      dat[*,rawspec-1:rawspec+ygap1+1] = 0.
      invv[*,2*rawspec+ygap1-1:2*rawspec+ygap1+ygap2+1] = -1.
      dat[*,2*rawspec+ygap1-1:2*rawspec+ygap1+ygap2+1] = 0.
      bad = where(abs(dat-1.) LT 1e-3, nbad)
      if nbad GT 0 then invv[bad] = -1.

;      if keyword_set(GSKLUDGE) then begin
;          invv[1172:1334,*] = -1
;      endif
      
      ;; Write
      mkhdr, hdr, dat, /extend
      for jj=0L,ncrds-1 do $
             sxaddpar, hdr, cards[jj], sxpar(fullhead,cards[jj])
      sxaddpar, hdr, 'MJD-OBS', sxpar(rawhead,'MJD-OBS')
      sxaddpar, hdr, 'XGMOS', 1
      
      mwrfits, dat, outfil, hdr, /create
      mwrfits, invv, outfil
      spawn, 'gzip -f '+outfil
  endfor

  return
end
