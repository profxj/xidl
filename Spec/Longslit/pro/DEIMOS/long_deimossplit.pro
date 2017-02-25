;+
; NAME:
;   long_deimossplit
;
; PURPOSE:
;   Split a DEIMOS FITS file into the chip of interest that covers
;   the longslit footprint
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
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
;   2006  Written by JXP 
;-
;------------------------------------------------------------------------------

pro long_deimossplit, root, chip

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'long_deimossplit, root, chip  [v1.0]'
      return
  endif 

  if not keyword_set(CHIP) then chip = [3L, 7L]

  ;; Frames
;  if keyword_set(all) then begin
;     files = findfile(ROOT+'*')
;     nfrm = n_elements(files)
;     frame = lonarr(nfrm)
;     for ii=0L,nfrm-1 do begin
;        ipos = strpos(files[ii], '.fits')
;        frame[ii] = strtrim(strmid(files[ii],ipos-4,4))
;     endfor
;  endif else if not keyword_set(FRAME) then return
  nfrm = n_elements(chip)

  fil = findfile(root+'*', count=nfil)
  if nfil EQ 0 then return

  for qq=0L,nfil-1 do begin
      for kk=0L,nfrm-1 do begin
          print, 'Reading ', fil[qq], chip[kk]
          ;; Read
          idat = xmrdfits(fil[qq],chip[kk],headi,/sile)
          if chip[kk] GT 4 then dat = rotate(idat,7) $ ;; Flipping to keep blue -> red
          else dat = idat
          sz = size(dat,/dime)

          datsec = sxpar(headi,'DATASEC')
          ;; More header
          head = xheadfits(fil[qq])
;          head = [head,head2]
;          gd = where(strlen(strtrim(head,2)) GT 0 and strlen(strtrim(head,2)) LT 80)
;          head =  head[gd]

          if chip[kk] LE 4 then instrum = 'DEIMOSb' else instrum='DEIMOSr'
          
          ;; Header
          sxaddpar, head, 'INSTRUME', instrum
          sxaddpar, head, 'SIMPLE', 'T'
          sxaddpar, head, 'NAXIS', 2
          sxaddpar, head, 'NAXIS1', sz[0]
          sxaddpar, head, 'NAXIS2', sz[1]  ;; Need these up top
          nhead = n_elements(head)
          newhead = [head[0:2], head[nhead-3:nhead-2], head[3:nhead-4],head[nhead-1]]
          sxaddpar, newhead, 'DATASEC', datsec
          sxdelpar, newhead, 'DATASUM'
          sxdelpar, newhead, 'CHECKSUM'

          ;; Outfil
          pos = strpos(fil[qq], '.fits')
          outfil = strmid(fil[qq],0,pos)+'_'+strtrim(chip[kk],2)+'.fits'
      
          ;; Write
          mwrfits, dat, outfil, newhead, /create
          spawn, 'gzip -f '+outfil

          print, 'long_deimssplit: Wrote ', outfil
      endfor
  endfor
  
  return
END
