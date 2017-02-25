;+ 
; NAME:
; esi_mkzero   
;     Version 1.1
;
; PURPOSE:
;    Create bias frames from zero exposures
;    Default is to create both IMG and SPEC bias frames
;
; CALLING SEQUENCE:
;   
;  esi_mkzero, esi, /NOSPEC, /NOIMG, MMEM=
;
; INPUTS:
;   esi   -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;   Creates bias frames in 'Bias/BiasS.fits' and 'Bias/BiasI.fits'
;    for spectroscopy and imaging respectively
;
; OPTIONAL KEYWORDS:
;    /NOSPEC -- Dont create spectroscopy bias
;    /NOIMG -- Dont create image bias
;    MMEM= -- Maximum memory to use in combining zero flats
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_mkzero, esi, /NOIMG
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_mkzero, esi, NOSPEC=nospec, NOIMG=noimg, MMEM=mmem, CLOBBER = CLOBBER

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_mkzero, esi, /NOSPEC, /NOIMG, MMEM= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  ;; Find all binnings
  gd = where(esi.flg_anly NE 0)
  bintot = esi[gd].cbin + 10*esi[gd].rbin
  bintyp = bintot[uniq(bintot, sort(bintot))]
  nbin = n_elements(bintyp)

  for qq=0L,nbin-1 do begin
      row = bintyp[qq]/10
      col = bintyp[qq] - row*10

      ;; SPEC First
      if not keyword_set( NOSPEC ) then begin
          print, 'esi_mkzero: Creating Spec Bias file!'
          ;;  Outfil
          outfil = esi_getfil('bias_fil', 1, cbin=col, rbin=row, /name)
          a = findfile(outfil, count=na)
          if na NE 0 and not keyword_set( CLOBBER ) then $
            print, 'esi_mkzero: Outfil ', outfil, ' already exists!' $
          else begin
              ;; Grab all the ZRO frames
              zro = where(esi.type EQ 'ZRO' AND esi.mode GT 0 AND $
                          esi.rbin EQ row AND esi.cbin EQ col AND $
                          esi.flg_anly NE 0, nzro)
              if nzro EQ 0 then begin
                  print, 'esi_mkzero: No ZRO frames for SPEC with binning ', $
                    col, row, ' Should use an old one!!' 
                  stop
              endif else begin
              ;; Combine
                  if nzro GT 1 then $
                    xcombine, esi[zro].rootpth+esi[zro].img_root, $
                    img_zro, head, FCOMB=0 $
                  else img_zro = xmrdfits(esi[zro].rootpth+esi[zro].img_root, $
                                          /fscale, /silent)
                  ;; Output
                  mwrfits, img_zro, outfil, /create, /silent
                  spawn, 'gzip -f '+outfil
                  print, 'esi_mkzro: Bias file created ', outfil
              endelse
          endelse
      endif

      ;; IMAGE
      if not keyword_set( NOIMG ) then begin
          flg_img = 0
          ;; Check binning
          if col NE 1 or row NE 1 then begin
              print, 'esi_mkzero:  Skipping IMG for non 1x1 imaging'
              flg_img = 1
          endif
          if flg_img EQ 1 then break
          print, 'esi_mkzero: Creating IMG Bias file!'
          ;;  Outfil
          outfil = esi_getfil('bias_fil', 0, cbin=col, rbin=row, /name)
          a = findfile(outfil, count=na)
          if na NE 0 and not keyword_set( CLOBBER ) then $
            print, 'esi_mkzero: Outfil ', outfil, ' already exists!' $
          else begin
              ;; Grab all the ZRO frames
              zro = where(esi.type EQ 'ZRO' AND esi.mode EQ 0 $
                          AND esi.rbin EQ row AND esi.cbin EQ col $
                          AND esi.flg_anly NE 0, nzro)
              if nzro EQ 0 then begin
                  print, 'esi_mkzero: No ZRO frames for IMG with binning: ', col,$
                    row, '! Should use an old one!!' 
                  stop
              endif else begin
                  ;; Combine
                  if nzro GT 1 then $
                    xcombine, esi[zro].rootpth+esi[zro].img_root, img_zro, $
                    head, FCOMB=0, MMEM=mmem $
                  else img_zro = xmrdfits(esi[zro].rootpth+esi[zro].img_root, $
                                          /fscale, /silent)
                  ;; Output
                  if esi[zro[0]].namp NE 2 then stop
                  mwrfits, img_zro[160:910,80:1425], outfil, /create, /silent
                  print, 'esi_mkzro: Bias file created ', outfil
              endelse
          endelse
      endif
  endfor

  return
end

  
