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

pro esi_mkzero, esi, NOSPEC=nospec, NOIMG=noimg, MMEM=mmem

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_mkzero, esi, /NOSPEC, /NOIMG, MMEM= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  

; SPEC First
  if not keyword_set( NOSPEC ) then begin
      print, 'esi_mkzero: Creating Spec Bias file!'
      ;;  Outfil
      outfil = 'Bias/BiasS.fits'
      a = findfile(outfil, count=na)
      if na NE 0 and not keyword_set( CLOBBER ) then $
        print, 'esi_mkzero: Outfil ', outfil, ' already exists!' $
      else begin
          ;; Grab all the ZRO frames
          zro = where(esi.type EQ 'ZRO' AND esi.mode GT 0 AND $
                      esi.flg_anly NE 0, nzro)
          if nzro EQ 0 then $
            print, 'esi_mkzero: No ZRO frames for SPEC! Should use an old one!!' $
          else begin
              ;; Combine
              if nzro GT 1 then $
                xcombine, esi[zro].rootpth+esi[zro].img_root, $
                img_zro, head, FCOMB=0 $
              else img_zro = xmrdfits(esi[zro].rootpth+esi[zro].img_root, $
                                     /fscale, /silent)
              ;; Output
              mwrfits, img_zro, outfil, /create, /silent
;              case esi[zro[0]].namp of
;                  1: mwrfits, img_zro[12:2059,*], outfil, /create, /silent
;                  2: mwrfits, img_zro[25:2072,*], outfil, /create, /silent
;                  else: stop
;              endcase
              spawn, 'gzip -f '+outfil
              print, 'esi_mkzro: Bias file created (Bias/BiasS.fits.gz)'
          endelse
      endelse
  endif

; IMAGE
  if not keyword_set( NOIMG ) then begin
      print, 'esi_mkzero: Creating IMG Bias file!'
      ;;  Outfil
      outfil = 'Bias/BiasI.fits'
      a = findfile(outfil, count=na)
      if na NE 0 and not keyword_set( CLOBBER ) then $
          print, 'esi_mkzero: Outfil ', outfil, ' already exists!' $
      else begin
          ;; Grab all the ZRO frames
          zro = where(esi.type EQ 'ZRO' AND esi.mode EQ 0 $
                      AND esi.flg_anly NE 0, nzro)
          if nzro EQ 0 then $
            print, 'esi_mkzero: No ZRO frames for IMG! Should use an old one!!' $
          else begin
              ;; Combine
              if nzro GT 1 then $
                xcombine, esi[zro].rootpth+esi[zro].img_root, img_zro, $
                head, FCOMB=0, MMEM=mmem $
              else img_zro = xmrdfits(esi[zro].rootpth+esi[zro].img_root, /fscale, $
                                     /silent)
              ;; Output
              if esi[zro[0]].namp NE 2 then stop
              mwrfits, img_zro[160:910,80:1425], outfil, /create, /silent
              print, 'esi_mkzro: Bias file created (Bias/BiasI.fits)'
          endelse
      endelse
  endif

  return
end

  
