;+ 
; NAME:
; proc_ech
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  proc_ech, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   proc_ech
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;   29-Aug-2002 Modified by JXP
;-
;------------------------------------------------------------------------------

pro proc_ech, esi, slit, MKALL=mkall, MKSTRCT=mkstrct, $
              MKZERO=mkzero, SETUP=setup, CLOBBER=clobber, $
              MKMAP=mkmap, MKFLAT=mkflat, IFLAT=iflat, CHK=chk, $
              NRMFLAT=nrmflat,FLATCHK=flatchk, ARCCHK=arcchk, $
              MKARC=mkarc, FARC_CHK=farc_chk, PROCSTD=procstd, $
              INTER=inter, STDCHK=stdchk, CLOBAIMG=clobaimg, PINTER=pinter,$
              FITARC=fitarc, TRCARC=trcarc, CUAR=cuar, GUESSARC=guessarc,$
              LINLIST=linlist, PIXSHFT=pixshft

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'proc_ech, esi, [slit], /IFLAT, /REDOOV, /MKALL, /CLOBBER, /MKMAP [v1.0]'
      print, '/CHK, /NRMFLAT, /FLATCHK, /ARCCHK, /MKARC, /FARC_CHK, /PROCSTD'
      print, '/INTER, /STDCHK, /CLOBAIMG, /PINTER, /FITARC, /TRCARC, /CUAR'
      print, 'GUESSARC=, LINLIST=, PIXSHFT= '
      return
  endif 

  ; MKALL
  if keyword_set( MKALL ) then begin
      MKFLAT = 1
      MKMAP = 1
      MKSLITS = 1
      MKAIMG = 1
      MKNFLAT = 1
      MKOBJ = 1
      FLGCR = 1
      CLOBBER = 1
      clobaimg = 1
  endif

  ;;;;;;;;;;;;;;
  ; Make the structure
  if keyword_set( MKSTRCT ) then begin
      esi_strct, esi, /mkdir, NOEDIT=noedit
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif
                 
  ;;;;;;;;;;;;;;
  ; Make the zero image
  if keyword_set( MKZERO ) then esi_mkzero, esi, /NOIMG  ;; 5 min

  ;;;;; SETUP ;;;;;;;;
  if keyword_set( SETUP ) then esi_echsetup, esi

  ;;;;; PINHOLE ;;;;;;;;
  if keyword_set( MKMAP ) then begin
      ;; 5-10min
      esi_echtrcholes, esi
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
      esi_echmkmap, esi
  endif
  ;; OR COPY ECH_map.fits.gz, img_hole.fits.gz, hole_fit.idl 

  ;;;;; FLAT ;;;;;;;;;;
  if keyword_set( MKFLAT ) then begin
      ;; 30s/image
      esi_echmkflat, esi, slit, IFLAT=iflat
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  if keyword_set( NRMFLAT ) then begin
      ;; 1min
      if not keyword_set( SLIT ) then begin
          print, 'proc_ech: Need to specifiy slit size for this step!'
          return
      endif
      esi_echfltsct, esi, slit, IFLAT=iflat, CHK=flatchk
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;  ARC  ;;;;;;;;;;
  if keyword_set( MKARC ) then begin
      ;; <2min
      if not keyword_set( SLIT ) then return
      esi_echmkarc, esi, slit, CLOBBER=clobber
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif
  if keyword_set( FITARC ) then begin
      ;; 2min
      esi_echfitarc, esi, slit, CHK=FARC_CHK, INTER=inter, PINTER=pinter, $
        CUAR=cuar, CLOBBER=clobber, GUESSARC=guessarc, LINLIST=linlist
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif
  if keyword_set( TRCARC ) then begin
      ;; 10min :: Human interaction
      if keyword_set( INTER) then AUTO=0 else AUTO=1
      esi_echtrcarc, esi, slit, AUTO=auto, CUAR=cuar, LINLIST=linlist, $
        GUESSARC=guessarc, PIXSHFT=pixshft
      ;; 5min
      esi_echmkaimg, esi, slit, CHK=arcchk, CUAR=cuar
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;; STD ;;;;;;

  if keyword_set( PROCSTD ) then begin
      ;; 10min (sky sub)
      esi_echtrcstd, esi, slit, CHK=stdchk, CUAR=cuar
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  return
end
