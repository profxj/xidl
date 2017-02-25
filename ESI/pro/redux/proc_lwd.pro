;+ 
; NAME:
; proclwd_11sep04   
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  esi_lwdfndobj, esi, obj_id
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
;   proclwd_11sep04
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;   29-Aug-2002 Modified by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro proc_lwd, esi, slit, obj_id, MKBIAS=mkbias, MKFLAT=mkflat, $
                     SETUP=setup, INTER=inter, MMEM=mmem, $
                     MKAIMG=mkaimg, STDCHK=stdchk, TRCSTD=trcstd, CLOBBER=clobber,$
                     STD_NOFND=std_nofnd, ARC_CHK=arc_chk, FLT_CHK=flt_chk

;;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'proclwd_11sep04, esi, slit, [obj_id], /MKBIAS, /MKFLAT, /SETUP'
      print, '     /INTER, MMEM=, /MKAIMG, /STDCHK, /TRCSTD [v1.0]'
      return
  endif 

  ;; Make Bias
  if keyword_set( MKBIAS ) then begin
      esi_mkzero, esi, /NOIMG, MMEM=mmem
  endif

  ;; SETUP
  if keyword_set( SETUP ) then begin
      esi_lwdsetup, esi
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Make Arc Image
  if keyword_set( MKAIMG ) then begin
      if not keyword_set(SLIT) then begin
          print, 'proclwd_11sep04: Slit not set!!'
          return
      endif
      esi_lwdmkarc, esi, slit, CLOBBER=clobber
      esi_lwdfitarc, esi, slit, INTER=inter
      esi_lwdtrcarc, esi, slit
      esi_lwdmkaimg, esi, slit, CHK=arc_chk
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif
      
  ;; Process Flat
  if keyword_set( MKFLAT ) then begin
      esi_lwdmkflat, esi, slit, CHK=flt_chk, INTER=inter
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Process Standard
  if keyword_set( TRCSTD ) then begin
      esi_lwdtrcstd, esi, slit, CHK=stdchk, NOFND=std_nofnd
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  return
end

