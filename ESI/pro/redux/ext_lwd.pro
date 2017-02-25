;+ 
; NAME:
; ext_lwd
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
;   extlwd_11sep04
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

pro ext_lwd, esi, obj_id, expsr, FNDOBJ=fndobj, EXTRCT=extrct,  $
             MKALL=mkall, PROCOBJ=procobj, FND_CLOB=fnd_clob, $
             SKYSUB=skysub, SKYCHK=skychk, SPECCHK=specchk, $
             CLOBBER=clobber, APER=aper, FLUX=flux, COMB=comb, $
             TRCCHK=trcchk, TRCOBJ=trcobj, FORCE=force, STD=std, $
             FILSTD=filstd

;;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'extlwdf_11sep04, esi, obj_id, /MKOBJ, /EXTRCT, /MKALL, /SKYSUB, '
      print, '     /TRCCHK, /TRCOBJ, /SPECCHK, /FORCE [v1.0]'
      return
  endif 

  ;; MKALL
  if keyword_set( MKALL ) then begin
      PROCOBJ=1
      FNDOBJ=1
      EXTRCT=1
      SKYSUB=1
      TRCOBJ=1
      FLUX=1
      COMB=1
  endif
  
  ;; Process Image
  if keyword_set( PROCOBJ ) then begin
      esi_lwdproc, esi, obj_id, CLOBBER=clobber
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Make Object structure
  if keyword_set( FNDOBJ ) then begin
      esi_lwdfndobj, esi, obj_id, expsr, CLOBBER=fnd_clob
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; SkySub
  if keyword_set( SKYSUB ) then begin
      esi_lwdskysub, esi, obj_id, CHK=skychk
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; TRACE
  if keyword_set( TRCOBJ ) then begin
      esi_lwdtrcobj, esi, obj_id, STD=std, CHK=trcchk, FILSTD=filstd
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Extract
  if keyword_set( EXTRCT ) then begin
      esi_lwdextobj, esi, obj_id, APER=aper
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Flux
  if keyword_set( FLUX ) then begin
      esi_lwdflux, esi, obj_id, FORCE=FORCE
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Combine
  if keyword_set( COMB ) then begin
      esi_lwdcomb, esi, obj_id
      ;; Write
      esi_wrstrct, esi, FITS='esi_11sep04.fits'
  endif

  ;; Check Spec
  if keyword_set( SPECCHK ) then begin
      esi_lwdpltobj, /fspec, /NOIMG
  endif

  return
end

