;+ 
; NAME:
; extract_ech   
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  procech_07apr00, esi, obj_id
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
;   extract_ech
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro extract_ech, esi, obj_id, exp_id, MKALL=mkall, CLOBBER=clobber, $
                 FAINT=faint, STPTRC=stptrc, PSTSKY=pstsky, MINOFF=minoff, $
                 PROCOBJ=procobj, TRCOBJ=trcobj, FNDOBJ=fndobj, $
                 SUBSCAT=subscat, SCICLM=sciclm, FOBJCLOB=fobjclob, $
                 EXTOBJ=extobj, SKYSUB=skysub, CHK=chk, TRCCHK=trcchk, $
                 INTER=inter, COMBOBJ=combobj, EXTCHK=extchk, FND_CHK=fnd_chk,$
                 REFORDR=refordr, REFWV=refwv, OBJCR=objcr, CUAR=cuar, $
                 FLUX=flux, COADD=coadd, PSKYSUB=pskysub, BORDR=bordr,$
                 ORDR=ordr

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'extract_ech, esi, obj_id, exp_id, /MKALL, /CLOBBER, /FAINT, /SUBSCAT '
      print, '    /PROCOBJ, /TRCOBJ, /FNDOBJ, /SUBSCAT, /EXTOBJ, /SKYSUB'
      print, '    /INTER, /TRCCHK, /COMBOBJ, /EXTCHK, REFORDR=, REFWV='
      print, '    /OBJCR, /FLUX, /COADD, /PSKYSUB, BORDR=, /CUAR, ORDR='
      print, '    [v1.1]'
      return
  endif 

  ; MKALL
  if keyword_set( MKALL ) then begin
      PROCOBJ=1
      FNDOBJ=1
      SKYSUB=1
      TRCOBJ=1
      EXTOBJ=1
      COMBOBJ=1
      FLUX=1
      COADD=1
  endif

  ; Stop at TRACE
  if keyword_set( STPTRC ) then begin
      PROCOBJ=1
      FNDOBJ=1
      SKYSUB=1
  endif
  
  if keyword_set( PSTSKY ) then begin
      TRCOBJ=1
      EXTOBJ=1
      COMBOBJ=1
      FLUX=1
      COADD=1
  endif


  ;; INDX
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)

  ;;;;;;;;;;;;;;
  ; Process
  if keyword_set( PROCOBJ ) then begin
      esi_echproc, esi, indx, CLOBBER=clobber, SUBSCAT=subscat
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; OBJ CR
  if keyword_set( OBJCR ) then begin
      esi_echobjcr, esi, obj_id
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Find OBJ
  if keyword_set( FNDOBJ ) then begin
      esi_echfndobj, esi, obj_id, exp_id, NOCLOB=fobjclob, INTER=inter, $
        REFORDR=refordr, REFWV=refwv, SCICLM=sciclm, CHK=fnd_chk, MAXOFF=maxoff
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif
                 
  ;;;;;;;;;;;;;;
  ; SkySub OBJ (STD)
  if keyword_set( SKYSUB ) then begin
      if keyword_set( CUAR ) then begin
          ordr = [0L, 8L]
          clobber = 1L
      endif
      esi_echskysub, esi, obj_id, exp_id, CHK=chk, BORDR=bordr, ORDR=ordr, $
        CLOBBER=clobber
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; SkySub OBJ (POLY)
  if keyword_set( PSKYSUB ) then begin
      esi_echpskysub, esi, obj_id, exp_id, CHK=chk
      esi_wrstrct, esi, FITS='esi_13jan02.fits', /ANONLY
  endif
                 
  ;;;;;;;;;;;;;;
  ; Trace OBJ
  if keyword_set( TRCOBJ ) then begin
      esi_echtrcobj, esi, obj_id, exp_id, FAINT=faint, /USESTD, CHK=trcchk
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Extract OBJ
  if keyword_set( EXTOBJ ) then begin
      esi_echextobj, esi, obj_id, exp_id, CHK=extchk
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Combine OBJ
  if keyword_set( COMBOBJ ) then begin
      esi_echcombspec, esi, obj_id
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; FLUX
  if keyword_set( FLUX ) then begin
      esi_echfluxfin, esi, obj_id, CLOBBER=clobber
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; COADD
  if keyword_set( COADD ) then begin
      esi_echcoaddfin, esi, obj_id
      esi_wrstrct, esi, FITS='esi_n1.fits', /ANONLY
  endif

;      esi_echspecplt, /fspec

  return
end
