;+ 
; NAME:
; extrct_ech   
;     Version 1.2
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  extrct_ech, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   obj_id  -   Object ID  (e.g. 0L, 1L, etc)
;   [exp_id]  -  Exposure frames (e.g. [0L, 1L])
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /STANDARD - use this if processing a standard here
;   /OPTIMAL - use if analyzing point source and want optimal
;              extraction rather than boxcar extraction
;   /OLD - Use this to use archived standard solution
;   
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   extract_ech, esi, 1
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   06-Sep-2002 Written by JXP
;   11-Mar-2008 Updated by MAR
;-
;------------------------------------------------------------------------------

pro extrct_ech, esi, obj_id, exp_id, MKALL=mkall, CLOBBER=clobber, $
                FAINT=faint, STPTRC=stptrc, PSTSKY=pstsky, MINOFF=minoff, $
                PROCOBJ=procobj, TRCOBJ=trcobj, FNDOBJ=fndobj, $
                SUBSCAT=subscat, SCICLM=sciclm, FOBJCLOB=fobjclob, $
                EXTOBJ=extobj, SKYSUB=skysub, CHK=chk, TRCCHK=trcchk, $
                INTER=inter, COMBOBJ=combobj, EXTCHK=extchk, FND_CHK=fnd_chk,$
                REFORDR=refordr, REFWV=refwv, OBJCR=objcr, CUAR=cuar, $
                FLUX=flux, COADD=coadd, PSKYSUB=pskysub, BORDR=bordr,$
                ORDR=ordr, SKY_FCHK=sky_fchk, USESTD=usestd, STDFIL=stdfil, $
                OPTIMAL=optimal, SEXTCHK=sextchk, OLD=old, flxfil=flxfil, standard=standard

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'extract_ech, esi, obj_id, exp_id, /MKALL, /CLOBBER, /FAINT, /SUBSCAT '
      print, '    /PROCOBJ, /TRCOBJ, /FNDOBJ, /SUBSCAT, /EXTOBJ, /SKYSUB'
      print, '    /INTER, /TRCCHK, /COMBOBJ, /EXTCHK, REFORDR=, REFWV='
      print, '    /OBJCR, /FLUX, /COADD, /PSKYSUB, BORDR=, /CUAR, ORDR='
      print, '    /old, /STANDARD [v1.2]'
      return
  endif 

  ; MKALL
  if keyword_set( MKALL ) then begin
      PROCOBJ=1
      FNDOBJ=1
      SKYSUB=1
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

  if keyword_set(standard) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'STD', nindx)
  endif

  if nindx EQ 0 then stop


  ;;;;;;;;;;;;;;
  ; Process
  if keyword_set( PROCOBJ ) then begin
      esi_echproc, esi, indx, CLOBBER=clobber, SUBSCAT=subscat
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; OBJ CR
  if keyword_set( OBJCR ) then begin
      esi_echobjcr, esi, obj_id, CHK=chk
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Find OBJ
  if keyword_set( FNDOBJ ) then begin
      esi_echfndobj, esi, obj_id, exp_id, NOCLOB=fobjclob, CHK=fnd_chk, std=standard
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif
                 
  ;;;;;;;;;;;;;;
  ; SkySub OBJ (STD)
  if keyword_set( SKYSUB ) then begin
      if keyword_set( CUAR ) then begin
          ordr = [0L, 8L]
          clobber = 1L
      endif
      esi_echskysub, esi, obj_id, exp_id, CHK=chk, BORDR=bordr, ORDR=ordr, $
        FCHK=sky_fchk, std=standard
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; SkySub OBJ (POLY)
  if keyword_set( PSKYSUB ) then begin
      esi_echpskysub, esi, obj_id, exp_id, CHK=chk
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif
                 
  ;;;;;;;;;;;;;;
  ; Trace OBJ
  ; This function is obsolete - here for backward compatibility only
  if keyword_set( TRCOBJ ) then begin
      esi_echtrcobj, esi, obj_id, exp_id, FAINT=faint, USESTD=usestd, CHK=trcchk,$
        STDFIL=stdfil
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Extract OBJ
  if keyword_set( EXTOBJ ) then begin
      esi_echextobj, esi, obj_id, exp_id, CHK=extchk, OPTIMAL=optimal, $
        SPCCHK=sextchk, std=standard
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; Combine OBJ
  if keyword_set( COMBOBJ ) then begin
      esi_echcombspec, esi, obj_id
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; FLUX
  if keyword_set( FLUX ) then begin
      esi_echfluxfin, esi, obj_id, CLOBBER=clobber, old=old, fluxfil=flxfil
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

  ;;;;;;;;;;;;;;
  ; COADD
  if keyword_set( COADD ) then begin
      esi_echcoaddfin, esi, obj_id
      esi_wrstrct, esi, FITS='esistrct.fits', /ANONLY
  endif

;      esi_echspecplt, /fspec

  return
end
