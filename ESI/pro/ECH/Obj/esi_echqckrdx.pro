;+ 
; NAME:
; esi_echqckrdx   
;     Version 1.0
;
; PURPOSE:
;    Perform a quick extraction from the data (a la makee)
;
; CALLING SEQUENCE:
;   
;  esi_echqckrdx, obj_fil
;
; INPUTS:
;   obj_fil   -  Object file name
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echqckrdx, esi, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echqckrdx, obj_fil, ARCFIL=arcfil, FLATFIL=flatfil, ORDRS=ordrs, $
                   CHK=chk, NOPROC=noproc, REFORDR=refordr, REFWV=refwv, $
                   NOFND=nofnd, SKYCHK=skychk, NOSKY=nosky, TRCCHK=trcchk, $
                   NOTRC=notrc, NOEXT=noext

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echqckrdx, obj_fil, ARCFIL=, FLATFIL=, /CHK, /NOPROC [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( ARCFIL ) then arcfil = $
    getenv('ESI_CALIBS')+'/ArcECH_50IMG.fits'
  if not keyword_set( FLATFIL ) then flatfil = $
    getenv('ESI_CALIBS')+'/FlatECH50N.fits'
  if not keyword_set( BIASFIL ) then biasfil = $
    getenv('ESI_CALIBS')+'/BiasS.fits'
  if not keyword_set( SEDG_FIL ) then sedg_fil = $
    getenv('ESI_CALIBS')+'/SEdg_ECH50.fits'
  if not keyword_set( FITFIL ) then fitfil = $
    getenv('ESI_CALIBS')+'/hole_fit.idl'
  if not keyword_set( STDFIL ) then stdfil = $
    getenv('ESI_CALIBS')+'/ECH_STDOBJ.fits'
  if not keyword_set(ORDRS) then ordrs=[0L,9L]

; Create the structure
  print, systime()
  print, 'esi_echqckrdx: Creating the structure'
  esi_strct, esi, IMG=[obj_fil], /NOEDIT, /NOFILE, /MKDIR

; Set tags
  esi.obj_id = 1L
  esi.arc_fil = arcfil
  esi.obj_fil = 'Extract/Obj_'+esi[0].img_root

; Process the image
  if not keyword_set( NOPROC ) then begin
      esi_echproc, esi, [0], /clobber, FLATFIL=flatfil, BIASFIL=biasfil
      if keyword_set( CHK ) then begin
          xatv, 'Final/f_'+esi[0].img_root, /block, min=-10, max=100.
          stop
      endif
  endif

; Identify the object
  if not keyword_set( NOFND ) then begin
      esi_echfndobj, esi, 1L, REFORDR=refordr, REFWV=refwv, CHK=chk, $
        SEDG_FIL=sedg_fil, FITFIL=fitfil
  endif

; Sky subtract
  if not keyword_set( NOSKY ) then begin
      esi_echskysub, esi, 1L, CHK=skychk, SEDG_FIL=sedg_fil, FITFIL=fitfil,$
        ORDR=ordrs
  endif

; Trace
  if not keyword_set( NOTRC ) then begin
      esi_echtrcobj, esi, 1L, FAINT=faint, STDFIL=stdfil, CHK=trcchk, $
        SEDG_FIL=sedg_fil
  endif

; Extract
  if not keyword_set( NOEXT ) then begin
      esi_echextobj, esi, 1L, CHK=extchk, SEDG_FIL=sedg_fil, ORDRS=ordrs
      print, 'esi_echqckrdx: Type [esi_echspecplt] to view the output'
  endif

; SN

  obj = xmrdfits(esi[0].obj_fil,1)
  obj.sig = sqrt(obj.var > 0.)
  snfil = 's2n_'+strmid(esi[0].img_root,0,strlen(esi[0].img_root)-5)+'.ps'
  ps_open, file=snfil, font=1
  for i=0L,n_elements(obj)-1 do begin
      bad= where(obj[i].sig EQ 0.)
      obj[i].sig[bad] = 1e30
      plot, obj[i].wave, obj[i].fx / obj[i].sig, psym=10
  endfor
  ps_close, /noprint, /noid
  print, 'esi_echqckrdx: s2n ps file is ', snfil

;  DONE
  print, 'esi_echqckrdx: All done! '
  print, systime()
  return
end
