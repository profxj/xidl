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
;  NFIND=   ;; Number of objects to extract (brightest first)
;  [default: 1]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echqckrdx, 'e140423_0064.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echqckrdx, obj_fil, arcs, flats, $
                   ARCFIL=arcfil, FLATFIL=flatfil, ORDRS=ordrs, OLD=old, $
                   CHK=chk, NOPROC=noproc, REFORDR=refordr, REFWV=refwv, $
                   TWOFND=twofnd, IFLAT=iflat, NOOSCAN=nooscan, $
                   NOFND=nofnd, SKYCHK=skychk, NOSKY=nosky, TRCCHK=trcchk, $
                   NOTRC=notrc, NOEXT=noext, SCICLM=sciclm, NOFLAT=noflat, APERV=aperv, $
                   FINTER=finter, NOARC=noarc, NFIND=nfind, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echqckrdx, obj_fil, ARCFIL=, FLATFIL=, SCICLM=, /CHK, /NOPROC, APERV= [v1.0]'
      return
  endif 
  if not keyword_set(NFIND) then nfind = 1L
  if not keyword_set(FWHM) then fwhm = 5.0

  flg_arc = 0
  flg_flat = 0
  if keyword_set(ARCS) then flg_arc = 1 ;; Process the arcs
  if keyword_set(FLATS) then flg_flat = 1 ;; Process the flats
;  Optional Keywords
  if not keyword_set( ARCFIL ) then arcfil = $
    getenv('ESI_CALIBS')+'/ArcECH75_1x1IMG.fits' 
  if not keyword_set( FLATFIL ) then flatfil = $
    getenv('ESI_CALIBS')+'/FlatECH75_1x1N.fits' 
  if not keyword_set( BIASFIL ) then biasfil = $
    getenv('ESI_CALIBS')+'/BiasS.fits'
  if not keyword_set( SEDG_FIL ) then sedg_fil = $
    getenv('ESI_CALIBS')+'/SEdgECH75_1x1.fits'
  if not keyword_set( FITFIL ) then fitfil = $
    getenv('ESI_CALIBS')+'/hole_fit.idl'
  if not keyword_set( STDFIL ) then stdfil = $
    getenv('ESI_CALIBS')+'/ECH_STDOBJ.fits'
;  if not keyword_set(ORDRS) then ordrs=[0L,9L]

; Create the structure
  print, systime()
  print, 'esi_echqckrdx: Creating the structure'
  files = [obj_fil]
  if flg_arc then files = [files,arcs]
  if flg_flat then files = [files,flats]

  esi_strct, esi, IMG=files, /NOEDIT, /NOFILE

; Set tags
  esi.obj_id = 1L
  esi.arc_fil = arcfil
  esi.obj_fil = 'Extract/Obj_'+esi[0].img_root
  esi.flat_fil = FLATFIL

  slit = esi[0].slit

  ;; Process the Flats
  if flg_flat and not keyword_set(NOFLAT) then begin
      esi_echmkflat, esi, slit, HOT_THRESH = 30, /CLOBBER, BIASFIL=biasfil, IFLAT=iflat
      esi_echtraceorders, esi, slit, chk = chk, SEDG_FIL=sedg_fil, IFLAT=iflat
  endif

  ;; Process the Arcs
  if flg_arc and not keyword_set(NOARC) then begin
      esi_echmkarc, esi, slit, /clobber, BIASFIL=biasfil
      esi_echfitarc, esi, slit, /clobber
      esi_fit2darc, esi, slit, /clobber
      esi_echtrcarc, esi, slit, /auto, pksig = 5.0
      esi_echmkaimg, esi, slit, chk = chk, OUTFIL=arcfil
      esi.arc_fil = arcfil
  endif

  ;; Back to the Flat (if necessary)
  if flg_flat and not keyword_set(NOFLAT) then begin
      esi_echfltsct, esi, slit, chk = chk, OUTFIL=flatfil, IFLAT=iflat
      esi.flat_fil = flatfil
  endif

  ;; Process the image
  if not keyword_set( NOPROC ) then begin
      esi_echproc, esi, [0], /CLOBBER, FLATFIL=flatfil, $
                   BIASFIL=biasfil, NOOSCAN=nooscan, _EXTRA=extra
      if keyword_set( CHK ) then begin
          xatv, 'Final/f_'+esi[0].img_root, /block, min=-10, max=100.
      endif
  endif

  obj_id = 1L
  if keyword_set(EXTRA) then begin
      if tag_exist(EXTRA,'STD') then begin
          flg_std = 1
          obj_id = 0L
      endif
  endif

  ;; Identify the object
  if not keyword_set( NOFND ) then begin
      esi_echfndobj, esi, obj_id, CHK=chk, SEDG_FIL=sedg_fil, _EXTRA=extra, $
                     nfind=nfind, FWHM=fwhm, /USESTD, FILSTD=stdfil
  endif

  ;; Sky subtract
  if not keyword_set( NOSKY ) then begin
      if not keyword_set(OLD) then begin
          esi_echskysub, esi, obj_id, FCHK=skychk, SEDG_FIL=sedg_fil, $
                         BORDR=3, FITFIL=fitfil, ORDR=ordrs, /NO_SKYLINE, _EXTRA=EXTRA
      endif else begin
          esi_echold_skysub, esi, obj_id, FCHK=skychk, SEDG_FIL=sedg_fil, $
                         ORDR=ordrs, _EXTRA=EXTRA
      endelse
  endif
      
  ;; Another trace
  if keyword_set( TWOFND ) then begin
      esi_echfndobj, esi, obj_id, CHK=chk, /SKYSUB, FWHM=FWHM, /USESTD, $
                     nfind=nfind, SEDG_FIL=sedg_fil, FILSTD=stdfil, _EXTRA=extra
  endif

  ;; Extract
  if not keyword_set( NOEXT ) then begin
      if keyword_set(OLD) then begin
          esi_echold_extract, esi, obj_id, /OPTIMAL, CHK=chk, SEDG_FIL=sedg_fil
      endif else begin
          esi_echextobj, esi, obj_id, CHK=chk, SEDG_FIL=sedg_fil, NITER=1, ORDRS=ordrs, $
                         _EXTRA=extra
      endelse
      print, 'esi_echqckrdx: Type [esi_echspecplt] to view the output'
  endif

; SN

  obj = xmrdfits(esi[0].obj_fil,1)
  obj.sig = sqrt(obj.var > 0.)
  snfil = 's2n_'+strmid(esi[0].img_root,0,strlen(esi[0].img_root)-5)+'.ps'
  x_psopen, snfil, /maxs
  for i=0L,n_elements(obj)-1 do begin
      bad= where(obj[i].sig EQ 0., complement=good)
      obj[i].sig[bad] = 1e30
      s2n = obj[i].fx / obj[i].sig
      gs2n = s2n[good]
      srt = sort(gs2n)
      nsrt = n_elements(good)
      yval = gs2n[srt[nsrt*0.95]]
      plot, obj[i].wave, s2n, psym=10, yrang=[0., yval*1.3]
  endfor
  x_psclose;, /noprint, /noid
  print, 'esi_echqckrdx: s2n ps file is ', snfil

;  DONE
  print, 'esi_echqckrdx: All done! '
  print, systime()
  return
end
