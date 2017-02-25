;+ 
; NAME:
; esi_echpipe   
;     Version 1.0
;
; PURPOSE:
;    Given a list of files, process them in full.  This code is
;    designed to receive calls from the Reduction package designed and
;    built by Brad Holden (UCO)
;
; CALLING SEQUENCE:
;   
;  esi_echpipe, obj_fil
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
;   esi_echpipe, esi, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echpipe, obj_fil, arcs, flats, $
                 ARCFIL=arcfil, FLATFIL=flatfil, ORDRS=ordrs, OLD=old, $
                 CHK=chk, NOPROC=noproc, REFORDR=refordr, REFWV=refwv, $
                 TWOFND=twofnd, IFLAT=iflat, FLUXFIL=fluxfil, $
                 NOFND=nofnd, SKYCHK=skychk, NOSKY=nosky, TRCCHK=trcchk, $
                 NOTRC=notrc, NOEXT=noext, SCICLM=sciclm, NOFLAT=noflat, $
                 APERV=aperv, STANDARD=standard, STD=std, $
                 FINTER=finter, NOARC=noarc, NFIND=nfind, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
            'esi_echpipe, obj_fil, ARCFIL=, FLATFIL=, ' + $
            'SCICLM=, /CHK, /NOPROC, APERV= [v1.0]'
      return
  endif 
  if not keyword_set(NFIND) then nfind = 1L
  if not keyword_set(FWHM) then fwhm = 5.0

  flg_arc = 0
  flg_flat = 0
  if keyword_set(ARCS) then flg_arc = 1 ;; Process the arcs
  if keyword_set(FLATS) then flg_flat = 1 ;; Process the flats
  if keyword_set(STANDARD) then flg_std = 1 else flg_std = 0 ;; Process a standard
;  Optional Keywords
  if not keyword_set( FLUXFIL ) then fluxfil = $
    getenv('ESI_CALIBS')+'/sens_2008may09.idl' 

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

  if not keyword_set( STDFIL ) and not keyword_set(FLG_STD) then stdfil = $
    getenv('ESI_CALIBS')+'/ECH_STDOBJ.fits'

  ;; Create the structure
  print, systime()
  print, 'esi_echpipe: Creating the structure'
  files = [obj_fil]
  n_inobj = n_elements(files)
  if flg_arc then files = [files,arcs]
  if flg_flat then files = [files,flats]
  if flg_std then begin
      files = [files,standard]
      nstd = n_elements(standard)
  endif else nstd = 0
  nfil = n_elements(files)

  esi_strct, esi, IMG=files, /NOEDIT, /NOFILE

  ;; Set tags
  esi.obj_id = 1L
  esi.arc_fil = arcfil
  esi.obj_fil = 'Extract/Obj_'+esi.img_root
  esi.flat_fil = FLATFIL
  slit = esi[0].slit

  esi[0:n_inobj-1].type = 'OBJ'  ;; Short exposures would have STD
  if flg_std then $
    esi[nfil-nstd:nfil-1].type = 'STD' 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Process the Flats
  if flg_flat and not keyword_set(NOFLAT) then begin
      esi_echmkflat, esi, slit, HOT_THRESH = 30, /CLOBBER, BIASFIL=biasfil, IFLAT=iflat
      esi_echtraceorders, esi, slit, chk = chk, SEDG_FIL=sedg_fil, IFLAT=iflat
  endif

  ;; Process the Arcs
  if flg_arc and not keyword_set(NOARC) then begin
      esi_echmkarc, esi, slit, /clobber, BIASFIL=biasfil
      esi_echfitarc, esi, slit, /clobber, SEDG_FIL=sedg_fil
      esi_fit2darc, esi, slit, /clobber
      esi_echtrcarc, esi, slit, /auto, pksig = 5.0, SEDG_FIL=sedg_fil
      esi_echmkaimg, esi, slit, chk = chk, OUTFIL=arcfil, SEDG_FIL=sedg_fil
      esi.arc_fil = arcfil
  endif 

  ;; Back to the Flat (if necessary)
  if flg_flat and not keyword_set(NOFLAT) then begin
      esi_echfltsct, esi, slit, chk = chk, OUTFIL=flatfil, IFLAT=iflat,$
                     ARC_FILE=arcfil
      esi.flat_fil = flatfil
  endif


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Standard star  (when an science exposure is present too)
  if flg_std then begin  
      ;; Only taking the first exposure
      objind = (where(esi.type EQ 'STD'))[0]
      print, 'esi_echpipe: Using standard star from exposure ', esi[objind].img_root

      esi_echtrcstd, esi, slit, chk = chk, fwhm = 5.0, /clobber, $
                   BIASFIL=biasfil, FLATFIL=flatfil, SEDG_FIL=sedg_fil
      
  endif


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Science Object

  ;; Query on multiple objects
  all_obj = where(esi.type EQ 'OBJ', nobj)
  uni_objnm = uniq(esi[all_obj].obj, sort(esi[all_obj].obj))
  if n_elements(uni_objnm) GT 1 then begin
      print, 'esi_echpipe:  WARNING -- Not expecting differing names...'
      print, 'esi_echpipe:  WARNING -- Proceeding to extract them all but beware!'
      print, esi[all_obj[uni_objnm]].obj
  endif

  obj_id = 1L
  if keyword_set(EXTRA) then begin
      if tag_exist(EXTRA,'STD') then begin
          stop ;; Not functional right now
          flg_std = 1
          obj_id = 0L
      endif
  endif

  ;; Process the image
  if not keyword_set( NOPROC ) then begin
      objind = where(esi.obj_id EQ obj_id AND $
                     (esi.type EQ 'OBJ') AND $
                     esi.mode EQ 2 and esi.flg_anly NE 0)
      esi_echproc, esi, objind, /CLOBBER, FLATFIL=flatfil, $
                   BIASFIL=biasfil, _EXTRA=extra
  endif


  ;; Cosmic rays
  if nobj GT 1 then esi_echobjcr, esi, obj_id

  ;; Identify the object
  if not keyword_set( NOFND ) then begin
      esi_echfndobj, esi, obj_id, CHK=chk, SEDG_FIL=sedg_fil, _EXTRA=extra, $
                     nfind=nfind, FWHM=fwhm, /USESTD, FILSTD=stdfil
  endif

  ;; Sky subtract
  if not keyword_set( NOSKY ) then begin
      if not keyword_set(OLD) then begin
         if keyword_set(STD) then obj_in=0L else obj_in = obj_id
          esi_echskysub, esi, obj_in, FCHK=skychk, SEDG_FIL=sedg_fil, $
                         BORDR=3, FITFIL=fitfil, ORDR=ordrs, STD=std, _EXTRA=EXTRA
      endif else begin
          esi_echold_skysub, esi, obj_id, FCHK=skychk, SEDG_FIL=sedg_fil, $
                         ORDR=ordrs, _EXTRA=EXTRA
      endelse
  endif
      
  ;; Another trace
  if not keyword_set( NOFND ) and not keyword_set(STD) then begin
      esi_echfndobj, esi, obj_id, CHK=chk, /SKYSUB, FWHM=FWHM, /USESTD, $
                     nfind=nfind, SEDG_FIL=sedg_fil, _EXTRA=extra, FILSTD=stdfil
  endif

  ;; Extract
  if not keyword_set( NOEXT ) then begin
      if keyword_set(OLD) then begin
          esi_echold_extract, esi, obj_id, /OPTIMAL, CHK=chk, SEDG_FIL=sedg_fil
      endif else begin
         if keyword_set(STD) then obj_in=0L else obj_in = obj_id
          esi_echextobj, esi, obj_in, CHK=chk, SEDG_FIL=sedg_fil, ORDRS=ordrs, $
                         STD=std, _EXTRA=extra
      endelse
;      print, 'esi_echpipe: Type [esi_echspecplt] to view the output'
  endif

  ;;;;;;;;;;;;;;;
  ;; Coadd
  esi_echcombspec, esi, obj_id, obj_nm = 'a'
  esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, obj_nm = 'a'
  esi_echcoaddfin, esi, obj_id, obj_nm = 'a', /SKY, /NOVAR

  ;; SN
  for qq=0L,nobj-1 do begin
      obj = xmrdfits(esi[all_obj[qq]].obj_fil,1)
      obj.sig = sqrt(obj.var > 0.)
      snfil = 's2n_'+strmid(esi[all_obj[qq]].img_root,0, $
                            strlen(esi[all_obj[qq]].img_root)-5)+'.ps'
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
      x_psclose                 ;, /noprint, /noid
      print, 'esi_echpipe: s2n ps file is ', snfil
  endfor

;  DONE
  print, 'esi_echpipe: All done! '
  print, systime()
  return
end
