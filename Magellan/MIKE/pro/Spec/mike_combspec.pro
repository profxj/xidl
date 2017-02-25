;+ 
; NAME:
; mike_combspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run even on an object with a single exposure
;
; CALLING SEQUENCE:
;   
;   mike_combspec, mike, setup, obj_id, side, exp_id
;
; INPUTS:
;   mike      - MIKE structure
;   obj_id   -  Object ID  (e.g. 0L, 1L, etc)
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   echfspec      -  ESI fspec structure (fits file; FSpec/name_ech.fits)
;
; OPTIONAL KEYWORDS:
;    /SILENT   - No text output
;    LIST=     - Use an input list to combine (required for multiple
;                night obs).  Output is in 'FSpec/tmp_ech.fits', unless OUTFIL flag set to something.
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file
;    /NOFLUX   - Fill echfspec structure with fx, rather than flux array 
;    ENDTRM=   - Trim order edges by this when calculating ratios
;    MINPIX1=  - Minimum 'good' pixels to calculate fitting profile
;    ORDNM=    - Fitting profile order (default: 1)
;    SNRMIN=   - Minimum S/N per pixel for consideration in fitting
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   mike_combspec, mike, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;   10-Jun-2004 Revisited by GEP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro mike_combspec, mike, setup, obj_id, side, exp_id, SILENT=silent, LIST=list, $
                   OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, SNRMIN=snrmin, $
                   FLUX=flux, KEYINDX=keyindx, MEDINDX=medindx, NOFLUX=noflux, $
                   MINPIX1=minpix1, ENDTRM=endtrm, ORDNM=ordnm, CHK=chk, $
                   SIGREJ=sigrej, MCHK=mchk, $
                   LSIDE=lside, OUTFIL=outfil, _extra=extra

;
  if  N_params() LT 4  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'mike_combspec, mike, setup, obj_id, side, [exp_id], /SILENT, ' + $
      'OBJ_NM=, /PATCH '
    print, '      LIST= [v1.1]'
    return
  endif 
  if keyword_set(LIST) and keyword_set(LSIDE) then side = lside
  if keyword_set(LIST) and not keyword_set(LSIDE) then begin
      print, 'Warning -- with LIST keyword the side needs to be explicitly set,' + $
              ' with use of LSIDE keyword, defaulting to blue!'
      side=1
  endif

;  Optional Keywords
  if not keyword_set(KEYINDX) then keyindx = 0L
  if not keyword_set(MEDINDX) then medindx = 100
  if not keyword_set(ENDTRM) then endtrm = 0
  if not keyword_set(ORDNM) then ordnm = 1
  if not keyword_set(SIGREJ) then sigrej=4.
  if not keyword_set(MINPIX1) then minpix1 = 100
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( CORDR ) then begin
      if side EQ 1 then CORDR = 85L else CORDR=55L
  endif
  if not keyword_set( SORDR ) then begin
      if side EQ 1 then SORDR = 70L else SORDR=30L
  endif
  if not keyword_set( SNRMIN ) then snrmin = 2.

  if keyword_set( LIST ) then $
    print, 'mike_echcombspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where((mike.type EQ 'OBJ' OR mike.type EQ 'STD') AND $
                          mike.flg_anly NE 0 AND $
                          mike.setup EQ setup AND mike.obj_id EQ obj_id $
                          AND mike.side EQ side)
          if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
          nexp = n_elements(exp)
      endif else begin  ;; STD
          exp = obj_id[0]
          nexp = 1
          obj_id = 99L
      endelse
  endif else begin
      obj_id = 99L
      readcol, list, files, FORMAT='A'
      nexp = n_elements(files)
  endelse


;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'mike_echcombspec: Loading up the files...'

; Open Obj files
  for q=0L,nexp-1 do begin
      ;; Read
      if not keyword_set( LIST ) then begin
          print, 'mike_combspec: Reading ', mike[exp[q]].obj_fil
          if x_chkfil(mike[exp[q]].obj_fil+'*') EQ 0 then begin
              print, 'mike_echcombspec: Obj file doesnt exist! Returning..', $
                mike[exp[q]].obj_fil
              return
          endif
          tmp = xmrdfits(mike[exp[q]].obj_fil, 1, /silent)
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'mike_echcombspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
          tmp = xmrdfits(files[q], 1, /silent)
      endelse
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm, ngd)
      ;; Add to total structure
      if q EQ 0 then mikeobj = tmp[gdobj] else mikeobj = [mikeobj, tmp[gdobj]]
   endfor

  ;; Kludge in case we are using old versions of the MIKE Obj structure
  if tag_exist(mikeobj,'NOSIG') EQ 0 then begin
     tmp2 = replicate({mikeobjstrct}, n_elements(mikeobj))
     copy_struct, mikeobj, tmp2
     mikeobj = tmp2
  endif

; Setup orders
  all_ordr = mikeobj.order
  all_ordr = all_ordr[uniq(all_ordr,sort(all_ordr))]

  ordrs = [min(all_ordr), max(all_ordr)]

  nspec = nexp

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(mike[exp[0]].Obj,2),/remove_all)+obj_nm
;      outfil = mike_getfil('fspec_fil',subfil=subfil, SIDE=side, /name)
      if side EQ 1 then clrc = '_b' else clrc = '_r'
      outfil = 'FSpec/'+subfil+clrc+'.fits'
  endif else begin
      if not keyword_set(OUTFIL) then outfil = 'FSpec/tmp_ech.fits' 
  endelse

  echfspec = { mikefspecstrct }

  ;; Copy
  echfspec.nexp = nspec
  ;; Center ordr
  all_cen = where(mikeobj.order EQ CORDR, ncen)
  if ncen NE nspec then stop

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = mikeobj[all_cen[i]].exp

  ;; Other tags
  copy_struct, mikeobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  if not keyword_set(LIST) then $
      for i=0L,nspec-1 do echfspec.obj_fil[i] = mike[exp[i]].obj_fil $
  else for i=0L,nspec-1 do echfspec.obj_fil[i] = files[i]

  x_echcombspec, mikeobj, echfspec, ordrs, keyindx, SIGREJ=sigrej, $
    ORDNM=ordnm, MCHK=mchk, NOFLUX=noflux, _EXTRA=extra

  ;;;; OUTPUT  ;;;;
  mike_wrfspec, echfspec, outfil 

;  close, 56

  print, 'mike_combspec:  Output is in ', outfil
  print, 'mike_combspec:  All done!'


  return
end
  

