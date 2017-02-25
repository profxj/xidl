;+ 
; NAME:
; uves_combspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run even on an object with a single exposure
;
; CALLING SEQUENCE:
;   
;   uves_combspec, uves, setup, obj_id, chip, exp_id
;
; INPUTS:
;   uves      - MIKE structure
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
;    REJSIG=   - Parameter passsed to x_combine [default = 4.]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   uves_combspec, uves, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;   10-Jun-2004 Revisited by GEP
;   Sep-2005 Revisited by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro uves_combspec, uves, setup, obj_id, side, exp_id, SILENT=silent, $
                    LIST=list, $
                    OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, $
                    SNRMIN=snrmin, $
                    FLUX=flux, KEYINDX=keyindx, MEDINDX=medindx, $
                    NOFLUX=noflux, $
                    MINPIX1=minpix1, ENDTRM=endtrm, ORDNM=ordnm, $
                    SIGREJ=sigrej, NOREJ=norej, $
                    LSIDE=lside, OUTFIL=outfil, _EXTRA=extra

;
  if  N_params() LT 4  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'uves_combspec, uves, setup, obj_id, side, [exp_id], /SILENT, ' + $
      'OBJ_NM=, /PATCH, /CHK, /MCHK, /NOREJ'
    print, '      LIST= [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(KEYINDX) then keyindx = 0L
  if not keyword_set(MEDINDX) then medindx = 100
  if size(ORDNM,/type) EQ 0 then ordnm = 1
  if not keyword_set(SIGREJ) then sigrej=4.
  if not keyword_set(REJSIG) then rejsig=4.
  if not keyword_set(MINPIX1) then minpix1 = 100
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( SNRMIN ) then snrmin = 2.

  if keyword_set(NOREJ) then REJSIG = -1.

  if keyword_set( LIST ) then $
    print, 'uves_echcombspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where((uves.type EQ 'OBJ' OR uves.type EQ 'STD') AND $
                          uves.flg_anly NE 0 AND $
                          uves.setup EQ setup AND uves.obj_id EQ obj_id $
                          AND uves.side EQ side)
          if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
          nexp = n_elements(exp)
      endif else begin  ;; STD
          exp = obj_id[0]
          nexp = 1
          obj_id = 99L
      endelse
  endif else begin
;      if not keyword_set(LCHIP) then begin
;          print, 'uves_combspec: Need to specify LCHIP'
;          return
;      endif
      obj_id = 99L
      readcol, list, files, FORMAT='A'
      nexp = n_elements(files)
  endelse


;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'uves_combspec: Loading up the files...'

; Open Obj files
  svexpt = fltarr(nexp)
  for q=0L,nexp-1 do begin
      ;; Read
      if not keyword_set( LIST ) then begin
          print, 'uves_combspec: Reading ', uves[exp[q]].obj_fil
          if x_chkfil(uves[exp[q]].obj_fil+'*') EQ 0 then begin
              print, 'uves_combspec: Obj file doesnt exist! Returning..', $
                uves[exp[q]].obj_fil
              return
          endif
          tmp = xmrdfits(uves[exp[q]].obj_fil, 1, STRUCTYP='uvesobjstrct', $
                         /silent)
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'uves_combspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
          tmp = xmrdfits(files[q], 1, STRUCTYP='uvesobjstrct', /silent)
      endelse
      if q EQ 0 then no1 = n_elements(tmp)
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm, ngd)
      ;; Add to total structure
      if q EQ 0 then uvesobj = tmp[gdobj] else uvesobj = [uvesobj, tmp[gdobj]]
      svexpt[q] = tmp[gdobj[0]].exp
  endfor


  ;; Setup orders
  all_ordr = uvesobj.order
  all_ordr = all_ordr[uniq(all_ordr,sort(all_ordr))]

  ordrs = [min(all_ordr), max(all_ordr)]

  nspec = nexp

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(uves[exp[0]].Obj,2),/remove_all)+obj_nm
      ;; Setup
      clrc = strtrim(round(uves[exp[0]].xdangl),2)
      outfil = 'FSpec/'+subfil+'_'+clrc+'.fits'
  endif else begin
      if not keyword_set(OUTFIL) then outfil = 'FSpec/tmp_ech.fits' 
  endelse

  echfspec = { uvesfspecstrct }

  ;; Copy
  echfspec.nexp = nspec

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = svexpt[i]

  ;; Other tags
  copy_struct, uvesobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  if not keyword_set(LIST) then $
      for i=0L,nspec-1 do echfspec.obj_fil[i] = uves[exp[i]].obj_fil $
  else for i=0L,nspec-1 do echfspec.obj_fil[i] = files[i]

  x_echcombspec, uvesobj, echfspec, ordrs, keyindx, SIGREJ=sigrej, $
    REJSIG=rejsig, ORDNM=ordnm, NOFLUX=noflux, _EXTRA=extra

  ;;;; OUTPUT  ;;;;
  uves_wrfspec, echfspec, outfil 

;  close, 56

  print, 'uves_combspec:  Output is in ', outfil
  print, 'uves_combspec:  All done!'


  return
end
  

