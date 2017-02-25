;+ 
; NAME:
; hamspec_combspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run with the pipeline, even on an object observed for
;    only a single exposure.  
;
;  The main call is to x_echcombspec
;
; CALLING SEQUENCE:
;   hamspec_combspec, hamspec, setup, obj_id, exp_id
;
; INPUTS:
;   hamspec    - Hamspec structure
;   setup
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   echfspec   -  Hamspec fspec structure (fits file; FSpec/name.fits)
;
; OPTIONAL KEYWORDS:
;    KEYINDX=  - Exposure to use as the 'fidcuial'.   [Default: 0]
;    LIST=  -- ASCII file containing a list of exposures to combine
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file (Good for
;                list)
;    /NOFLUX   - Do not use the fluxed array
;    /BOXCAR   - Use the boxcar extraction
;    MINPIX1=  - Minimum 'good' pixels to calculate fitting profile
;    ORDNM=    - Fitting profile order [default: 1]
;    SNRMIN=   - Minimum S/N per pixel for stats [Default: 2]
;    REJSIG=   - Parameter passsed to x_combine [default = 4.]
;    MEDINDX=  - Value for median smoothing [default: 100]
;    /NOREJ    - Do not search for and reject bad pixels (e.g. cosmic
;                rays)
;    ORDRS     - Orders to combine [default: all]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hamspec_combspec, hamspec, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_echcombspec
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;   10-Jun-2004 Revisited by GEP
;   Sep-2005 Revisited by JXP
;   24-May-2013 Slightly modified from hires_combspec by ENK
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hamspec_combspec, hamspec, setup, obj_id, exp_id, SILENT=silent, $
                    LIST=list, OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, $
                    SNRMIN=snrmin, KEYINDX=keyindx, MEDINDX=medindx, $
                    NOFLUX=noflux,  ORDRS=ordrs, BOXCAR=boxcar, $
                    MINPIX1=minpix1, ORDNM=ordnm, $
                    SIGREJ=sigrej, NOREJ=norej, OUTFIL=outfil, _EXTRA=extra

;
  if  N_params() LT 3  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'hamspec_combspec, hamspec, setup, obj_id, [exp_id], /SILENT, ' + $
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
    print, 'hamspec_combspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where((hamspec.type EQ 'OBJ' OR hamspec.type EQ 'STD') AND $
                          hamspec.flg_anly NE 0 AND $
                          hamspec.setup EQ setup AND hamspec.obj_id EQ obj_id)
          if n_elements(exp_id) GT 0 then exp = allexp[exp_id] else exp=allexp
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
    print, 'hamspec_combspec: Loading up the files...'

; Open Obj files
  svexpt = fltarr(nexp)
  for q=0L,nexp-1 do begin
      ;; Read
      if not keyword_set( LIST ) then begin
          print, 'hamspec_combspec: Reading ', hamspec[exp[q]].obj_fil
          if x_chkfil(hamspec[exp[q]].obj_fil+'*') EQ 0 then begin
              print, 'hamspec_combspec: Obj file doesnt exist! Returning..', $
                hamspec[exp[q]].obj_fil
              return
          endif
;          tmp = xmrdfits(hamspec[exp[q]].obj_fil, 1, STRUCTYP='hamspecobjstrct', $
          tmp = xmrdfits(hamspec[exp[q]].obj_fil, 1, /silent)
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'hamspec_combspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
;          tmp = xmrdfits(files[q], 1, STRUCTYP='hamspecobjstrct', /silent)
          print, 'hamspec_combspec: Reading ', files[q]
          tmp = xmrdfits(files[q], 1, /silent)
      endelse
      if q EQ 0 then no1 = n_elements(tmp)
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm, ngd)
      msko = bytarr(ngd)
      ;; Deal with updates to the Hamspec obj structure
      tmp_str = replicate({hamspecobjstrct},ngd)
      for ss=0L,ngd-1 do begin
          dums = {hamspecobjstrct}
          ;; Check for variance 
          if max(tmp[gdobj[ss]].var) GT 0. then msko[ss] = 1B
          ;;
          copy_struct, tmp[gdobj[ss]], dums
          tmp_str[ss] = dums
       endfor

      ;; Add to total structure
      if q EQ 0 then hamspecobj = tmp_str else hamspecobj = [hamspecobj, tmp_str[where(msko)]]
      svexpt[q] = tmp[gdobj[0]].exp
  endfor


  ;; Setup orders
  all_ordr = hamspecobj.order
  all_ordr = all_ordr[uniq(all_ordr,sort(all_ordr))]

  if not keyword_set(ORDRS) then ordrs = [min(all_ordr), max(all_ordr)]
  nspec = nexp

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(hamspec[exp[0]].Obj,2),/remove_all)+obj_nm
      outfil = 'FSpec/'+subfil+'.fits'
  endif else begin
      if not keyword_set(OUTFIL) then outfil = 'FSpec/tmp_ech.fits' 
  endelse

  echfspec = { hamspecfspecstrct }

  ;; Copy
  echfspec.nexp = nspec

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = svexpt[i]

  ;; Other tags
  copy_struct, hamspecobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  if not keyword_set(LIST) then $
      for i=0L,nspec-1 do echfspec.obj_fil[i] = hamspec[exp[i]].obj_fil $
  else for i=0L,nspec-1 do echfspec.obj_fil[i] = files[i]

  hamspecobj.sig = sqrt(hamspecobj.var)
  hamspecobj.flux = hamspecobj.fx
  x_echcombspec, hamspecobj, echfspec, ordrs, keyindx, SIGREJ=sigrej, $
    REJSIG=rejsig, ORDNM=ordnm, BOXCAR=boxcar, NOFLUX=noflux, _EXTRA=extra

  ;;;; OUTPUT  ;;;;
  hamspec_wrfspec, echfspec, outfil 

;  close, 56

  print, 'hamspec_combspec:  Output is in ', outfil
  print, 'hamspec_combspec:  All done!'


  return
end
  

