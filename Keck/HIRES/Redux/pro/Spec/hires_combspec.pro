;+ 
; NAME:
; hires_combspec
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
;   hires_combspec, hires, setup, obj_id, chip, exp_id
;
; INPUTS:
;   hires    - HIRES structure
;   setup
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   echfspec   -  HIRES fspec structure (fits file; FSpec/name.fits)
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
;   hires_combspec, hires, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_echcombspec
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;   10-Jun-2004 Revisited by GEP
;   Sep-2005 Revisited by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hires_combspec, hires, setup, obj_id, chip, exp_id, SILENT=silent, $
                    LIST=list, OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, $
                    SNRMIN=snrmin, KEYINDX=keyindx, MEDINDX=medindx, $
                    NOFLUX=noflux,  ORDRS=ordrs, BOXCAR=boxcar, $
                    MINPIX1=minpix1, ORDNM=ordnm, $
                    SIGREJ=sigrej, NOREJ=norej, OUTFIL=outfil, _EXTRA=extra

;
  if  N_params() LT 4  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'hires_combspec, hires, setup, obj_id, chip, [exp_id], /SILENT, ' + $
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
    print, 'hires_echcombspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') AND $
                          hires.flg_anly NE 0 AND $
                          hires.setup EQ setup AND hires.obj_id EQ obj_id $
                          AND hires.chip EQ chip)
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


;  if not keyword_set( CORDR ) then begin
;      case chip of 
;          1: CORDR = 95
;          2: CORDR = 85
;          3: CORDR = 70
;      endcase
;  endif

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'hires_echcombspec: Loading up the files...'

; Open Obj files
  svexpt = fltarr(nexp)
  for q=0L,nexp-1 do begin
      ;; Read
      if not keyword_set( LIST ) then begin
          print, 'hires_combspec: Reading ', hires[exp[q]].obj_fil
          if x_chkfil(hires[exp[q]].obj_fil+'*') EQ 0 then begin
              print, 'hires_echcombspec: Obj file doesnt exist! Returning..', $
                hires[exp[q]].obj_fil
              return
          endif
;          tmp = xmrdfits(hires[exp[q]].obj_fil, 1, STRUCTYP='hiresobjstrct', $
          tmp = xmrdfits(hires[exp[q]].obj_fil, 1, /silent)
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'hires_echcombspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
;          tmp = xmrdfits(files[q], 1, STRUCTYP='hiresobjstrct', /silent)
          print, 'hires_combspec: Reading ', files[q]
          tmp = xmrdfits(files[q], 1, /silent)
      endelse
      if q EQ 0 then no1 = n_elements(tmp)
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm, ngd)
      msko = bytarr(ngd)
      ;; Deal with updates to the HIRES obj structure
      tmp_str = replicate({hiresobjstrct},ngd)
      for ss=0L,ngd-1 do begin
          dums = {hiresobjstrct}
          ;; Check for variance 
          if max(tmp[gdobj[ss]].var) GT 0. then msko[ss] = 1B
          ;;
          copy_struct, tmp[gdobj[ss]], dums
          tmp_str[ss] = dums
       endfor

      ;; Add to total structure
      if q EQ 0 then hiresobj = tmp_str else hiresobj = [hiresobj, tmp_str[where(msko)]]
      svexpt[q] = tmp[gdobj[0]].exp
  endfor


  ;; Setup orders
  all_ordr = hiresobj.order
  all_ordr = all_ordr[uniq(all_ordr,sort(all_ordr))]

  if not keyword_set(ORDRS) then ordrs = [min(all_ordr), max(all_ordr)]
  nspec = nexp

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then subfil = outnm+obj_nm $
      else subfil = strcompress(strtrim(hires[exp[0]].Obj,2),/remove_all)+obj_nm
      case chip of 
          -1: clrc = '_S' 
          1: clrc = '_B' 
          2: clrc = '_G'
          3: clrc = '_R'
      endcase
      outfil = 'FSpec/'+subfil+clrc+'.fits'
  endif else begin
      if not keyword_set(OUTFIL) then outfil = 'FSpec/tmp_ech.fits' 
  endelse

  echfspec = { hiresfspecstrct }

  ;; Copy
  echfspec.nexp = nspec

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = svexpt[i]

  ;; Other tags
  copy_struct, hiresobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  if not keyword_set(LIST) then $
      for i=0L,nspec-1 do echfspec.obj_fil[i] = hires[exp[i]].obj_fil $
  else for i=0L,nspec-1 do echfspec.obj_fil[i] = files[i]

  x_echcombspec, hiresobj, echfspec, ordrs, keyindx, SIGREJ=sigrej, $
    REJSIG=rejsig, ORDNM=ordnm, BOXCAR=boxcar, NOFLUX=noflux, _EXTRA=extra

  ;;;; OUTPUT  ;;;;
  hires_wrfspec, echfspec, outfil 

;  close, 56

  print, 'hires_combspec:  Output is in ', outfil
  print, 'hires_combspec:  All done!'


  return
end
  

