;+ 
; NAME:
; esi_echcombspec
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run even on an object with a single exposure
;
; CALLING SEQUENCE:
;   
;   esi_echcombspec, esi, obj_id, exp_id
;
; INPUTS:
;   esi      - ESI structure
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
;                night obs).  Output is in 'FSpec/tmp_ech.fits'.
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file
;    ORDRS=    - Orders to combine (default: lindgen(10L)+6 )
;    CHK       - Various QA for call to long_combspec for combining
;                spectra. 
;    NPOLY=  Order polynomial for fit.  Best to let long_combspec
;    decide if you can
;    /USE_OLD      - Old coadd.  Best for very high S/N obs
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echcombspec, esi, obj_id, exp_id
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Aug-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------

pro esi_echcombspec_out, echfspec

  ;; Output to ASCII
  printf, 56, echfspec.field
  for j=0L,echfspec.nexp-1 do begin
      printf, 56, FORMAT='(10x,f7.1,1x,2f10.3,1x,a25)',$
        echfspec.texp[j], $
        echfspec.wvmnx[j,0], $
        echfspec.wvmnx[j,1], $
        echfspec.obj_fil[j]
  endfor
  return
end

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro esi_echcombspec, esi, obj_id1, exp_id, SILENT=silent, LIST=list $
                     , OBJ_NM = OBJ_NM1, STD = std, OUTNM = outnm $
                     , ORDRS = ordrs, REFO = refo, CHK = CHK $
                     , NOSHIFT = NOSHIFT, SHIFT_SN = SHIFT_SN1 $
                     , SN_MIN_MEDSCALE = SN_MIN_MEDSCALE1 $
                     , SN_MAX_MEDSCALE = SN_MAX_MEDSCALE1 $
                     , IREF = IREF, npoly = npoly, USE_OLD = use_old $
                     , MXSHIFT = MXSHIFT, DEFRINGE = DEFRINGE $
                     , HIZQSO = HIZQSO, LAM_MASK_MIN = LAM_MASK_MIN1 $
                     , USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS1
  
 
;
  if  N_params() LT 2  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'esi_echcombspec, esi, obj_id, [exp_id], /SILENT, OBJ_NM=, /PATCH '
    print, '      LIST=, REFO=, ORDRS= [v1.1]'
    return
  endif 


;;Bug fix added by JFH, prevents return of obj_id =99 value which breaks 
;; later calls to esi routines. 
IF n_elements(obj_id1) GT 0 THEN obj_id = obj_id1 
;  Optional Keywords
  if not keyword_set(ORDRS) then ordrs=[0L,9L]
  if not keyword_set(REFO) then refo = 0L
  nord = ordrs[1]-ordrs[0] + 1L

  IF n_elements(SHIFT_SN1) GT 0 THEN BEGIN
     IF n_elements(SHIFT_SN1) EQ 1 THEN SHIFT_SN = replicate(SHIFT_SN1, nord) $
     ELSE IF n_elements(SHIFT_SN1) EQ nord THEN SHIFT_SN = reverse(SHIFT_SN1) $
     ELSE message, 'Problem with size of SHIFT_SN'
  ENDIF ELSE SHIFT_SN = replicate(0, nord)

  IF n_elements(SN_MIN_MEDSCALE1) GT 0 THEN BEGIN
     IF n_elements(SN_MIN_MEDSCALE1) EQ 1 THEN SN_MIN_MEDSCALE = replicate(SN_MIN_MEDSCALE1, nord) $
     ELSE IF n_elements(SN_MIN_MEDSCALE1) EQ nord THEN SN_MIN_MEDSCALE = reverse(SN_MIN_MEDSCALE1) $
     ELSE message, 'Problem with size of SN_MIN_MEDSCALE'
  ENDIF ELSE SN_MIN_MEDSCALE = replicate(0, nord)

  IF n_elements(SN_MAX_MEDSCALE1) GT 0 THEN BEGIN
     IF n_elements(SN_MAX_MEDSCALE1) EQ 1 THEN SN_MAX_MEDSCALE = replicate(SN_MAX_MEDSCALE1, nord) $
     ELSE IF n_elements(SN_MAX_MEDSCALE1) EQ nord THEN SN_MAX_MEDSCALE = reverse(SN_MAX_MEDSCALE1) $
     ELSE message, 'Problem with size of SN_MAX_MEDSCALE'
  ENDIF ELSE SN_MAX_MEDSCALE = replicate(0, nord)

  IF n_elements(USE_AVG_SN_WEIGHTS1) EQ 0 THEN USE_AVG_SN_WEIGHTS1 = 0
  
  if keyword_set( LIST ) then $
    print, 'esi_echcombspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
                         esi.mode EQ 2 AND esi.obj_id EQ obj_id)
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

  IF KEYWORD_SET(OBJ_NM1) THEN BEGIN
     IF n_elements(obj_nm1) EQ 1 THEN obj_nm =  replicate(obj_nm1, nexp) $
     ELSE IF n_elements(obj_nm1) EQ nexp THEN obj_nm = obj_nm1 $
     ELSE message, 'Invalid size for obj_nm'
  ENDIF ELSE obj_nm = replicate('a', nexp)



;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) AND NOT KEYWORD_SET(LIST) then BEGIN
     print, 'esi_echcombspec: Loading up the files...'
     IF KEYWORD_SET(DEFRINGE) THEN obj_files = repstr(esi[exp].obj_fil, 'Obj', 'FringeObj') $
     ELSE obj_files = esi[exp].obj_fil
     forprint, obj_files, textout = 1
  ENDIF

; Open Obj files
  svobj = strarr(nexp)
  for q=0L,nexp-1 do begin
      ;; Read
     if not keyword_set( LIST ) then begin
        IF KEYWORD_SET(DEFRINGE) THEN obj_fil = repstr(esi[exp[q]].obj_fil, 'Obj', 'FringeObj') $
        ELSE obj_fil = esi[exp[q]].obj_fil
        if x_chkfil(obj_fil+'*') EQ 0 then begin
           print, 'esi_echcombspec: Obj file doesnt exist! Returning..', obj_fil 
           return
        endif
        tmp = xmrdfits(obj_fil, 1, /silent)
        svobj[q] = obj_fil
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'esi_echcombspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
          tmp = xmrdfits(files[q], 1, /silent)
          svobj[q] = files[q]  ;; This might not work (JXP 8/16/05)
      endelse
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm[q])
      ;; deal with new tags being a problem for combining with old data (MR 5/13/14)
      testtags = where(tag_names(tmp) eq 'HAND_SUB', ntest)
      if ntest gt 0 then remove_tags, tmp, ['HAND_SUB', 'HAND_MAXX', 'HAND_MINX'], tmp_new
      tmp = tmp_new
      ;; Add to total structure
      if q EQ 0 then esiobj = tmp[gdobj] else BEGIN
          tags = tag_names(esiobj[0])
          tmp_now =  struct_selecttags(tmp[gdobj], select_tags = tags)
          esiobj = [esiobj, tmp_now]
      ENDELSE
  endfor

  nspec = n_elements(esiobj) / 10  ; 10 orders per spec
  if nspec NE nexp then stop
  

;; ASCII OUTPUT ;;;

;  if obj_id LT 10L then txtfil = 'FSpec/txt/ECHF_0'+$
;    string(obj_id,format='(i1)')+'.txt' $
;  else txtfil = 'FSpec/txt/ECHF_'+string(obj_id,format='(i2)')+obj_nm+'.txt'
; Change made by JFH 9/03/06
  txtfil = 'FSpec/txt/ECHF_' + strcompress(string(obj_id, format = '(I2.2)') $
                                           , /rem)
  close, 56
  openw, 56, txtfil

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then $
         outfil = 'FSpec/'+outnm+obj_nm[0]+'_ech.fits' $
      else outfil = 'FSpec/'+strtrim(esi[exp[0]].Obj, 2)+obj_nm[0]+'_ech.fits' 
   endif else BEGIN
      IF KEYWORD_SET(OUTNM)THEN  outfil = 'FSpec/'+outnm+obj_nm[0]+'_ech.fits' $
     ELSE outfil = 'FSpec/tmp_ech.fits' 
  ENDELSE

  echfspec = { echfspecstrct }

  ;; Copy
  echfspec.nexp = nspec
  ;; Zero ordr
  all_zero = where(esiobj.order EQ 0, nzero)
  if nzero NE nspec then stop

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = esiobj[all_zero[i]].exp

  ;; Other tags
  copy_struct, esiobj[0], echfspec $
               , EXCEPT = ["wave", "fx", "var", "novar", "sky", "npix"]
  for i=0L,nspec-1 do begin
;      ipos = strpos(esiobj[all_zero[i]].spec2d_fil, 'esi')
;      obj_fil = 'Extract/Obj_'+strmid(esiobj[all_zero[i]].spec2d_fil, ipos)
      ipos = strpos(svobj[i], 'esi')
      IF KEYWORD_SET(DEFRINGE) THEN  $
         obj_fil = 'Extract/FringeObj_'+strmid(svobj[i], ipos) $
      ELSE obj_fil = 'Extract/Obj_'+strmid(svobj[i], ipos)
      echfspec.obj_fil[i] = obj_fil
  endfor

  svmedian = fltarr(nspec)
  svsnr = fltarr(nspec)
  IF KEYWORD_SET(HIZQSO) THEN BEGIN
     ;; Find the first order moving from the red that contains Ly-a,
     ;; i.e. LAM_MASK_MIN
     ;; Loop on ORDER (reverse)
     FOR qq = ordrs[1], ordrs[0], -1 DO BEGIN
        all_ordr = where(esiobj.order EQ qq)
        if nspec EQ 1 then npix = esiobj[all_ordr].npix $
        ELSE npix = 5000L
        inloglam = alog10(esiobj[all_ordr[REFO]].wave[0:npix-1])
        inside = WHERE(LAM_MASK_MIN1 GE min(10.0^(inloglam)) AND LAM_MASK_MIN1 LE max(10.0^(inloglam)), nin)
        IF nin GT 0 THEN BEGIN
           QQ_MIN = qq
           BREAK
        ENDIF
     ENDFOR
  ENDIF
  
  ;; Loop on ORDER (reverse)
  for qq=ordrs[1],ordrs[0],-1 do begin
      all_ordr = where(esiobj.order EQ qq)
      ;; Coadd
      if nspec EQ 1 then begin
          npix = esiobj[all_ordr].npix
          echfspec.wave[0:npix-1,qq] = esiobj[all_ordr].wave[0:npix-1]
          echfspec.fx[0:npix-1,qq] = esiobj[all_ordr].fx[0:npix-1]
          echfspec.var[0:npix-1,qq] = esiobj[all_ordr].var[0:npix-1]
          echfspec.novar[0:npix-1, qq] = esiobj[all_ordr].novar[0:npix-1]
          echfspec.sky[0:npix-1, qq] = esiobj[all_ordr].sky[0:npix-1]
      endif else begin
          npix = 5000L
          ;; Find median ratio of flux and SNR
          med_rto = fltarr(nspec)
          snr_rto = fltarr(nspec)
          med_rto[REFO] = 1.
          snr_rto[REFO] = 1.
          gdpix = where(esiobj[all_ordr[REFO]].var GT 0.)
          med0 = median(esiobj[all_ordr[REFO]].fx[gdpix])
          if med0 LT 10 then begin
              if qq EQ ordrs[1] then begin
                  ;; Use exposure time ratio
                  for kk=1L,nspec-1 do begin
                      med_rto[kk] = esiobj[all_ordr[REFO]].exp/ $
                        esiobj[all_ordr[kk]].exp
                      snr_rto[kk] = sqrt(esiobj[all_ordr[REFO]].exp/ $
                        esiobj[all_ordr[kk]].exp)
                  endfor
                  svmedian[*] = med_rto
                  svsnr[*] = snr_rto
              endif else begin
                  med_rto[*] = svmedian
                  snr_rto[*] = svsnr
              endelse
          endif else begin
              ;; Take ratio
              for kk=0L,nspec-1 do begin
                  med_rto[kk] = median(esiobj[all_ordr[REFO]].fx[gdpix]/$
                                       esiobj[all_ordr[kk]].fx[gdpix])
                  snr_rto[kk] = median(esiobj[all_ordr[kk]].fx[gdpix]/ $
                                       sqrt(esiobj[all_ordr[kk]].var[gdpix])/ $
                                       esiobj[all_ordr[REFO]].fx[gdpix] * $
                                       sqrt(esiobj[all_ordr[REFO]].var[gdpix]))
              endfor
              ;; Save
              svmedian[*] = med_rto
              svsnr[*] = snr_rto
          endelse

          ;; Set bad pix to -1

;          for kk=0L,nspec-1 do begin
;              a = where(esiobj[all_ordr[kk]].var[0:npix-1] LE 0., na)
;              if na GT 0 then var[a,kk] = -1.
;          endfor
          print, 'esi_echcombspec: Order ', 15-qq, med0
;          for kk=0L,nspec-1 do $
;            print, 'esi_echcombspec: Medflux ', kk, med_rto[kk], snr_rto[kk], $
;            FORMAT='(a26,1x,i2,1x,2f7.2)'
          ;; Change made by JFH now uses more sophisticated combining 
          ;; routine
          var   = double(esiobj[all_ordr].var[0:npix-1])
          if keyword_set(USE_OLD) then begin
             x_combspec, esiobj[all_ordr].fx[0:npix-1], var $
                         , fflux, fvar $
                         , novar = esiobj[all_ordr].novar[0:npix-1] $
                         , sky   =  esiobj[all_ordr].sky[0:npix-1]  $
                         , WAVE = esiobj[all_ordr[REFO]].wave[0:npix-1] $
                         , SCALE = med_rto, SNR = snr_rto $
                         , fnovar = fnovar, FSKY = FSKY
          endif else begin
             novar =  double(esiobj[all_ordr].novar[0:npix-1])
             mask = (var GT 0.0)
             inivar  = mask/(var + (var LE 0.0))
             innivar = mask/(novar + (novar LE 0.0))
             inloglam = alog10(esiobj[all_ordr[REFO]].wave[0:npix-1])
             influx = esiobj[all_ordr].fx[0:npix-1]
             insky  =  esiobj[all_ordr].sky[0:npix-1]
             if qq EQ 0 then NOREJ = 1 else NOREJ = 0 ; Order 15 (partial data)
             IF KEYWORD_SET(HIZQSO) THEN BEGIN
                IF qq EQ QQ_MIN THEN BEGIN
                   LAM_MASK_MIN = LAM_MASK_MIN1
                   USE_AVG_SN_WEIGHTS = 1
                ENDIF ELSE IF qq LT QQ_MIN THEN BEGIN
                   INPUT_WEIGHTS1 = mean_SN2_Lya
                   LAM_MASK_MIN = 0.0
                   USE_AVG_SN_WEIGHTS = 0
                ENDIF ELSE BEGIN
                   LAM_MASK_MIN = 0.0
                   USE_AVG_SN_WEIGHTS = 0
                ENDELSE
             ENDIF ELSE BEGIN
                LAM_MASK_MIN = 0.0
                USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS1
             ENDELSE
             long_combspec, influx, inivar, inloglam, insky = insky $
                            , innivar = innivar $
                            , newflux = fflux, newivar = newivar $
                            , in_npoly = npoly $
                            , newnivar = newnivar, newsky = fsky $
                            , newmask = newmask, CHECK = CHK $
                            , MEDSCALE = MEDSCALE $
                            , SN_MIN_MEDSCALE = SN_MIN_MEDSCALE[qq] $
                            , SN_MAX_MEDSCALE = SN_MAX_MEDSCALE[qq] $
                            , NOSHIFT = NOSHIFT, SHIFT_SN = SHIFT_SN[qq] $  
                            , IREF = IREF, NOREJ = norej, MXSHIFT = MXSHIFT $
                            , XCORR_SKY = (qq GT 2) $ ;; For orders < order 13 use sky for xcorr
                            , MEAN_SN2 = MEAN_SN2 $
                            , LAM_MASK_MIN = LAM_MASK_MIN, INPUT_WEIGHTS = INPUT_WEIGHTS1 $
                            , USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS
             IF KEYWORD_SET(HIZQSO) THEN BEGIN
                IF qq EQ QQ_MIN THEN mean_SN2_Lya = mean_SN2
             ENDIF
             
             fvar   = newmask/(newivar + (newivar LE 0.0))
             fnovar = newmask/(newnivar + (newnivar LE 0.0))
             bad = where(newmask EQ 0, nbad)
             IF nbad NE 0 THEN BEGIN
                fvar[bad] = -1.0D
                fnovar[bad] = -1.0D
             ENDIF
          endelse
          echfspec.wave[0:npix-1, qq] = esiobj[all_ordr[REFO]].wave[0:npix-1]
          echfspec.fx[0:npix-1, qq] = fflux[0:npix-1]
          echfspec.var[0:npix-1, qq] = fvar[0:npix-1]
          echfspec.novar[0:npix-1, qq] = fnovar[0:npix-1]
          echfspec.sky[0:npix-1, qq] = fsky[0:npix-1]
          ;stop
;;
          ;x_splot, echfspec.wave[*,qq], esiobj[qq].fx, /block, $
          ;  YTWO=echfspec.fx[*,qq], PSYM_Y2=10
      endelse
      ;; Output
;      esi_echcomb_out, echfspec
  endfor
  

;;;; OUTPUT  ;;;;

  x_wrechfspec, echfspec, outfil 

  close, 56

  print, 'esi_echcombspec:  Output is in ', outfil
  print, 'esi_echcombspec:  All done!'


  return
end
  

