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
;    ORDRS=    - Orders to combine (default: [0L,9L])
;
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

pro esi_echcombspec, esi, obj_id, exp_id, SILENT=silent, LIST=list, $
                    OBJ_NM=OBJ_NM, STD=std, OUTNM=outnm, ORDRS=ordrs

;
  if  N_params() LT 2  and not keyword_set(LIST) then begin 
    print,'Syntax - ' + $
      'esi_echcombspec, esi, obj_id, [exp_id], /SILENT, OBJ_NM=, /PATCH '
    print, '      LIST= [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(ORDRS) then ordrs=[0L,9L]

  if keyword_set( LIST ) then $
    print, 'esi_echcombspec: Combining Obj files from this list -- ', list

; Set exp
  if not keyword_set( LIST ) then begin
      if not keyword_set( STD ) then begin
          allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
                         esi.mode EQ 2 AND esi.obj_id EQ obj_id)
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
    print, 'esi_echcombspec: Loading up the files...'

; Open Obj files
  for q=0L,nexp-1 do begin
      ;; Read
      if not keyword_set( LIST ) then begin
          if x_chkfil(esi[exp[q]].obj_fil+'*') EQ 0 then begin
              print, 'esi_echcombspec: Obj file doesnt exist! Returning..', $
                esi[exp[q]].obj_fil
              return
          endif
          tmp = xmrdfits(esi[exp[q]].obj_fil, 1, STRUCTYP='dblsobjstrct', $
                         /silent)
      endif else begin
          if x_chkfil(files[q]+'*') EQ 0 then begin
              print, 'esi_echcombspec: Obj file doesnt exist! Returning..', $
                files[q]
              return
          endif
          tmp = xmrdfits(files[q], 1, STRUCTYP='dblsobjstrct', /silent)
      endelse
      ;; Grab good obj
      gdobj = where(tmp.obj_id EQ obj_nm)
      ;; Add to total structure
      if q EQ 0 then esiobj = tmp[gdobj] else esiobj = [esiobj, tmp[gdobj]]
  endfor

  nspec = n_elements(esiobj) / 10  ; 10 orders per spec


;; ASCII OUTPUT ;;;

  if obj_id LT 10L then txtfil = 'FSpec/txt/ECHF_0'+$
    string(obj_id,format='(i1)')+'.txt' $
  else txtfil = 'FSpec/txt/ECHF_'+string(obj_id,format='(i2)')+obj_nm+'.txt'
  close, 56
  openw, 56, txtfil

;;; CREATE FINAL 2D ;;;

  if not keyword_set(LIST) then begin
      if keyword_set( OUTNM ) then $
        outfil = 'FSpec/'+outnm+obj_nm+'_ech.fits' $
      else outfil = 'FSpec/'+strtrim(esi[exp[0]].Obj,2)+obj_nm+'_ech.fits' 
  endif else outfil = 'FSpec/tmp_ech.fits' 

  echfspec = { echfspecstrct }

  ;; Copy
  echfspec.nexp = nspec
  ;; Zero ordr
  all_zero = where(esiobj.slit_id EQ 0, nzero)
  if nzero NE nspec then stop

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = esiobj[all_zero[i]].exp

  ;; Other tags
  copy_struct, esiobj[0], echfspec, EXCEPT=["wave","fx","var","npix"]
  for i=0L,nspec-1 do begin
      ipos = strpos(esiobj[all_zero[i]].spec2d_fil, 'esi')
      obj_fil = 'FSpec/Obj_'+strmid(esiobj[all_zero[i]].spec2d_fil, $
                                      ipos)
      echfspec.obj_fil[i] = obj_fil
  endfor

  svmedian = fltarr(nspec)
  svsnr = fltarr(nspec)
  ;; Loop on ORDER (reverse)
  for qq=ordrs[1],ordrs[0],-1 do begin
      all_ordr = where(esiobj.slit_id EQ qq)
      ;; Coadd
      if nspec EQ 1 then begin
          npix = esiobj[all_ordr].npix
          echfspec.wave[0:npix-1,qq] = esiobj[all_ordr].wave[0:npix-1]
          echfspec.fx[0:npix-1,qq] = esiobj[all_ordr].fx[0:npix-1]
          echfspec.var[0:npix-1,qq] = esiobj[all_ordr].var[0:npix-1]
      endif else begin
          npix = 5000L
          ;; Find median ratio of flux and SNR
          med_rto = fltarr(nspec)
          snr_rto = fltarr(nspec)
          med_rto[0] = 1.
          snr_rto[0] = 1.
          gdpix = where(esiobj[all_ordr[0]].var GT 0.)
          med0 = median(esiobj[all_ordr[0]].fx[gdpix])
          if med0 LT 10 then begin
              if qq EQ ordrs[1] then begin
                  ;; Use exposure time ratio
                  for kk=1L,nspec-1 do begin
                      med_rto[kk] = esiobj[all_ordr[0]].exp/ $
                        esiobj[all_ordr[kk]].exp
                      snr_rto[kk] = sqrt(esiobj[all_ordr[0]].exp/ $
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
              for kk=1L,nspec-1 do begin
                  med_rto[kk] = median(esiobj[all_ordr[0]].fx[gdpix]/$
                                       esiobj[all_ordr[kk]].fx[gdpix])
                  snr_rto[kk] = median(esiobj[all_ordr[kk]].fx[gdpix]/ $
                                       sqrt(esiobj[all_ordr[kk]].var[gdpix])/ $
                                       esiobj[all_ordr[0]].fx[gdpix] * $
                                       sqrt(esiobj[all_ordr[0]].var[gdpix]))
              endfor
              ;; Save
              svmedian[*] = med_rto
              svsnr[*] = snr_rto
          endelse

          ;; Set bad pix to -1
          var = double(esiobj[all_ordr].var[0:npix-1])
          for kk=0L,nspec-1 do begin
              a = where(esiobj[all_ordr[kk]].var[0:npix-1] LE 0., na)
              if na GT 0 then var[a,kk] = -1.
          endfor
          print, 'esi_echcombspec: Order ', 15-qq, med0
          for kk=0L,nspec-1 do $
            print, 'esi_echcombspec: Medflux ', kk, med_rto[kk], snr_rto[kk], $
            FORMAT='(a22,1x,i2,1x,2f7.2)'
          x_combspec, esiobj[all_ordr].fx[0:npix-1], var, $
            fflux, fvar, WAVE=esiobj[all_ordr[0]].wave[0:npix-1], $
            SCALE=med_rto, SNR=snr_rto
          echfspec.wave[0:npix-1,qq] = esiobj[all_ordr[0]].wave[0:npix-1]
          echfspec.fx[0:npix-1,qq] = temporary(fflux[0:npix-1])
          echfspec.var[0:npix-1,qq] = temporary(fvar[0:npix-1])
          ;;
;          x_splot, echfspec.wave[*,qq], esiobj[qq].fx, /block, $
;            YTWO=echfspec.fx[*,qq], PSYM_Y2=10
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
  

