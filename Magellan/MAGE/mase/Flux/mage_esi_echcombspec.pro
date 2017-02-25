;+ 
; NAME:
; mage_echcombspec
;    Version 1.0
;
; PURPOSE:
;   Combines multiple exposures of the same obj
;    Must be run even on an object with a single exposure
;
; CALLING SEQUENCE:
;   
;   mage_echcombspec, allframes, fspec, orders, keyindex
;
;INPUTS:
;   allobj   -- Array of echspec structures
;   ordrs    -- Orders to coadd 
;   keyindx  -- Exposure to serve as the fiducial for stats
;
; RETURNS:
;
; OUTPUTS:
;   finspec  --  echfspec structure containing the combined spectrum
;;
; OPTIONAL KEYWORDS:
;    /SILENT   - No text output
;    LIST=     - Use an input list to combine (required for multiple
;                night obs).  Output is in 'FSpec/tmp_ech.fits'.
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file
;    ORDRS=    - Orders to combine (default: [0L,9L])
;    CHK       - Various QA for call to long_combspec for combining
;                spectra. 
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
;   05-Mar-2009 Ported to MAGE by JJBq
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

pro mage_esi_echcombspec, allobj, echfspec, SILENT=silent $
                     , LIST=list $
                     , OBJ_NM = OBJ_NM1, STD = std, OUTNM = outnm $
                     , ORDRS = ordrs , REFO = refo, CHK = CHK $
                     , NOSHIFT = NOSHIFT, IREF = IREF

;
  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
      'x_echcombspec, allobj, finspec, ordrs, /SILENT, ' + $
      'OBJ_NM=, /PATCH, /CHK, /MCHK  [v1.1]'
    return
  endif 
  
;;I know all the objects already, and the final output structure has
;;been made already


;;Bug fix added by JFH, prevents return of obj_id =99 value which breaks 
;; later calls to esi routines. 
IF n_elements(obj_id1) GT 0 THEN obj_id = obj_id1 
;  Optional Keywords
  if not keyword_set(ORDRS) then ordrs=[6L,20L]  ;removed the orders from the call
  if not keyword_set(REFO) then refo=0L

  if keyword_set( LIST ) then $
    print, 'esi_echcombspec: Combining Obj files from this list -- ', list

; Set exp
;;   if not keyword_set( LIST ) then begin
;;       if not keyword_set( STD ) then begin
;;           allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
;;                          esi.mode EQ 2 AND esi.obj_id EQ obj_id)
;;           if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
;;           nexp = n_elements(exp)
;;       endif else begin  ;; STD
;;           exp = obj_id[0]
;;           nexp = 1
;;           obj_id = 99L
;;       endelse
;;   endif else begin
;;       obj_id = 99L
;;       readcol, list, files, FORMAT='A'
;;       nexp = n_elements(files)
;;    endelse

  nexp = echfspec.nexp
  nspec = echfspec.nexp*15.

  IF KEYWORD_SET(OBJ_NM1) THEN BEGIN
     IF n_elements(obj_nm1) EQ 1 THEN obj_nm =  replicate(obj_nm1, nexp) $
     ELSE IF n_elements(obj_nm1) EQ nexp THEN obj_nm = obj_nm1 $
     ELSE message, 'Invalid size for obj_nm'
  ENDIF ELSE obj_nm = replicate('a', nexp)

;;;;;;;;;;
; Open Files

  ;; if not keyword_set( SILENT ) AND NOT KEYWORD_SET(LIST) then BEGIN
;;       print, 'esi_echcombspec: Loading up the files...'
;;       forprint, esi[exp].OBJ_FIL, textout = 1
;;   ENDIF

;;; CREATE FINAL 2D ;;;

 ;;  if not keyword_set(LIST) then begin
;;       if keyword_set( OUTNM ) then $
;;          outfil = 'FSpec/'+outnm+obj_nm[0]+'_ech.fits' $
;;       else outfil = 'FSpec/'+strtrim(esi[exp[0]].Obj, 2)+obj_nm[0]+'_ech.fits' 
;;    endif else BEGIN
;;       IF KEYWORD_SET(OUTNM)THEN  outfil = 'FSpec/'+outnm+obj_nm[0]+'_ech.fits' $
;;      ELSE outfil = 'FSpec/tmp_ech.fits' 
;;   ENDELSE

 ; echfspec = { echfspecstrct }

  ;; Copy
  ;; echfspec.nexp = nspec
  
  nspec=nexp
  
  ;; Zero ordr
  all_zero = where(allobj.order EQ 6, nzero)
  if nzero NE nspec then stop

  ;; Set texp
  for i=0L,nspec-1 do echfspec.texp[i] = allobj[all_zero[i]].exp

  ;; Other tags
  copy_struct, allobj[0], echfspec $
               , EXCEPT = ["wave", "fx", "var", "novar", "sky", "npix"]
  ;; Loop on ORDER (reverse)
  for qq=ordrs[1],ordrs[0],-1 do begin
      all_ordr = where(allobj.order EQ 20-qq)
      ;fill in order numbers
      echfspec.phys_ordr[qq] = qq
      npix = allobj[all_ordr[REFO]].npix
      ;; Coadd
      if nspec EQ 1 then begin 
         echfspec.wave[0:npix-1,qq] = allobj[all_ordr].wave[0:npix-1]
          echfspec.fx[0:npix-1,qq] = allobj[all_ordr].flux[0:npix-1]
          echfspec.var[0:npix-1,qq] = allobj[all_ordr].sig[0:npix-1]^2
          echfspec.novar[0:npix-1, qq] = allobj[all_ordr].nosig[0:npix-1]^2
          echfspec.sky[0:npix-1, qq] = allobj[all_ordr].sky[0:npix-1]
       endif else begin
         print, 'mage_echcombspec: Order ', qq
         var   = double(allobj[all_ordr].sig[0:npix-1]^2)
         novar =  double(allobj[all_ordr].nosig[0:npix-1]^2)
         mask = (var GT 0.0)
         inivar  = mask/(var + (var EQ 0.0))
         innivar = mask/(novar + (novar EQ 0.0))
         inloglam = alog10(allobj[all_ordr[REFO]].wave[0:npix-1])
         influx = allobj[all_ordr].flux[0:npix-1]
         insky  =  allobj[all_ordr].sky[0:npix-1]
         ;;??????? BUG  ????
         ;; Currently we are not allowing for any shifts between
         ;; exposures. But once wavelengths are fixed we may want to 
         ;; allow shifts so that we can corect flexure. 
         NOSHIFT=1
         
         long_combspec, influx, inivar, inloglam, insky = insky $
                         , innivar = innivar $
                         , newflux = fflux, newivar = newivar $
                         , newnivar = newnivar, newsky = fsky $
                         , newmask = newmask, CHECK = CHK $
                         , NOSHIFT = NOSHIFT, IREF = IREF
          fvar   = newmask/(newivar + (newivar EQ 0.0))
          fnovar = newmask/(newnivar + (newnivar EQ 0.0))
          bad = where(newmask EQ 0, nbad)
          IF nbad NE 0 THEN BEGIN
              fvar[bad] = 0.0D
              fnovar[bad] = 0.0D
          ENDIF
          echfspec.wave[0:npix-1, qq] = allobj[all_ordr[REFO]].wave[0:npix-1]
          echfspec.fx[0:npix-1, qq] = fflux[0:npix-1]
          echfspec.var[0:npix-1, qq] = fvar[0:npix-1]
          echfspec.novar[0:npix-1, qq] = fnovar[0:npix-1]
          echfspec.sky[0:npix-1, qq] = fsky[0:npix-1]
      endelse
   endfor
  ;;;; OUTPUT  ;;;;

  ;x_wrechfspec, echfspec, outfil 
  ;close, 56

  print, 'mage_esi_echcombspec:  All done!'


  return
end
  

