;+ 
; NAME:
; esi_echproc   
;     Version 1.2
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;   
;  esi_echproc, esi, indx, /IFLAT, /REDDOV
;
; INPUTS:
;   esi     -  ESI structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /SVOV    - Save OV files
;  /REDOOV  - Redo OV subtraction
;  /CLOBBER - Clobber existing image
;  /SUBSCAT - Subtract scattered light from data image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echproc, esi, [20L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echproc, esi, indx, REDOOV=redoov, CLOBBER=clobber, SUBSCAT=subscat, $
                 SVOV=svov, FLATFIL=flatfil, BIASFIL=biasfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echproc, esi, indx, /SVOV, /REDOOV, /CLOBBER, /SUBSCAT [v1.2]'
      return
  endif 
  
;  Optional Keywords

; NIMG
  nimg = n_elements(indx)

  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      for qq=0L,nimg-1 do begin
          ;; Chk for outfil
          outfil = 'Final/f_'+esi[indx[qq]].img_root
          a = findfile(outfil+'*', count=na)
          if na NE 0 then begin
              print, 'esi_echproc: File ', outfil, ' exists, continuing...'
              esi[indx[qq]].img_final = outfil
              continue
          endif
          flg_proc = 1
      endfor
  endif else flg_proc = 1

  if flg_proc EQ 0 then begin
      print, 'esi_echproc: No files to process. Returning!'
      return
  endif

; Slit
  slit = esi[indx[0]].slit  ; Takes first value so only 1 slit width at a time
  c_s = esi_slitnm(slit)

; BIAS sub if necessary
  if not keyword_set( REDOOV ) then bias = where(esi[indx].flg_ov EQ 0, nbias) $
  else begin
      nbias = n_elements(indx)
      bias = lindgen(nbias)
  endelse
  if nbias NE 0 then begin
      print, 'esi_echproc: Bias subtracting'
      esi_subbias, esi, indx[bias], FORCE=REDOOV, BIASFIL=biasfil
  endif

; Scatter subtract (or just do gain)
  print, 'esi_echproc: Fixing gain'
  if keyword_set(SUBSCAT ) then print, 'esi_echproc: Scatter subtracting'
  esi_echsubscat, esi, indx, GAINONLY=(not keyword_set(SUBSCAT)), $
    FLATFIL=flatfil

; Open Flat

  if not keyword_set( FLATFIL ) then $
    flatfil = 'Flats/FlatECH'+c_s+'N.fits'
  if x_chkfil(flatfil+'*') EQ 0 then begin
      print, 'esi_echproc: Flat file not found', flatfil
      return
  endif
  print, 'esi_echproc: Loading the flat: ', flatfil
  img_flat = xmrdfits(flatfil, 0, fhead, /silent)
  esi[indx].flat_fil = flatfil ; Kludge

  ;; Zero out bad pix
  esi_echbadpix, img_flat
  gdflt = where(img_flat GT 0.)
  scatt = sxpar(fhead,'SCATTER')
  norm = sxpar(fhead,'NORM')
  if scatt NE 1 OR norm NE 1 then begin
      print, 'esi_echproc: Flat not processed! Returning...'
      return
  endif
          

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  for qq=0L,nimg-1 do begin
      ;; Outfil
      outfil = 'Final/f_'+esi[indx[qq]].img_root
      ;; Read IMG
      img_fil = 'OV/ov_'+esi[indx[qq]].img_root
      if x_chkfil(img_fil+'*') EQ 0 then begin
          print, 'esi_echproc: Image file does not exist!', img_fil
          return
      endif
      img = xmrdfits(img_fil, 0, head, /silent)
      sz_img = size(img, /dimensions)
  
      ;; Create VAR image
      var = ((img>1.)*esi[indx[qq]].gain + esi[indx[qq]].readno^2) 

      ;; Flatten
      print, 'esi_echproc: Flattening'
      img_nrm = x_normspec(img, img_flat, var, var_nrm, PIX=gdflt)
      delvarx, img, var

      ;; Add in gain
      fin_img = temporary(img_nrm)*esi[indx[qq]].gain
      
      ;; OUTPUT
      print, 'esi_echproc: Writing... ', outfil
      sxaddpar, head, 'FLAT', 'T', flatfil
      esi[indx[qq]].img_final = outfil
      mwrfits, fin_img, outfil, head, /create, /silent
      mwrfits, var_nrm, outfil, /silent
      ;; COMPRESS
      print, 'esi_echproc: Compressing... '
      spawn, 'gzip -f '+outfil
      ;; Release
      delvarx, fin_img, var_nrm, img, var, img_nrm

  endfor

  ;; DEL OV
  if not keyword_set( SVOV ) then esi_delov, esi, indx

  return
end


