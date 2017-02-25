;+ 
; NAME:
; esi_lwdtrcobj   
;     Version 1.0 
;
; PURPOSE:
;    Trace the obj [science only for now]
;
; CALLING SEQUENCE:
;   
;  esi_lwdtrcobj, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with trcobjered light removed
;
; OPTIONAL KEYWORDS:
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdtrcobj, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdtrcobj, esi, obj_id, STD=std, CHK=chk, NCOLL=ncoll, FAINT=FAINT, $
                   NMED=nmed, FILSTD=filstd


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdtrcobj, esi, obj_id, /STD, /CHK, /FAINT, NCOLL=, /STDONLY [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then begin
      if keyword_set( FAINT ) then ncoll = 40L else ncoll = 30L
  endif
  if not keyword_set( YERRMX ) then yerrmx = 1.0
  if not keyword_set( NMED ) then nmed = 15L

; INDX
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
               esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
;;;;;;
;  Find standard star
  if keyword_set( STD ) then begin
      istd = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
                   esi.slit EQ esi[indx[0]].slit AND $
                   strtrim(esi.type,2) EQ 'STD', nstd)
      case nstd of 
          0 : begin
              print, 'esi_lwdtrcobj: No standard star images! Returning..'
              return
          end 
          1 : print, 'esi_lwdtrcobj: Using standard star ', $
            esi[istd].img_root
          else : begin
              print, 'esi_lwdtrcobj: Warning -- Multiple standard star images'
              istd = istd[0]
              print, 'esi_lwdtrcobj: Using first one -- ', esi[istd].img_root
          end
      endcase
  endif else begin
      if not keyword_set( FILSTD ) then stop
      stdfil = filstd
  endelse
          
;;;;;;;;;
; Open Stuff

  ;; Open Slit file
  c_s = esi_slitnm(esi[indx[0]].slit)

  ;; OBJ
  objfil = esi[indx[0]].obj_fil
  objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
  if strtrim(objstr[0].obj_id,2) NE 'a' then stop
  
  ;; STD OBJ
  if not keyword_set( STDFIL ) then stdfil = esi[istd].obj_fil
  stdstr = xmrdfits(stdfil, 1, STRUCTYP='specobjstrct', /silent)
  if strtrim(stdstr[0].obj_id,2) NE 'a' then stop

  ;; Open Image, Variance
  print, 'esi_lwdtrcobj: Opening [combining] images...'
  nimg = n_elements(indx)

  ;; COMBINE as needed
  img = esi_lwdcombimg(esi, indx, VAR=var, IMGINDX=2L, /SKY) 
  sz_img = size(img, /dimensions)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  subimg = fltarr(ncoll,nmed*2+1)
  subvar = fltarr(ncoll,nmed*2+1)
  nfit = sz_img[0]/ncoll + 10
  xfit = dblarr(nfit)
  ysig = dblarr(nfit)
  yfit = dblarr(nfit)
  gdfit = 0L
  
  ;; Zero out subimg
  subimg[*] = 0.
  svoff = 0.

  ;; Trace
  if keyword_set( STD ) or keyword_set( FILSTD) then begin ; STANDARD STAR AS GUIDE
      yoff = objstr[0].ycen - stdstr[0].ycen
      std_trc = stdstr[0].trace[0:sz_img[0]-1] + yoff
      rnd_trc = round(std_trc)
  endif else begin              ; PIN HOLE AS GUIDE
      rnd_trc = round(objstr[0].trace[0:sz_img[0]-1])
  endelse

;  if keyword_set( CHK ) then begin
;      tmp = img
;      trc_msk = lindgen(sz_img[0]) + rnd_trc*sz_img[0]
;      tmp[trc_msk] = -1000
;      xatv, tmp, /block
;  endif

  ;;;;;;;;;;;;;;;;;;;;
  ;; LOOPS

  ;; MIDDLE to END
  istrt = 1700L ;; First 400 columns are crummy
  iend = sz_img[0]-1
  for i=istrt,iend,ncoll do begin
      ;; Create Subimg
      ip = 0L
      for isub = (i-ncoll),i-1 do begin
          subimg[ip,*] = img[isub,rnd_trc[isub]-nmed:rnd_trc[isub]+nmed]
          subvar[ip,*] = var[isub,rnd_trc[isub]-nmed:rnd_trc[isub]+nmed]
          ip = ip+1
      endfor
      ;; Median
      mdn = djs_median(subimg,1)
      mdn_var = djs_median(subvar,1)
      ivar = 1./(mdn_var > 0.)
      ycen = double(nmed) + svoff
      ;; Find Centroid
      for k=0,19 do $
        ycen = trace_fweight(mdn, ycen, 0L, radius=3., xerr=yerr, invvar=ivar)
      ;; Keep the good points
      if yerr LT YERRMX AND yerr GT 0.00000001 then begin
          ysig[gdfit] = yerr
          xfit[gdfit] = (i-ncoll/2)
          yfit[gdfit] = rnd_trc[xfit[gdfit]] + ycen - double(nmed)
          svoff = ycen-double(nmed)
          gdfit = gdfit + 1L
      endif else begin
          ;; Use STD as a guide
          if keyword_set( STD ) then begin
              ysig[gdfit] = 1.0
              xfit[gdfit] = (i-ncoll/2)
              yfit[gdfit] = std_trc[xfit[gdfit]] 
              gdfit = gdfit + 1L
          endif
      endelse
  endfor

  svoff = 0.
  ;; MIDDLE to START
  istrt = 1700L ;; First 400 columns are crummy
  iend = ncoll
  for i=istrt,iend,-ncoll do begin
      ;; Create Subimg
      ip = 0L
      for isub = (i-ncoll),i-1 do begin
          subimg[ip,*] = img[isub,rnd_trc[isub]-nmed:rnd_trc[isub]+nmed]
          subvar[ip,*] = var[isub,rnd_trc[isub]-nmed:rnd_trc[isub]+nmed]
          ip = ip+1
      endfor
      ;; Median
      mdn = djs_median(subimg,1)
      mdn_var = djs_median(subvar,1)
      ivar = 1./(mdn_var > 0.)
      ycen = double(nmed) + svoff
      ;; Find Centroid
      for k=0,19 do $
        ycen = trace_fweight(mdn, ycen, 0L, radius=3., xerr=yerr, invvar=ivar)
      ;; Keep the good points
;      plot, mdn
;      oplot, [ycen,ycen], [-10000,10000]
;      print, ycen, yerr
;      print, i-ncoll/2
;      stop
      if yerr LT YERRMX AND yerr GT 0.00000001 then begin
          ysig[gdfit] = yerr
          xfit[gdfit] = (i-ncoll/2)
          yfit[gdfit] = std_trc[xfit[gdfit]] + ycen - double(nmed)
          svoff = ycen-double(nmed)
          gdfit = gdfit + 1L
      endif else begin
          ;; Use STD as a guide
          if keyword_set( STD ) then begin
              ysig[gdfit] = 1.0
              xfit[gdfit] = (i-ncoll/2)
              yfit[gdfit] = rnd_trc[xfit[gdfit]] 
              gdfit = gdfit + 1L
          endif
      endelse
  endfor

  ;; FIT
  trc_fit = x_setfitstrct(NITER=2L, NORD=4L, FLGREJ=1L, HSIG=3., LSIG=3., $
                          FUNC='POLY')
  new_trc = x_fitrej(xfit[0:gdfit-1], yfit[0:gdfit-1], SIG=ysig[0:gdfit-1], $
                     FITSTR=trc_fit)
  
  print, 'esi_lwdtrcobj: RMS = ', trc_fit.rms
  ;; objstr (Take trace to iend)
  objstr[0].trace[0:sz_img[0]-1] = x_calcfit(findgen(sz_img[0]-1+1), FITSTR=trc_fit)
  if keyword_set( CHK ) then begin
      x_splot, xfit[0:gdfit-1], yfit[0:gdfit-1], PSYM1=1, $
        YTWO=new_trc, YTHR=ysig, PSYM3=2, /block
  endif

  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      trc_msk = lindgen(sz_img[0]) + $
        round(objstr[0].trace[0:sz_img[0]-1])*sz_img[0]
      tmp[trc_msk] = -1000
      xatv, tmp, /block
  endif

; OUTPUT

  mwrfits, objstr, objfil, /create, /silent

  ;; For multiple images
  for ii=1L,nimg-1 do begin
      objfil = esi[indx[ii]].obj_fil
      tmpstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      tmpstr[0].trace = objstr[0].trace
      mwrfits, tmpstr, objfil, /create, /silent
  endfor
  
  print, 'esi_lwdtrcobj: All done!'

  return
end
              
      
      
