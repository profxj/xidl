;+ 
; NAME:
; esi_lwdtrcstd   
;     Version 1.0
;
; PURPOSE:
;    Trace a standard star in each ordrer
;
; CALLING SEQUENCE:
;   
;  esi_lwdtrcstd, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with trcstdered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdtrcstd, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdtrcstd, esi, slit, NOFND=nofnd, NOSKY=nosky, CHK=chk 

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdtrcstd, esi, slit, /NOFND, /NOSKY, /CHK [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then ncoll = 5L
  if not keyword_set( NMED ) then nmed = 30L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'

;;;;;;
;  Find standard star
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
               esi.slit EQ slit AND esi.type EQ 'STD', nindx)
  case nindx of 
      0 : begin
          print, 'esi_lwdtrcstrd: No standard star images! Returning..'
          return
      end 
      1 : print, 'esi_lwdtrcstd: Tracing standard star image --- ', $
        esi[indx].img_root
      else : begin
          print, 'esi_lwdtrcstd: Warning -- Multiple standard star images'
          indx = indx[0]
          print, 'esi_lwdtrcstd: Taking first one -- ', esi[indx].img_root
      end
  endcase
          
;;;;;;;;;
; Process
  esi_lwdproc, esi, indx, /INDEX

;;;;;;;;;
; FIND Obj + Create Obj Structure
  if not keyword_set( NOFND ) then $
    esi_lwdfndobj, esi, indx, /STD, /AUTO, /CLOBBER

;;;;;;;;;
; Sky subtract
  if not keyword_set( NOSKY ) then $ 
    esi_lwdskysub, esi, indx, /STD

  
;;;;;;;;;
; Open Stuff

  ;; Open Image, Variance
  print, 'esi_lwdtrcstd: Opening images...'
  img = xmrdfits(esi[indx].img_final, 2, /silent)
  sz_img = size(img, /dimensions)
  var = xmrdfits(esi[indx].img_final, 1, /silent)

  ;; OBJ
  objfil = esi[indx].obj_fil
  objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  subimg = fltarr(ncoll,nmed*2+1)
  subvar = fltarr(ncoll,nmed*2+1)
  nfit = sz_img[0]/ncoll
  xfit = dblarr(nfit)
  ysig = dblarr(nfit)
  yfit = dblarr(nfit)
  gdfit = 0L
  
  ;; Zero out subimg
  subimg[*] = 0.
  svoff = 0.
  ;; Trace
  rnd_trc = round(objstr[0].trace[0:sz_img[0]-1])
  ;; CHK
;  if keyword_set( CHK ) then begin
;      trc_msk = lindgen(sz_img[0]) + rnd_trc*sz_img[0]
;      tmp = img
;      tmp[trc_msk] = -1
;      xatv, tmp, /block
;      stop
;  endif
  ;;;;;;;;;;;;;;;;;;;;
  ;; LOOP
  istrt = ncoll + 410L  ;; First 400 columns are crummy
  iend = sz_img[0]
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
        ycen = trace_fweight(mdn, ycen, 0L, radius=7., xerr=yerr, invvar=ivar)
      ;; Keep the good points
      if yerr LT 0.2 AND yerr GT 0.00000001 then begin
          ysig[gdfit] = yerr
          xfit[gdfit] = (i-ncoll/2)
          yfit[gdfit] = rnd_trc[xfit[gdfit]] + ycen - double(nmed)
          svoff = ycen-double(nmed)
          gdfit = gdfit + 1L
      endif
;      plot, mdn
;      oplot, [ycen,ycen], [-10000,10000]
;      print, ycen, yerr
  endfor

  ;; FIT
  trc_fit = x_setfitstrct(NITER=2L, NORD=4L, FLGREJ=1L, HSIG=3., LSIG=3., $
                          FUNC='POLY')
  new_trc = x_fitrej(xfit[0:gdfit-1], yfit[0:gdfit-1], SIG=ysig[0:gdfit-1], $
                     FITSTR=trc_fit)
  
  ;; CHK
  if keyword_set( CHK ) then $
    x_splot, xfit[0:gdfit-1], yfit[0:gdfit-1], YTWO=new_trc, /block

  print, 'esi_lwdtrcstd: RMS = ', trc_fit.rms
  ;; objstr
  objstr[0].trace[0:sz_img[0]-1] = x_calcfit(findgen(sz_img[0]), $
                                                  FITSTR=trc_fit)
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

  c_s = esi_slitnm(esi[indx[0]].slit)
  std_trc = 'Extract/STD_LWD'+c_s+'_TRC.fits'
  mwrfits, objstr[0].trace[0:sz_img[0]-1], std_trc, /create, /silent
  
  print, 'esi_lwdtrcstd: All done!'

  return
end
