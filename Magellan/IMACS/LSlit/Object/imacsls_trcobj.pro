;+ 
; NAME:
; imacsls_trcobj   
;     Version 1.1
;
; PURPOSE:
;    Trace the object along the CCD.
;
; CALLING SEQUENCE:
;  imacsls_trcobj, imacsls, setup, obj_id, side, [exp_id], /STD, /FAINT, /CHK, NCOLL=
;
; INPUTS:
;   imacsls  -  IMACS long slit structure
;   setup    -  Setup ID value
;   side     -  Specific CCD (1='blue', 2='red')
;   obj_id   -  Object ID value
;   [exp]    -  Expsoure index values
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK    - Show final trace
;   /STD    - Use standard star as a crutch
;   /FAINT  - Faint object; sum up more rows (40) to search for flux 
;   NCOLL=  - Set number of rows to sum by hand (default: 25)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_trcobj, imacsls, 1L, [0L], /CHK, /STD
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Added imacsls_trcobj_sngl
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_trcobj, imacsls, setup, obj_id, side, exp, indx=indx, $
                    CHK=chk, NCOLL=ncoll, FAINT=FAINT, $
                    STDFIL=stdfil, GUIDE=guide, DEBUG=debug, STAR=star
;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'imacsls_trcobj, imacsls, setup, obj_id, side, [exp] ' + $
        '/CHK, /FAINT, NCOLL=, /GUIDE /DEBUG [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then begin
      if keyword_set( FAINT ) then ncoll = 30L else ncoll = 15L
  endif
  if keyword_set( STAR ) then ncoll = 5L
  if not keyword_set( STAR ) then type = 'OBJ' else type = 'STD'
  if not keyword_set( XERRMX ) then xerrmx = 0.8

; INDX
  if NOT keyword_set(indx) then $
    indx = where(imacsls.flg_anly NE 0 AND imacsls.side EQ side AND $
                 imacsls.obj_id EQ obj_id AND imacsls.setup EQ setup AND $
                 strtrim(imacsls.type,2) EQ type, nindx) $
  else nindx = n_elements(indx)

  if nindx GT 1 then print, 'imacsls_trcobj: Combining multiple frames. '

  if keyword_set( exp_id ) then indx = indx[exp_id]

          
;;;;;;;;;
; Open Stuff

  ;; OBJ
  objfil = imacsls[indx[0]].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'imacsls_trcobj: Obj file does not exist or obj_fil tag not ' + $
        'set! Check it..'
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
  
  ;; Open Image, Variance
  print, 'imacsls_trcobj: Opening [combining] images...'
  nimg = n_elements(indx)

  ;; COMBINE as needed
  if keyword_set( STAR ) then begin
      img = xmrdfits(imacsls[indx[0]].img_final, /silent)
      ivar = xmrdfits(imacsls[indx[0]].img_final, 1, /silent)
  endif else begin
      img = imacsls_combimg(imacsls, indx, IVAR=ivar, IMGINDX=2L, /SKY) 
  endelse
  sz_img = size(img, /dimensions)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  subimg = fltarr(21,ncoll)
  subivar = fltarr(21,ncoll)
  nfit = sz_img[1]/ncoll + 1
  xfit = dblarr(nfit)
  xsig = dblarr(nfit)
  yfit = dblarr(nfit)
  
  ;; Index
  ;; Zero out subimg
  subimg[*] = 0.
  svoff = 0.
  ;; Crude guide Trace
  crude_trc = objstr.trace[0:sz_img[1]-1]
  rnd_trc = round(crude_trc)
  flg_guide = 1
  ;; Set looping range
  jmid = sz_img[1]/2
  jstrt = ncoll
  jend = sz_img[1]-1
  gdfit = 0L

  ;;;;;;;;;;;;;;;;;;;;
  ;; LOOP
  ;; MIDDLE to TOP
  svoff = 0.
  for j=jmid,jend,ncoll do begin
      ;; Create Subimg
      jp = 0L
      for jsub = (j-ncoll),j-1 do begin
          subimg[*,jp] = 0.
          subivar[*,jp] = 0.
          i1 = (rnd_trc[jsub] - 10L) > 0L
          i2 = 0L > (rnd_trc[jsub] + 10L) < (sz_img[0] - 1)
          if i2 GT i1 then begin
              s1 = (10L - rnd_trc[jsub]) > 0L
              s2 = i2-i1+s1
              subimg[s1:s2,jp] = img[i1:i2,jsub]
              subivar[s1:s2,jp] = ivar[i1:i2,jsub]
          endif
          jp = jp+1
      endfor
      ;; Median
      mdn = djs_median(subimg,2)
      mdn_var = djs_median(subivar,2)
      mivar = mdn_var > 0.
      xcen = 10.d + svoff
      ;; Find Centroid
      for k=0,19 do $
        xcen = trace_fweight(mdn, xcen, 0L, radius=3., xerr=xerr, $
                             invvar=mivar)
      ;; Keep the good points
      if xerr LT XERRMX AND xerr GT 0.000001 then begin
          xsig[gdfit] = xerr
          yfit[gdfit] = (j-ncoll/2)
          xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - 10.
          svoff = xcen-10.
          gdfit = gdfit + 1L
      endif else begin
          ;; GUIDE
          if keyword_set(flg_guide) then begin
              xsig[gdfit] = 1.0
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] 
              gdfit = gdfit + 1L
          endif 
      endelse
;          plot, mdn, psym=10
;          oplot, [xcen,xcen], [-1000,1000]
;          print, xcen, xerr
;          stop
  endfor

  svoff = 0.
  ;; MIDDLE to START
  for j=jmid,jstrt,-ncoll do begin
      ;; Create Subimg
      jp = 0L
      for jsub = (j-ncoll),j-1 do begin
          subimg[*,jp] = 0.
          subivar[*,jp] = 0.
          i1 = (rnd_trc[jsub] - 10L) > 0L
          i2 = (rnd_trc[jsub] + 10L) < (sz_img[0] - 1)
          if i2 GT i1 then begin
              s1 = (10L - rnd_trc[jsub]) > 0L
              s2 = i2-i1+s1
              subimg[s1:s2,jp] = img[i1:i2,jsub]
              subivar[s1:s2,jp] = ivar[i1:i2,jsub]
          endif
          jp = jp+1
      endfor
      ;; Median
      mdn = djs_median(subimg,2)
      mdn_var = djs_median(subivar,2)
      mivar = mdn_var > 0.
      xcen = 10.d + svoff
      ;; Find Centroid
      for k=0,19 do $
        xcen = trace_fweight(mdn, xcen, 0L, radius=3., xerr=xerr, $
                             invvar=mivar)
      ;; Keep the good points
      if xerr LT XERRMX AND xerr GT 0.000001 then begin
          xsig[gdfit] = xerr
          yfit[gdfit] = (j-ncoll/2)
          xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - 10.
          svoff = xcen-10.
          gdfit = gdfit + 1L
      endif else begin
          ;; GUIDE
          if keyword_set(flg_guide) then begin
              xsig[gdfit] = 1.0
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] 
              gdfit = gdfit + 1L
          endif
      endelse
;          plot, mdn, psym=10
;          oplot, [xcen,xcen], [-1000,1000]
;          print, xcen, xerr
;          stop
  endfor

  ;; FIT
  trc_fit = x_setfitstrct(NITER=2L, NORD=9L, FLGREJ=1L, HSIG=5., LSIG=5., $
                          FUNC='POLY')
  new_trc = x_fitrej(yfit[0:gdfit-1], xfit[0:gdfit-1], SIG=xsig[0:gdfit-1], $
                     FITSTR=trc_fit)
  
  print, 'imacsls_trcobj: RMS = ', trc_fit.rms

  ;; objstr (Take trace to jend)
  objstr.trace[0:jend] = x_calcfit(findgen(jend+1), FITSTR=trc_fit)
  if keyword_set( DEBUG ) then begin
      srt = sort(yfit[0:gdfit-1])
      x_splot, yfit[srt], xfit[srt], PSYM1=1, $
            YTWO=new_trc[srt], YTHR=xsig[srt], PSYM3=2, /block
  endif
  

  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      trc_msk = round(objstr.trace[0:sz_img[1]-1]) $
        + lindgen(sz_img[1])*sz_img[0]
      tmp[trc_msk] = -1000
      xatv, tmp, /block, min=-80, max=80
  endif

; OUTPUT

  print, 'imacsls_trcobj: Updating trace in ', objfil
  mwrfits, objstr, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
  ;; For multiple images
  for ii=1L,nimg-1 do begin
      objfil = imacsls[indx[ii]].obj_fil
      print, 'imacsls_trcobj: Updating trace in ', objfil
      tmpstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
      tmpstr.trace = objstr.trace
      mwrfits, tmpstr, objfil, /silent, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'imacsls_trcobj: All done!'

  return
end
              
