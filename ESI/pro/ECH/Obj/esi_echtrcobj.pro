;+ 
; NAME:
; esi_echtrcobj   
;     Version 1.1
;
; PURPOSE:
;    Trace the object through each order of the image.  The program
;    uses a standard star as a crutch (/STD) through regions where the object
;    has very low flux.  The trace is written into the object
;    structure (e.g. Extract/Obj_esi0024.fits)
;
; CALLING SEQUENCE:
;   
;  esi_echtrcobj, esi, obj_id, [exp_id], /STD, /FAINT, /CHK, NCOLL=
;
; INPUTS:
;   esi     -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp_id]   -  Exposure frames (e.g. [0L, 1L])
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
;   esi_echtrcobj, esi, 1L, [0L], /CHK, /STD
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtrcobj, esi, obj_id, exp_id, ORDRS=ordrs, $
                   USESTD=usestd, CHK=chk, NCOLL=ncoll, FAINT=FAINT, $
                   STDFIL=stdfil, SEDG_FIL=sedg_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtrcobj, esi, indx, /USESTD, /CHK, /FAINT, NCOLL=, [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then begin
      if keyword_set( FAINT ) then ncoll = 40L else ncoll = 25L
  endif
  if not keyword_set( XERRMX ) then xerrmx = 0.8
  if not keyword_set(ORDRS) then ordrs=[0L,9L]
  if not keyword_set(USESTD) and not keyword_set(STDFIL) then stop

; INDX
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
  if indx[0] EQ -1 then begin
      print, 'esi_echtrcobj: Found no obj!', obj_id
      return
  endif

  if keyword_set( exp_id ) then indx = indx[exp_id]

;;;;;;
;  Find standard star
  if keyword_set( USESTD ) then begin
      istd = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.slit EQ esi[indx[0]].slit AND $
                   strtrim(esi.type,2) EQ 'STD', nstd)
      case nstd of 
          0 : begin
              print, 'esi_echtrcobj: No standard star images! Returning..'
              return
          end 
          1 : print, 'esi_echtrcobj: Tracing standard star image ', $
            esi[istd].img_root
          else : begin
              print, 'esi_echtrcobj: Warning -- Multiple standard star images'
              istd = istd[0]
              print, 'esi_echtrcobj: Taking first one -- ', esi[istd].img_root
          end
      endcase
  endif
          
;;;;;;;;;
; Open Stuff

  ;; Open Slit file
  c_s = esi_slitnm(esi[indx[0]].slit)
  if not keyword_set( SEDG_FIL ) then $
    sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
  if x_chkfil(sedg_fil+'*') EQ 0 then begin
      print, 'esi_echtrcobj: Slit edge file doesnt exist: ', sedg_fil
      return
  endif
  print, 'esi_echtrcobj: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)

  ;; OBJ
  objfil = esi[indx[0]].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'esi_echtrcobj: Obj file does not exist or obj_fil tag not ' + $
        'set! Check it..'
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
  
  ;; STD OBJ
  if not keyword_set( STDFIL ) then begin
      if keyword_set(USESTD) then stdfil = esi[istd].obj_fil else stop
  endif
  if x_chkfil(stdfil+'*') EQ 0 then begin
      print, 'esi_echtrcobj: STD Obj file does not exist or obj_fil tag not ' + $
        'set! Check it..', stdfil
      return
  endif
  stdstr = xmrdfits(stdfil, 1, STRUCTYP='dblsobjstrct', /silent)

  ;; Open Image, Variance
  print, 'esi_echtrcobj: Opening [combining] images...'
  nimg = n_elements(indx)
  ;; COMBINE as needed
  img = esi_echcombimg(esi, indx, VAR=var, IMGINDX=2L) 
  sz_img = size(img, /dimensions)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  subimg = fltarr(21,ncoll)
  subvar = fltarr(21,ncoll)
  nfit = sz_img[1]/ncoll + 1
  xfit = dblarr(nfit)
  xsig = dblarr(nfit)
  yfit = dblarr(nfit)
  
  for qq=ordrs[0],ordrs[1] do begin
      print, 'esi_echtrcobj: Tracing order ', string(15L-qq, FORMAT='(i3)')
      ;; Zero out subimg
      subimg[*] = 0.
      svoff = 0.
      ;; Trace
      if keyword_set( STDFIL ) then begin  ; STANDARD STAR AS GUIDE
          jdum = objstr[qq].xcen
;          std_frac = (stdstr[qq].trace[jdum]-slit_edg[jdum,qq,0]) $
;            /(slit_edg[jdum,qq,1]-slit_edg[jdum,qq,0])
;          obj_frac = (objstr[qq].ycen-slit_edg[jdum,qq,0]) $
;            /(slit_edg[jdum,qq,1]-slit_edg[jdum,qq,0])
;          big_frac = obj_frac/std_frac
          str_off = objstr[qq].ycen - stdstr[qq].trace[jdum] 
;          rnd_trc = round(big_frac*stdstr[qq].trace[0:sz_img[1]-1] + $
;                          (1.-big_frac)*slit_edg[0:sz_img[1]-1,qq,0])
          rnd_trc = round(stdstr[qq].trace[0:sz_img[1]-1] + str_off)
      endif else begin  ; PIN HOLE AS GUIDE
          rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
          stop
      endelse
      ;; Order 14 offset!
;      if qq EQ 1 AND not keyword_set( NOFIX14 ) $
;        AND not keyword_set( STD ) then begin
;          print, 'esi_echtrcobj: Applying the Fix to Order 14'
;          rnd_trc[0:500] = rnd_trc[0:500] - 7L
;      endif
      ;; CHK
;      if keyword_set( CHK ) then begin
;          rnd_trc2 = round(stdstr[qq].trace[0:sz_img[1]-1])
;          trc_msk = rnd_trc2 + lindgen(sz_img[1])*sz_img[0]
;          tmp = img
;          tmp[trc_msk] = -10000
;          xatv, tmp, /block
;          stop
;      endif
      ;; Set ooping range
      jmid = 2048L
      case qq of 
          0: begin
              jstrt =1500L
              jend = 3800L
          end
          9: begin
              jstrt = ncoll + 3L
              jend = 2170L
          end
          else: begin
              jstrt = ncoll + 3L  ;; First 3 rows are crummy
              jend = sz_img[1]-1
          end
      endcase
      gdfit = 0L

      ;;;;;;;;;;;;;;;;;;;;
      ;; LOOP
      ;; MIDDLE to TOP
      for j=jmid,jend,ncoll do begin
          ;; Create Subimg
          jp = 0L
          for jsub = (j-ncoll),j-1 do begin
              subimg[*,jp] = 0.
              subvar[*,jp] = 0.
              i1 = (rnd_trc[jsub] - 10L) > 0L
              i2 = (rnd_trc[jsub] + 10L) < (sz_img[0] - 1)
              s1 = (10L - rnd_trc[jsub]) > 0L
              s2 = i2-i1+s1
              subimg[s1:s2,jp] = img[i1:i2,jsub]
              subvar[s1:s2,jp] = var[i1:i2,jsub]
              jp = jp+1
          endfor
          ;; Median
          mdn = djs_median(subimg,2)
          mdn_var = djs_median(subvar,2)
          ivar = 1./(mdn_var > 0.)
          xcen = 10.d + svoff
          ;; Find Centroid
          for k=0,19 do $
            xcen = trace_fweight(mdn, xcen, 0L, radius=3., xerr=xerr, $
                                 invvar=ivar)
          ;; Keep the good points
          if xerr LT XERRMX AND xerr GT 0.000001 then begin
              xsig[gdfit] = xerr
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - 10.
              svoff = xcen-10.
              gdfit = gdfit + 1L
          endif else begin
              ;; Use STD as a guide
              if keyword_set( STDFIL ) then begin
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
              subvar[*,jp] = 0.
              i1 = (rnd_trc[jsub] - 10L) > 0L
              i2 = (rnd_trc[jsub] + 10L) < (sz_img[0] - 1)
              s1 = (10L - rnd_trc[jsub]) > 0L
              s2 = i2-i1+s1
              subimg[s1:s2,jp] = img[i1:i2,jsub]
              subvar[s1:s2,jp] = var[i1:i2,jsub]
              jp = jp+1
          endfor
          ;; Median
          mdn = djs_median(subimg,2)
          mdn_var = djs_median(subvar,2)
          ivar = 1./(mdn_var > 0.)
          xcen = 10.d + svoff
          ;; Find Centroid
          for k=0,19 do $
            xcen = trace_fweight(mdn, xcen, 0L, radius=3., xerr=xerr, $
                                 invvar=ivar)
          ;; Keep the good points
          if xerr LT XERRMX AND xerr GT 0.000001 then begin
              xsig[gdfit] = xerr
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - 10.
              svoff = xcen-10.
              gdfit = gdfit + 1L
          endif else begin
              ;; Use STD as a guide
              if keyword_set( STDFIl ) then begin
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
      
      print, 'esi_echtrcobj: RMS = ', trc_fit.rms
      ;; objstr (Take trace to jend)
      objstr[qq].trace[0:jend] = x_calcfit(findgen(jend+1), FITSTR=trc_fit)
;      if keyword_set( CHK ) then begin
;          srt = sort(yfit[0:gdfit-1])
;          x_splot, yfit[srt], xfit[srt], PSYM1=1, $
;            YTWO=new_trc[srt], YTHR=xsig[srt], PSYM_Y3=2, /block
;      endif

  endfor

  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      for qq=ordrs[0],ordrs[1] do begin
          trc_msk = round(objstr[qq].trace[0:sz_img[1]-1]) $
            + lindgen(sz_img[1])*sz_img[0]
          tmp[trc_msk] = -1000
      endfor
      xatv, tmp, /block, min=-80, max=80
  endif

; OUTPUT

  print, 'esi_echtrcobj: Updating trace in ', objfil
  mwrfits, objstr, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
  ;; For multiple images
  for ii=1L,nimg-1 do begin
      objfil = esi[indx[ii]].obj_fil
      print, 'esi_echtrcobj: Updating trace in ', objfil
      tmpstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
      tmpstr[0:9].trace = objstr[0:9].trace
      mwrfits, tmpstr, objfil, /silent, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'esi_echtrcobj: All done!'

  return
end
              
      
      
