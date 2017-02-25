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
;  esi_echtrcobj, esi, obj_id, [exp_id],PEAKTHRESH=, /STD, /FAINT, /CHK, NCOLL=
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
;   PEAKTHRESH = Fraction of brightest object on the slit to trace. 
;   EXFIL=  - Kludge for standard star tracing

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

pro esi_echtrcobj, esi, obj_id, exp_id, ORDRS=ordrs $
                   , USESTD = usestd, CHK = chk, NCOLL = ncoll, FAINT = FAINT $
                   , STDFIL = stdfil, SEDG_FIL = sedg_fil $
                   , CBIN = cbin, RBIN = rbin $
                   , RADIUS = RADIUS, OLDTRC = OLDTRC, FLG_APPEND = FLG_APPEND $
                   , PEAKTHRESH = PEAKTHRESH, OBJSEP = OBJSEP, NFIND = NFIND $
                   , FWHM = FWHM, EXFIL=exfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtrcobj, esi, indx, /USESTD, /CHK, /FAINT, NCOLL=, [v1.1]'
      return
  endif 

  message, 'esi_echtrcobj: This routine is obsolete.  Use esi_echfndobj ' + $
           'to find and trace.'
;  Optional Keywords
  alphabet = ['a', 'b', 'c', 'd', 'e', 'f']

  IF NOT KEYWORD_SET(RADIUS) THEN RADIUS = 3.0D
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( NCOLL ) then begin
      if keyword_set( FAINT ) then ncoll = 40L else ncoll = 25L
  endif
  if not keyword_set( XERRMX ) then xerrmx = 0.8
  if not keyword_set(ORDRS) then ordrs=[0L,9L]
  if not keyword_set(USESTD) and not keyword_set(STDFIL) then stop

  if not keyword_set(OBJSEP) then objsep = 0.05

; INDX
  if not keyword_set(EXFIL) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
      if indx[0] EQ -1 then begin
          print, 'esi_echtrcobj: Found no obj!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
      ncoll = 10L
  endelse
  if size( exp_id, /type ) NE 0 then indx = indx[exp_id]
  if keyword_set(EXFIL) then STDFIL = exfil

;;;;;;
;  Find standard star
  if keyword_set( USESTD ) then begin
      istd = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
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

  if not keyword_set( SEDG_FIL ) then begin
      slit_edg = esi_getfil('sedg_fil', SLIT=esi[indx[0]].slit, $
                            cbin=cbin, rbin=rbin)
  endif else begin
      slit_edg = xmrdfits(strtrim(SEDG_FIL,2), 0, /silent)
  endelse
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)

  ;; OBJ
  objfil = esi[indx[0]].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'esi_echtrcobj: Obj file does not exist or obj_fil tag not ' + $
        'set! Check it..'
      return
  endif
  objstr = xmrdfits(objfil, 1, /silent)
  objfinal = objstr
  nobj = n_elements(objstr)/10L

  ;; STD OBJ
  if not keyword_set(STDFIL) then begin
      if keyword_set(USESTD) then stdfil = esi[istd].obj_fil else stop
  endif
  if x_chkfil(stdfil+'*') EQ 0 then begin
      print, 'esi_echtrcobj: STD Obj file does not exist or obj_fil tag not ' + $
             'set! Check it..', stdfil
      return
  endif
  stdstr = xmrdfits(stdfil, 1, /silent)

  ;; Open Image, Variance
  print, 'esi_echtrcobj: Opening [combining] images...'
  nimg = n_elements(indx)
  ;; COMBINE as needed
  img = esi_echcombimg(esi, indx, VAR=var, IMGINDX=2L) 
  sz_img = size(img, /dimensions)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Trace

  if keyword_set(OLDTRC) then begin
      subimg = fltarr(21,ncoll)
      subvar = fltarr(21,ncoll)
      nfit = sz_img[1]/ncoll + 1
      xfit = dblarr(nfit)
      xsig = dblarr(nfit)
      yfit = dblarr(nfit)
      
      FOR iobj = 0L, nobj-1L DO BEGIN
          indo1 = iobj*10L
          indo2 = (iobj+1L)*10L-1L
          objnow = objstr[indo1:indo2]
          for qq=ordrs[0],ordrs[1] do begin
              print, 'esi_echtrcobj: Tracing order ', string(15L-qq, FORMAT='(i3)')
              ;; Zero out subimg
              subimg[*] = 0.
              svoff = 0.
              ;; Trace
              if keyword_set( STDFIL ) then begin ; STANDARD STAR AS GUIDE
                  jdum = objnow[qq].xcen
;          std_frac = (stdstr[qq].trace[jdum]-slit_edg[jdum,qq,0]) $
;            /(slit_edg[jdum,qq,1]-slit_edg[jdum,qq,0])
;          obj_frac = (objstr[qq].ycen-slit_edg[jdum,qq,0]) $
;            /(slit_edg[jdum,qq,1]-slit_edg[jdum,qq,0])
;          big_frac = obj_frac/std_frac
                  str_off = objnow[qq].ycen - stdstr[qq].trace[jdum] 
;          rnd_trc = round(big_frac*stdstr[qq].trace[0:sz_img[1]-1] + $
;                          (1.-big_frac)*slit_edg[0:sz_img[1]-1,qq,0])
                  rnd_trc = round(stdstr[qq].trace[0:sz_img[1]-1] + str_off)
              endif else begin  ; PIN HOLE AS GUIDE
                  rnd_trc = round(objnow[qq].trace[0:sz_img[1]-1])
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
                      jstrt =1500L/rbin
                      jend = 3800L/rbin
                  end
                  9: begin
                      jstrt = ncoll + 3L
                      jend = 2170L/rbin
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
                      IF (i2 GT i1) AND i2 GT 0 THEN BEGIN
                          subimg[s1:s2, jp] = img[i1:i2, jsub]
                          subvar[s1:s2, jp] = var[i1:i2, jsub]
                      ENDIF
                      jp = jp+1
                  endfor
                  ;; Median
                  mdn = djs_median(subimg,2)
                  mdn_var = djs_median(subvar,2)
                  ivar = 1./(mdn_var > 0.)
                  xcen = 10.d + svoff
                  ;; Find Centroid
                  for k=0,19 do $
                        xcen = trace_fweight(mdn, xcen, 0L, $
                                             radius = radius, xerr = xerr, $
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
                  mdn = djs_median(subimg, 2)
                  mdn_var = djs_median(subvar,2)
                  ivar = 1./(mdn_var > 0.)
                  xcen = 10.d + svoff
                  ;; Find Centroid
                  for k=0,19 do $
                        xcen = trace_fweight(mdn, xcen, 0L, $
                                             radius = radius, xerr = xerr, $
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
        
              trc_fit = x_setfitstrct(NITER=2L, NORD=9L, $
                                      FLGREJ=1L, HSIG=5., LSIG=5., $
                                      FUNC='POLY')
              new_trc = x_fitrej(yfit[0:gdfit-1], xfit[0:gdfit-1], $
                                 SIG=xsig[0:gdfit-1], $
                                 FITSTR=trc_fit)
      
              print, 'esi_echtrcobj: RMS(pix) = ', trc_fit.rms
              ;; objstr (Take trace to jend)
              objnow[qq].trace[0:jend] = x_calcfit(findgen(jend+1), FITSTR=trc_fit)
;      if keyword_set( CHK ) then begin
;          srt = sort(yfit[0:gdfit-1])
;          x_splot, yfit[srt], xfit[srt], PSYM1=1, $
;            YTWO=new_trc[srt], YTHR=xsig[srt], PSYM_Y3=2, /block
;      endif
              
              ;; CHK
              if keyword_set( CHK ) then begin
                  tmp = img
                  for qq=ordrs[0],ordrs[1] do begin
                      trc_msk = round(objnow[qq].trace[0:sz_img[1]-1]) $
                                + lindgen(sz_img[1])*sz_img[0]
                      tmp[trc_msk] = -1000
                  endfor
                  xatv, tmp, /block, min=-80, max=80
              endif
          ENDFOR
          objfinal[indo1:indo2] = objnow
      ENDFOR
  endif else begin ;; Joe trace

      gdtrc = stdstr.trace

      ;; Trace
      objstruct = esi_echjoefind(img, var, gdtrc, slit_edg, cbin $
                                 , GFRAC = gfrac, FWHM = FWHM $
                                 , PEAKTHRESH = PEAKTHRESH)
      nobjs = n_elements(objstruct)/10
      
      ;; Labels
      srtg = sort(abs(gfrac-0.5))
      lbl = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h','i']

      ;; Existing traces
      nold = n_elements(objfinal)/10

      mskg = intarr(nobjs)
      for kk=0L,nold-1 do begin
          ;; Find closest (use 4th order since it is in the middle)
          mtg = where(abs(objfinal[kk*10 + 4].trace[4999L] - gfrac) LT OBJSEP $
                      , nmtg)
          if nmtg NE 1 then message $
            , 'Could not find one of the slits within the OBJSEP tolerance'
          mskg[mtg] = 1
          
          i1 = kk*10
          ;; Fill it up
          for jj=0L,9 do begin
              mt = where(abs(objstruct.xfracpos-gfrac[mtg[0]]) LT OBJSEP AND $
                         objstruct.slitid EQ (jj+1), nmt)
              if nmt NE 1 then stop
              ;; Trace
              objfinal[i1+jj].trace[0:sz_img[1]-1] = objstruct[mt].xpos
              objfinal[i1+jj].obj_id = lbl[srtg[mtg[0]]]
          endfor
          ;; CHK
          if keyword_set( CHK ) then begin
              print, 'Object: ', lbl[srtg[mtg[0]]]
              tmp = img
              for qq=ordrs[0],ordrs[1] do begin
                  trc_msk = round(objfinal[i1+qq].trace[0:sz_img[1]-1]) $
                            + lindgen(sz_img[1])*sz_img[0]
                  tmp[trc_msk] = -1000
              endfor
              xatv, tmp, /block, min=-80, max=80
          endif
              
      endfor

      ;; New obj?
      FLG_APPEND = 0L
      if nobjs GT nold then begin
          stop 
          FLG_APPEND = 1L
          ;; JXP -- The following could be buggy but could be ok  (Sep 2006)
          newg = where(mskg EQ 0)
          lbl = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h','i']
          medj = sz_img[1]/2

          for kk=0L,(nobjs-nold)-1 do begin
              tmp = { esiobjstrct }
              objstr = replicate(tmp, 10L)
              objstr.img_fil = ' '
              objstr.UT = ' '
              objstr.field = ' '
              objstr.exp = objfinal[0].exp

              objstr[*].xcen = medj

              objstr[*].obj_id = lbl[kk+nold]

              ;; Trace
              for jj=0L,9 do begin
                  mt = where(abs(objstruct.xfracpos-gfrac[newg[kk]]) LT OBJSEP AND $
                             objstruct.slitid EQ (jj+1), nmt)
                  if nmt NE 1 then stop
                  ;; xcen, ycen
                  objstr[i1+jj].xcen = medj
                  objstr[i1+jj].ycen = objstruct[mt].xpos[medj]
                  ;; Trace
                  objstr[i1+jj].trace[0:sz_img[1]-1] = objstruct[mt].xpos
              endfor

              ;; Fill it upt
              objfinal = [objfinal, objstr]
          endfor
      endif
      
  endelse
  
  if keyword_set(CHK) then begin
      nobj = n_elements(objstr)/10L
      tmp = img
      FOR iobj = 0L, nobj-1L DO BEGIN
          FOR qq = ordrs[1], ordrs[0], -1 do begin
              this = where(objfinal.obj_id EQ alphabet[iobj] $
                           AND objfinal.order EQ qq $
                           , COMPLEMENT = b, NCOMPLEMENT = nb)
              rnd_trc = round(objfinal[this].trace[0:sz_img[1]-1])
              trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
              tmp[trc_msk] = -10000
          ENDFOR
      ENDFOR
      xatv, tmp, /block, min = -70, max = 70
  endif
      
  

; OUTPUT
  print, 'esi_echtrcobj: Updating trace in ', objfil
  mwrfits, objfinal, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
  ;; For multiple images
  for ii=1L,nimg-1 do begin
      objfil = esi[indx[ii]].obj_fil
      print, 'esi_echtrcobj: Updating trace in ', objfil
      tmpstr = xmrdfits(objfil, 1, /silent)
      if keyword_set(OLDTRC) or flg_append EQ 0 then begin
          IF n_elements(tmpstr) NE n_elements(objfinal) THEN $
            message, 'The number of objects found changed. Need to fix bookeeping'
          tmpstr.trace = objfinal.trace
      endif else begin
          stop
          ;; Number of objects changed so need to work out bookeeping
      endelse
;      tmpstr[0:9].trace = objstr[0:9].trace
      mwrfits, tmpstr, objfil, /silent, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'esi_echtrcobj: All done!'

  return
end
              
      
      
