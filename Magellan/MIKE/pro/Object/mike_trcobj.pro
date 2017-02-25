;+ 
; NAME:
; mike_trcobj   
;     Version 1.0
;
; PURPOSE:
;    Trace the object through each order of the image.  The program
;    uses a standard star as a crutch (/STD) through regions where the object
;    has very low flux.  The trace is written into the object
;    structure (e.g. Extract/Obj_mike0024.fits)
;
; CALLING SEQUENCE:
;   
;  mike_trcobj, mike, obj_id, [exp_id], /STD, /FAINT, /CHK, NCOLL=
;
; INPUTS:
;   mike     -  ESI structure
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
;   mike_trcobj, mike, 1L, [0L], /CHK, /STD
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Sep-2003 Added mike_trcobj_sngl
;-
;------------------------------------------------------------------------------

pro mike_trcobj_sngl, img, ivar, setup, side, order, trace, $
                      NCOLL=ncoll, GUIDE=guide, FRAD=frad, $
                      FAINT=faint, XYOFF=xyoff, DEBUG=debug, MXORDR=mxordr, $
                      SUBSKY=subsky, FORCGD=forcgd

  if  N_params() LT 6  then begin 
      print,'Syntax - ' + $
        'mike_trcstd, img, ivar, setup, side, order, trace ' + $
        '/USESTD, /CHK, /FAINT, NCOLL=, /GUIDE, /SUBSKY'
      print,   '/DEBUG [v1.0]'
      return
  endif

  ;; NCOLL
  if keyword_set(FAINT) and not keyword_set( NCOLL ) then ncoll = 20L
  if not keyword_set( NCOLL ) then ncoll = 10L
  if not keyword_set( XYOFF ) then xyoff = 0.
  if not keyword_set( XERRMX ) then xerrmx = 0.8
  if not keyword_set( FRAD ) then frad = 3.
  if not keyword_set( MXORDR ) then mxordr = 50L

  ;; Order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)

  ;; Offset
  ordr_str.lhedg = ordr_str.lhedg + xyoff
  ordr_str.rhedg = ordr_str.rhedg + xyoff

  ;; Image size
  sz_img = size(img, /dimensions)

  ;; Setup stuff
  winsz = 30L*sz_img[0]/2048L
  subimg = fltarr(2*winsz+1,ncoll)
  subivar = fltarr(2*winsz+1,ncoll)
  nfit = sz_img[1]/ncoll + 1
  xfit = dblarr(nfit)
  xsig = dblarr(nfit)
  yfit = dblarr(nfit)

  
  ;; Set looping range
  jmid = sz_img[1]/2
  jstrt = ncoll
  jend = sz_img[1]-1
  gdfit = 0L

  ;; Set rnd_trc
  if keyword_set( GUIDE ) then rnd_trc = round(guide) $
  else rnd_trc = (ordr_str[order].lhedg+ordr_str[order].rhedg)/2.

  for qq=0L,1 do begin
      flg_edge = 0
      svoff = 0.
      if  qq EQ 0 then begin
          ;; MIDDLE to TOP
          j1 = jmid
          j2 = jend
          jstep = ncoll
      endif else begin
          ;; MIDDLE to BOTTOM
          j1 = jmid
          j2 = jstrt
          jstep = -ncoll
      endelse

      ;; LOOP
      for j=j1,j2,jstep do begin
      
          ;; Check order edges
          if ordr_str[order].lhedg[j] LT 0 OR $
            ordr_str[order].rhedg[j] GT (sz_img[0]-1) then begin
              flg_edge = 1
              continue
          endif
          ;; Create Subimg
          jp = 0L
          for jsub = (j-ncoll),j-1 do begin
              subimg[*,jp] = 0.
              subivar[*,jp] = 0.
              i1 = (rnd_trc[jsub] - winsz) > 0L
              i2 = 0L > (rnd_trc[jsub] + winsz) < (sz_img[0] - 1)
              if i2 GT i1 then begin
                  s1 = (winsz - rnd_trc[jsub]) > 0L
                  s2 = i2-i1+s1
                  subimg[s1:s2,jp] = img[i1:i2,jsub]
                  subivar[s1:s2,jp] = ivar[i1:i2,jsub]
              endif
              jp = jp+1
          endfor
          ;; Median
          mdn = djs_median(subimg,2)
          if keyword_set(SUBSKY) then mdn = mdn - median(mdn) ; !!
          mdn_var = djs_median(subivar,2)
          mivar = mdn_var > 0.
          xcen = double(winsz) + svoff
          ;; Find Centroid
          for k=0,19 do $
            xcen = trace_fweight(mdn, xcen, 0L, radius=frad, xerr=xerr, $
                                 invvar=mivar)
          ;; Keep the good points
          if xerr LT XERRMX AND xerr GT 0.000001 then begin
              xsig[gdfit] = xerr
              yfit[gdfit] = (j-ncoll/2)
              xfit[gdfit] = rnd_trc[yfit[gdfit]] + xcen - double(winsz)
              if keyword_set(FORCGD) then begin
                  if xfit[gdfit]-guide[yfit[gdfit]] GT 0.5 then begin
                      print, 'mike_trcobj_sngl: Forcing the guide!', j
                      xsig[gdfit] = 0.1
                      yfit[gdfit] = (j-ncoll/2)
                      xfit[gdfit] = guide[yfit[gdfit]] 
                  endif else svoff = xcen-float(winsz)
              endif else svoff = xcen-float(winsz)
              ;; Update gdfit
              gdfit = gdfit + 1L
          endif else begin
              ;; GUIDE
              if keyword_set(GUIDE) then begin
;                  print, 'mike_trcobj_sngl: Guiding..', j
                  xsig[gdfit] = 0.1
                  yfit[gdfit] = (j-ncoll/2)
                  xfit[gdfit] = guide[yfit[gdfit]] 
                  gdfit = gdfit + 1L
              endif 
          endelse
      endfor
  endfor
              

  ;; FIT frac
  frac = (xfit[0:gdfit-1]-ordr_str[order].lhedg[round(yfit[0:gdfit-1])])/$
    (ordr_str[order].rhedg[round(yfit[0:gdfit-1])] $
     -ordr_str[order].lhedg[round(yfit[0:gdfit-1])])
  frac_fit = x_setfitstrct(NITER=4L, NORD=9L, FLGREJ=1L, HSIG=2.2, $
                           LSIG=2.2, FUNC='POLY')

  ;; Change fit order in some cases
  if side EQ 2 AND ordr_str[order].order LT 48L then frac_fit.nord = 13
  if flg_edge EQ 1 then begin
      frac_fit.nord = 5  ;; Edges
  endif

  ;; Set MAX
  frac_fit.nord = frac_fit.nord < MXORDR
  
  ;; FIT
;  print, 'gdfit', gdfit
  fit = x1dfit(yfit[0:gdfit-1], frac, FITSTR=frac_fit)
  print, 'mike_trcobj_sngl: FRAC RMS = ', frac_fit.rms

  ;; DEBUG
  if keyword_set( DEBUG ) then begin
      srt = sort(yfit[0:gdfit-1])
      x_splot, yfit[srt], frac[srt], psym1=1, ytwo=fit[srt], /block
  endif
      
  ;; CALCULATE TRACE
  trace = x_calcfit(findgen(jend+1),FITSTR=frac_fit)* $
    (ordr_str[order].rhedg-ordr_str[order].lhedg) + ordr_str[order].lhedg
  
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_trcobj, mike, setup, obj_id, side, exp, indx=indx, $
                   USESTD=usestd, CHK=chk, NCOLL=ncoll, FAINT=FAINT, $
                   STDFIL=stdfil, SEDG_FIL=sedg_fil, GUIDE=guide, $
                 DEBUG=debug 

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_trcobj, mike, setup, obj_id, side, [exp] ' + $
        '/USESTD, /CHK, /FAINT, NCOLL=, /GUIDE'
      print,   '/DEBUG [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( NCOLL ) then begin
      if keyword_set( FAINT ) then ncoll = 30L else ncoll = 15L
  endif
  if not keyword_set( XERRMX ) then xerrmx = 0.8
  if not keyword_set(USESTD) and not keyword_set(STDFIL) then flg_trc = 1 
  if keyword_set( USESTD ) or keyword_set( GUIDE ) then flg_guide = 1

; INDX
  if NOT keyword_set(indx) then $
  indx = where(mike.flg_anly NE 0 AND mike.side EQ side AND $
               mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
               strtrim(mike.type,2) EQ 'OBJ', nindx) $
  else nindx = n_elements(indx)

  if nindx NE 1 then $
    print, 'mike_trcobj: Combining multiple frames. Beware of thermal expansion '+ $
    'issues'

  if keyword_set( exp_id ) then indx = indx[exp_id]

;  Setup and setup
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  case side of 
      1: begin 
          nm = 'B'
          ostep = 1
      end
      2: begin
          nm = 'R'
          ostep = -1
      end
      else: stop
  endcase

;;;;;;
;  Find standard star
  if keyword_set( USESTD ) then begin
      istd = where(mike.flg_anly NE 0 AND mike.side EQ side AND $
                   mike.setup EQ setup AND $
                   strtrim(mike.type,2) EQ 'STD', nstd)
      case nstd of 
          0 : begin
              print, 'mike_trcobj: No standard star images! Returning..'
              return
          end 
          1 : print, 'mike_trcobj: Tracing standard star image ', $
            mike[istd].img_root
          else : begin
              print, 'mike_trcobj: Warning -- Multiple standard star images'
              istd = istd[0]
              print, 'mike_trcobj: Taking first one -- ', mike[istd].img_root
          end
      endcase
  endif
          
;;;;;;;;;
; Open Stuff

;  Read in order structure
  ordr_fil = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
  if x_chkfil(ordr_fil+'*',/silent) EQ 0 then begin
      print, 'mike_mkaimg: Order structre doesnt exist. ' + $
        'Run mike_fndobj first!'
      return
  endif
  ordr_str = xmrdfits(ordr_fil,1,/silent)
  nordr = n_elements(ordr_str)

  ;; Offset order edges due to thermal expansion
  ordr_str.lhedg = ordr_str.lhedg + mike[indx[0]].arc_xyoff[0]
  ordr_str.rhedg = ordr_str.rhedg + mike[indx[0]].arc_xyoff[0]

  ordr_cen = round((ordr_str.lhedg + ordr_str.rhedg)/2.)

; Orders
  ordrs = [ordr_str[0].order,ordr_str[nordr-1].order]

  ;; OBJ
  objfil = mike[indx[0]].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'mike_trcobj: Obj file does not exist or obj_fil tag not ' + $
        'set! Check it..'
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
  
  ;; STD OBJ
  if flg_trc NE 1 then begin
      if not keyword_set( STDFIL ) then begin
          if keyword_set(USESTD) then stdfil = mike[istd].obj_fil 
      endif
      if x_chkfil(stdfil+'*') EQ 0 then begin
          print, 'mike_trcobj: STD Obj file does not exist or obj_fil tag not ' + $
            'set! Check it..', stdfil
          return
      endif
      stdstr = xmrdfits(stdfil, 1, STRUCTYP='mikeobjstrct', /silent)
  endif

  ;; Open Image, Variance
  print, 'mike_trcobj: Opening [combining] images...'
  nimg = n_elements(indx)

  ;; COMBINE as needed
  img = mike_combimg(mike, indx, IVAR=ivar, IMGINDX=2L, /SKY) 
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
  
  for qq=ordrs[0],ordrs[1],ostep do begin
      ;; Index
      mm = where(ordr_str.order EQ qq, nmm)
      if nmm EQ 0 then stop else mm = mm[0]
      print, 'mike_trcobj: Tracing order ', string(qq,2)
      ;; Zero out subimg
      subimg[*] = 0.
      svoff = 0.
      ;; Trace
      if keyword_set( STDFIL ) then begin  ; STANDARD STAR AS GUIDE
          jdum = objstr[qq].xcen
          ;; Use fraction
          std_frac = (stdstr[mm].trace[jdum]-ordr_str[mm].lhedg[jdum]) $
            /(ordr_str[mm].rhedg[jdum]-ordr_str[mm].lhedg[jdum])
          obj_frac = (objstr[mm].ycen-ordr_str[mm].lhedg[jdum]) $
            /(ordr_str[mm].rhedg[jdum]-ordr_str[mm].lhedg[jdum])
          big_frac = obj_frac/std_frac
          rnd_trc = round(big_frac*stdstr[mm].trace[0:sz_img[1]-1] + $
                          (1.-big_frac)*ordr_str[mm].lhedg[0:sz_img[1]-1])
;          str_off = objstr[qq].ycen - stdstr[qq].trace[jdum] 
;          rnd_trc = round(stdstr[qq].trace[0:sz_img[1]-1] + str_off)
      endif else begin  ; USE CRUDE FIT AS GUIDE
          crude_trc = objstr[mm].trace[0:sz_img[1]-1]
          rnd_trc = round(crude_trc)
      endelse
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
      
      print, 'mike_trcobj: RMS = ', trc_fit.rms
      ;; objstr (Take trace to jend)
      objstr[mm].trace[0:jend] = x_calcfit(findgen(jend+1), FITSTR=trc_fit)
      if keyword_set( DEBUG ) then begin
         srt = sort(yfit[0:gdfit-1])
          x_splot, yfit[srt], xfit[srt], PSYM1=1, $
            YTWO=new_trc[srt], YTHR=xsig[srt], PSYM3=2, /block
      endif

  endfor

  ;; CHK
  if keyword_set( CHK ) then begin
      tmp = img
      for qq=ordrs[0],ordrs[1],ostep do begin
          mm = where(ordr_str.order EQ qq, nmm)
          if nmm EQ 0 then stop else mm = mm[0]
          trc_msk = round(objstr[mm].trace[0:sz_img[1]-1]) $
            + lindgen(sz_img[1])*sz_img[0]
          tmp[trc_msk] = -1000
      endfor
      xatv, tmp, /block, min=-80, max=80
  endif

; OUTPUT

  print, 'mike_trcobj: Updating trace in ', objfil
  mwrfits, objstr, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
  ;; For multiple images
  for ii=1L,nimg-1 do begin
      objfil = mike[indx[ii]].obj_fil
      print, 'mike_trcobj: Updating trace in ', objfil
      tmpstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
      tmpstr[0:nordr-1].trace = objstr[0:nordr-1].trace
      mwrfits, tmpstr, objfil, /silent, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'mike_trcobj: All done!'

  return
end
              
