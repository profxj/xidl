;+ 
; NAME:
; esi_lwdfndobj   
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  esi_lwdfndobj, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdfndobj, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdfndobj, esi, obj_id, exp, REFWV=refwv, CHK=chk, CLOBBER=clobber, $
                   SCICLM=sciclm, STD=std, AUTO=auto

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdfndobj, esi, obj_id, [exp]  [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(REFWV) then refwv = [6500., 6700.]
  if not keyword_set(SCICLM) then sciclm = 1709L

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
                   esi.obj_id EQ obj_id AND esi.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_lwdproc: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = obj_id[0]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = esi[indx[exp[q]]].obj_fil
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOBBER) then begin
          print, 'esi_lwdfndobj: Using Obj structure -- ', objfil
          objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          print, 'esi_lwdfndobj: Creating Obj structure -- ', objfil
          objfil = 'Extract/Obj_'+esi[indx[exp[q]]].img_root
          esi[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { specobjstrct }
          objstr = replicate(tmp, 10)
          objstr.slit_fil = ' '
          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = esi[indx[exp[q]]].exp
          nobj = 0L
      endelse

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+esi[indx[exp[q]]].img_root
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'esi_lwdfndobj: No Final file ', imgfil
          stop
          continue
      endif
      objstr.spec2d_fil = imgfil
      img = xmrdfits(imgfil, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      flg_arc = 0
      objstr.slit_fil = strtrim(esi[indx[exp[q]]].arc_fil,2)
      if q EQ 0 then begin
          img_arc = xmrdfits(strtrim(esi[indx[exp[q]]].arc_fil,2), /silent) 
      endif else begin
          if esi[indx[exp[q]]].arc_fil NE esi[indx[exp[q-1]]].arc_fil then $
            img_arc = xmrdfits(esi[indx[exp[q]]].arc_fil, /silent) $
          else flg_arc=1
      endelse
            
      ;; Mask
      if flg_arc NE 1 then begin
          msk = bytarr(sz_img[0], sz_img[1])
          sva = where(img_arc GT refwv[0] AND img_arc LT refwv[1])
          msk[sva] = 1
      endif
      
      ;; Smash
      smsh = fltarr(sz_img[1])
      print, 'esi_lwdfndobj: Smashing'
      for jj=0L,sz_img[1]-1 do begin
          b = where(msk[*,jj] EQ 1, nb)
          if nb NE 0 then smsh[jj] = median(img[b,jj])
      endfor

      ;; Peaks
      print, 'esi_lwdfndobj: Finding obj'
      gdsmsh = where(smsh GT 0.)
      x_fndobj, smsh[gdsmsh], center, NSIG=3., PEAK=peak, $
        PKWDTH=4L, EDGES=edges 
      gdpeak = gdsmsh[peak]
      npk = n_elements(peak)

      ;; Chk
      if keyword_set( CHK ) then $
        x_splot, smsh, /block, XTWO=gdpeak, $
        YTWO=fltarr(n_elements(center)), psym_y2=1

      ;; Identify SCI
      mn = min(abs(gdpeak-sciclm), sci)
      off = gdpeak[sci]-center[sci]
      xmsk = where(msk[*,gdpeak[sci]] EQ 1)

      ;; Setup Obj
      if flg_objstr EQ 0 then begin
          ;; Science obj
          objstr[0].xcen = median(xmsk)
          objstr[0].ycen = center[sci]+off
          objstr[0].obj_id = 'a'
          if keyword_set( STD ) then objstr[0].aper[*] = 50. $
          else objstr[0].aper[*] = 15.
          nobj = nobj + 1

          if not keyword_set( STD ) then begin
              ;; Other obj
              for jj=0L,npk-1 do begin
                  ;; Skip sci
                  if jj EQ sci then continue
                  ;; Set
                  off = float(gdpeak[jj]-peak[jj])
                  xmsk = where(msk[*,gdpeak[jj]] EQ 1)
                  objstr[nobj].xcen = median(xmsk)
                  objstr[nobj].ycen = center[jj]+off
                  objstr[nobj].obj_id = x_objnumid(nobj)
                  objstr[nobj].aper[*] = 15.
                  nobj = nobj + 1
              endfor
          endif
          ;; Reset
          objstr = objstr[0:nobj-1]
      endif
      
      ;; AUTO
      if not keyword_set( AUTO ) then begin
          dum = x_setapergui( smsh, OBJSTR=objstr )
          nobj = n_elements(objstr)
          ;; Set xcen for new obj
          newobj = where(objstr.xcen EQ 0., nnew)
          if nnew NE 0 then begin
              for i=0L,nnew-1 do begin
                  peak = round(objstr[newobj[i]].ycen)
                  xmsk = where(msk[*,peak] EQ 1)
                  objstr[newobj[i]].xcen = median(xmsk)
                  objstr[newobj[i]].aper[*] = 15.
              endfor
          endif
      endif

      ;; Set Trace 'by hand'
      trc_fil = getenv('XIDL_DIR')+'/ESI/CALIBS/LWD_TRC.fits'
      gdtrc = xmrdfits(trc_fil, /silent)

      for jj=0L,nobj-1 do begin
          ;; Offset
          cen_gd = gdtrc[objstr[jj].xcen]
          off = objstr[jj].ycen - cen_gd
          ;; Trace
          trc = gdtrc + off
          ;; Save
          objstr[jj].trace[0:sz_img[0]-1] = trc
      endfor
          

      ;; Write Obj structure
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'esi_lwdfndobj: All done! '
  return
end
