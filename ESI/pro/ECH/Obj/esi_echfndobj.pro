;+ 
; NAME:
; esi_echfndobj   
;     Version 1.1
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  esi_echfndobj, esi, obj_id, [exp], REFWV=, SCICLM=,
;  REFORDR=, /INTER, /STD, /CHK, MAXOFF=, /NOCLOB
;
; INPUTS:
;   esi     -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Find object for standard star
;  /CHK     - Show overall trace
;  FITFIL=  - Map of pinholes (default: Maps/hole_fit.idl )
;  REFWV=   - Reference wavelength region (default: [5600., 5610.] )
;  REFORDR= - Reference order  (default: 4L)
;  /INTER   - Interactively identify objects and set apertures
;             (recommended)
;  MAXOFF=  - Minimum offset between max obj and center of slit
;             (default: 20.)
;  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echfndobj, esi, 0L, [0L, 1L], /CHK, /INTER, 
;      REFWV=[6500., 6520.], REFORDR=5L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Jul-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------

function esi_echfndobj_off, slit_edg, img, gdtrc, qq, MEDJ=medj, SCICLM=sciclm

  ;; cen
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  ;; SIZE
  sz_img = size(img, /dimensions)
  ;; Mask
  msk = bytarr(sz_img[0], sz_img[1])
  for j=2000L,2100L do begin
      imn = slit_edg[j,qq,0]
      imx = slit_edg[j,qq,1]
      msk[imn:imx,j] = 1
  endfor
  ;; Smash
  smsh = fltarr(sz_img[0])
  svb = [0]
  svi = [0]
  for ii=0L,sz_img[0]-1 do begin
      b = where(msk[ii,*] EQ 1, nb)
      if nb NE 0 then begin
          svb = [svb, b]        ; Grab j value
          svi = [svi, ii]
          smsh[ii] = median(img[ii,b])
      endif
  endfor
  medj = round(median(float(svb)))
  ;; Peaks
  gdsmsh = where(smsh[round(slit_edg[medj,qq,0]): $
                      round(slit_edg[medj,qq,1])] NE 0.)
  gdsmsh = gdsmsh + round(slit_edg[medj,qq,0])
  x_fndobj, smsh[gdsmsh], center, PEAK=peak, NSIG=3., PKWDTH=4L, $
    EDGES=edges, NORDB=2L, /silent
  if center[0] EQ -1. AND n_elements(center) EQ 1 then begin
      print, 'esi_echfndobj_off: No obj found!!', qq
      return, -1.
  end
      
  ;; Parse peaks
  if keyword_set( INTER ) then $
    x_prspeaks, smsh[gdsmsh], center, pk_msk, /block $
  else pk_msk = bytarr(n_elements(peak)) + 1B
  ;; Reset peaks
  peak = peak[where(pk_msk EQ 1)]
  center = center[where(pk_msk EQ 1)] + gdsmsh[0]
  gdpeak = gdsmsh[peak]
  npk = n_elements(peak)
  ;; Set Science obj as max + within 20 pixels of center?
  cen_ordr = slit_cen[medj,qq]
  if keyword_set( SCICLM ) then cen_ordr = cen_ordr + SCICLM

  nearobj = where(abs(gdpeak-cen_ordr) LT 10., nnear)
  if nnear EQ 0 then begin
      print, 'esi_echfndobj_off: No obj near center of order!!', qq
      print, 'esi_echfndobj_off: Using last off_sci'
      return, -1
  endif
  ;; Find Max near the center
  mx = max(smsh[gdpeak[nearobj]], sci)
  sci = nearobj[sci]
  off = gdpeak[sci]-gdtrc[medj,qq]
  return, off
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfndobj, esi, obj_id, exp, REFWV=refwv, FITFIL=fitfil, $
                   SCICLM=sciclm, REFORDR=refordr, INTER=inter, STD=std, $
                   CHK=chk, MAXOFF=maxoff, NOCLOB=noclob, SEDG_FIL=sedg_fil, $
                   NSIG=nsig

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echfndobj, esi, obj_id, [exp], /INTER, /STD [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(OBJ_ID) then obj_id = 0L
  if not keyword_set(MAXOFF) then maxoff=20.
  if not keyword_set(REFORDR) then refordr = 4L
  if not keyword_set(REFWV) then refwv = [5600., 5610.]
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( NSIG ) then nsig = 3.

;  Find all relevant obj images
  if not keyword_set( STD) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', $
                   nindx)
      if nindx EQ 0 then begin
          print, 'esi_echfndobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = obj_id[0]
      nindx = 1L
  endelse
  

; Open Slit file
  c_s = esi_slitnm(esi[indx[0]].slit)
  if not keyword_set( SEDG_FIL ) then $
    sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
  if x_chkfil(sedg_fil+'*') EQ 0 then begin
      print, 'esi_echfndobj: Slit edge file ', sedg_fil, ' does not exist!'
      return
  endif
  print, 'esi_echfndobj: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)

; Open Trace file
  trc_fil = getenv('XIDL_DIR')+'/ESI/CALIBS/ECH_TRC.fits'
  if x_chkfil(trc_fil+'*') EQ 0 then begin
      print, 'esi_echfndobj: Slit edge file ', trc_fil, ' does not exist!'
      return
  endif
  gdtrc = xmrdfits(trc_fil, /silent)

; Hole trace
  restore, fitfil  ; fin_fit

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = esi[indx[exp[q]]].obj_fil
      ;; REFORDR
      esi[indx[exp[q]]].refordr = refordr
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'esi_echfndobj: Using Obj structure -- ', objfil
          if x_chkfil(objfil+'*') EQ 0 then begin
              print, 'esi_echfndobj: File ', objfil, $
                'doesnt exist! Turn off NOCLOB'
              return
          endif
          objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          objfil = 'Extract/Obj_'+esi[indx[exp[q]]].img_root
          esi[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { dblsobjstrct }
          objstr = replicate(tmp, 50)
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
          print, 'esi_echfndobj: No Final file ', imgfil
          stop
          continue
      endif
      objstr.spec2d_fil = imgfil
      print, 'esi_echfndobj: Opening files... ', imgfil
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
      smsh = fltarr(sz_img[0])
      print, 'esi_echfndobj: Smashing...'
      svb = [0]
      svi = [0]
      for ii=0L,sz_img[0]-1 do begin
          b = where(msk[ii,*] EQ 1, nb)
          if nb NE 0 then begin
              svb = [svb, b]  ; Grab j value
              svi = [svi, ii]
              smsh[ii] = median(img[ii,b])
          endif
      endfor
      medj = round(median(float(svb)))

      ;; Peaks
      gdsmsh = where(smsh[round(slit_edg[medj,refordr,0]): $
                          round(slit_edg[medj,refordr,1])] NE 0.)
      gdsmsh = gdsmsh + round(slit_edg[medj,refordr,0])

      ;; Interactive
      if keyword_set( INTER ) then begin
;          tmp = { specobjstrct }
;          tmp.obj_id
          if keyword_set( tmpobj ) then delvarx, tmpobj
          dum = x_setapergui( smsh[gdsmsh], OBJSTR=tmpobj )
          npk = n_elements(tmpobj)
          sciobj = where(tmpobj.obj_id EQ 'a', nsci)
          sci = sciobj[0]
          if nsci EQ 0 then begin
              print, 'esi_echfndobj: You failed to identify the sci obj!'+$
                '  Returning!'
              return
          endif
          off_sci = gdsmsh[round(tmpobj[sciobj].ycen)]-gdtrc[medj,refordr]
          center = gdsmsh[round(tmpobj.ycen)]
          gdpeak = round(center)
      endif else begin ;; AUTO
          x_fndobj, smsh[gdsmsh], center, PEAK=peak, NSIG=nsig, PKWDTH=4L, $
            EDGES=edges, NORDB=2L
          ;; Check for center
          if center[0] EQ -1. AND n_elements(center) EQ 1 then begin
              print, 'esi_echfndobj: No obj found!!'
              stop
              return
          end
          pk_msk = bytarr(n_elements(peak)) + 1B
          ;; Reset peaks
          peak = peak[where(pk_msk EQ 1)]
          center = center[where(pk_msk EQ 1)] + gdsmsh[0]
          gdpeak = gdsmsh[peak]
          npk = n_elements(peak)
          ;; OFFSET center of slit
          cen_ordr = slit_cen[medj,refordr]
          if keyword_set( SCICLM ) then cen_ordr = cen_ordr + SCICLM
          ;; Set Science obj as max within 20 pixels of center
          nearobj = where(abs(gdpeak-cen_ordr) LT maxoff, nnear)
          if nnear EQ 0 then begin
              print, 'esi_echfndobj: No obj near center of order!!  Punting'
              stop
              return
          endif
          ;; Find Max near the center
          mx = max(smsh[gdpeak[nearobj]], sci)
          sci = nearobj[sci]
          off_sci = gdpeak[sci]-gdtrc[medj,refordr]
      endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Setup objstr (default aper to 11)
      if flg_objstr EQ 0 then begin
          print, 'esi_echfndobj: Tracing..'

          ;; Science first
          if not keyword_set( INTER ) then begin
              ;; AUTO
              if keyword_set( STD ) then objstr[0:9].aper[*] = 25. $
              else objstr[0:9].aper[*] = 10.
          endif else begin
              ;; INTER
              objstr[0:9].aper = tmpobj[sci].aper 
          endelse
          ;; xcen
          objstr[0:9].obj_id = 'a'
          objstr[0:9].slit_id = lindgen(10)

          for qq=REFORDR,9L do begin ;; GOING UP
              ;; TRACE
              if qq EQ REFORDR then begin
                  off = off_sci 
                  off_sv = off_sci
                  medx = medj
              endif else begin
                  off = esi_echfndobj_off(slit_edg, img, gdtrc, qq, MEDJ=medx,$
                                         SCICLM=sciclm)
                  if off EQ -1. or abs(off-off_sv) GT 3. then off = off_sv
                  off_sv = off
              endelse
              ;; xcen, ycen
              objstr[qq].xcen = medx
              objstr[qq].ycen = gdtrc[medx,qq] + off
              ;; TRACE
              objstr[qq].trace[0:sz_img[1]-1] = gdtrc[0:sz_img[1]-1,qq] + off
              ;; Fix 9 centroids
              if qq EQ 9 then begin
                  objstr[qq].xcen = 1000L
                  objstr[qq].ycen = objstr[9].trace[1000L] 
              endif 
          endfor
          ;; Going down
          off_sv = off_sci
          for qq=REFORDR-1,0L,-1 do begin
              if not keyword_set( INTER ) then begin
                  ;; AUTO
                  if keyword_set( STD ) then objstr[qq].aper[*] = 25. $
                  else objstr[qq].aper[*] = 10.
              endif else begin
                  ;; INTER
                  objstr[qq].aper[*] = tmpobj[sci].aper 
              endelse
              ;; Off
              off = esi_echfndobj_off(slit_edg, img, gdtrc, qq, MEDJ=medx, $
                                     SCICLM=sciclm)
              if off EQ -1. or abs(off-off_sv) GT 3. then off = off_sv
              off_sv = off

              ;; ycen
              objstr[qq].xcen = medx
              objstr[qq].ycen = gdtrc[medx,qq] + off
              ;; TRACE
              objstr[qq].trace[0:sz_img[1]-1] = gdtrc[0:sz_img[1]-1,qq] + off
              ;; Fix 9 centroids
              if qq EQ 9 then begin
                  objstr[qq].xcen = 1000L
                  objstr[qq].ycen = objstr[9].trace[1000L] 
              endif 
          endfor
          nobj = nobj + 1
          ;; CHK
          if keyword_set( CHK ) then begin
              tmp = img
              for qq=0L,9 do begin
                  rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
                  trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
                  tmp[trc_msk] = -10000
              endfor
              print, objstr[refordr].xcen, objstr[refordr].ycen
              xatv, tmp, /block, min=-70, max=770
          endif

          ;; Set other obj
          if not keyword_set( STD ) then begin
              for jj=0L,npk-1 do begin
                  ;; Skip sci
                  if jj EQ sci then continue
                  off_obj = gdpeak[jj]-gdtrc[medj,refordr]
                  for qq=0L,9 do begin
                      objstr[nobj*10+qq].xcen = medj
                      objstr[nobj*10+qq].ycen = center[jj]+off_obj
                      objstr[nobj*10+qq].aper[0] = 10.
                      objstr[nobj*10+qq].aper[1] = 10.
                      objstr[nobj*10+qq].obj_id = x_objnumid(nobj)
                      objstr[nobj*10+qq].slit_id = qq
                      ;; TRACE (not ideal!)
                      objstr[nobj*10+qq].trace[0:sz_img[1]-1] = $
                        gdtrc[0:sz_img[1]-1,qq] + off_obj
                      ;; Fix 9 centroids
                      if qq EQ 9 then begin
                          objstr[nobj*10+qq].xcen = 1000L
                          objstr[nobj*10+qq].ycen = objstr[9].trace[1000L] 
                      endif 
                  endfor
                  ;; Update nobj
                  nobj = nobj + 1
                  ;; CHK
                  if keyword_set( CHK ) then begin
                      tmp = img
                      for qq=0L,9 do begin
                          rnd_trc = round(objstr[nobj*10+qq].trace[0:sz_img[1]-1])
                          trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
                          tmp[trc_msk] = -10000
                      endfor
                      print, objstr[refordr].xcen, objstr[refordr].ycen
                      xatv, tmp, /block, min=-70, max=70
                  endif
              endfor
          endif
          ;; Take only the relevant structures
          objstr = objstr[0:nobj*10L-1]
      endif
      

      ;; Write Obj structure
      print, 'esi_echfndobj: Creating ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'esi_echfndobj: All done! '
  return
end
