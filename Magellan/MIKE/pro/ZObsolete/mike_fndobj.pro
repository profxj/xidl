;+ 
; NAME:
; mike_fndobj   (Retired)
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  mike_fndobj, mike, obj_id, [exp], REFWV=, SCICLM=,
;  REFORDR=, /INTER, /STD, /CHK, MAXOFF=, /NOCLOB
;
; INPUTS:
;   mike     -  ESI structure
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
;   mike_fndobj, mike, 0L, [0L, 1L], /CHK, /INTER, 
;      REFWV=[6500., 6520.], REFORDR=5L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function mike_fndobj_off, ordr_str, img, gdtrc, qq, MEDJ=medj, SCICLM=sciclm,$
                          RBIN=rbin

  ;; Set region for obj
  ymin = 2200/rbin
  ymax = 2300/rbin


  ;; SIZE
  sz_img = size(img, /dimensions)

  ;; Mask
  msk = bytarr(sz_img[0], sz_img[1])
  for j=ymin,ymax do begin
      imn = round(ordr_str[qq].lhedg[j]) > 0
      imx = 0 > round(ordr_str[qq].rhedg[j]) < (sz_img[0]-1)
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

  ;; Check edges
  if round(ordr_str[qq].lhedg[medj]) LT 0 OR $
    round(ordr_str[qq].rhedg[medj]) GE sz_img[0] then begin
      print, 'mike_fndobj_off: Off edge of CCD!!', ordr_str[qq].order
      return, -1
  endif
    

  ;; Peaks
  gdsmsh = where(smsh[round(ordr_str[qq].lhedg[medj]): $
                      round(ordr_str[qq].rhedg[medj])] NE 0.)
  gdsmsh = gdsmsh + round(ordr_str[qq].lhedg[medj])
  x_fndobj, smsh[gdsmsh], center, PEAK=peak, NSIG=3., PKWDTH=4L, $
    EDGES=edges, NORDB=2L, /silent
  if center[0] EQ -1. AND n_elements(center) EQ 1 then begin
      print, 'mike_fndobj_off: No obj found!!', qq
      print, "giving up on X's crappy fndobj, let's use find_npeaks,  Hehe, just joking X..."
      center = find_npeaks(smsh[gdsmsh], nfind=1, width=3L, ypeak=peak)
      print, 'found a peak with height', peak[0], ' at pixel ', center[0]

      return, -1.
  end
      
  ;; Parse peaks
  if keyword_set( INTER ) then $
    x_prspeaks, smsh[gdsmsh], center, pk_msk, /block $
  else pk_msk = bytarr(n_elements(peak)) + 1B
  ;; Reset peaks
  peak = peak[where(pk_msk EQ 1)]
  gdpeak = gdsmsh[center[where(pk_msk EQ 1)]]
  center = center[where(pk_msk EQ 1)] + gdsmsh[0]
  npk = n_elements(peak)

  ;; Set Science obj as max + within 20 pixels of center?
  cen_ordr = round((ordr_str[qq].lhedg[medj] + ordr_str[qq].rhedg[medj])/2.)
  if keyword_set( SCICLM ) then cen_ordr = cen_ordr + SCICLM

  nearobj = where(abs(gdpeak-cen_ordr) LT 10., nnear)
  if nnear EQ 0 then begin
      print, 'mike_fndobj_off: No obj near center of order!!', qq
      print, 'mike_fndobj_off: Using last off_sci'
      return, -1
  endif
  ;; Find Max near the center
  mx = max(smsh[gdpeak[nearobj]], sci)
  sci = nearobj[sci]
  off = gdpeak[sci]-gdtrc[medj,qq]
  ;; Fraction
  frac = (gdpeak[sci] - ordr_str[qq].lhedg[medj]) / $
    (ordr_str[qq].rhedg[medj]-ordr_str[qq].lhedg[medj])
  return, frac
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_fndobj, mike, setup, obj_id, side, exp, REFWV=refwv, FITFIL=fitfil, $
                   SCICLM=sciclm, REFORDR=refordr, INTER=inter, STD=std, $
                   CHK=chk, MAXOFF=maxoff, NOCLOB=noclob, SEDG_FIL=sedg_fil, $
                   NSIG=nsig

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_fndobj, mike, setup, obj_id, side, [exp], /INTER, /STD, /CHK [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(MAXOFF) then maxoff=5.
  if not keyword_set( NSIG ) then nsig = 3.

;; REFWV
  if not keyword_set(REFWV) then begin
      if side EQ 1 then refwv = alog10([4347., 4349.d]) $
      else refwv = alog10([6630., 6633.d]) 
  endif
  if not keyword_set(REFORDR) then begin
      if side EQ 1 then refordr = 82L else refordr = 52L
  endif

;  Setup and setup
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  case side of 
      1: nm = 'B'
      2: nm = 'R'
      else: stop
  endcase

;  Find all relevant obj images
  if not keyword_set( STD) then begin
      indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
                   mike.side EQ side  AND $
                   mike.obj_id EQ obj_id AND strtrim(mike.type,2) EQ 'OBJ', $
                   nindx)
      if nindx EQ 0 then begin
          print, 'mike_fndobj: No images to find obj for!', obj_id, ' and side ',$
            side
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = obj_id[0]
      nindx = 1L
  endelse
  
; APERTURES
  if not keyword_set( STDAPER ) then stdaper = [14.,14.]/mike[indx[0]].colbin
  if not keyword_set( OBJAPER ) then objaper = [10.,10.]/mike[indx[0]].colbin

;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  nordr = n_elements(ordr_str)

; Set mm
  mm = where(ordr_str.order EQ refordr, nmm)
  mm = mm[0]
  if nmm EQ 0 then stop

; Cen ref
  cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2.
  stop  ;; JXP (9/22) -- Need to deal with shifting image!
 
;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = mike[indx[exp[q]]].obj_fil
      ;; xyoff
      xyoff = mike[indx[exp[q]]].arc_xyoff[0]
      ;; REFORDR
;      mike[indx[exp[q]]].refordr = refordr
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'mike_fndobj: Using Obj structure -- ', objfil
          if x_chkfil(objfil+'*') EQ 0 then begin
              print, 'mike_fndobj: File ', objfil, $
                'doesnt exist! Turn off NOCLOB'
              return
          endif
          objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          objfil = 'Extract/Obj_'+mike[indx[exp[q]]].img_root
          mike[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { mikeobjstrct }
          objstr = replicate(tmp, nordr)
          objstr.slit_fil = ' '
          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = mike[indx[exp[q]]].exp
          nobj = 0L
      endelse

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+mike[indx[exp[q]]].img_root
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'mike_fndobj: No Final file ', imgfil
          stop
          continue
      endif
      objstr.spec2d_fil = imgfil
      print, 'mike_fndobj: Opening files... ', imgfil
      img = xmrdfits(imgfil, /silent)
;      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      flg_arc = 0
      objstr.slit_fil = strtrim(mike[indx[exp[q]]].arc_img,2)
      if q EQ 0 then begin
          img_arc = xmrdfits(strtrim(mike[indx[exp[q]]].arc_img,2), /silent) 
      endif else begin
          if mike[indx[exp[q]]].arc_img NE mike[indx[exp[q-1]]].arc_img then $
            img_arc = xmrdfits(mike[indx[exp[q]]].arc_img, /silent) $
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
      print, 'mike_fndobj: Smashing...'
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
      gdsmsh = where(smsh[round(ordr_str[mm].lhedg[medj]+xyoff): $
                          round(ordr_str[mm].rhedg[medj]+xyoff)] NE 0.)
      gdsmsh = gdsmsh + round(ordr_str[mm].lhedg[medj])

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
              print, 'mike_fndobj: You failed to identify the sci obj!'+$
                '  Returning!'
              return
          endif
          off_sci = gdsmsh[round(tmpobj[sciobj].ycen)]-cen_ref[medj,mm]
          center = gdsmsh[round(tmpobj.ycen)]
          gdpeak = round(center)
      endif else begin ;; AUTO
          x_fndobj, smsh[gdsmsh], center, PEAK=peak, NSIG=nsig, PKWDTH=3L, $
            EDGES=edges, NORDB=2L
          ;; Check for center
          if center[0] EQ -1. AND n_elements(center) EQ 1 then begin
              print, 'mike_fndobj: No obj found!!'
              print, "giving up on X's crappy fndobj, let's use find_npeaks"
              center = find_npeaks(smsh[gdsmsh], nfind=1, width=3L, ypeak=peak)
              print, 'found a peak with height', peak[0], ' at pixel ', center[0]
          end
          pk_msk = bytarr(n_elements(peak)) + 1B
          ;; Reset peaks
          peak = peak[where(pk_msk EQ 1)]
          gdpeak = gdsmsh[center[where(pk_msk EQ 1)]]
          center = center[where(pk_msk EQ 1)] + gdsmsh[0]
          npk = n_elements(peak)
          ;; OFFSET center of slit
          cen_ordr = cen_ref[medj,mm]
          if keyword_set( SCICLM ) then cen_ordr = cen_ordr + SCICLM
          ;; Set Science obj as max within 20 pixels of center
          nearobj = where(abs(gdpeak-cen_ordr) LT maxoff, nnear)
          if nnear EQ 0 then begin
              print, 'mike_fndobj: No obj near center of order!!  Punting'
              print, 'mike_fndobj: Try running with interactive or choosing' + $
                ' another region of the chip for centroiding.' 

              stop
              return
          endif
          ;; Find Max near the center
          mx = max(smsh[gdpeak[nearobj]], sci)
          sci = nearobj[sci]
          off_sci = gdpeak[sci]-cen_ref[medj,mm]
          ;; Fraction
          sci_frac = (gdpeak[sci] - ordr_str[mm].lhedg[medj]) / $
            (ordr_str[mm].rhedg[medj]-ordr_str[mm].lhedg[medj])
      endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Setup objstr (default aper to 11)
      if flg_objstr EQ 0 then begin
          print, 'mike_fndobj: Tracing..'

          ;; Science first
          if not keyword_set( INTER ) then begin
              ;; AUTO
              if keyword_set( STD ) then objstr[0:9].aper = stdaper $
              else objstr[0:nordr-1].aper = objaper
          endif else begin
              ;; INTER
              objstr[0:nordr-1].aper = tmpobj[sci].aper 
          endelse
          ;; xcen
          objstr[0:nordr-1].obj_id = 'a'
          objstr[0:nordr-1].slit_id = ordr_str.order
print, mm
          for qq=mm,nordr-1 do begin ;; GOING UP
              ;; TRACE
              if qq EQ mm then begin
                  off = off_sci 
                  frac_sv = sci_frac
                  medx = medj
              endif else begin
                  frac = mike_fndobj_off(ordr_str, img, cen_ref, $
                                        qq, MEDJ=medx,$
                                        SCICLM=sciclm, $
                                        RBIN=mike[indx[exp[q]]].rowbin)
                  if frac EQ -1. or abs(frac-frac_sv) GT 0.1 then frac = frac_sv
                  frac_sv = frac
              endelse
              ;; xcen, ycen
              objstr[qq].xcen = medx
              objstr[qq].ycen = cen_ref[medx,qq] + off
              ;; TRACE
              objstr[qq].trace[0:sz_img[1]-1] = $
                (ordr_str[qq].rhedg-ordr_str[qq].lhedg)*frac_sv + $
                ordr_str[qq].lhedg
          endfor
          ;; Going down
          frac_sv = sci_frac
          for qq=mm-1,0L,-1 do begin
              if not keyword_set( INTER ) then begin
                  ;; AUTO
                  if keyword_set( STD ) then objstr[qq].aper = stdaper $
                  else objstr[qq].aper = objaper
              endif else begin
                  ;; INTER
                  objstr[qq].aper[*] = tmpobj[sci].aper 
              endelse
              ;; TRACE
              frac = mike_fndobj_off(ordr_str, img, cen_ref, $
                                     qq, MEDJ=medx,$
                                     SCICLM=sciclm, $
                                     RBIN=mike[indx[exp[q]]].rowbin)
              if frac EQ -1. or abs(frac-frac_sv) GT 0.1 then frac = frac_sv
              frac_sv = frac
              ;; xcen, ycen
              objstr[qq].xcen = medx
              objstr[qq].ycen = cen_ref[medx,qq] + off
              ;; TRACE
              objstr[qq].trace[0:sz_img[1]-1] = $
                (ordr_str[qq].rhedg-ordr_str[qq].lhedg)*frac_sv + $
                ordr_str[qq].lhedg
          endfor
          nobj = nobj + 1
          ;; CHK
          if keyword_set( CHK ) then begin
              tmp = img
              for qq=0L,nordr-1 do begin
                  rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
                  trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
                  tmp[trc_msk] = -10000
              endfor
              print, objstr[mm].xcen, objstr[mm].ycen
              xatv, tmp, /block, min=-70, max=770
          endif

      endif
      

      ;; Write Obj structure
      print, 'mike_fndobj: Creating ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'mike_fndobj: All done! '
  return
end
