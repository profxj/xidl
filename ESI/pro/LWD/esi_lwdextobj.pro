;+ 
; NAME:
; esi_lwdextobj   
;     Version 1.0
;
; PURPOSE:
;    Sky Subtract image
;
; CALLING SEQUENCE:
;   
;  esi_lwdextobj, esi, obj_id
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdextobj, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Jul-2002 Written by JXP
;   28-Aug-2002 Revised (uses trace)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdextobj, esi, obj_id, expsr, EXTREG=extreg, REFWV=refwv, APER=aper, $
                   DEBUG=debug, STD=std, CHK=chk

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdextobj, esi, obj_id, [exspr], EXTREG=, APER=, /CHK, /STD [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( EXTREG ) then extreg = [100., 100.]
  if not keyword_set(REFWV) then refwv = [6550., 6650.]
  if not keyword_set(WVMNX) then wvmnx = [3900., 10000.]
  if not keyword_set(WVSOL) then $
    wvsol = getenv('XIDL_DIR')+'/ESI/pro/LWD/Arcs/LWD_wvsol.fits'

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
                   esi.obj_id EQ obj_id AND esi.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_lwdextobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Wavelength solution
  newwv = xmrdfits(wvsol, /dscale, /silent)

;  Loop

  for q=0L,n_elements(exp)-1 do begin

      print, 'esi_lwdextobj: Reading files...'

      ;; Open Obj file
      objfil = esi[indx[q]].obj_fil
      a = findfile(objfil, count=na)
      if na EQ 0 then begin
          print, 'esi_lwdextobj: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)

      ;; SKY SUB Fil 
      imgfil = objstr[0].spec2d_fil
      img = xmrdfits(imgfil, 2, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      ivar = 1./(var>1)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      img_arc = xmrdfits(strtrim(esi[indx[exp[q]]].arc_fil,2), /silent) 


      ;;  HELIO correction
      helio = x_keckhelio(esi[indx[exp[q]]].RA, esi[indx[exp[q]]].dec, $
                          esi[indx[exp[q]]].equinox, jd=esi[indx[exp[q]]].date)
      hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
      img_arc = img_arc * hel_corr
      head = xheadfits(imgfil)
      sxaddpar, head, 'HELIO', helio
;      modfits, imgfil, 0, head
      print, 'esi_echextobj: Helio correction applied -- ', helio, hel_corr

; Parse out the Object
      sci = where(objstr.obj_id EQ 'a', COMPLEMENT=b, NCOMPLEMENT=nb)
      ;; Reset extreg
      if nb NE 0 then begin
          mnsep = min( objstr[b].trace[1600L] - objstr[sci].trace[1600L], $
                       max=mxsep)
          if mnsep LT 0. then begin
              if mxsep GT 0. then extreg[0] = extreg[0] < abs(mnsep) $
              else extreg[0] = extreg[0] < abs(mxsep) 
          endif
          if mxsep GT 0. then begin
              if mnsep GT 0. then extreg[1] = extreg[1] < abs(mnsep) $
              else extreg[1] = extreg[1] < abs(mxsep) 
          endif
      endif

      ;; Mask
      jmn = 0 > min(objstr[sci].trace[0:sz_img[0]-1] - extreg[0])
      jmx = (sz_img[1]-1) < max(objstr[sci].trace[0:sz_img[0]-1] + extreg[1])
      mnj = round(jmn)
      mxj = round(jmx)

      ;; Run x_extobjbox
      print, 'esi_lwdextobj: Extracting...'
      x_extobjbox, img[*,mnj:mxj], img_arc[*,mnj:mxj], $
        [objstr[sci].xcen, objstr[sci].ycen-mnj], fin_spec, $
        VAR=var[*,mnj:mxj], WVMNX=wvmnx, $
        NEWWV=newwv, APER=aper, DEBUG=debug, COLLMNX=[4400.,7000.], $
        TOT_TRC=objstr[sci].trace[0:sz_img[0]-1]-mnj, /REBINC, /REJ_CR

      ;; Write to structure
      if fin_spec.npix NE 0 then begin
          objstr[sci].npix = fin_spec.npix
          objstr[sci].wave[0:fin_spec.npix-1] = fin_spec.wv
          objstr[sci].fx[0:fin_spec.npix-1] = fin_spec.fx
          objstr[sci].var[0:fin_spec.npix-1] = fin_spec.var
          ;; Trace
          ntrc = n_elements(fin_spec.trc)
          objstr[sci].trace[0:ntrc-1] = $
            fin_spec.trc+replicate(mnj, ntrc) ; Offset!
          ;; Aper
          objstr[sci].aper = fin_spec.aper
          ;; Flag
          objstr[sci].flg_anly = 1

          ;; CHK
          if keyword_set( CHK ) then x_splot, fin_spec.wv, fin_spec.fx, $
            YTWO=sqrt(fin_spec.var>0.), /block
      endif else objstr[objslit[obj]].flg_anly = 0

      ;; Ouptut Spectra
      print, 'esi_lwdextobj: Output spectrum in -- ', objfil
      mwrfits, objstr, objfil, /create, /silent
  endfor
  
;  DONE
  print, 'esi_lwdextobj: All done! '
  return
end
