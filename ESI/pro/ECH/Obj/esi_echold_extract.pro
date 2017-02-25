;+ 
; NAME:
; esi_echold_extract
;     Version 1.1
;
; PURPOSE:
;    Extract flux from 2D image to create ten 1D spectra (1 per order)
;    Output is written to the object structure (e.g. Extract/Obj_esi0024.fits)
;    The code only does boxcar extraction for now.
;
; CALLING SEQUENCE:
;   
;  esi_ech_oldextract, esi, obj_id, [exp], /DEBUG, /CHK, /STD, APER=,
;  RADIUS=
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
;   /CHK    - Show intermidate diagnostics (e.g. profile)
;   /SPCCHK - Show extracted spectra
;   /STD    - Extraction should be set for a standard star
;   /DEBUG  - Stop within extraction routine to check stuff
;   APER=   - Set aperture by hand (e.g. [5., 7.] )
;   RADIUS= - Size of window for setting aperture size (default: 20L)
;   ORDRS=  - Orders to extract (default: [0L,9L])
;   /OPTIMAL -  Perform an optimal extraction
;   INPGAU  - Gaussian sigam to use for optimal [default: Code
;             calculates the best value]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  1)  The program begins extracting in order 9L (physical 6) and will
;  automatically calculate an aperture for that order.  If there is
;  insufficient flux in the following orders, it will adopt the value
;  from the next higher order. 
;
; EXAMPLES:
;   esi_echextobj, esi, 1L, [0L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Aug-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echold_extract, esi, obj_id, exp, DEBUG=debug, CHK=chk $
                   , STD = std, APER = aper, RADIUS = radius, ORDRS = ordrs1  $
                   , SEDG_FIL = sedg_fil, OPTIMAL = optimal, INPGAU = inpgau $
                   , FIT_APER = FIT_APER $
                   , SPCCHK = spcchk, MAXGAU = maxgau, _EXTRA = extra $
                   , CBIN = cbin, RBIN = rbin $
                   , LSLITE = LSLITE, RSLITE = RSLITE $
                   , OBJ_NM = OBJ_NM, OUTGAU = OUTGAU

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_ech_oldextract, esi, obj_id, [exspr], RADIUS=, APER=, /DEBUG, /CHK'
      print, '          /STD, ORDRS=, /SPCCHK, /OPTIMAL [v1.1]'
      return
  endif 

  ;; Added by JFH to prevent rejection from standard stars
  IF KEYWORD_SET(STD) THEN BEGIN
      REJ_CR = 0
      NOCRREJ = 1
  ENDIF ELSE BEGIN
      REJ_CR = 1
      NOCRREJ = 0
  ENDELSE
  ;;  Optional Keywords
  ;;  if not keyword_set( EXTREG ) then extreg = [25., 25.]
  if not keyword_set( RADIUS ) then radius = 20L
  if n_elements(ORDRS1) EQ 0 THEN ordrs = lindgen(10) $
  ELSE ordrs = ordrs1
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( LSLITE ) then lslite = round(17./cbin)
  if not keyword_set( RSLITE ) then rslite = round(9./cbin)
  alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
  plate_scale = reverse([0.168, 0.163, 0.158, 0.153, 0.149, 0.144, 0.137 $
                         , 0.134, 0.127, 0.120])
  wave_ordr   = reverse([9800.0, 8750.0, 7800.0, 6850.0, 6200.0, 5650.0 $
                         , 5200.0, 4800.0, 4500.0, 4200.0])
  ;; make sure orders is sorted
  nord = n_elements(ordrs)
  outgau = fltarr(nord)
  icheck = WHERE(sort(ordrs) - lindgen(nord) NE 0, nbad)
  IF nbad GT 0 THEN message, 'The numbers in ordrs must be sorted'
  IF KEYWORD_SET(INPGAU) THEN BEGIN
     IF n_elements(inpgau) EQ nord THEN inpgau_vec = inpgau $
     ELSE IF n_elements(inpgau) EQ 1 THEN inpgau_vec = replicate(inpgau, nord) $
     ELSE message, 'invalid nubmer of elements in inpgau'
 ENDIF
;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_echextobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STD star
      indx = obj_id[0]
      nindx = 1L
;      extreg = [45., 45.]
      radius = 40L
  endelse

;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

;  Wavelength solution (put in a file)
  cdelt = 0.00001447624d
  tot_wave = 10^(alog10(3900.d) + dindgen(33000L)*cdelt)
  all_crval1 = dblarr(10)
  all_crval1[0] = 3900.
  all_crval1[1] = 4200.
  all_crval1[2] = 4450.
  all_crval1[3] = 4675.
  all_crval1[4] = 5100.
  all_crval1[5] = 10^3.74818803
  all_crval1[6] = 10^3.79448805
  all_crval1[7] = 7000.
  all_crval1[8] = 8000.
  all_crval1[9] = 9300.
  for qq=0L,9 do begin
      mn = min(abs(all_crval1[qq]-tot_wave),imn)
      all_crval1[qq] = alog10(tot_wave[imn])
  endfor

  all_coll = dblarr(10,2)
  all_coll[0,*] = [4100., 4250.]
  all_coll[1,*] = [4300., 4490.]
  all_coll[2,*] = [4600., 4850.]
  all_coll[3,*] = [5100., 5450.]
  all_coll[4,*] = [5500., 5850.]
  all_coll[5,*] = [5900., 6290.]
  all_coll[6,*] = [6600., 6850.]
  all_coll[7,*] = [7390., 7650.]
  all_coll[8,*] = [8650., 8900.]
  all_coll[9,*] = [9500., 9900.]
  all_mnxwv = dblarr(10,2)
  all_mnxwv[0,*] = [3900., 4380.]
  all_mnxwv[1,*] = [4200., 4700.]
  all_mnxwv[2,*] = [4450., 5060.]
  all_mnxwv[3,*] = [4660., 5500.]
  all_mnxwv[4,*] = [5055., 5980.]
  all_mnxwv[5,*] = [5618., 6572.]
  all_mnxwv[6,*] = [6230., 7300.]
  all_mnxwv[7,*] = [7000., 8210.]
  all_mnxwv[8,*] = [8000., 9380.]
  all_mnxwv[9,*] = [9330., 10200.]

; Open Slit file
  if not keyword_set( SEDG_FIL ) then $
    sedg_fil = esi_getfil('sedg_fil', SLIT=esi[indx[0]].slit, cbin=cbin, $
                          rbin=rbin, /name)
  if x_chkfil(sedg_fil+'*') EQ 0 then begin
      print, 'esi_echextobj: Slit edge file doesnt exist: ', sedg_fil
      return
  endif
  print, 'esi_echechextobj: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop

  for q=0L,n_elements(exp)-1 do begin

      print, 'esi_echextobj: Reading files...'

      ;;;;;;;;;;;;;;
      ;; Open Obj file
      objfil = esi[indx[exp[q]]].obj_fil
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'esi_echextobj: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
;      if not keyword_set(NOBJ) then 
      IF KEYWORD_SET(OBJ_NM) THEN BEGIN
         nobj = n_elements(obj_nm)
         alpha_redux = obj_nm
      ENDIF ELSE BEGIN
         nobj = n_elements(objstr)/10L
         alpha_redux = alphabet[lindgen(nobj)]
      ENDELSE
      ;;;;;;;;;;;;;;
      ;; SKY SUB Fil 
      imgfil = esi_getfil('fin_fil', subfil=esi[indx[exp[q]]].img_root,/name)
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'esi_echextobj: No Image file!  Returning...'
          return
      endif
      print, 'esi_echextobj: Image file -- ', imgfil
      head = xheadfits(imgfil)
      ;; IF TRULY SUPER OLD, THEN THE FOLLOWING IS OK
      img = xmrdfits(imgfil, 2, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      arc_fil = strtrim(esi[indx[exp[q]]].arc_fil,2)
      print, 'esi_echextobj: Arc -- ', arc_fil
      if x_chkfil(arc_fil+'*') EQ 0 then begin
          print, 'esi_echextobj: No Arc file!  Returning...', arc_fil
          return
      endif
      img_arc = xmrdfits(arc_fil, /silent) 
      
      ;; Read Sky if OPTIMAL
      if keyword_set( OPTIMAL ) then begin
          proc = xmrdfits(imgfil, 0, /silent)
          sky = proc - img
          delvarx, proc
      endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  HELIO correction
; MR fixed to work with helio update based on JXP 2009 fix
      x_radec, esi[indx[exp[q]]].ra, esi[indx[exp[q]]].dec, radeg, decdeg
      helio = (-1.)*x_keckhelio(radeg, decdeg $
                          , esi[indx[exp[q]]].equinox $
                          , jd = esi[indx[exp[q]]].date)
      hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
      img_arc = img_arc * hel_corr
      sxaddpar, head, 'HELIO', helio
      print, 'esi_echextobj: Helio correction applied -- ', helio, hel_corr
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     LOOP ON ORDERS
      print, systime()
    FOR iobj = 0L, nobj-1L DO BEGIN
        for iorder = 0L, nord-1L DO BEGIN
        qq = ordrs[nord-1L-iorder]
;        qq = ordrs[1], ordrs[0], -1 do begin
        ;;
          print, '-----------------------------------------------'
          print, 'esi_echextobj: Extracting order ' $
                 , string(15L-qq, FORMAT = '(i3)') +  $
                 ' with  wave_cen = ' +  $
                 string(wave_ordr[qq], FORMAT = '(I4)')
          print, '               Object = '  + alpha_redux[iobj]
          ;; Grab sci obj
          sci = where(objstr.obj_id EQ alpha_redux[iobj] $
                      AND objstr.order EQ qq $
                      , COMPLEMENT = b, NCOMPLEMENT = nb)
          ;; Mask Sky
          msk = lonarr(sz_img[0],sz_img[1])
          lhs = (rnd_edg[*, qq, 0]+LSLITE) > 0L ;LHS has funny edge (0.5" only?)
          rhs = (rnd_edg[*, qq, 1]-RSLITE) < (sz_img[0]-1)
;          lhs = (rnd_edg[*,qq,0]+17L) > 0L  ;; LHS has funny edge (0.5" only?)
;          rhs = (rnd_edg[*,qq,1]-9L) < (sz_img[0]-1)
          mni = min(round(lhs))
          mxi = max(round(rhs))
          a = where(lhs LE rhs, na)
          for ii=0L,na-1 do begin
              j = a[ii]
              msk[lhs[j]:rhs[j],j] = 1
          endfor
          
          ;; Transpose
          timg = transpose(img[mni:mxi,*]*msk[mni:mxi,*])
          tarc = transpose(img_arc[mni:mxi,*]*msk[mni:mxi,*])
          ;;          tarc = transpose(10.^img_arc[mni:mxi,*]*msk[mni:mxi,*])
          tvar = transpose(var[mni:mxi,*]*msk[mni:mxi,*])
          if keyword_set( OPTIMAL ) then tsky = $
            transpose(sky[mni:mxi,*]*msk[mni:mxi,*])
          
          ;; Run x_extobjbox
          print, 'esi_echextobj: Extracting...'
          if qq EQ 9 then rejsig = 15. else rejsig = 7.
          if not keyword_set( OPTIMAL ) then begin
              ;; Standard
              if keyword_set(SKY) then aper = objstr[sci].aper

              ;; Boxcar
              x_extobjbox, timg, tarc, $
                [objstr[sci].xcen, objstr[sci].ycen-mni], fin_spec, $
                VAR=tvar, WVMNX=all_mnxwv[qq,*], $
                APER=aper, DEBUG=debug, COLLMNX=all_coll[qq,*], $
                CRVAL1=all_crval1[qq], CDELT=cdelt, NPIX=5000L, $
                REJ_CR = REJ_CR, TOT_TRC=objstr[sci].trace[0:sz_img[1]-1]-mni, $
                REJSIG=rejsig, /REBINC, RADIUS=radius, BKAPER=svaper, CHK=chk, $
                NOCRREJ = NOCRREJ
              ;; Save the aperture
              svaper=fin_spec.aper
          endif else begin
              ;; Optimal 
              IF KEYWORD_SET(INPGAU) THEN INPGAU1 = INPGAU_VEC[nord-1L-iorder]
              x_extobjopt, timg, tarc, tvar, tsky, $
                [objstr[sci].xcen, objstr[sci].ycen-mni], fin_spec, $
                APER=aper, DEBUG=debug, COLLMNX=all_coll[qq,*], $
                CRVAL1=all_crval1[qq], CDELT=cdelt, NPIX=5000L, $
                TOT_TRC=objstr[sci].trace[0:sz_img[1]-1]-mni, $
                REJSIG=rejsig, /REBINC, RADIUS=radius, $
                READNO=esi[indx[exp[q]]].readno, DEFGAU=defgau, $
                FIT_APER = FIT_APER, USEGAU = INPGAU1, MAXGAU = maxgau, $
                CHK = chk
                outgau[nord-1L-iorder] = fin_spec.gauss_sig
              print, '             Seeing FWHM = ' $
                     + string(fin_spec.gauss_sig*plate_scale[qq]*2.3548D $
                              , FORMAT = '(F5.2)')
              ;; Help with Gaussian
              if not keyword_set(DEFGAU) then defgau = inpgau1
              if iorder EQ 0 then svgau = defgau $
;              if qq EQ ordrs[1] then svgau = defgau $
              else begin
                  maxgau = (defgau/svgau + 0.2)*defgau
                  svgau = defgau
              endelse
          endelse
          delvarx, timg, tarc, tvar


          ;; CHK
          if keyword_set( SPCCHK ) AND fin_spec.npix NE 0 then begin
              x_splot, fin_spec.wv, fin_spec.fx, YTWO=sqrt(fin_spec.var > 0.), $
                /block
          endif

          ;; Write to structure
          if fin_spec.npix NE 0 then begin
              objstr[sci].npix = fin_spec.npix
              objstr[sci].wave[0:fin_spec.npix-1] = fin_spec.wv
              objstr[sci].fx[0:fin_spec.npix-1] = fin_spec.fx
              objstr[sci].var[0:fin_spec.npix-1] = fin_spec.var
              objstr[sci].novar[0:fin_spec.npix-1] = fin_spec.novar
              objstr[sci].sky[0:fin_spec.npix-1] = fin_spec.sky
              ;; Box
              if not keyword_set(OPTIMAL) then begin
                  objstr[sci].box_wv[0:fin_spec.npix-1] = fin_spec.wv
                  objstr[sci].box_fx[0:fin_spec.npix-1] = fin_spec.fx
                  objstr[sci].box_var[0:fin_spec.npix-1] = fin_spec.var
              endif
              ;; Aper
              objstr[sci].aper = fin_spec.aper
              objstr[sci].gauss_sig = fin_spec.gauss_sig
              ;; Flag
              objstr[sci].flg_anly = 1
          endif else objstr[sci].flg_anly = 0
      endfor
  endfor
  print, systime()
  ;; Ouptut Spectra
  print, 'esi_echextobj: Output spectrum in -- ', objfil
  mwrfits, objstr, objfil, /create, /silent
  spawn, 'gzip -f '+objfil
endfor


  
;  DONE
  print, 'esi_echextobj: All done! '
  return
end
