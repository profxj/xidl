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
;  esi_echfndobj, esi, obj_id, [exp], REFWV=, SCICLM=, PEAKTHRESH=,
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
;  MAXOFF=  - Minimum offset between max obj and center of slit
;             (default: 20.)
;  APERV=   - Aperture value (for masking sky)  [default: 10]
;  SCICLM=  - Guess at the position of the object relative to center
;  PEAKTHRESH = Fraction of brightest object on the slit to trace. 
;  ABSTHRESH  = Absolute threshold in units of peakflux to trace. 
;  MINPEAK  -- Absolute threshold required to be a prospective object
;              (default = 0.2)

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfndobj, esi, obj_id, exp, STD = std, CHK = chk $
                   , NOCLOB = noclob, SEDG_FIL = sedg_fil $
                   , CBIN = cbin, RBIN = rbin, APERV = aperv, APERS = apers $
                   , NFIND = NFIND, PEAKTHRESH = PEAKTHRESH, SKYSUB = SKYSUB $
                   , USESTD = USESTD, FILSTD = FILSTD, FWHM = FWHM $
                   , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN, FOFSEP = FOFSEP
; old keywords
;, FITFIL = fitfil $,refwv=refwv, SCICLM = sciclm, REFORDR = refordr
;, INTER = inter , FRACPK = FRACPK, MAXOFF = maxoff, NSIG = nsig, NEDG = NEDG, 


lbl = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

if  N_params() LT 2  then begin 
    print, 'Syntax - ' + $
           'esi_echfndobj, esi, obj_id, [exp]'
    return
endif 

;  Optional Keywords
  if not keyword_set(ORDRS) then ordrs = [0L, 9L]
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set(OBJ_ID) then obj_id = 0L
  if not keyword_set(REFORDR) then refordr = 4L
  if not keyword_set(REFWV) then refwv = [5600., 5610.]
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
;  if not keyword_set( NSIG ) then nsig = 3.
  if not keyword_set( APERV ) then aperv = 10. / float(cbin)
  if not keyword_set( APERS ) then apers = 25. / float(cbin)
  
;  Find all relevant obj images
  if not keyword_set(STD) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   esi.obj_id EQ obj_id AND strtrim(esi.type, 2) EQ 'OBJ', $
                   nindx)
      if nindx EQ 0 then begin
          print, 'esi_echfndobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin              ; STANDARD STAR
      indx = [obj_id[0]]
      nindx = 1L
  endelse
  alphabet = ['a', 'b', 'c', 'd', 'e', 'f']
  
; Open Slit file
  if not keyword_set( SEDG_FIL ) then $
     slit_edg_fil = esi_getfil('sedg_fil', SLIT = esi[indx[0]].slit $
                               , cbin = cbin, rbin = rbin, /name) $
  ELSE slit_edg_fil = strtrim(SEDG_FIL, 2)
  slit_edg = xmrdfits(slit_edg_fil, 0, /silent)
  tset_slits = xmrdfits(slit_edg_fil, 1, /silent)
  ordermask = esi_echordermask(tset_slits)
  ;; Open trace file to guide object finding and tracing. If USESTD is set
  ;; use the STD star already traced. Otherwise use an archived standard trace. 
  IF KEYWORD_SET(USESTD) THEN BEGIN
      IF NOT keyword_set(FILSTD) then BEGIN
          istd = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                       esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                       esi.slit EQ esi[indx[0]].slit AND $
                       strtrim(esi.type, 2) EQ 'STD', nstd)
          case nstd of 
              0 : begin
                  print, 'esi_echfndobj: No standard star images! Returning..'
                  return
              end 
              1 : print, 'esi_echfndobj: Tracing standard star image ', $
                esi[istd].img_root
              else : begin
                  print $ 
                    , 'esi_echfndobj: Warning -- Multiple standard star images'
                  istd = istd[0]
                  print, 'esi_echfndobj: Taking first one -- ' $
                         , esi[istd].img_root
              end
          endcase
          filstd = esi[istd].obj_fil
      ENDIF
      IF x_chkfil(filstd+'*') EQ 0 then begin
          print, 'esi_echfndobj: STD Obj file does not exist or' + $
                 ' obj_fil tag not set! Check it..', filstd
          stop
          return
      ENDIF ELSE stdstr = xmrdfits(filstd, 1, /silent)
      gdtrc = stdstr.trace
  ENDIF ELSE BEGIN
      if rbin NE 1 then stop
      ;trc_fil = getenv('XIDL_DIR')+'/ESI/CALIBS/ECH_TRC.fits'
      ;if x_chkfil(trc_fil+'*') EQ 0 then begin
      ;    print, 'esi_echfndobj: Slit edge file ', trc_fil, ' does not exist!'
      ;    trc_fil = getenv('ESI_CALIBS')+'/ECH_TRC.fits'
      ;    if x_chkfil(trc_fil+'*') EQ 0 then return
      ;    print, 'esi_echfndobj: Using... ', trc_fil
      ;endif
      ;;gdtrc = xmrdfits(trc_fil, /silent)
      ;; Binning
      ;;gdtrc = gdtrc/ float(cbin)
      dim = size(slit_edg, /dim)
      gdtrc = fltarr(5000, dim[1])
      gdtrc[0:dim[0]-1L, *] = (slit_edg[*, *, 0] + slit_edg[*, *, 1])/2.0d
  ENDELSE
  
; Hole trace
;      restore, fitfil           
  
;  Exposures
  if size(exp, /type) EQ 0 then exp = lindgen(nindx)
  nexp = n_elements(exp)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  FOR q = 0L, nexp-1 do begin
      ;; Read IMG+VAR
      imgfil = 'Final/f_'+esi[indx[exp[q]]].img_root
      a = findfile(imgfil+'*', count = na)
      if na EQ 0 then begin
          print, 'esi_echfndobj: No Final file ', imgfil
          stop
          continue
      endif
      print, 'esi_echfndobj: Opening files... ', imgfil
      IF KEYWORD_SET(SKYSUB) THEN BEGIN
          img1 = xmrdfits(imgfil, 0, /silent)
          sky = xmrdfits(imgfil, 2, /silent)
          img = img1-sky
      ENDIF ELSE img = xmrdfits(imgfil, /silent) 
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)
      var = var*(ordermask GT 0.0) ;; explicit masking of bad regions
      ;; Read ARC
      flg_arc = 0
      if q EQ 0 then begin
          img_arc = xmrdfits(strtrim(esi[indx[exp[q]]].arc_fil, 2), /silent) 
      endif else begin
          if esi[indx[exp[q]]].arc_fil NE esi[indx[exp[q-1]]].arc_fil then $
            img_arc = xmrdfits(esi[indx[exp[q]]].arc_fil, /silent) $
          else flg_arc = 1
      endelse
      
      ;; Mask
      if flg_arc NE 1 then begin
          msk = bytarr(sz_img[0], sz_img[1])
          sva = where(img_arc GT refwv[0] AND img_arc LT refwv[1])
          msk[sva] = 1
       endif
      ;; JFH Sep 11, 2016
      ;;objstruct = esi_echjoefind(img, var, gdtrc, slit_edg, cbin $
      ;;                           , GFRAC = gfrac, FWHM = FWHM $
      ;;                           , PEAKTHRESH = PEAKTHRESH, NFIND = NFIND $
      ;;                           , ABSTHRESH = ABSTHRESH,
      ;;                           MINPEAK=minpeak)
      ;print, "FOFSEP, MIN_SN", FOFSEP
      objstruct = mage_echjoefind(img, var, img_arc, gdtrc, slit_edg $
                                  , FITSTD = FITSTD $
                                  , CBIN = CBIN, GFRAC = gfrac, FWHM = FWHM $
                                  , PEAKTHRESH = PEAKTHRESH, NFIND = NFIND $
                                  , TRCMASK = TRCMASK $
                                  , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN $
                                  , FOFSEP = FOFSEP, /ESI)
      
      nobj = n_elements(objstruct)/10
      IF KEYWORD_SET(NFIND) THEN BEGIN
         peakflux = fltarr(nobj)
         FOR kk = 0L, nobj-1L DO BEGIN
            i1 = kk*10
            i2 = 9+kk*10
            ;; Avg flux across 10 orders
            peakflux[kk] = total(objstruct[i1:i2].PEAKFLUX)/10
         ENDFOR
         ;; Sort by average flux and keep nfind brightest objects
         peakind = reverse(sort(peakflux))
         keepind = lonarr(nfind*10)
         FOR j = 0L, nfind-1L DO BEGIN
            i1 = j*10
            i2 = 9+j*10
            keepind[i1:i2] = peakind[j]*10 + lindgen(10)
         ENDFOR
         objstruct = objstruct[keepind]
         nobj = nfind
     ENDIF
     
     ;; Now sort by fractional slit position
     fracpos = fltarr(nobj)
     FOR kk = 0L, nobj-1L DO BEGIN
         i1 = kk*10
         i2 = 9+kk*10
         ;; Avg flux across 10 orders
         fracpos[kk] = total(objstruct[i1:i2].XFRACPOS)/10
     ENDFOR
     ;; Sort by average flux and keep nfind brightest objects
     fracind = sort(fracpos)
     keepind = lonarr(nobj*10)
     FOR j = 0L, nobj-1L DO BEGIN
         i1 = j*10
         i2 = 9+j*10
         keepind[i1:i2] = fracind[j]*10 + lindgen(10)
     ENDFOR
     objstruct = objstruct[keepind]
     
     objfil = 'Extract/Obj_'+esi[indx[exp[q]]].img_root
      esi[indx[exp[q]]].obj_fil = objfil
      ;; Create objects
      tmp = { esiobjstrct }
      tmp2 = objstruct[0]
      dum = {junk:0.0}
      struct_assign, dum, tmp2
      proto = struct_addtags(tmp, tmp2)
      objstr = replicate(proto, nobj*10)
      objstr.img_fil = ' '
      objstr.UT = ' '
      objstr.field = ' '
      objstr.exp = esi[indx[exp[q]]].exp
      objstr.arc_fil = strtrim(esi[indx[exp[q]]].arc_fil, 2)
      ;; Set aperture
      if keyword_set(STD) then objstr.aper[*] = APERS $
      else objstr.aper[*] = APERV
      
      for kk = 0L, nobj-1 do begin
          i1 = kk*10
          i2 = 9+kk*10
          objstr[i1:i2].obj_id = lbl[kk]
;          objstr[i1:i2].slit_id = lindgen(10)
          objstr[i1:i2].order = lindgen(10)
          objstr[i1:i2].trace[4999L] = gfrac[kk]
          medj = sz_img[1]/2
          for jj = 0L, 9 do begin
              ;mt = where(abs(objstruct.xfracpos-gfrac[kk]) LT 0.03 AND $
              ;           objstruct.slitid EQ (jj+1), nmt)
              ;if nmt NE 1 then stop
              ;; xcen, ycen
              objstr[i1+jj].xcen = medj
              objstr[i1+jj].ycen = objstruct[i1+jj].xpos[medj]
              ;; Trace
              objstr[i1+jj].trace[0:sz_img[1]-1] = objstruct[i1+jj].xpos
              ;; copy over other relevant tags for these indices
              copy_struct_inx, objstruct, objstr $
                               , index_from = (i1+jj), index_to = (i1+jj)
          endfor
      endfor
      
      ;; Take only the relevant structures
      objstr = objstr[0:nobj*10L-1]
      
      if keyword_set(CHK) then begin
          nobj = n_elements(objstr)/10L
          tmp = img*(ordermask GT 0.0)
          FOR iobj = 0L, nobj-1L DO BEGIN
              FOR qq = ordrs[1], ordrs[0], -1 do begin
                  this = where(objstr.obj_id EQ alphabet[iobj] $
                               AND objstr.order EQ qq $
                               , COMPLEMENT = b, NCOMPLEMENT = nb)
                  igood = WHERE(trcmask[*, qq], ngood)
                  ;;rnd_trc = round(objstr[this].trace[0:sz_img[1]-1])
                  ;;trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
                  rnd_trc = round(objstr[this].trace[igood])
                  trc_msk = rnd_trc + igood*sz_img[0]
                  tmp[trc_msk] = -10000
              endfor
          ENDFOR
          xatv, tmp, /block, min = -70, max = 70
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

