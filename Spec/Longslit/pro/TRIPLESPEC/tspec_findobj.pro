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

FUNCTION tspec_findobj, img_minsky,sciivar,waveimg,tset_slits $
                        , FILSTD = FILSTD, CHK = chk $
                        , APERV = aperv, APERS = apers $
                        , NFIND = NFIND, PEAKTHRESH = PEAKTHRESH $
                        , FWHM = FWHM $
                        , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN
  
;;  Optional Keywords
  if not keyword_set( APERV ) then aperv = 10. 
  if not keyword_set( APERS ) then apers = 25. 
  IF NOT KEYWORD_SET(NFIND) THEN NFIND=1L
  cbin=1.0d

  alphabet = ['a', 'b', 'c', 'd', 'e', 'f','g','h','i']
  
  dim = size(img_minsky, /dim)
  nx = dim[0]
  ny = dim[1] 
  dimt = size(tset_slits.coeff, /dimen)
  norders = dimt[1]
  
  slit_edg = fltarr(ny, norders, 2)
  traceset2xy, tset_slits[0], rows, left_edge
  traceset2xy, tset_slits[1], rows, right_edge
  slit_edg[*, *, 0] = left_edge
  slit_edg[*, *, 1] = right_edge
  ordermask = tspec_ordermask(tset_slits)
  IF KEYWORD_SET(FILSTD) THEN BEGIN
     IF x_chkfil(filstd+'*') EQ 0 then begin
        print, 'tspec_findobj: STD Obj file does not exist or' + $
               ' obj_fil tag not set! Check it..', filstd
        return,0
     ENDIF ELSE stdstr = xmrdfits(filstd, 1, /silent)
     gdtrc = stdstr.trace[0:ny-1L,*]
  ENDIF ELSE gdtrc = (slit_edg[*, *, 0] + slit_edg[*, *, 1])/2.0d
  ;stop
  objstruct = tspec_echjoefind(img_minsky, sciivar*(ordermask GT 0) $
                               , waveimg, gdtrc, slit_edg $
                               , GFRAC = gfrac, FWHM = FWHM $
                               , PEAKTHRESH = PEAKTHRESH, NFIND = NFIND $
                               , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN)
  if size(objstruct,/type) EQ 2 then begin
    print, 'tspec_findobj:  No object.  Punting'
    return, -1
  endif
  nobj = n_elements(objstruct)/norders

  IF KEYWORD_SET(NFIND) THEN BEGIN
     peakflux = fltarr(nobj)
     FOR kk = 0L, nobj-1L DO BEGIN
        i1 = kk*norders
        i2 = (norders-1L) +kk*norders
        ;; Avg flux across 15 orders
        peakflux[kk] = total(objstruct[i1:i2].PEAKFLUX)/norders
     ENDFOR
     ;; Sort by average flux and keep nfind brightest objects
     peakind = reverse(sort(peakflux))
     keepind = lonarr(nfind*norders)
     FOR j = 0L, nfind-1L DO BEGIN
        i1 = j*norders
        i2 = (norders-1L) +j*norders
        keepind[i1:i2] = peakind[j]*norders + lindgen(norders)
     ENDFOR
     objstruct = objstruct[keepind]
     nobj = nfind
  ENDIF
  ;; Now sort by fractional slit position
  fracpos = fltarr(nobj)
  FOR kk = 0L, nobj-1L DO BEGIN
     i1 = kk*norders
     i2 = (norders-1L)+kk*norders
     ;; Avg flux across 15 orders
     fracpos[kk] = total(objstruct[i1:i2].XFRACPOS)/norders
  ENDFOR
  ;; Sort by average flux and keep nfind brightest objects
  fracind = sort(fracpos)
  keepind = lonarr(nobj*norders)
  FOR j = 0L, nobj-1L DO BEGIN
     i1 = j*norders
     i2 = (norders-1L)+j*norders
     keepind[i1:i2] = fracind[j]*norders + lindgen(norders)
  ENDFOR
  objstruct = objstruct[keepind]
  ;; Create objects
  tmp = {esiobjstrct}
  tmp2 = objstruct[0]
  dum = {junk:0.0}
  struct_assign, dum, tmp2
  proto = struct_addtags(tmp, tmp2)
  objstr = replicate(proto, nobj*norders)
  objstr.img_fil = ' '
  objstr.UT = ' '
  objstr.field = ' '
  ;;objstr.exp = esi[indx[exp[q]]].exp
  ;;objstr.arc_fil = strtrim(esi[indx[exp[q]]].arc_fil, 2)
  ;; Set aperture
  if keyword_set(STD) then objstr.aper[*] = APERS $
  else objstr.aper[*] = APERV
  
  for kk = 0L, nobj-1 do begin
     i1 = kk*norders
     i2 = (norders-1L)+kk*norders
     objstr[i1:i2].obj_id = alphabet[kk]
     objstr[i1:i2].order = lindgen(norders) 
     ;;objstr[i1:i2].trace[4999L] = gfrac[kk]
     medj = ny/2
     for jj = 0L,norders-1L do begin
        ;; xcen, ycen
        objstr[i1+jj].xcen = medj
        objstr[i1+jj].ycen = objstruct[i1+jj].xpos[medj]
        ;; Trace
        objstr[i1+jj].trace[0:ny-1L] = objstruct[i1+jj].xpos
        ;; copy over other relevant tags for these indices
        copy_struct_inx, objstruct, objstr $
                         , index_from = (i1+jj), index_to = (i1+jj)
     endfor
  endfor
  ;; Take only the relevant structures
  objstr = objstr[0:nobj*norders-1]
      
  if keyword_set(CHK) then begin
     nobj = n_elements(objstr)/norders
     tmp = img_minsky*(ordermask GT 0.0)
     FOR iobj = 0L, nobj-1L DO BEGIN
        FOR qq = 0, norders-1L DO BEGIN
           this = where(objstr.obj_id EQ alphabet[iobj] $
                        AND objstr.order EQ qq $
                        , COMPLEMENT = b, NCOMPLEMENT = nb)
           rnd_trc = round(objstr[this].trace[0:ny-1L])
           trc_msk = rnd_trc + lindgen(ny)*nx
           tmp[trc_msk] = -10000
        ENDFOR
     ENDFOR
     xatv, tmp, /block, min = -70, max = 70
  endif
  
  ;; Write Obj structure
  ;;print, 'esi_echfndobj: Creating ', objfil
  ;;mwrfits, objstr, objfil, /create, /silent
  ;;spawn, 'gzip -f '+objfil
  
  
;  DONE
  print, 'mage_findobj: All done! '
  return,objstr
end
