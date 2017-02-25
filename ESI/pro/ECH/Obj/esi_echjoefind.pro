;+ 
; NAME:
; esi_echjoefind
;     Version 1.1
;
; PURPOSE:
;    Uses Joe+Schlegel routine to find+trace objects in
;      each order (independently)
;
; CALLING SEQUENCE:
;   objstruct = esi_echjoefind(img, var, gdtrc, slit_edg, cbin, FITSTD =
;                      , GFRAC = , FWHM = , PEAKTHRESH = , NFIND =  
;                      , ABSTHRESH =, )
;
; INPUTS:
;   img      -  Image to search for objects
;   var      -  Variance image  
;   gdtrc    -  Guide trace (standard)
;   slit_edg -  Araay with slit edges
;   CBIN     -  Column binning. 
;
; RETURNS:
;  objstruct - Object structure
;
; OPTIONAL KEYWORDS:
;
;  FWHM       - FWHM of profiles parameter for long_objfind
;  PEAKTHRESH - only reduce objects within this fraction of the brightest 
;               on each slit
;  ABSTHRESH  - only reduce objects above this absolute threshold (overrides 
;               peakthresh if both are given. 
;  NFIND      - Find only this number of objects on each slit. 
;
; OPTIONAL OUTPUTS:
;  FITSTD     - xy2traceset fit with of the standard trace
;  GFRAC      - [nobj] array of fractional positions of each object on the slit 
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-Sep-2006 Written by JH
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function esi_echjoefind, img, var, gdtrc, slit_edg, cbin, FITSTD = fitstd $
  , GFRAC = gfrac, FWHM = FWHM, PEAKTHRESH = PEAKTHRESH $
  , NFIND = NFIND, ABSTHRESH = ABSTHRESH, MINPEAK = minpeak
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' + $
             'ostr = esi_echjoefind( img, var, gdtrc, slit_edg, cbin )  [v1.0]'
      return, -1
  endif 

  if not keyword_set(PEAKTHRESH) then peakthresh = 0
  if not keyword_set(OBJSEP) then objsep = 0.07
  if not keyword_set(MINPEAK) then MINPEAK = 0.2
  IF NOT KEYWORD_SET(FWHM) THEN FWHM = 4.0/CBIN
;; ESI pixels are 0.15". This presumes seeing is 0.6"
  
  ymin_ord15 = 1500
  ymax_ord15 = 3700
  ymax_ord06 = 2170
;          obj = mrdfits(imgfile, 0)
;          var = mrdfits(imgfile, 1)
  ivar = (var GT 0.0)/(var + 3.0*(var EQ 0.0))
;          skysub = mrdfits(imgfile, 2)
;          skysubmask = (skysub NE 0.0)
;          std = mrdfits(stdfile, 1)
  dims = size(img, /dimen)
  ny = dims[1]
  
  ;; Set standard
  x_std = gdtrc[0:ny-1L, *]
  y_std = rebin(findgen(ny), ny, 10)

;; Mask the bad tracee (and off edge of CCD)
trcmask = lonarr(ny, 10) + 1L
trcmask[0:ymin_ord15, 0] = 0
trcmask[ymax_ord15:*, 0] = 0
trcmask[ymax_ord06:*, 9] = 0

zeroind = WHERE(x_std EQ 0.0, nbad)
IF nbad NE 0 THEN trcmask[zeroind] = 0

edg_arr = slit_edg
dums = size(edg_arr, /dim)
nslit = dums[1]

func = 'legendre'
ncoeff = 6
xy2traceset, y_std, x_std, stdset $
             , func = func, ncoeff = ncoeff $
             , invvar = float(trcmask), yfit = fitstd
tset_slits = [stdset, stdset]

medslit = fltarr(nslit, 2)
FOR iord = 0L, nslit-1L DO BEGIN
    djs_iterstat, fitstd[*, iord]-edg_arr[*, iord, 0], sigrej = 2.0 $
                  , median = left_med, mask = (trcmask[*,iord] EQ 1)
    medslit[iord, 0] = left_med
    djs_iterstat, edg_arr[*, iord, 1] - fitstd[*, iord], sigrej = 2.0 $
                  , median = right_med, mask = (trcmask[*,iord] EQ 1)
    medslit[iord, 1] = right_med
    tset_slits[0].coeff[0, iord] = tset_slits[0].coeff[0, iord] $
      - left_med
    tset_slits[1].coeff[0, iord] = tset_slits[1].coeff[0, iord] $
      + right_med
ENDFOR
tset_slits = struct_addtags(tset_slits $
                            , replicate(create_struct('DIMS', dims) $
                                        , size(tset_slits, /dimens)))

objstruct = long_objfind(img, tset_slits = tset_slits, invvar = ivar $
                         , ncoeff = ncoeff, skymask = skymask $
                         , objmask = objmask $
                         , PEAKTHRESH = peakthresh, FWHM = FWHM $
                         , NPERSLIT = NFIND, ABSTHRESH = ABSTHRESH $
                         , CLIP_SLIT=[10])

;; Error checking
nobj = n_elements(objstruct)
msk = intarr(nobj)
for jj = 0L, nobj-1 do begin
    ;; Slit
    is = objstruct[jj].slitid - 1
    
    ;; Check RHS/LHS
    bd = where(objstruct[jj].xpos LT slit_edg[*, is, 0] OR $
               objstruct[jj].xpos GT slit_edg[*, is, 1] AND $
               trcmask[*, is] EQ 1, nbad)
    if nbad NE 0 then begin
          print, 'esi_echjoefind: Bad trace in slit ', is
          msk[jj] = 1 
      endif
      ;; Expunge very faint objects (for multiple detections)
      sobj = where(objstruct[jj].slitid EQ objstruct.slitid, nslit)
      if nslit GT 1 then begin
          bad = where(objstruct[sobj].peakflux LT MINPEAK, nbad)
          if nbad NE 0 then begin msk[sobj[bad]] = 1
              print, 'esi_echjoefind: Too faint in slit', is, $
                     objstruct[sobj[bad]].peakflux, minpeak
          endif
      endif
   endfor

  objstruct = objstruct[where(msk EQ 0)]
  ;; Book-keeping. Identify unique objects by running friends of friends 
  ;; on their fractional positions 
  lbl = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

  x = objstruct.xfracpos/1000.0D
  y = 0.0*x
  group = spheregroup(x, y, OBJSEP/1000.0D, mult = mult)
  ;; Require that objects occur on at least two orders to be considered. 
  ;; If you are interested in emission line sources with no detectable 
  ;; continuum this should be changed. 
  keep = WHERE(mult[group] GE 2) 
  objstruct = objstruct[keep]
  group = group[keep]
  objstruct.objid = lindgen(n_elements(objstruct)) + 1
  nobjs = n_elements(objstruct)

  uni_ind = uniq(group, sort(group))
  uni_group = group[uni_ind]
  nobj = n_elements(uni_group)
  gfrac = fltarr(nobjs)
  FOR j = 0L, nobj-1L DO BEGIN
      ind = WHERE(group EQ uni_group[j])
      this_frac = djs_median(objstruct[ind].XFRACPOS)
      gfrac[ind] = this_frac
  ENDFOR
  uni_frac = gfrac[uni_ind]
  
;   for kk = 0L, nobjs-1 do begin
;       if kk EQ 0 then begin
;           gfrac = [objstruct[kk].xfracpos]
;           nobj = 1
;       endif else begin
;           mt = where(abs(gfrac-objstruct[kk].xfracpos) LT OBJSEP, nmt)
;           if nmt EQ 0 then begin
;               gfrac = [gfrac, objstruct[kk].xfracpos]
;               nobj = nobj + 1
;           endif
;       endelse
;   endfor
;   gfrac = gfrac[sort(gfrac)]
;   for jj = 0L, nobj-1 do begin
;       mt = where(abs(gfrac[jj]-objstruct.xfracpos) LT OBJSEP, nmt)
;       gfrac[jj] = djs_median(objstruct[mt].xfracpos)
;   endfor
  ;; Require that objects occur on at least two orders to be considered. 
  ;; If you are interested in emission line sources with no detectable 
  ;; continuum this should be changed. 
 ;  msk = intarr(nobjs)
;   gmsk = intarr(nobj)
;   FOR jj = 0L, nobj-1L DO BEGIN
;       mult = where(abs(gfrac[jj]-objstruct.xfracpos) LT OBJSEP, nmult)
;       IF nmult EQ 0 THEN stop $
;       ELSE IF nmult LT 2 THEN BEGIN
;           msk[mult] = 1
;           gmsk[jj] = 1
;       ENDIF
;   ENDFOR
  
  ;; Fill in traces
  for kk=10,1,-1 do begin
      for ss=0L,nobj-1 do begin
          ;; Is it there?
          mt = where(group EQ uni_group[ss] AND objstruct.slitid EQ kk, nmt)
          if nmt EQ 0 then begin ;; Nope
              print, 'esi_echjoefind: Extrapolating object ', lbl[ss], ' in '+$
                     'order ',kk-1
              ;; Find closest
              this_group = WHERE(group EQ uni_group[ss])
              mn = min(abs(objstruct[this_group].slitid-kk), imn)
              cslit = objstruct[this_group[imn]].slitid-1 ;; Closest slit
              ;; Calculate frac
              ledg = fitstd[*,cslit]-medslit[cslit,0]
              redg = fitstd[*,cslit]+medslit[cslit,1]
              frac = (objstruct[this_group[imn]].xpos-ledg)/(redg-ledg)
              ;; Find new pos
              ledg = fitstd[*,kk-1]-medslit[kk-1,0]
              redg = fitstd[*,kk-1]+medslit[kk-1,1]
              newpos = ledg + (redg-ledg)*frac
              ;; Fill up
              tmp = objstruct[0]
              tmp.objid = n_elements(objstruct)+1
              tmp.slitid = kk
              tmp.xfracpos = objstruct[this_group[imn]].xfracpos
              tmp.xpos = newpos
              ;;
              objstruct = [objstruct, tmp]
              group = [group, uni_group[ss]]
              gfrac = [gfrac, uni_frac[ss]]
          endif
      endfor
  endfor
  ;; Some code to make sure the objstruct is sorted by fractional 
  ;; position and by slit, respectively. 
  ind_sort = sort(uni_frac)
  FOR k = 0L, nobj-1L DO BEGIN
      this_frac = WHERE(gfrac EQ uni_frac[ind_sort[k]])
      this_obj = objstruct[this_frac]
      this_obj = this_obj[sort(this_obj.SLITID)]
      IF NOT KEYWORD_SET(OBJOUT) THEN objout = [this_obj] $
      ELSE objout = [objout, this_obj]
  ENDFOR

  ;; Did we find all the objects the user desired?
  if keyword_set(NFIND) then begin
      if NOBJ LT NFIND then begin
          print, 'esi_echjoefind:  Found fewer objects than you desired.', $
                 nfind, nobj
          print, 'esi_echjoefind:  Consider setting MINPEAK to a value lower'
          print, 'esi_echjoefind:  than ', MINPEAK
          nfind = nobj
      endif
  endif
  
  return, objout
end
