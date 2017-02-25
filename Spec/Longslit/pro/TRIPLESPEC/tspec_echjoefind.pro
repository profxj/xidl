;+ 
; NAME:
; tspec_echjoefind
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
function tspec_echjoefind, img, sciivar, waveimg,gdtrc, slit_edg $
                          , FITSTD = fitstd $
                          , GFRAC = gfrac, FWHM = FWHM $
                          , PEAKTHRESH = PEAKTHRESH $
                          , NFIND = NFIND, ABSTHRESH = ABSTHRESH $
                          , MIN_SN = MIN_SN
  if  N_params() LT 4  then begin 
      print, 'Syntax - ' + $
             'ostr = esi_echjoefind( img, var, gdtrc, slit_edg, cbin )  [v1.0]'
      return, -1
  endif 

  if not keyword_set(PEAKTHRESH) then peakthresh = 0
  if not keyword_set(STDSEP) then STDSEP = 1.0d
  if not keyword_set(MIN_SN) then MIN_SN = 0.5
  IF NOT KEYWORD_SET(FWHM) THEN FWHM = 2.0
   ;; MagE pixels are 0.3". This FWHM presumes seeing is 0.6"
  cbin=1.0
  if not keyword_set( LSLITE ) then lslite = round(3./cbin)
  if not keyword_set( RSLITE ) then rslite = round(3./cbin)
  IF NOT KEYWORD_SET(BOX_RAD1) THEN box_rad = 10L/cbin
  SN_BOX_RAD=BOX_RAD/2.0 
  OBJ_SHIFT=BOX_RAD*0.8
  ;; SN_BOX_RAD used only for estimating SN ratio core only 1.5"
  lbl = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
  
  dims = size(img, /dimen)
  ny = dims[1]
  dimg = size(gdtrc, /dim)
  norders = dimg[1]
  
  ;; Set standard
  x_std = gdtrc[0:ny-1L, *]
  y_std = rebin(findgen(ny), ny, norders)
  
  ;; Mask the bad trace (and off edge of CCD)
  trcmask = lonarr(ny, norders) + 1L
  ;; Good regions for each order

  ymin_ord7 = 1100L
  ymin_ord6 = 150L
  ymin_ord5 = 1300L
  ymax_ord5 = 1540L
  ymax_ord4 = 1940L

  trcmask[0:ymin_ord7, 0] = 0
  trcmask[0:ymin_ord6, 1] = 0
  trcmask[ymin_ord5:ymax_ord5, 2] = 0
  trcmask[ymax_ord4:*, 3] = 0
  ;stop
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
  objstruct = long_objfind(img, tset_slits = tset_slits, invvar = sciivar $
                           , ncoeff = ncoeff, skymask = skymask $
                           , objmask = objmask $
                           , PEAKTHRESH = peakthresh, FWHM = FWHM $
                           , NPERSLIT = NFIND, ABSTHRESH = ABSTHRESH)
  traceset2xy, tset_slits[0], yy1, ledg
  traceset2xy, tset_slits[1], yy2, redg
  ;; Error checking
  nobj_str = n_elements(objstruct)
  msk = intarr(nobj_str)
  sn_arr=fltarr(nobj_str)
  FOR iobj = 0L, nobj_str-1 DO BEGIN
     ;; Slit
     qq = objstruct[iobj].slitid - 1
     ;; Check RHS/LHS
     bd = where(objstruct[iobj].xpos LT slit_edg[*, qq, 0] OR $
                objstruct[iobj].xpos GT slit_edg[*, qq, 1] AND $
                trcmask[*, qq] EQ 1, nbad)
     if nbad NE 0 then begin
        print, 'esi_echjoefind: Bad trace in slit ';, is
        msk[iobj] = 1 
     endif
  ENDFOR
  rnd_edg = round(slit_edg)
  var = (sciivar GT 0.0)/(sciivar + (sciivar EQ 0.0))
  ordermask = tspec_ordermask(tset_slits, order_vec = order_vec) 
  ;; Compute S/N per object and expunge very low S/N detections
  FOR iobj = 0L, nobj_str-1L DO BEGIN
     qq = objstruct[iobj].slitid - 1
     ;; extract asymmetric boxcars for S/N calculation
     left_edge  = (objstruct[iobj].xpos - sn_box_rad/cbin) $
                  > (rnd_edg[*, qq, 0]+LSLITE)
     right_edge = (objstruct[iobj].xpos + sn_box_rad/cbin) $
                  < (rnd_edg[*, qq, 1]-RSLITE)
     left_shift = (objstruct[iobj].xpos + obj_shift - sn_box_rad/cbin) $
                  > (rnd_edg[*, qq, 0]+LSLITE)
     right_shift = (objstruct[iobj].xpos + obj_shift + sn_box_rad/cbin) $
                   < (rnd_edg[*, qq, 1]-RSLITE)
     snmask = ordermask EQ  order_vec[qq] AND sciivar GT 0.0 AND waveimg GT 0.0
     ;; boxcar extraction
     fx    = extract_asymbox2(img*snmask, left_edge, right_edge)
     fvar  = extract_asymbox2(var*snmask, left_edge, right_edge)
     pixtot = extract_asymbox2(float(var*0 + 1), left_edge, right_edge)
     mask_box = extract_asymbox2(float(snmask EQ 0) $
                                 , left_edge, right_edge) NE pixtot
     fi = mask_box/(fvar + (fvar EQ 0))
     box_denom = extract_asymbox2((waveimg GT 0.0) $
                                  , left_edge, right_edge)
     wave  = extract_asymbox2(waveimg, left_edge, right_edge) $
             /(box_denom + (box_denom EQ 0))
     ;; boxcar shifted
     fx_shf    = extract_asymbox2(img*snmask, left_shift, right_shift)
     fvar_shf  = extract_asymbox2(var*snmask, left_shift, right_shift)
     pixtot_shf = extract_asymbox2(float(var*0 + 1), left_shift, right_shift)
     mask_box_shf = extract_asymbox2(float(snmask EQ 0) $
                                     , left_shift, right_shift) NE pixtot_shf
     fi_shf = mask_box_shf/(fvar_shf + (fvar_shf EQ 0))
     ;; evaluate sn 
     ;; median filter the resulting fx and fi to remove hot pixels
     flux_sm = djs_median(fx, width = 5, boundary = 'reflect')
     fluxivar_sm =  djs_median(fi, width = 5, boundary = 'reflect')
     fluxivar_sm = fluxivar_sm*(fi GT 0.0)
     indsp = WHERE(wave GT 0.94000 AND wave LT 2.5 $
                   AND finite(flux_sm) AND flux_sm LT 5.0d5 $
                   AND flux_sm GT -1000.0d $
                   AND fluxivar_sm GT 0.0, nsp)
     sn = flux_sm[indsp]*sqrt(fluxivar_sm[indsp] >  0) ;;*bmask2
     djs_iterstat, sn,median=med_sn,mean=mean_sn
     ;;sn_mode = 3.0*med_sn - 2.0*mean_sn
     ;; evaulate sn for shifted aperture
     flux_sm_shf = djs_median(fx_shf, width = 5, boundary = 'reflect')
     fluxivar_sm_shf =  djs_median(fi_shf, width = 5, boundary = 'reflect')
     fluxivar_sm_shf = fluxivar_sm_shf*(fi_shf GT 0.0)
     indsp_shf = WHERE(wave GT 0.9400 AND wave LT 2.5 $
                       AND finite(flux_sm_shf) AND flux_sm_shf LT 5.0d5 $
                       AND flux_sm_shf GT -1000.0d $
                       AND fluxivar_sm_shf GT 0.0, nsp_shf)
     if (nsp_shf GT 0) then begin
        sn_shf = flux_sm_shf[indsp_shf]*sqrt(fluxivar_sm_shf[indsp_shf] >  0) 
        djs_iterstat, sn_shf,median=med_sn_shf,mean=mean_sn_shf
        ;;sn_mode_shf = 3.0*med_sn_shf - 2.0*mean_sn_shf
        sn_arr[iobj]=med_sn - med_sn_shf
     endif else begin
        sn_arr[iobj]=0
     endelse
  ENDFOR
  ;; Expunge very low S/N ratio objects
  bad = where(sn_arr LT MIN_SN,nbad_sn)
  IF nbad_sn NE 0 THEN BEGIN
     msk[bad] = 1
     print, 'tspec_echjoefind: S/N ratio too low in orders:' 
     forprint, '                                 ' + $
               string(order_vec[objstruct[bad].slitid - 1]), textout = 1
  ENDIF
  gdobj = where(msk EQ 0, ngdo)
  if ngdo EQ 0 then begin
     print, 'tspec_echjoefind:  No objects found in any order'
     print, 'tspec_echjoefind:  Giving up on this exposure.'
     return, -1
  endif

  objstruct = objstruct[gdobj]
  nobj_str = n_elements(objstruct)
  ;; Identify unique objects by running FOF on the relative 
  ;; distance from the standard star trace. 
  std_dist=fltarr(nobj_str)
  FOR iobj=0L,nobj_str-1L DO BEGIN
     qq=objstruct[iobj].slitid-1L
     igood=where(trcmask[*, qq] EQ 1)
     std_dist[iobj]=djs_median(objstruct[iobj].xpos[igood] - $
                               x_std[igood,qq])*0.3D
  ENDFOR
  x = std_dist/1000.0D
  y = 0.0*x
  group = spheregroup(x, y, STDSEP/1000.0D, mult = mult)
  ;; Require that objects occur on at least two orders to be considered. 
  ;; If you are interested in emission line sources with no detectable 
  ;; continuum this should be changed. 
  ;;keep = WHERE(mult[group] GE 2) 
  ;;objstruct = objstruct[keep]
  ;;group = group[keep]
  objstruct.objid = lindgen(nobj_str) + 1
  uni_ind = uniq(group, sort(group))
  uni_group = group[uni_ind]
  nobj = n_elements(uni_group)
  gdist = fltarr(nobj_str)
  FOR j = 0L, nobj-1L DO BEGIN
     ind = WHERE(group EQ uni_group[j])
     this_dist = djs_median(std_dist[ind])
     gdist[ind] = this_dist
  ENDFOR
  uni_dist = gdist[uni_ind]
  ;; Fill in traces
  for kk=norders,1,-1 do begin
      for ss=0L,nobj-1 do begin
          ;; Is it there?
          mt = where(group EQ uni_group[ss] AND objstruct.slitid EQ kk, nmt)
          if nmt EQ 0 then begin ;; Nope
             print, 'tspec_echjoefind: Extrapolating object ', lbl[ss], ' in '+$
                    'order ',20-(kk-1)
             newpos = x_std[*,kk-1] + uni_dist[ss]/0.3D
             ;; Fill up new objstruct
             tmp = objstruct[0]
             tmp.objid = n_elements(objstruct)+1
             tmp.slitid = kk
             tmp.xfracpos = djs_median((newpos - ledg[*,kk-1L])/ $
                                       (redg[*,kk-1L]-ledg[*,kk-1L]))
             tmp.xpos = newpos
             ;; find closest slit for FWHM
             this_group = WHERE(group EQ uni_group[ss])
             mn = min(abs(objstruct[this_group].slitid-kk), imn)
             fwhm_now=objstruct[this_group[imn]].FWHM
             tmp.fwhm=fwhm_now
             tmp.maskwidth=3.*fwhm_now
             objstruct = [objstruct, tmp]
             group = [group, uni_group[ss]]
             gdist = [gdist, uni_dist[ss]]
          endif
      endfor
  endfor
  ;; Some code to make sure the objstruct is sorted by fractional 
  ;; position and by slit, respectively. 
  ind_sort = sort(uni_frac)
  FOR k = 0L, nobj-1L DO BEGIN
     this_dist = WHERE(gdist EQ uni_dist[k])
     this_obj = objstruct[this_dist]
     this_obj = this_obj[sort(this_obj.SLITID)]
     IF NOT KEYWORD_SET(OBJOUT) THEN objout = [this_obj] $
     ELSE objout = [objout, this_obj]
  ENDFOR

  ;; Did we find all the objects the user desired?
  if keyword_set(NFIND) then begin
      if NOBJ LT NFIND then begin
          print, 'tspec_echjoefind:  Found fewer objects than you desired.', $
                 nfind, nobj
          print, 'tspec_echjoefind:  Consider setting MINPEAK to a value lower'
          print, 'tspec_echjoefind:  than ', MINPEAK
          nfind = nobj
      endif
  endif
  
  return, objout
end
