;+ 
; NAME:
; esi_echfltsct   
;     Version 1.1
;
; PURPOSE:
;    Subtracts a scattered light model from the flat
;    Normalizes the Image and Outputs to same file with
;    updated headers
;
; CALLING SEQUENCE:
;   
;  esi_echfltsct, esi, /IFLAT, GAIN_RTO=, /CHK
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  Image with scattered light removed (e.g. Flats/FlatECH##_D[I].fits)
;
; OPTIONAL KEYWORDS:
;   /IFLAT    - Use internal flats 
;   GAIN_RTO= - Value for gain ratio for two amp mode (default:
;      calculate from image)
;   /CHK      - Display final image and a few other steps
;   MINWID=   - Minimum width between orders for measuring scattered light
;   /PINH     - Indicates pinhole slit (used in early 2000)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning -- Many 'hard wired' numbers
;
; EXAMPLES:
;   esi_echfltsct, esi, 0.75, /CHK
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfltsct, esi, slit, IFLAT=iflat, GAIN_RTO=gain_rto, CHK=chk, $
                   MINWID=minwid, PINH=pinh

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echfltsct, esi, slit, /IFLAT, /CHK, GAIN_RTO=, MINWID= '
      print, '       /PINH [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if keyword_set(IFLAT) then ftype = 'I' else ftype = 'D'
  if not keyword_set(SMSHROW) then smshrow = 2150L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( MINWID ) then minwid = 6.
  if not keyword_set( NRMMED ) then NRMMED = 20L
  if not keyword_set( MED_SMTH ) then med_smth = 5L
  if not keyword_set( XSEP ) then xsep = 2L

; Slit

  c_s = esi_slitnm(slit)

; Open Flat

  flatfil = 'Flats/FlatECH'+c_s+'_'+ftype+'.fits'
          
  print, 'esi_echfltsct: Opening flat ', flatfil
  if x_chkfil(flatfil+'*') EQ 0 then begin
      print, 'esi_echfltsct: File ', flatfil, ' does not exist.  Returning!'
      return
  endif
  flat = xmrdfits(flatfil, 0, head, /silent)
  sz_flat = size(flat, /dimensions)

; Saw

  print, 'esi_echfltsct: Finding slit edges...'
  saw = shift(flat,1) - shift(flat,-1)
  smsh = djs_median(saw[*,smshrow-5:smshrow+5],2)

;;;;;;;;
;  Find Slit edges
  guess_edg =  [ [129., 308], $
                 [401, 569], $
                 [641, 802], $
                 [856, 1010], $
                 [1049, 1196], $
                 [1223, 1365], $
                 [1382, 1520], $
                 [1531, 1665], $
                 [1675, 1806], $
                 [1825, 1952] $
               ]
  if keyword_set(PINH) then begin
      for i=0,9 do guess_edg[0,i] = guess_edg[0,i] + 13.
  endif
  fin_edg = guess_edg

  ;; Right edges
  x_fndpeaks, smsh, cen_pos, NSIG=5., /silent
  for i=0,9 do begin
      mn = where(abs(cen_pos-guess_edg[1,i]) LT 15, imn)
      if imn NE 0 then begin
          mx = max(smsh(cen_pos[mn]), imx)
          fin_edg[1,i] = cen_pos[mn[imx]]
      endif else fin_edg[1,i] = guess_edg[1,i]
  endfor
;  x_splot, smsh, xtwo=fin_edg[1,*], ytwo=fltarr(10)+100., /block, psym_y2=1

  ;; Left edges
  x_fndpeaks, -(smsh), cen_neg, NSIG=5., /silent
  for i=0,9 do begin
      mn = where(abs(cen_neg-guess_edg[0,i]) LT 15, imn)
      if imn NE 0 then begin
          mx = max(-smsh(cen_neg[mn]), imx)
          fin_edg[0,i] = cen_neg[mn[imx]]
      endif else fin_edg[0,i] = guess_edg[0,i]
  endfor
;  x_splot, -smsh, xtwo=fin_edg[0,*], ytwo=fltarr(10)+100., /block, psym_y2=1

  ;; Trace out the slit edges
  restore, fitfil
  slit_edg = fltarr(sz_flat[1],10,2)
  for i=0,9 do begin
      for j=0,1 do begin
          slit_edg[*,i,j] = esi_echgettrc(fin_edg[j,i],float(smshrow), $
                                          fin_fit)
      endfor
  endfor

  ;; CHK (Plot slit edges)
  if keyword_set( CHK ) then begin
      tmp = flat
      for i=0,9 do begin
          for j=0,1 do begin
              rnd_trc2 = round(slit_edg[*,i,j])
              trc_msk = rnd_trc2 + lindgen(sz_flat[1])*sz_flat[0]
              tmp[trc_msk] = -10000
          endfor
      endfor
      xatv, tmp, /block
  endif

  sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
  print, 'esi_echfltsct: Writing slit edge info to: ', sedg_fil
  mwrfits, slit_edg, sedg_fil, /create


  print, 'esi_echfltsct: Normalizing the Gain'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Normalize Gain
  if not keyword_set( GAIN_RTO ) then begin
      ;; Centers
      cen_4 = (slit_edg[*,4,0]+slit_edg[*,4,1])/2.
      ;; Running median
      lft4 = round(cen_4)-7L
      rgt4 = round(cen_4)+7L
      arr4 = fltarr(sz_flat[1],15)
      for i=0L,sz_flat[1]-1 do arr4[i,*] = flat[lft4[i]:rgt4[i],i]
      ;; Smash
      smsh4 = djs_median(arr4,2)
      ;; Fit RHS
      gainr_fit = x_setfitstrct(NITER=1L, NORD=2L, FLGREJ=1L, HSIG=3., LSIG=3., $
                                FUNC='POLY')
      gainl_fit = gainr_fit
      
      
      fitrhs = x_fitrej(2700L+findgen(221),smsh4[2700L:2920L], FITSTR=gainr_fit)
      fitlhs = x_fitrej(3040L+findgen(361),smsh4[3040L:3400L], FITSTR=gainl_fit)
      mn = min(abs(cen_4[1500L:sz_flat[1]-1]-1022.5), imn)
      
      ;; CORRECT RHS
      val_rhs = x_calcfit(float(imn+1500L), FITSTR=gainr_fit)
      val_lhs = x_calcfit(float(imn+1500L), FITSTR=gainl_fit)
      gain_rto = val_lhs/val_rhs
      print, 'esi_echfltsct: Gain Ratio LHS/RHS = ', gain_rto
  endif
  sxaddpar, head, 'GAINFIX', gain_rto
  flat[1023:sz_flat[0]-1,*] = flat[1023:sz_flat[0]-1,*]*gain_rto

;;;;;;;;;;;;;;;;;;;;
; Mask Good pixels
  print, 'esi_echfltsct: Mask+Median'
  fit_val = fltarr(4096, 11, 2)
  msk = bytarr(4096, 11)
  
  ;; PINHOLE (circ early 2000)
  if keyword_set(PINH) then off = 15. else off = 0L

  ;; Left edge
  lft = round(slit_edg[*,0,0] - 10.) > 1L
  rgt = round(slit_edg[*,0,0] - 5.)
  a = where(rgt GE 1L, na)
  msk[a,0] = 1
  fit_val[a,0,0] = (lft[a]+rgt[a])/2.
  for i=0L,na-1 do fit_val[a[i],0,1] = median(flat[lft[a[i]]:rgt[a[i]],a[i]])
  ;; Right edge
  lft = round(slit_edg[*,9,1] + 5. + off) 
  rgt = round(slit_edg[*,9,1] + 10. + off) < (sz_flat[0]-1)
  a = where(lft LE (sz_flat[0]-1), na)
  msk[a,10] = 1
  fit_val[a,10,0] = (lft[a]+rgt[a])/2.
  for i=0L,na-1 do fit_val[a[i],10,1] = median(flat[lft[a[i]]:rgt[a[i]],a[i]])

  ;; Between slits
  for q=0L,8 do begin
      sep = round(slit_edg[*,q+1,0] - slit_edg[*,q,1] - off)
      cen = (slit_edg[*,q+1,0] + slit_edg[*,q,1] + off)/2.
      fit_val[*,q+1,0] = cen
      ;; Mask
      a = where(sep GT MINWID, na)
      msk[a,q+1] = 1
      ;; Median
      lft = round(slit_edg[a,q,1]+minwid/2.+off) > round(cen[a]-3.)
      rgt = round(slit_edg[a,q+1,0]-minwid/2.) < round(cen[a]+3.)
      for i=0L,na-1 do fit_val[a[i],q+1,1] = median(flat[lft[i]:rgt[i],a[i]])
  endfor


  ;; Additional mask stuff (e.g. by hand)
  msk[3800:3935L,1] = 0   ; Bad spot
  msk[3460:3680L,2] = 0   ; Bad spot
  msk[3850:4000L,5] = 0   ; Bad columns
  msk[3000:sz_flat[1]-1,6] = 0
  msk[2200:sz_flat[1]-1,7:8] = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PATCH ENDS 
  fit_ptch = x_setfitstrct(NITER=1L, NORD=1L, FLGREJ=1L, HSIG=3., LSIG=3., $
                           FUNC='POLY')
  ;; Match 5 to 4 (linear fit)
  fitmsk = lonarr(491) + 1L
  fitmsk[250:400]=0L
  fit45 = x_fitrej(3600+findgen(250), fit_val[3600:3849,5,1]/fit_val[3600:3849,4,1],$
                   MSK=fitmsk, FITSTR=fit_ptch)
  a = where(msk[*,5] EQ 0, na)
  if na NE 0 then fit_val[a,5,1] = fit_val[a,4,1]*x_calcfit(float(a),FITSTR=fit_ptch)
  ;; Match 6 to 5 (linear fit)
  fit_ptch = x_setfitstrct(NITER=1L, NORD=1L, FLGREJ=1L, HSIG=3., LSIG=3., $
                           FUNC='POLY')
  fit56 = x_fitrej(2600+findgen(400), fit_val[2600:2999,6,1]/fit_val[2600:2999,5,1],$
                   FITSTR=fit_ptch)
  a = where(msk[*,6] EQ 0, na)
  if na NE 0 then fit_val[a,6,1] = fit_val[a,5,1]*x_calcfit(float(a),FITSTR=fit_ptch)
  msk[a,6] = 1
  ;; Match 7 to 6 (assume constant)
  rto_67 = median( fit_val[1800:2199,7,1] / fit_val[1800:2199,6,1] )
  a = where(msk[*,7] EQ 0, na)
  msk[a,7] = 1
  if na NE 0 then fit_val[a,7,1] = fit_val[a,6,1]*rto_67
  ;; Match 8 to 7 (linear fit)
  fit_ptch = x_setfitstrct(NITER=1L, NORD=1L, FLGREJ=1L, HSIG=3., LSIG=3., $
                           FUNC='POLY')
  fit87 = x_fitrej(1900+findgen(300), fit_val[1900:2199,8,1]/fit_val[1900:2199,7,1],$
                   FITSTR=fit_ptch)
  a = where(msk[*,8] EQ 0, na)
  if na NE 0 then fit_val[a,8,1] = fit_val[a,7,1]*x_calcfit(float(a),FITSTR=fit_ptch)
  msk[a,8] = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FIT (e.g. SMOOTH)
  for q=0,10L do begin
      fit_smth = x_setfitstrct(NITER=1L, NORD=23L, FLGREJ=1L, HSIG=3., LSIG=3., $
                               FUNC='BSPLIN')
      a = where(msk[*,q] NE 0)
      fit = x_fitrej(float(a), fit_val[a,q,1], FITSTR=fit_smth)
      if keyword_set( CHK ) then begin
          x_splot, float(a), fit_val[a,q,1], $
            YTWO=x_calcfit(float(a),fitstr=fit_smth), /block, PSYM1=1
      endif
      fit_val[a,q,1] = fit
      print, 'SMOOTH RMS ', fit_smth.rms
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP + SPLINE
  print, 'esi_echfltsct: Splining and subtracting scattered light...' 
  clr = getcolor(/load)
  img_scatt = findgen(sz_flat[0], sz_flat[1])
  dumx = findgen(sz_flat[0])
  for j=0L,sz_flat[1]-1 do begin
      a = where(msk[j,*] NE 0)
      splin = spl_init(fit_val[j,a,0], fit_val[j,a,1], /double)
      img_scatt[*,j] = spl_interp(fit_val[j,a,0], fit_val[j,a,1], $
                                  splin, dumx, /double)
  endfor

  if keyword_set( CHK) then xatv, img_scatt, /block
  newflt = flat - img_scatt
  sxaddpar, head, 'SCATTER', 'T'
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NORMALIZE

  print, 'esi_echfltsct: Normalizing.. ' 
  ;; Stuff
  img_nrm = fltarr(sz_flat[0],sz_flat[1]) 
  med_img = fltarr(sz_flat[1], NRMMED*2+1L) 
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  fit_nrm = x_setfitstrct(NITER=1L, FLGREJ=1L, HSIG=3., LSIG=3.)
                           

  ;; LOOP to FIT along the ORDER
  for q=0L,9 do begin
      ;; Deal with off the edge order
      if q EQ 0 then begin
          b = where(msk[*,0] NE 0)
          frstj = min(b)
          lstj = 3700L
          fit_nrm.func = 'POLY'
          fit_nrm.nord = 6L
      endif else begin
          frstj = 0L
          lstj = sz_flat[1]-1
          fit_nrm.func = 'BSPLIN'
          fit_nrm.nord = 13L
      endelse
      ;; Take median along slit centers
      for j=frstj,lstj do $
        med_img[j,*] = newflt[(slit_cen[j,q]-NRMMED)>0L:slit_cen[j,q]+NRMMED,j]
      med_nrm = djs_median(med_img, 2)
      ;; FIT
      fit = x_fitrej(frstj+findgen(lstj-frstj+1),med_nrm[frstj:lstj], $
                     FITSTR=fit_nrm)
      ;; LEFT and RIGHT
      if q EQ 0 then lft = (slit_edg[*,q,0]-xsep) > 0L $
      else lft = (slit_edg[*,q-1,1]+xsep) > (slit_edg[*,q,0]-xsep) < (sz_flat[0]-1)
      if q EQ 9 then rgt = lft > (slit_edg[*,q,1]+xsep) < (sz_flat[0]-1) $
      else rgt = lft > (slit_edg[*,q,1]+xsep) < (slit_edg[*,q+1,0]-xsep)
      ;; MAP
      for j=frstj,lstj do img_nrm[lft[j]:rgt[j],j] = fit[j-frstj]
  endfor

  ;; DIVIDE
  a = where(img_nrm NE 0)
  nrm_flat = fltarr(sz_flat[0],sz_flat[1])
  nrm_flat[a] = newflt[a]/img_nrm[a]

  delvarx, img_nrm
  img_nrm = fltarr(sz_flat[0],sz_flat[1]) 
  ;; LOOP to FIT 'spatial' variation in Internal flat
  fit_spa = x_setfitstrct(NITER=1L, FLGREJ=1L, HSIG=3., LSIG=3.)
  fit_spa.func = 'POLY'
  fit_spa.nord = 3L
  for q=0L,9 do begin
      coll_img = fltarr(round(slit_edg[2400L,q,1])- $
                        round(slit_edg[2400L,q,0])+1, 201L)
      sz_coll = size(coll_img, /dimensions)
      for j=2300L,2500L do begin
          i1 = round(slit_edg[j,q,0])
          i2 = i1 + sz_coll[0] - 1
          coll_img[*,j-2300L] = nrm_flat[i1:i2,j]
      endfor

      ;; Take median along slit centers
      med_coll = djs_median(coll_img, 2)
      ;; FIT
      fit = x_fitrej(29.+findgen(sz_coll[0]-5-29), $
                     med_coll[30:sz_coll[0]-5], FITSTR=fit_spa)
      ;; Calculate the FIT at a large range and NOT at all for q>6
      if q LT 7 then $
        dumf = x_calcfit(findgen(350), FITSTR=fit_spa)

      ;; MAP
      if q EQ 0 then begin
          b = where(msk[*,0] NE 0)
          frstj = min(b)
          lstj = 3800L
      endif else begin
          frstj = 0L
          lstj = sz_flat[1]-1
      endelse

      lft = round(slit_edg[*,q,0]) > 0L
      rgt = round(slit_edg[*,q,1]) > 0L
      rnd_lft = round(slit_edg[*,q,0])
      for j=frstj,lstj do begin
          img_nrm[lft[j]:rgt[j],j] = dumf[lft[j]-rnd_lft[j]:rgt[j]-rnd_lft[j]]
      endfor
  endfor

  ;; FINAL DIVIDE
  a = where(img_nrm NE 0)
  nrm_flat[a] = nrm_flat[a]/img_nrm[a]
      
  if keyword_set( CHK ) then xatv, nrm_flat, min=0.9, max=1.1, /block

  sxaddpar, head, 'NORM', 'T'
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; OUTPUT
  outfil = 'Flats/FlatECH'+c_s+'N.fits'
  print, 'esi_echfltsct: Creating Normalized Flat ', outfil
  mwrfits, nrm_flat, outfil, head, /create, /silent
  ;; COMPRESS
  spawn, 'gzip -f '+outfil
  ;; Update header cards
  objstd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                 esi.slit EQ slit AND $
                 (strtrim(esi.type,2) EQ 'STD' OR $
                  strtrim(esi.type,2) EQ 'OBJ' OR $ 
                  strtrim(esi.type,2) EQ 'ARC' ), nobj)
  if nobj NE 0 then esi[objstd].flat_fil = outfil

  print, 'esi_echfltsct: All done!'

  return 
end

