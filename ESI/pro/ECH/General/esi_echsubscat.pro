;+ 
; NAME:
; esi_echsubscat   
;     Version 1.0
;
; PURPOSE:
;    Subtracts a scattered light model from the image
;
; CALLING SEQUENCE:
;   
;  esi_echsubscat, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with scattered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Not well tested
;
; EXAMPLES:
;   esi_echsubscat, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echsubscat, esi, indx, DFLAT=dflat, GAIN_RTO=gain_rto, CHK=chk, $
                    REDOOV=redoov, GAINONLY=gainonly, FLATFIL=flatfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echsubscat, esi, indx, /DFLAT, /CHK, GAIN_RTO= [v1.0]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(SMSHROW) then smshrow = 2150L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( MINWID ) then minwid = 6.
  if not keyword_set( MED_SMTH ) then med_smth = 3L
  if not keyword_set( XSEP ) then xsep = 2L

; NIMG
  nimg = n_elements(indx)

; Slit
  slit = esi[indx[0]].slit  ; Takes first value so only 1 slit width at a time
  c_s = esi_slitnm(slit)

; Open Flat
  if not keyword_set(FLATFIL) then $
    flatfil = strtrim(esi[indx[0]].flat_fil,2)
  if strlen(flatfil) EQ 0 then begin
     print, '    *********************************'
     print, 'esi_echsubscat: Flat file not set!'
;     print, 'esi_echsubscat: Assuming a dummy name!!'
     print, '    *********************************'
;     flatfil = 'Dummy_filename'
;      stop
;      return
  endif
  IF KEYWORD_SET(FLATFIL) AND file_test(flatfil+'*') THEN  $
    headflt = xheadfits(flatfil, /silent) $
  ELSE BEGIN
      ;; JFH added 04/08. Kludge the header of the flat for the 
      ;; case where the flat doesn't exist yet (i.e. for arcs). 
      headflt = strarr(1)
      sxaddpar, headflt, 'SCATTER', 'T'
      sxaddpar, headflt, 'GAINFIX', 0.947916
  ENDELSE
        
  scatt = sxpar(headflt, 'SCATTER', COUNT = nval)
  if scatt NE 1 then begin
      print, 'esi_echsubscat: Flat not created! Returning...'
      return
  endif
          
; Slit Edges

  if not keyword_set( GAINONLY ) then begin
      sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
      print, 'esi_echsubscat: Grabbing slit edges from: ', sedg_fil
      slit_edg = xmrdfits(sedg_fil, /silent)
  endif


; Gain Fix
  if not keyword_set( GAIN_RTO ) then begin
      ;; Grab value from Flat header
      gain_rto = sxpar(headflt,'GAINFIX',count=nval)
      if nval EQ 0 then begin
          print, 'esi_echsubscat: GAINFIX not measured! Returning..'
          return
      endif
  endif

; BIAS sub if necessary
  if not keyword_set( REDOOV ) then bias = where(esi[indx].flg_ov EQ 0, nbias) $
  else begin
      nbias = n_elements(indx)
      bias = lindgen(nbias)
  endelse
  if nbias NE 0 then esi_subbias, esi, indx[bias], FORCE=REDOOV

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  for qq=0L,nimg-1 do begin
      ;; Read IMG
      img_fil = 'OV/ov_'+esi[indx[qq]].img_root
      img = xmrdfits(img_fil, 0, head, /silent)
      scatt = sxpar(head,'SCATTER',COUNT=nval)
      if scatt EQ 1 then begin
          print, 'esi_echsubscat: IMG ', img_fil, ' already subtracted! Continuing'
          continue
      endif
      sz_img = size(img, /dimensions)
  
      ;;;;;;;;;
      ;; Normalize Gain
      print, 'esi_echsubscat: Normalizing the Gain'
      img[1023:sz_img[0]-1,*] = img[1023:sz_img[0]-1,*]*gain_rto
      if keyword_set( GAINONLY ) then begin
          esi[indx[qq]].flg_ov = 1
          mwrfits, img, img_fil, head, /create
          continue
      endif

      ;;;;;;;;;;;;;;;;;;;;
      ;; ID PIX Between Slits + Median
      print, 'esi_echsubscat: Mask+Median'
      fit_val = fltarr(4096, 11, 2)
      msk = bytarr(4096, 11)
      
      ;; Left edge
      lft = round(slit_edg[*,0,0] - 10.) > 1L
      rgt = round(slit_edg[*,0,0] - 5.)
      a = where(rgt GE 1L, na)
      msk[a,0] = 1
      fit_val[a,0,0] = (lft[a]+rgt[a])/2.
      for i=0L,na-1 do fit_val[a[i],0,1] = median(img[lft[a[i]]:rgt[a[i]],a[i]])
      ;; Right edge
      lft = round(slit_edg[*,9,1] + 15.) 
      rgt = round(slit_edg[*,9,1] + 20.) < (sz_img[0]-1)
      a = where(lft LE (sz_img[0]-1), na)
      msk[a,10] = 1
      fit_val[a,10,0] = (lft[a]+rgt[a])/2.
      for i=0L,na-1 do fit_val[a[i],10,1] = median(img[lft[a[i]]:rgt[a[i]],a[i]])
      
      ;; Between slits
      for q=0L,8 do begin
          sep = round(slit_edg[*,q+1,0] - slit_edg[*,q,1])
          cen = (slit_edg[*,q+1,0] + slit_edg[*,q,1])/2.
          fit_val[*,q+1,0] = cen
          ;; Mask
          a = where(sep GT MINWID, na)
          msk[a,q+1] = 1
          ;; Median
          lft = round(slit_edg[a,q,1]+minwid/2.) > round(cen[a]-3.)
          rgt = round(slit_edg[a,q+1,0]-minwid/2.) < round(cen[a]+3.)
          for i=0L,na-1 do fit_val[a[i],q+1,1] = $
            median(img[lft[a[i]]:rgt[a[i]],a[i]])
      endfor

      ;; Additional mask stuff (e.g. by hand)
      msk[3200:sz_img[1]-1,6] = 0
      msk[2400:sz_img[1]-1,7:8] = 0
      msk[3765:3995,1] = 0  ; Bad square
      msk[3490:3589,2] = 0  ; Bad columns

      ;;;;;;;;;;;;;
      ;; PATCH ENDS 
      fit_ptch = x_setfitstrct(NITER=1L, NORD=1L, FLGREJ=1L, HSIG=3., LSIG=3., $
                               FUNC='POLY')
      ;; Match 6 to 5 (linear fit)
      fit56 = x_fitrej(2400+findgen(801), $
                       fit_val[2400:3200,6,1]/fit_val[2400:3200,5,1],$
                       FITSTR=fit_ptch)
      a = where(msk[*,6] EQ 0, na)
      if na NE 0 then fit_val[a,6,1] = $
        fit_val[a,5,1]*x_calcfit(float(a),FITSTR=fit_ptch)
      msk[a,6] = 1
      ;; Match 7 to 6 (assume constant)
      rto_67 = median( fit_val[1800:2300,7,1] / fit_val[1800:2300,6,1] )
      a = where(msk[*,7] EQ 0, na)
      msk[a,7] = 1
      if na NE 0 then fit_val[a,7,1] = fit_val[a,6,1]*rto_67
      ;; Match 8 to 7 (linear fit)
      fit87 = x_fitrej(1800+findgen(601), $
                       fit_val[1800:2400,8,1]/fit_val[1800:2400,7,1],$
                       FITSTR=fit_ptch)
      a = where(msk[*,8] EQ 0, na)
      if na NE 0 then fit_val[a,8,1] = $
        fit_val[a,7,1]*x_calcfit(float(a),FITSTR=fit_ptch)
      msk[a,8] = 1


      ;; CHK
      if keyword_set( CHK ) then begin
          xatv, fit_val[*,*,1], /block
          stop
      endif

      ;;;;;;;;;;;;;;
      ;; FIT (e.g. SMOOTH)
      fit_smth = x_setfitstrct(NITER=1L, NORD=7L, FLGREJ=1L, HSIG=3., LSIG=3., $
                               FUNC='BSPLIN')
      for q=0,10L do begin
          a = where(msk[*,q] NE 0)
          fit = x_fitrej(float(a), fit_val[a,q,1], FITSTR=fit_smth)
          if q EQ 0 then fit_val[a,q,1] = fit else $
            fit_val[*,q,1] = x_calcfit(findgen(sz_img[1]), FITSTR=fit_smth) 
      endfor

      ;;;;;;;;;;;;;;;;;;;;
      ;; LOOP + SPLINE
      print, 'esi_echsubscat: Splining and subtracting scattered light...' 
      img_scatt = findgen(sz_img[0], sz_img[1])
      dumx = findgen(sz_img[0])
      for j=0L,sz_img[1]-1 do begin
          a = where(msk[j,*] NE 0)
          splin = spl_init(fit_val[j,a,0], fit_val[j,a,1], /double)
          img_scatt[*,j] = spl_interp(fit_val[j,a,0], fit_val[j,a,1], $
                                      splin, dumx, /double)
      endfor
      if keyword_set( CHK) then xatv, img_scatt, /block

      ;; OUTPUT
      print, 'esi_echsubscat: Overwriting ', img_fil
      sxaddpar, head, 'SCATTER', 'T'
      esi[indx[qq]].flg_ov = 2
      mwrfits, img-img_scatt, img_fil, head, /create
      
  endfor

  return
end



;  newflt = flat - img_scatt
;  sxaddpar, head, 'SCATTER', 'T'
  
